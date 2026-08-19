"""Render direct JADE JSS-RTP magnetic field and full-band Waves.

All science values shown are direct fields from archived files:
  * Magnetic field: JADE Level-5 MAG_VECTOR_JSSRTP at MAG_UTC
  * JADE: official Level-5 electron-moment SC_POS_R records
  * Waves: one calibrated E-field PSD source record nearest each 10-min
    display-bin center (no average or interpolation)

No magnetic-coordinate transformation is performed.  At the user's explicit
request, |B| is calculated from the direct Br, Btheta and Bphi samples.

This renderer mirrors plot_juno_orbit7_fullorbit_layout.m.  It is retained as
a fallback for MATLAB R2025a+ systems whose WebGL helper cannot initialize
with a virtual display adapter.
"""

from __future__ import annotations

import csv
import math
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


DATABASE_ROOT = Path(r"Z:\SPART-WORK\Data\Juno")
SCRIPT_DIR = Path(__file__).resolve().parent
TIME_START = datetime(2017, 6, 27, 10, 0, 0, tzinfo=timezone.utc)
TIME_END = datetime(2017, 7, 9, 1, 54, 5, tzinfo=timezone.utc)
DISPLAY_SECONDS = 600
FREQUENCY_RANGE_HZ = (50.0, 41.0e6)

DATE_TICKS = [
    datetime(2017, 6, 29, tzinfo=timezone.utc) + timedelta(days=2 * i)
    for i in range(6)
]


def doy_code(day: datetime) -> str:
    return f"{day.year:04d}{day.timetuple().tm_yday:03d}"


def days_in_interval() -> list[datetime]:
    first = datetime(TIME_START.year, TIME_START.month, TIME_START.day,
                     tzinfo=timezone.utc)
    last = datetime(TIME_END.year, TIME_END.month, TIME_END.day,
                    tzinfo=timezone.utc)
    days: list[datetime] = []
    day = first
    while day <= last:
        days.append(day)
        day += timedelta(days=1)
    return days


def highest_version(folder: Path, pattern: str) -> Path:
    candidates = list(folder.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"No matching file in {folder}: {pattern}")

    def version(path: Path) -> int:
        match = re.search(r"_V(\d+)\.", path.name, re.IGNORECASE)
        if match is None:
            raise ValueError(f"No version in filename: {path}")
        return int(match.group(1))

    return max(candidates, key=version)


def discover_files() -> tuple[list[Path], list[Path], list[Path]]:
    days = days_in_interval()
    wave_folder = DATABASE_ROOT / "Waves_data" / "2017165_ORBIT_07"
    wave_files: list[Path] = []
    jade_rtp_files: list[Path] = []
    moment_files: list[Path] = []
    for day in days:
        code = doy_code(day)
        wave_files.append(highest_version(
            wave_folder, f"WAV_{code}T*_E_V*.CSV"))
        moment_folder = (DATABASE_ROOT / "JADE_data" /
                         f"{day.year:04d}" / code / "ELECTRONS")
        jade_rtp_files.append(highest_version(
            moment_folder, f"JAD_L50_LRS_ELC_ANY_DEF_{code}_V*.DAT"))
        moment_files.append(highest_version(
            moment_folder,
            f"JAD_L50_HLS_ELC_MOM_ISO_2D_ELECTRONS_{code}_V*.CSV"))
    return wave_files, jade_rtp_files, moment_files


def label_integer(label_text: str, keyword: str) -> int:
    match = re.search(rf"^\s*{re.escape(keyword)}\s*=\s*(\d+)",
                      label_text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"Missing {keyword} in PDS label")
    return int(match.group(1))


def label_column_start(label_text: str, column_name: str) -> int:
    pattern = (r"OBJECT\s*=\s*COLUMN\s+NAME\s*=\s*\"?" +
               re.escape(column_name) +
               r"\"?.*?START_BYTE\s*=\s*(\d+).*?" +
               r"END_OBJECT\s*=\s*COLUMN")
    match = re.search(pattern, label_text, flags=re.DOTALL)
    if match is None:
        raise ValueError(f"Missing PDS column {column_name}")
    return int(match.group(1))


def read_jade_rtp(files: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    times: list[datetime] = []
    values: list[tuple[float, float, float]] = []
    for index, path in enumerate(files, start=1):
        print(f"JADE RTP {index:2d}/{len(files):2d}: {path}")
        label_path = path.with_suffix(".LBL")
        label_text = label_path.read_text(encoding="ascii")
        record_bytes = label_integer(label_text, "RECORD_BYTES")
        record_count = label_integer(label_text, "FILE_RECORDS")
        utc_offset = label_column_start(label_text, "MAG_UTC") - 1
        rtp_offset = label_column_start(
            label_text, "MAG_VECTOR_JSSRTP") - 1
        if utc_offset + 21 > record_bytes or rtp_offset + 12 > record_bytes:
            raise ValueError(f"JADE field lies outside record: {label_path}")
        if path.stat().st_size != record_bytes * record_count:
            raise ValueError(f"JADE DAT size disagrees with label: {path}")
        dtype = np.dtype({
            "names": ("mag_utc", "mag_rtp"),
            "formats": ("S21", ("<f4", (3,))),
            "offsets": (utc_offset, rtp_offset),
            "itemsize": record_bytes,
        })
        records = np.memmap(path, mode="r", dtype=dtype,
                            shape=(record_count,))
        for raw_time, raw_rtp in zip(records["mag_utc"],
                                     records["mag_rtp"]):
            try:
                timestamp = parse_doy_utc(raw_time.decode("ascii"))
            except (UnicodeDecodeError, ValueError):
                continue
            if not (TIME_START <= timestamp < TIME_END):
                continue
            rtp = np.asarray(raw_rtp, dtype=float)
            if (not np.all(np.isfinite(rtp)) or
                    np.any(rtp == 9990000.0) or
                    np.any(np.abs(rtp) > 1.6e6)):
                rtp[:] = np.nan
            times.append(timestamp)
            values.append((float(rtp[0]), float(rtp[1]), float(rtp[2])))
    order = np.argsort(np.array([t.timestamp() for t in times]))
    sorted_times = np.asarray(times, dtype=object)[order]
    sorted_values = np.asarray(values, dtype=float)[order]
    epoch = np.asarray([item.timestamp() for item in sorted_times])
    _, reverse_unique = np.unique(epoch[::-1], return_index=True)
    keep = np.sort(len(epoch) - 1 - reverse_unique)
    return sorted_times[keep], sorted_values[keep]


def parse_doy_utc(text: str) -> datetime:
    # Example: 2017-179T14:16:31.960
    year = int(text[0:4])
    doy = int(text[5:8])
    hour = int(text[9:11])
    minute = int(text[12:14])
    second = float(text[15:])
    return (datetime(year, 1, 1, tzinfo=timezone.utc) +
            timedelta(days=doy - 1, hours=hour, minutes=minute,
                      seconds=second))


def read_jade_position(
        files: list[Path],
) -> tuple[np.ndarray, np.ndarray]:
    times: list[datetime] = []
    radius_rj: list[float] = []
    for index, path in enumerate(files, start=1):
        print(f"JADE {index:2d}/{len(files):2d}: {path}")
        with path.open("r", encoding="ascii", newline="") as handle:
            for fields in csv.reader(handle, skipinitialspace=True):
                if len(fields) != 24:
                    continue
                timestamp = parse_doy_utc(fields[0].strip())
                if not (TIME_START <= timestamp < TIME_END):
                    continue
                numeric = [float(item) for item in fields[1:]]
                radius = numeric[8]       # Direct SC_POS_R, R_J
                if radius < 0.0 or radius > 1000.0:
                    radius = math.nan
                times.append(timestamp)
                radius_rj.append(radius)
    order = np.argsort(np.array([t.timestamp() for t in times]))
    return (
        np.asarray(times, dtype=object)[order],
        np.asarray(radius_rj)[order],
    )


def parse_scet_seconds(raw: bytes, year_start_cache: dict[int, float]) -> float:
    # Fixed-width SCET example: b'2017-179T00:00:00.993'
    year = ((raw[0] - 48) * 1000 + (raw[1] - 48) * 100 +
            (raw[2] - 48) * 10 + raw[3] - 48)
    doy = (raw[5] - 48) * 100 + (raw[6] - 48) * 10 + raw[7] - 48
    hour = (raw[9] - 48) * 10 + raw[10] - 48
    minute = (raw[12] - 48) * 10 + raw[13] - 48
    second = (raw[15] - 48) * 10 + raw[16] - 48
    fraction = 0.0
    scale = 0.1
    for value in raw[18:]:
        if 48 <= value <= 57:
            fraction += (value - 48) * scale
            scale *= 0.1
    if year not in year_start_cache:
        year_start_cache[year] = datetime(
            year, 1, 1, tzinfo=timezone.utc).timestamp()
    return (year_start_cache[year] + (doy - 1) * 86400.0 +
            hour * 3600.0 + minute * 60.0 + second + fraction)


def quality_flag_set(text: str, number: int) -> bool:
    parts = text.strip('"').split(":")
    return len(parts) >= number and parts[number - 1] == "1"


def parse_wave_record(raw_line: bytes) -> np.ndarray:
    fields = raw_line.decode("ascii").rstrip("\r\n").split(",")
    if len(fields) != 154:
        raise ValueError(f"Incomplete selected Waves record: {len(fields)} fields")
    row = np.fromiter(
        (float(item) if item.strip() else math.nan
         for item in fields[28:154]),
        dtype=float, count=126,
    )
    if quality_flag_set(fields[3], 2):
        row[0:43] = np.nan
    if quality_flag_set(fields[3], 3):
        row[43:61] = np.nan
    row[row <= 0.0] = np.nan
    return row


def read_waves_direct(
        files: list[Path],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    start_seconds = TIME_START.timestamp()
    end_seconds = TIME_END.timestamp()
    n_bins = int(math.ceil((end_seconds - start_seconds) / DISPLAY_SECONDS))
    time_edges_seconds = start_seconds + np.arange(n_bins + 1) * DISPLAY_SECONDS
    bin_centers_seconds = 0.5 * (
        time_edges_seconds[:-1] + time_edges_seconds[1:])
    source_seconds = np.full(n_bins, np.nan)
    source_distance = np.full(n_bins, np.inf)
    psd = np.full((n_bins, 126), np.nan)
    frequencies: np.ndarray | None = None
    year_start_cache: dict[int, float] = {}

    for file_index, path in enumerate(files, start=1):
        print(f"Waves {file_index:2d}/{len(files):2d}: {path}")
        with path.open("rb", buffering=1024 * 1024) as handle:
            headers = [handle.readline() for _ in range(5)]
            header_fields = headers[3].decode("ascii").rstrip().split(",")
            if len(header_fields) < 154:
                raise ValueError(f"Incomplete Waves header: {path}")
            this_frequency = np.asarray(
                [float(item.strip('"')) for item in header_fields[28:154]])
            if frequencies is None:
                frequencies = this_frequency
            elif not np.allclose(this_frequency, frequencies, rtol=1e-10,
                                 atol=1e-9):
                raise ValueError(f"Waves frequency grid changed: {path}")

            active_bin = -1
            best_distance = math.inf
            best_seconds = math.nan
            best_line: bytes | None = None

            def flush_best() -> None:
                nonlocal best_line
                if active_bin < 0 or best_line is None:
                    return
                if best_distance < source_distance[active_bin]:
                    psd[active_bin, :] = parse_wave_record(best_line)
                    source_seconds[active_bin] = best_seconds
                    source_distance[active_bin] = best_distance

            for raw_line in handle:
                first_comma = raw_line.find(b",")
                if first_comma < 0:
                    continue
                second_comma = raw_line.find(b",", first_comma + 1)
                if second_comma < 0:
                    continue
                scet = raw_line[first_comma + 1:second_comma]
                if len(scet) < 21:
                    continue
                timestamp = parse_scet_seconds(scet, year_start_cache)
                if timestamp < start_seconds or timestamp >= end_seconds:
                    continue
                bin_index = int((timestamp - start_seconds) // DISPLAY_SECONDS)
                if not 0 <= bin_index < n_bins:
                    continue
                if bin_index != active_bin:
                    flush_best()
                    active_bin = bin_index
                    best_distance = math.inf
                    best_seconds = math.nan
                    best_line = None
                distance = abs(timestamp - bin_centers_seconds[bin_index])
                if distance < best_distance:
                    best_distance = distance
                    best_seconds = timestamp
                    best_line = raw_line
            flush_best()

    if frequencies is None:
        raise RuntimeError("No Waves frequency grid was read")
    receiver_channels = (range(0, 43), range(43, 61),
                         range(61, 88), range(88, 126))
    receiver_boundaries = []
    for lower, upper in zip(receiver_channels[:-1], receiver_channels[1:]):
        lower_maximum = np.max(frequencies[list(lower)])
        upper_centers = frequencies[list(upper)]
        next_upper = np.min(upper_centers[upper_centers > lower_maximum])
        receiver_boundaries.append(math.sqrt(lower_maximum * next_upper))
    order = np.argsort(frequencies)
    frequencies = frequencies[order]
    psd = psd[:, order]
    time_edges = np.asarray([
        datetime.fromtimestamp(value, tz=timezone.utc)
        for value in time_edges_seconds
    ], dtype=object)
    source_time = np.asarray([
        datetime.fromtimestamp(value, tz=timezone.utc)
        if np.isfinite(value) else None for value in source_seconds
    ], dtype=object)
    return (time_edges, source_time, psd, frequencies,
            np.asarray(receiver_boundaries))


def nearest_direct_radius(
        times: np.ndarray, radii: np.ndarray, max_offset_hours: float = 8.0,
) -> tuple[list[datetime | None], list[float]]:
    source_seconds = np.asarray([time.timestamp() for time in times])
    outputs_time: list[datetime | None] = []
    outputs_radius: list[float] = []
    for tick in DATE_TICKS:
        distance = np.abs(source_seconds - tick.timestamp())
        valid = np.isfinite(radii)
        if not np.any(valid):
            outputs_time.append(None)
            outputs_radius.append(math.nan)
            continue
        valid_indices = np.flatnonzero(valid)
        index = valid_indices[np.argmin(distance[valid])]
        if distance[index] <= max_offset_hours * 3600.0:
            outputs_time.append(times[index])
            outputs_radius.append(float(radii[index]))
        else:
            outputs_time.append(None)
            outputs_radius.append(math.nan)
    return outputs_time, outputs_radius


def frequency_edges(frequencies: np.ndarray) -> np.ndarray:
    middle = np.sqrt(frequencies[:-1] * frequencies[1:])
    return np.concatenate((
        [frequencies[0] ** 2 / middle[0]], middle,
        [frequencies[-1] ** 2 / middle[-1]],
    ))


def symmetric_log(values: np.ndarray, threshold: float = 1.0) -> np.ndarray:
    return np.sign(values) * np.log10(1.0 + np.abs(values) / threshold)


def robust_color_limits(wave_psd: np.ndarray) -> tuple[float, float]:
    values = np.log10(wave_psd[np.isfinite(wave_psd) & (wave_psd > 0.0)])
    values.sort()
    lower_index = max(0, min(len(values) - 1,
                             round(0.005 * (len(values) - 1))))
    upper_index = max(0, min(len(values) - 1,
                             round(0.997 * (len(values) - 1))))
    lower = math.floor(2.0 * values[lower_index]) / 2.0
    upper = math.ceil(2.0 * values[upper_index]) / 2.0
    if upper - lower < 5.0:
        center = 0.5 * (lower + upper)
        lower, upper = center - 2.5, center + 2.5
    return 10.0 ** lower, 10.0 ** upper


def render_figure(
        mag_time: np.ndarray, mag_b: np.ndarray,
        mag_total: np.ndarray,
        radius_tick_rj: list[float], wave_edges: np.ndarray,
        wave_psd: np.ndarray, frequencies: np.ndarray,
        receiver_boundaries: np.ndarray,
) -> tuple[Path, Path]:
    plt.rcParams.update({
        "font.size": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.facecolor": "white",
    })
    figure = plt.figure(figsize=(18.0, 12.6), facecolor="white")
    grid = figure.add_gridspec(
        2, 2, width_ratios=(1.0, 0.018), height_ratios=(1.0, 2.0),
        left=0.075, right=0.945, bottom=0.13, top=0.96,
        hspace=0.045, wspace=0.035,
    )
    ax_mag = figure.add_subplot(grid[0, 0])
    ax_wave = figure.add_subplot(grid[1, 0], sharex=ax_mag)
    ax_cbar = figure.add_subplot(grid[1, 1])

    for axis in (ax_mag, ax_wave):
        axis.grid(False)
        axis.set_xlim(TIME_START, TIME_END)
        axis.set_xticks(DATE_TICKS)
        axis.tick_params(which="both", direction="out")

    transformed_b = symmetric_log(mag_b)
    transformed_total = symmetric_log(mag_total)
    ax_mag.plot(mag_time, transformed_b[:, 0], color=(0.19, 0.39, 0.78),
                linewidth=0.72, label=r"$B_r$", zorder=2)
    ax_mag.plot(mag_time, transformed_b[:, 1], color=(0.95, 0.48, 0.12),
                linewidth=0.72, label=r"$B_\theta$", zorder=2)
    ax_mag.plot(mag_time, transformed_b[:, 2], color=(0.14, 0.55, 0.25),
                linewidth=0.72, label=r"$B_\phi$", zorder=2)
    ax_mag.plot(mag_time, transformed_total, color="black",
                linewidth=0.9, label=r"$|B|$", zorder=3)
    max_mag = float(np.nanmax(np.concatenate(
        (np.abs(mag_b).ravel(), mag_total.ravel()))))
    plot_limit = 1.04 * float(symmetric_log(np.asarray(max_mag)))
    ax_mag.set_ylim(-plot_limit, plot_limit)
    physical_ticks = 10.0 ** np.arange(0, max(0, math.floor(math.log10(max_mag))) + 1)
    tick_values = np.concatenate((-physical_ticks[::-1], [0.0], physical_ticks))
    ax_mag.set_yticks(symmetric_log(tick_values))
    ax_mag.set_yticklabels([f"{value:g}" for value in tick_values])
    ax_mag.set_ylabel(r"$B^{JSS}_{r,\theta,\phi}$ (nT)" +
                      "\nsymmetric log")
    ax_mag.legend(loc="upper right", frameon=False, ncols=4, fontsize=9,
                  handlelength=2.0, columnspacing=1.1)

    selected_frequency = frequencies
    selected_psd = wave_psd.T
    color_limits = robust_color_limits(wave_psd)
    cmap = matplotlib.colormaps["turbo"].copy()
    cmap.set_bad("#e0e0e0")
    mesh = ax_wave.pcolormesh(
        mdates.date2num(wave_edges), frequency_edges(selected_frequency),
        np.ma.masked_invalid(selected_psd), shading="flat", cmap=cmap,
        norm=LogNorm(vmin=color_limits[0], vmax=color_limits[1]),
        rasterized=True,
    )
    ax_wave.set_yscale("log")
    ax_wave.set_ylim(*FREQUENCY_RANGE_HZ)
    ax_wave.set_yticks((50.0, 1.0e2, 1.0e3, 1.0e4, 1.0e5,
                        1.0e6, 1.0e7, 4.1e7))
    ax_wave.set_yticklabels(("50 Hz", "100 Hz", "1 kHz", "10 kHz",
                             "100 kHz", "1 MHz", "10 MHz", "41 MHz"))
    ax_wave.set_ylabel("Frequency")
    for boundary in receiver_boundaries:
        ax_wave.axhline(boundary, color="black", linestyle="--",
                        linewidth=0.9)
    colorbar = figure.colorbar(mesh, cax=ax_cbar)
    colorbar.set_label(r"E PSD  [(V m$^{-1}$)$^2$ Hz$^{-1}$]")

    ax_mag.tick_params(labelbottom=False)
    ax_wave.set_xticklabels([])
    start_number = mdates.date2num(TIME_START)
    end_number = mdates.date2num(TIME_END)
    for index, (tick, radius) in enumerate(zip(DATE_TICKS, radius_tick_rj)):
        x_normalized = ((mdates.date2num(tick) - start_number) /
                        (end_number - start_number))
        alignment = "center"
        if index == 0:
            alignment = "left"
        elif index == len(DATE_TICKS) - 1:
            alignment = "right"
        ax_wave.text(x_normalized, -0.042, tick.strftime("%m/%d"),
                     transform=ax_wave.transAxes, ha=alignment, va="top",
                     fontsize=8.5, clip_on=False)
        radius_text = f"R/R$_J$ {radius:.1f}" if math.isfinite(radius) else "R/R$_J$ --"
        ax_wave.text(x_normalized, -0.088, radius_text,
                     transform=ax_wave.transAxes, ha=alignment, va="top",
                     fontsize=8.5, clip_on=False)
    ax_wave.set_xlabel("")

    base_name = ("juno_orbit7_fullorbit_layout_jade_jss_"
                 "rtheta_phi_direct_btotal_waves_fullband_"
                 "20170627_20170709")
    png_path = SCRIPT_DIR / f"{base_name}.png"
    pdf_path = SCRIPT_DIR / f"{base_name}.pdf"
    figure.savefig(png_path, dpi=260)
    figure.savefig(pdf_path, dpi=260)
    plt.close(figure)
    return png_path, pdf_path


def main() -> None:
    if not DATABASE_ROOT.is_dir():
        raise FileNotFoundError(f"Juno database is unavailable: {DATABASE_ROOT}")
    wave_files, jade_rtp_files, moment_files = discover_files()
    if not (len(wave_files) == len(jade_rtp_files) ==
            len(moment_files) == 13):
        raise RuntimeError("Expected 13 files for each source product")
    print("Reading official direct-source records ...")
    mag_time, mag_spherical = read_jade_rtp(jade_rtp_files)
    mag_total = np.sqrt(np.sum(mag_spherical ** 2, axis=1))
    position_time, position_radius = read_jade_position(moment_files)
    wave_edges, wave_source_time, wave_psd, frequencies, receiver_boundaries = \
        read_waves_direct(wave_files)
    radius_source_time, radius_tick_rj = nearest_direct_radius(
        position_time, position_radius)
    png_path, pdf_path = render_figure(
        mag_time, mag_spherical, mag_total, radius_tick_rj, wave_edges, wave_psd,
        frequencies, receiver_boundaries)

    selected_wave_bins = int(np.count_nonzero(
        np.any(np.isfinite(wave_psd), axis=1)))
    offsets = []
    centers = [TIME_START + timedelta(seconds=DISPLAY_SECONDS * (i + 0.5))
               for i in range(wave_psd.shape[0])]
    for source, center in zip(wave_source_time, centers):
        if source is not None:
            offsets.append(abs((source - center).total_seconds()))
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")
    print(f"JADE direct MAG_VECTOR_JSSRTP records: {len(mag_time)}")
    print(f"Direct RTP records with derived |B|: "
          f"{np.count_nonzero(np.isfinite(mag_total))}")
    print(f"JADE SC_POS_R direct records: {len(position_radius)}")
    print(f"Waves direct 10-min display records: {selected_wave_bins}")
    print(f"Waves maximum source-center offset: {max(offsets):.3f} s")
    for tick, source_time, radius in zip(
            DATE_TICKS, radius_source_time, radius_tick_rj):
        if source_time is None:
            print(f"{tick:%Y-%m-%d}: no SC_POS_R within 8 h")
        else:
            offset_hours = (source_time - tick).total_seconds() / 3600.0
            print(f"{tick:%Y-%m-%d}: R={radius:.3f} R_J, source "
                  f"{source_time:%Y-%m-%d %H:%M:%S} UTC, "
                  f"offset={offset_hours:+.3f} h")


if __name__ == "__main__":
    main()
