"""Plot the visible Orbit 7 interval in the supplied Schok et al. crop.

All science values shown are direct fields from archived files:
  * MAG: official SS/JSO r60s BX, BY, BZ records
  * JADE: official Level-5 electron-moment N_CC and SC_POS_R records
  * Waves: one calibrated E-field PSD source record nearest each 5-min
    display-bin center (no average or interpolation)

This renderer mirrors plot_juno_orbit7_schok_interval.m.  It is retained as
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
TIME_START = datetime(2017, 6, 27, 8, 0, 0, tzinfo=timezone.utc)
TIME_END = datetime(2017, 7, 9, 12, 0, 0, tzinfo=timezone.utc)
DISPLAY_SECONDS = 300
FREQUENCY_RANGE_HZ = (73.244, 14929.0)
COLOR_LIMITS = (1.0e-14, 1.0e-10)

DATE_TICKS = [
    datetime(2017, 6, 29, tzinfo=timezone.utc) + timedelta(days=2 * i)
    for i in range(6)
]

# Intervals copied from the author-released Figure2.pkl object.
# Codes: 1=intermediate, 2=outer flux pileup, 3=plasmadisc.
PAPER_REGIONS = [
    (datetime(2017, 6, 19, 10, tzinfo=timezone.utc),
     datetime(2017, 6, 27, 10, tzinfo=timezone.utc), 1),
    (datetime(2017, 6, 29, 15, 10, tzinfo=timezone.utc),
     datetime(2017, 7, 3, 15, 3, tzinfo=timezone.utc), 2),
    (datetime(2017, 7, 4, 0, 1, tzinfo=timezone.utc),
     datetime(2017, 7, 5, 19, 16, tzinfo=timezone.utc), 1),
    (datetime(2017, 7, 5, 23, 46, tzinfo=timezone.utc),
     datetime(2017, 7, 9, 0, 47, tzinfo=timezone.utc), 3),
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
    mag_folder = (DATABASE_ROOT / "MAG" / "2017" / "JUPITER" / "SS" /
                  "PERI-07")
    wave_files: list[Path] = []
    mag_files: list[Path] = []
    moment_files: list[Path] = []
    for day in days:
        code = doy_code(day)
        wave_files.append(highest_version(
            wave_folder, f"WAV_{code}T*_E_V*.CSV"))
        mag_files.append(highest_version(
            mag_folder, f"fgm_jno_l3_{code}ss_r60s_v*.sts"))
        moment_folder = (DATABASE_ROOT / "JADE_data" / f"{day.year:04d}" /
                         code / "ELECTRONS")
        moment_files.append(highest_version(
            moment_folder,
            f"JAD_L50_HLS_ELC_MOM_ISO_2D_ELECTRONS_{code}_V*.CSV"))
    return wave_files, mag_files, moment_files


def decimal_doy_to_datetime(year: int, decimal_doy: float) -> datetime:
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(
        days=decimal_doy - 1.0)


def read_mag(files: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    times: list[datetime] = []
    values: list[tuple[float, float, float]] = []
    for index, path in enumerate(files, start=1):
        print(f"MAG {index:2d}/{len(files):2d}: {path}")
        with path.open("r", encoding="ascii", errors="strict") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) != 14 or len(fields[0]) != 4:
                    continue
                try:
                    year = int(fields[0])
                    decimal_doy = float(fields[6])
                    bx, by, bz = map(float, fields[7:10])
                except ValueError:
                    continue
                timestamp = decimal_doy_to_datetime(year, decimal_doy)
                if TIME_START <= timestamp < TIME_END:
                    if max(abs(bx), abs(by), abs(bz)) > 1.0e5:
                        bx = by = bz = math.nan
                    times.append(timestamp)
                    values.append((bx, by, bz))
    order = np.argsort(np.array([t.timestamp() for t in times]))
    return np.asarray(times, dtype=object)[order], np.asarray(values)[order]


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


def read_electron_moments(
        files: list[Path],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    times: list[datetime] = []
    density: list[float] = []
    sigma: list[float] = []
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
                ne = numeric[16]          # Direct N_CC, cm^-3
                ne_sigma = numeric[17]    # Direct N_SIGMA_CC, cm^-3
                if radius < 0.0 or radius > 1000.0:
                    radius = math.nan
                if ne <= 0.0 or ne > 1.0e8:
                    ne = math.nan
                if ne_sigma < 0.0 or ne_sigma > 1.0e8:
                    ne_sigma = math.nan
                times.append(timestamp)
                density.append(ne)
                sigma.append(ne_sigma)
                radius_rj.append(radius)
    order = np.argsort(np.array([t.timestamp() for t in times]))
    return (
        np.asarray(times, dtype=object)[order],
        np.asarray(density)[order],
        np.asarray(sigma)[order],
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    start_seconds = TIME_START.timestamp()
    end_seconds = TIME_END.timestamp()
    n_bins = int(round((end_seconds - start_seconds) / DISPLAY_SECONDS))
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
    return time_edges, source_time, psd, frequencies


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


def shade_regions(axis: plt.Axes) -> None:
    colors = {1: "#d8cbe2", 2: "#b7b7b7", 3: "#f5b8ae"}
    for start, end, code in PAPER_REGIONS:
        if end <= TIME_START or start >= TIME_END:
            continue
        axis.axvspan(max(start, TIME_START), min(end, TIME_END),
                     color=colors[code], alpha=0.48, linewidth=0, zorder=0)


def frequency_edges(frequencies: np.ndarray) -> np.ndarray:
    middle = np.sqrt(frequencies[:-1] * frequencies[1:])
    return np.concatenate((
        [frequencies[0] ** 2 / middle[0]], middle,
        [frequencies[-1] ** 2 / middle[-1]],
    ))


def render_figure(
        mag_time: np.ndarray, mag_b: np.ndarray,
        density_time: np.ndarray, density: np.ndarray,
        radius_tick_rj: list[float], wave_edges: np.ndarray,
        wave_psd: np.ndarray, frequencies: np.ndarray,
) -> tuple[Path, Path]:
    plt.rcParams.update({
        "font.size": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.facecolor": "white",
    })
    figure = plt.figure(figsize=(16.0, 8.8), facecolor="white")
    grid = figure.add_gridspec(
        3, 2, width_ratios=(1.0, 0.018), height_ratios=(1.05, 0.95, 1.55),
        left=0.075, right=0.945, bottom=0.11, top=0.91,
        hspace=0.13, wspace=0.045,
    )
    ax_mag = figure.add_subplot(grid[0, 0])
    ax_density = figure.add_subplot(grid[1, 0], sharex=ax_mag)
    ax_wave = figure.add_subplot(grid[2, 0], sharex=ax_mag)
    ax_cbar = figure.add_subplot(grid[2, 1])

    for axis in (ax_mag, ax_density, ax_wave):
        axis.grid(False)
        axis.set_xlim(TIME_START, TIME_END)
        axis.set_xticks(DATE_TICKS)
        axis.tick_params(which="both", direction="out")

    shade_regions(ax_mag)
    shade_regions(ax_density)
    ax_mag.plot(mag_time, mag_b[:, 0], color=(0.0, 0.35, 0.85),
                linewidth=0.65, label=r"$B_{x,JSO}$", zorder=2)
    ax_mag.plot(mag_time, mag_b[:, 1], color=(0.05, 0.58, 0.18),
                linewidth=0.65, label=r"$B_{y,JSO}$", zorder=2)
    ax_mag.plot(mag_time, mag_b[:, 2], color=(0.85, 0.12, 0.12),
                linewidth=0.65, label=r"$B_{z,JSO}$", zorder=2)
    ax_mag.set_ylim(-35.0, 35.0)
    ax_mag.set_ylabel(r"$B_{JSO}$ (nT)")
    ax_mag.legend(loc="upper left", frameon=False, ncols=3, fontsize=9,
                  handlelength=2.0, columnspacing=1.1)

    ax_density.plot(density_time, density, linestyle="none", marker=".",
                    markersize=2.6, color=(0.05, 0.35, 0.78), zorder=2)
    ax_density.set_yscale("log")
    ax_density.set_ylim(1.0e-4, 1.0)
    ax_density.set_ylabel(r"$n_e$ (cm$^{-3}$)")

    use_frequency = ((frequencies >= FREQUENCY_RANGE_HZ[0] / 1.2) &
                     (frequencies <= FREQUENCY_RANGE_HZ[1] * 1.2))
    selected_frequency = frequencies[use_frequency]
    selected_psd = wave_psd[:, use_frequency].T
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad("#e0e0e0")
    mesh = ax_wave.pcolormesh(
        mdates.date2num(wave_edges), frequency_edges(selected_frequency),
        np.ma.masked_invalid(selected_psd), shading="flat", cmap=cmap,
        norm=LogNorm(vmin=COLOR_LIMITS[0], vmax=COLOR_LIMITS[1]),
        rasterized=True,
    )
    ax_wave.set_yscale("log")
    ax_wave.set_ylim(*FREQUENCY_RANGE_HZ)
    ax_wave.set_yticks((1.0e2, 1.0e3, 1.0e4))
    ax_wave.set_yticklabels(("100 Hz", "1 kHz", "10 kHz"))
    ax_wave.set_ylabel("Frequency")
    colorbar = figure.colorbar(mesh, cax=ax_cbar)
    colorbar.set_ticks(10.0 ** np.arange(-14.0, -9.0))
    colorbar.set_label(r"E PSD  [(V m$^{-1}$)$^2$ Hz$^{-1}$]")

    ax_radius = ax_mag.twiny()
    ax_radius.set_xlim(TIME_START, TIME_END)
    ax_radius.set_xticks(DATE_TICKS)
    ax_radius.set_xticklabels([
        f"{value:.1f}" if math.isfinite(value) else "--"
        for value in radius_tick_rj
    ])
    ax_radius.set_xlabel(r"Radial distance ($R_J$)", labelpad=7)
    ax_radius.tick_params(direction="out")

    ax_mag.tick_params(labelbottom=False)
    ax_density.tick_params(labelbottom=False)
    ax_wave.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d",
                                                           tz=timezone.utc))
    for label in ax_wave.get_xticklabels():
        label.set_rotation(18)
        label.set_horizontalalignment("right")
    ax_wave.set_xlabel("")

    png_path = SCRIPT_DIR / "juno_orbit7_schok_crop_mag_ne_waves_20170627_20170709.png"
    pdf_path = SCRIPT_DIR / "juno_orbit7_schok_crop_mag_ne_waves_20170627_20170709.pdf"
    figure.savefig(png_path, dpi=260)
    figure.savefig(pdf_path, dpi=260)
    plt.close(figure)
    return png_path, pdf_path


def main() -> None:
    if not DATABASE_ROOT.is_dir():
        raise FileNotFoundError(f"Juno database is unavailable: {DATABASE_ROOT}")
    wave_files, mag_files, moment_files = discover_files()
    if not (len(wave_files) == len(mag_files) == len(moment_files) == 13):
        raise RuntimeError("Expected 13 files for each source product")
    print("Reading official direct-source records ...")
    mag_time, mag_b = read_mag(mag_files)
    density_time, density, density_sigma, density_radius = \
        read_electron_moments(moment_files)
    wave_edges, wave_source_time, wave_psd, frequencies = \
        read_waves_direct(wave_files)
    radius_source_time, radius_tick_rj = nearest_direct_radius(
        density_time, density_radius)
    png_path, pdf_path = render_figure(
        mag_time, mag_b, density_time, density, radius_tick_rj,
        wave_edges, wave_psd, frequencies)

    selected_wave_bins = int(np.count_nonzero(
        np.any(np.isfinite(wave_psd), axis=1)))
    finite_density = density[np.isfinite(density)]
    finite_sigma = density_sigma[np.isfinite(density_sigma)]
    offsets = []
    centers = [TIME_START + timedelta(seconds=DISPLAY_SECONDS * (i + 0.5))
               for i in range(wave_psd.shape[0])]
    for source, center in zip(wave_source_time, centers):
        if source is not None:
            offsets.append(abs((source - center).total_seconds()))
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")
    print(f"MAG direct records: {len(mag_time)}")
    print(f"JADE N_CC direct records: {len(finite_density)}")
    print(f"Waves direct 5-min display records: {selected_wave_bins}")
    print(f"Waves maximum source-center offset: {max(offsets):.3f} s")
    print(f"Electron density range: {finite_density.min():.6g} to "
          f"{finite_density.max():.6g} cm^-3")
    print(f"N_CC sigma range: {finite_sigma.min():.6g} to "
          f"{finite_sigma.max():.6g} cm^-3")
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
