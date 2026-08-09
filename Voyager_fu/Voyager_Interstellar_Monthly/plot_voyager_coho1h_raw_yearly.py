#!/usr/bin/env python3
"""Plot annual Voyager COHO one-hour records without interpolation."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parent
LOCAL_PACKAGES = PROGRAM_ROOT / "python_packages"
MATLAB_PACKAGES = Path(
    r"C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu"
    r"\Voyager_Interstellar_Monthly\python_packages"
)
for package_dir in (LOCAL_PACKAGES, MATLAB_PACKAGES):
    if package_dir.is_dir():
        sys.path.insert(0, str(package_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from cdflib import CDF, cdfepoch


STANDARD_VARIABLES = (
    "heliocentricDistance",
    "ABS_B",
    "F",
    "BR",
    "BT",
    "BN",
    "V",
    "protonDensity",
    "protonTemp",
)
PARTICLE_PATTERN = re.compile(r"^protonFlux(\d+)_(LECP|CRS)$")
FILE_PATTERN = re.compile(r"_(\d{8})_v(\d+)\.cdf$", re.IGNORECASE)
HELIOPAUSE = {
    1: np.datetime64("2012-08-25T00:00:00", "ns"),
    2: np.datetime64("2018-11-05T00:00:00", "ns"),
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", default=r"Z:\SPART-WORK\Data\Voyager")
    parser.add_argument(
        "--output",
        default=str(PROGRAM_ROOT / "UniformLowPrecision_Yearly_Raw"),
    )
    parser.add_argument("--spacecraft", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--gap-hours", type=float, default=1.5)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument(
        "--nasa-v2-mag",
        default=str(PROGRAM_ROOT / "NASA_COHO_Voyager2_MAG_1h_2021_2024.html"),
        help="NASA COHOWeb Voyager 2 one-hour MAG listing for 2021-2024.",
    )
    return parser.parse_args()


def monthly_files(data_root: Path, spacecraft: int) -> dict[tuple[int, int], Path]:
    root = (
        data_root
        / f"voyager{spacecraft}"
        / "coho"
        / "1hr"
        / "l2"
        / "merged_mag_plasma"
    )
    chosen: dict[tuple[int, int], tuple[int, Path]] = {}
    for path in root.rglob("*.cdf"):
        match = FILE_PATTERN.search(path.name)
        if not match:
            continue
        stamp, version_text = match.groups()
        key = (int(stamp[:4]), int(stamp[4:6]))
        version = int(version_text)
        if key not in chosen or version >= chosen[key][0]:
            chosen[key] = (version, path)
    return {key: value[1] for key, value in sorted(chosen.items())}


def number_attribute(attributes: dict, name: str) -> float | None:
    value = attributes.get(name)
    if value is None:
        return None
    try:
        return float(np.asarray(value).reshape(-1)[0])
    except (TypeError, ValueError, OverflowError):
        return None


def clean_numeric(raw: np.ndarray, attributes: dict) -> np.ndarray:
    values = np.asarray(raw, dtype=np.float64).reshape(-1)
    fill = number_attribute(attributes, "FILLVAL")
    minimum = number_attribute(attributes, "VALIDMIN")
    maximum = number_attribute(attributes, "VALIDMAX")
    if fill is not None:
        values[values == fill] = np.nan
    if minimum is not None:
        values[values < minimum] = np.nan
    if maximum is not None:
        values[values > maximum] = np.nan
    values[~np.isfinite(values) | (np.abs(values) >= 1.0e30)] = np.nan
    return values


def attribute_label(attributes: dict, fallback: str) -> str:
    for key in ("FIELDNAM", "LABLAXIS", "CATDESC"):
        value = attributes.get(key)
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        if isinstance(value, str) and value.strip():
            return value.strip()
    return fallback


def read_month(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, str]]:
    cdf = CDF(str(path))
    available = list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)
    raw_epoch = np.asarray(cdf.varget("Epoch"))
    epoch = np.asarray(cdfepoch.to_datetime(raw_epoch)).astype("datetime64[ns]")
    variables: dict[str, np.ndarray] = {}
    labels: dict[str, str] = {}
    requested = [name for name in STANDARD_VARIABLES if name in available]
    requested.extend(name for name in available if PARTICLE_PATTERN.match(name))
    for name in requested:
        attributes = cdf.varattsget(name)
        variables[name] = clean_numeric(cdf.varget(name), attributes)
        labels[name] = attribute_label(attributes, name)
    return epoch.reshape(-1), variables, labels


def load_year(paths: list[Path], spacecraft: int, year: int):
    chunks: list[tuple[np.ndarray, dict[str, np.ndarray]]] = []
    labels: dict[str, str] = {}
    all_names: set[str] = set()
    for path in paths:
        epoch, variables, current_labels = read_month(path)
        chunks.append((epoch, variables))
        labels.update(current_labels)
        all_names.update(variables)

    epoch_all: list[np.ndarray] = []
    value_all: dict[str, list[np.ndarray]] = {name: [] for name in all_names}
    for epoch, variables in chunks:
        epoch_all.append(epoch)
        for name in all_names:
            value_all[name].append(
                variables.get(name, np.full(epoch.size, np.nan, dtype=np.float64))
            )
    epoch = np.concatenate(epoch_all)
    values = {name: np.concatenate(parts) for name, parts in value_all.items()}

    start = np.datetime64(f"{year:04d}-01-01T00:00:00", "ns")
    end = np.datetime64(f"{year + 1:04d}-01-01T00:00:00", "ns")
    start = max(start, HELIOPAUSE[spacecraft])
    keep = (epoch >= start) & (epoch < end)
    epoch = epoch[keep]
    values = {name: data[keep] for name, data in values.items()}
    order = np.argsort(epoch, kind="stable")
    epoch = epoch[order]
    values = {name: data[order] for name, data in values.items()}
    if epoch.size:
        _, unique_index = np.unique(epoch, return_index=True)
        unique_index.sort()
        epoch = epoch[unique_index]
        values = {name: data[unique_index] for name, data in values.items()}
    return epoch, values, labels, start, end


def read_nasa_v2_mag(path: Path):
    """Read the numeric listing returned by NASA COHOWeb's nx1.cgi."""
    epoch: list[np.datetime64] = []
    distance: list[float] = []
    magnitude: list[float] = []
    br: list[float] = []
    bt: list[float] = []
    bn: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 8 or not re.fullmatch(r"20\d{2}", parts[0]):
            continue
        try:
            year, day_of_year, hour = map(int, parts[:3])
            numeric = np.asarray([float(item) for item in parts[3:8]], dtype=float)
        except ValueError:
            continue
        timestamp = (
            np.datetime64(f"{year:04d}-01-01T00:00:00", "ns")
            + np.timedelta64(day_of_year - 1, "D")
            + np.timedelta64(hour, "h")
        )
        numeric[np.abs(numeric) >= 900.0] = np.nan
        epoch.append(timestamp)
        distance.append(numeric[0])
        magnitude.append(numeric[1])
        br.append(numeric[2])
        bt.append(numeric[3])
        bn.append(numeric[4])
    values = {
        "heliocentricDistance": np.asarray(distance, dtype=float),
        "ABS_B": np.asarray(magnitude, dtype=float),
        "BR": np.asarray(br, dtype=float),
        "BT": np.asarray(bt, dtype=float),
        "BN": np.asarray(bn, dtype=float),
    }
    return np.asarray(epoch, dtype="datetime64[ns]"), values


def load_nasa_v2_year(
    all_epoch: np.ndarray, all_values: dict[str, np.ndarray], year: int
):
    start = np.datetime64(f"{year:04d}-01-01T00:00:00", "ns")
    end = np.datetime64(f"{year + 1:04d}-01-01T00:00:00", "ns")
    keep = (all_epoch >= start) & (all_epoch < end)
    return (
        all_epoch[keep],
        {name: data[keep] for name, data in all_values.items()},
        {},
        start,
        end,
    )


def particle_names(values: dict[str, np.ndarray], instrument: str) -> list[str]:
    selected: list[tuple[int, str]] = []
    for name in values:
        match = PARTICLE_PATTERN.match(name)
        if match and match.group(2) == instrument:
            selected.append((int(match.group(1)), name))
    return [name for _, name in sorted(selected)]


def plot_recorded(
    axis,
    epoch: np.ndarray,
    values: np.ndarray,
    color,
    gap_hours: float,
    *,
    positive: bool = False,
    label: str | None = None,
    linewidth: float = 0.55,
) -> bool:
    mask = ~np.isnat(epoch) & np.isfinite(values)
    if positive:
        mask &= values > 0
    time = epoch[mask]
    data = values[mask]
    if not time.size:
        return False
    breaks = np.flatnonzero(np.diff(time).astype("timedelta64[s]").astype(np.int64) > gap_hours * 3600) + 1
    time = np.insert(time, breaks, np.datetime64("NaT", "ns"))
    data = np.insert(data, breaks, np.nan)
    axis.plot(
        time,
        data,
        color=color,
        linewidth=linewidth,
        marker=".",
        markersize=1.4,
        label=label,
    )
    return True


def empty_panel(axis, message: str) -> None:
    axis.text(
        0.5,
        0.5,
        message,
        transform=axis.transAxes,
        ha="center",
        va="center",
        color="0.45",
        style="italic",
        fontsize=8,
    )


def first_valid(values: dict[str, np.ndarray], candidates: tuple[str, ...]):
    for name in candidates:
        if name in values and np.isfinite(values[name]).any():
            return name
    return None


def annual_figure(
    spacecraft: int,
    year: int,
    epoch: np.ndarray,
    values: dict[str, np.ndarray],
    labels: dict[str, str],
    output: Path,
    gap_hours: float,
    dpi: int,
    source_name: str,
) -> None:
    figure, axes = plt.subplots(5, 1, figsize=(12, 9), sharex=True)
    figure.patch.set_facecolor("white")
    for axis in axes:
        axis.grid(True, color="0.88", linewidth=0.45)
        axis.tick_params(labelsize=8)

    name = first_valid(values, ("ABS_B", "F"))
    if name:
        plot_recorded(axes[0], epoch, values[name], "0.08", gap_hours, linewidth=0.5)
    else:
        empty_panel(axes[0], "No recorded |B| value")
    axes[0].set_ylabel("Magnetic field |B| (nT)", fontsize=9)

    component_colors = ("#0072BD", "#D95319", "#EDB120")
    has_components = False
    for name, color, label in zip(("BR", "BT", "BN"), component_colors, ("B_R", "B_T", "B_N")):
        if name in values:
            has_components |= plot_recorded(
                axes[1], epoch, values[name], color, gap_hours, label=label
            )
    if has_components:
        axes[1].legend(loc="center left", bbox_to_anchor=(1.005, 0.5), fontsize=8)
    else:
        empty_panel(axes[1], "No recorded RTN magnetic components")
    axes[1].set_ylabel("RTN field $B_R$, $B_T$, $B_N$ (nT)", fontsize=9)

    for axis, instrument in zip(axes[2:4], ("LECP", "CRS")):
        names = particle_names(values, instrument)
        colors = plt.cm.turbo(np.linspace(0.05, 0.95, max(len(names), 2)))
        present = False
        for index, name in enumerate(names):
            present |= plot_recorded(
                axis,
                epoch,
                values[name],
                colors[index],
                gap_hours,
                positive=True,
                label=labels.get(name, name),
                linewidth=0.45,
            )
        if present:
            axis.set_yscale("log")
            axis.legend(
                loc="center left",
                bbox_to_anchor=(1.005, 0.5),
                fontsize=7.0,
                frameon=True,
            )
        else:
            empty_panel(axis, f"No recorded {instrument} proton flux")
        axis.set_ylabel(
            f"{instrument} proton flux\n(cm$^{{-2}}$ s$^{{-1}}$ sr$^{{-1}}$ MeV$^{{-1}}$)",
            fontsize=9,
        )

    name = "heliocentricDistance"
    if name in values and plot_recorded(
        axes[4], epoch, values[name], "#6236a8", gap_hours, linewidth=0.6
    ):
        pass
    else:
        empty_panel(axes[4], "No recorded heliocentric distance")
    axes[4].set_ylabel("Heliocentric distance (AU)", fontsize=9)
    axes[4].set_xlabel("Month (UTC)", fontsize=9)

    lower = datetime(year, 1, 1, tzinfo=timezone.utc)
    upper = datetime(year + 1, 1, 1, tzinfo=timezone.utc)
    axes[4].set_xlim(lower, upper)
    axes[4].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    figure.suptitle(
        f"Voyager {spacecraft} annual one-hour records — {year:04d} UTC — {source_name}",
        fontsize=12,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.008,
        "Source records are used as listed; this program performs no averaging, gap filling, interpolation, "
        f"resampling, or generated values. Dots mark records; lines break above {gap_hours:.1f} h.",
        ha="center",
        fontsize=7.5,
    )
    figure.subplots_adjust(left=0.105, right=0.78, top=0.945, bottom=0.055, hspace=0.28)
    figure.savefig(output, dpi=dpi, facecolor="white")
    plt.close(figure)


def finite_record_count(values: dict[str, np.ndarray], names: tuple[str, ...]) -> int:
    candidates = [values[name] for name in names if name in values]
    if not candidates:
        return 0
    return int(np.any(np.column_stack([np.isfinite(item) for item in candidates]), axis=1).sum())


def main() -> int:
    options = arguments()
    output = Path(options.output)
    output.mkdir(parents=True, exist_ok=True)
    nasa_path = Path(options.nasa_v2_mag)
    nasa_epoch = np.array([], dtype="datetime64[ns]")
    nasa_values: dict[str, np.ndarray] = {}
    if nasa_path.is_file():
        nasa_epoch, nasa_values = read_nasa_v2_mag(nasa_path)
    rows: list[dict[str, object]] = []
    for spacecraft in options.spacecraft:
        files = monthly_files(Path(options.data_root), spacecraft)
        year_set = {
            year for year, _ in files if year >= int(str(HELIOPAUSE[spacecraft])[:4])
        }
        if spacecraft == 2 and nasa_epoch.size:
            year_set.update(int(str(item)[:4]) for item in np.unique(nasa_epoch.astype("datetime64[Y]")))
        years = sorted(year_set)
        for year in years:
            paths = [path for (file_year, _), path in files.items() if file_year == year]
            use_nasa_mag = spacecraft == 2 and year > 2020 and nasa_epoch.size > 0
            source_name = (
                "NASA COHOWeb one-hour MAG"
                if use_nasa_mag
                else "COHO one-hour merged instrument product"
            )
            source_files = [nasa_path] if use_nasa_mag else paths
            print(f"[Annual one-hour raw] Voyager {spacecraft} {year}", flush=True)
            status = "ok"
            notes = ""
            epoch = np.array([], dtype="datetime64[ns]")
            values: dict[str, np.ndarray] = {}
            labels: dict[str, str] = {}
            data_start = max(
                np.datetime64(f"{year:04d}-01-01T00:00:00", "ns"),
                HELIOPAUSE[spacecraft],
            )
            data_end = np.datetime64(f"{year + 1:04d}-01-01T00:00:00", "ns")
            if use_nasa_mag:
                figure_file = output / f"V{spacecraft}_COHOweb1h_MAG_raw_{year:04d}.png"
            else:
                figure_file = output / f"V{spacecraft}_COHO1h_raw_{year:04d}.png"
            try:
                if use_nasa_mag:
                    epoch, values, labels, data_start, data_end = load_nasa_v2_year(
                        nasa_epoch, nasa_values, year
                    )
                else:
                    epoch, values, labels, data_start, data_end = load_year(
                        paths, spacecraft, year
                    )
                annual_figure(
                    spacecraft,
                    year,
                    epoch,
                    values,
                    labels,
                    figure_file,
                    options.gap_hours,
                    options.dpi,
                    source_name,
                )
            except Exception as error:
                status = "error"
                notes = f"{type(error).__name__}: {error}"
                print(f"ERROR V{spacecraft} {year}: {notes}", file=sys.stderr, flush=True)

            lecp = particle_names(values, "LECP")
            crs = particle_names(values, "CRS")
            rows.append(
                {
                    "Spacecraft": spacecraft,
                    "YearUTC": f"{year:04d}-01-01T00:00:00Z",
                    "DataStartUTC": f"{str(data_start)}Z",
                    "DataEndUTC": f"{str(data_end)}Z",
                    "SourceCDFCount": len(paths),
                    "SourceCDFs": ";".join(str(path) for path in paths),
                    "SourceProduct": source_name,
                    "SourceFiles": ";".join(str(path) for path in source_files),
                    "EpochRecords": int(epoch.size),
                    "MagneticRecords": finite_record_count(values, ("ABS_B", "F", "BR", "BT", "BN")),
                    "PlasmaRecords": finite_record_count(values, ("V", "protonDensity", "protonTemp")),
                    "LECPChannels": len(lecp),
                    "LECPValidValues": int(sum(np.isfinite(values[name]).sum() for name in lecp)),
                    "CRSChannels": len(crs),
                    "CRSValidValues": int(sum(np.isfinite(values[name]).sum() for name in crs)),
                    "OriginalInstrumentValuesOnly": True,
                    "GeneratedOrInterpolatedValues": 0,
                    "FigureFile": str(figure_file),
                    "Status": status,
                    "Notes": notes,
                }
            )

    fieldnames = list(rows[0]) if rows else []
    manifest = output / "UniformLowPrecision_Yearly_Raw_manifest.csv"
    with manifest.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Manifest: {manifest}", flush=True)
    print(f"Figures: {sum(row['Status'] == 'ok' for row in rows)}", flush=True)
    print(f"Errors: {sum(row['Status'] == 'error' for row in rows)}", flush=True)
    return 1 if any(row["Status"] == "error" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
