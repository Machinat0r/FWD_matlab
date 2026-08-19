"""Gap-safe MMS four-spacecraft FOTE and FOTE-V analysis.

This is a numerical port of ``kh_fote_event.m`` for workstations where the
IRFU EpochTT/CDF interface is unavailable.  It reads the same local MMS CDF
files directly, uses no smoothing by default, and applies the same affine
four-point gradient and topology classification.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

LOCAL_PACKAGES = Path(__file__).resolve().parents[1] / "python_packages"
if LOCAL_PACKAGES.exists():
    sys.path.insert(0, str(LOCAL_PACKAGES))

import numpy as np
import pandas as pd
from cdflib import CDF, cdfepoch
from matplotlib import dates as mdates
from matplotlib import pyplot as plt
from scipy.io import savemat
from scipy.ndimage import median_filter


@dataclass
class Analysis:
    gradient: np.ndarray
    eigenvalues: np.ndarray
    divergence: np.ndarray
    curl_magnitude: np.ndarray
    xi: np.ndarray
    distance_km: np.ndarray
    distance_l: np.ndarray
    scale_l: np.ndarray
    gradient_condition: np.ndarray
    geometry_condition: np.ndarray
    topology: np.ndarray
    screw_sense: np.ndarray
    finite: np.ndarray
    metric: np.ndarray | None = None
    good: np.ndarray | None = None
    near: np.ndarray | None = None


def parse_utc(text: str) -> datetime:
    value = datetime.fromisoformat(text.replace("Z", "+00:00"))
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def to_tt2000(dt: datetime) -> np.int64:
    dt = dt.astimezone(timezone.utc)
    millis, remainder = divmod(dt.microsecond, 1000)
    return np.int64(
        cdfepoch.compute_tt2000(
            [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, millis, remainder, 0]
        )
    )


def event_dates(start: datetime, stop: datetime) -> list[datetime]:
    day = datetime(start.year, start.month, start.day, tzinfo=timezone.utc)
    last = datetime(stop.year, stop.month, stop.day, tzinfo=timezone.utc)
    result = []
    while day <= last:
        result.append(day)
        day += timedelta(days=1)
    return result


def version_key(path: Path) -> tuple[int, ...]:
    match = re.search(r"_v(\d+(?:\.\d+)*)\.cdf$", path.name)
    return tuple(int(part) for part in match.group(1).split(".")) if match else (0,)


def file_start(path: Path) -> datetime | None:
    match = re.search(r"_(\d{14})_v", path.name)
    if not match:
        return None
    return datetime.strptime(match.group(1), "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)


def discover_files(
    root: Path, spacecraft: int, product: str, start: datetime, stop: datetime
) -> list[Path]:
    candidates: list[Path] = []
    for day in event_dates(start, stop):
        if product == "fgm":
            folder = (
                root
                / f"mms{spacecraft}"
                / "fgm"
                / "brst"
                / "l2"
                / f"{day.year:04d}"
                / f"{day.month:02d}"
                / f"{day.day:02d}"
            )
        elif product == "dis-moms":
            folder = (
                root
                / f"mms{spacecraft}"
                / "fpi"
                / "brst"
                / "l2"
                / "dis-moms"
                / f"{day.year:04d}"
                / f"{day.month:02d}"
                / f"{day.day:02d}"
            )
        else:
            raise ValueError(f"Unknown product: {product}")
        if folder.exists():
            candidates.extend(folder.glob("*.cdf"))

    # Burst files are short.  A two-hour look-back safely includes the file
    # containing the first requested sample without reading a whole day.
    chosen = []
    for path in candidates:
        stamp = file_start(path)
        if stamp is not None and start - timedelta(hours=2) <= stamp <= stop:
            chosen.append(path)

    # Keep the newest semantic version for each file start/product.
    latest: dict[str, Path] = {}
    for path in chosen:
        key = re.sub(r"_v\d+(?:\.\d+)*\.cdf$", "", path.name)
        if key not in latest or version_key(path) > version_key(latest[key]):
            latest[key] = path
    result = sorted(latest.values(), key=lambda p: file_start(p) or start)
    if not result:
        raise FileNotFoundError(
            f"No {product} CDF files for MMS{spacecraft} in {start.isoformat()}--{stop.isoformat()}"
        )
    return result


def clean_values(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values[np.abs(values) > 1.0e20] = np.nan
    return values


def read_concat(
    files: list[Path], epoch_name: str, variable_name: str, start_tt: np.int64, stop_tt: np.int64
) -> tuple[np.ndarray, np.ndarray]:
    times: list[np.ndarray] = []
    values: list[np.ndarray] = []
    margin_ns = np.int64(120 * 1_000_000_000)
    for path in files:
        cdf = CDF(str(path))
        epoch = np.asarray(cdf.varget(epoch_name), dtype=np.int64).reshape(-1)
        data = clean_values(cdf.varget(variable_name))
        if data.ndim == 1:
            data = data[:, None]
        if data.shape[0] != epoch.size:
            raise ValueError(f"Epoch/data length mismatch in {path.name}: {epoch.size} vs {data.shape}")
        mask = (epoch >= start_tt - margin_ns) & (epoch <= stop_tt + margin_ns)
        if np.any(mask):
            times.append(epoch[mask])
            values.append(data[mask])
    if not times:
        raise RuntimeError(f"No samples found for {variable_name}")
    t = np.concatenate(times)
    y = np.concatenate(values, axis=0)
    order = np.argsort(t, kind="stable")
    t, y = t[order], y[order]
    unique = np.concatenate(([True], np.diff(t) != 0))
    return t[unique], y[unique]


def load_spacecraft(
    root: Path, spacecraft: int, start: datetime, stop: datetime, start_tt: np.int64, stop_tt: np.int64
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    fgm_files = discover_files(root, spacecraft, "fgm", start, stop)
    fpi_files = discover_files(root, spacecraft, "dis-moms", start, stop)
    prefix = f"mms{spacecraft}"
    b = read_concat(
        fgm_files, "Epoch", f"{prefix}_fgm_b_gse_brst_l2", start_tt, stop_tt
    )
    r = read_concat(
        fgm_files, "Epoch_state", f"{prefix}_fgm_r_gse_brst_l2", start_tt, stop_tt
    )
    v = read_concat(
        fpi_files, "Epoch", f"{prefix}_dis_bulkv_gse_brst", start_tt, stop_tt
    )
    n = read_concat(
        fpi_files, "Epoch", f"{prefix}_dis_numberdensity_brst", start_tt, stop_tt
    )
    return {
        "B": (b[0], b[1][:, :3]),
        "R": (r[0], r[1][:, :3]),
        "V": (v[0], v[1][:, :3]),
        "N": (n[0], n[1][:, :1]),
    }


def gap_aware_resample(source_t: np.ndarray, source_y: np.ndarray, target_t: np.ndarray) -> np.ndarray:
    source_t = np.asarray(source_t, dtype=np.int64)
    source_y = np.asarray(source_y, dtype=float)
    target_t = np.asarray(target_t, dtype=np.int64)
    result = np.full((target_t.size, source_y.shape[1]), np.nan)
    if source_t.size < 2:
        return result
    dt = np.diff(source_t).astype(float) / 1e9
    positive = dt[np.isfinite(dt) & (dt > 0)]
    if positive.size == 0:
        return result
    gap_limit_ns = int(max(10 * np.median(positive), 1.0) * 1e9)
    split = np.flatnonzero(np.diff(source_t) > gap_limit_ns) + 1
    starts = np.r_[0, split]
    stops = np.r_[split, source_t.size]
    for first, last in zip(starts, stops, strict=True):
        if last - first < 2:
            continue
        segment_mask = (target_t >= source_t[first]) & (target_t <= source_t[last - 1])
        if not np.any(segment_mask):
            continue
        tx = (source_t[first:last] - source_t[first]).astype(float) / 1e9
        for column in range(source_y.shape[1]):
            valid = np.isfinite(source_y[first:last, column])
            if np.count_nonzero(valid) >= 2:
                valid_rows = np.flatnonzero(valid)
                # Restrict interpolation to the finite source support.  This
                # prevents constant edge extrapolation into the half-window
                # NaNs deliberately created beside a burst gap.
                column_mask = (
                    segment_mask
                    & (target_t >= source_t[first + valid_rows[0]])
                    & (target_t <= source_t[first + valid_rows[-1]])
                )
                tq = (target_t[column_mask] - source_t[first]).astype(float) / 1e9
                result[column_mask, column] = np.interp(
                    tq, tx[valid], source_y[first:last, column][valid]
                )
    return result


def align_four_spacecraft(
    raw: list[dict[str, tuple[np.ndarray, np.ndarray]]], start_tt: np.int64, stop_tt: np.int64
) -> dict[str, np.ndarray]:
    target = raw[0]["V"][0]
    target = target[(target >= start_tt) & (target <= stop_tt)]
    target = np.unique(target)
    n_samples = target.size
    b = np.full((n_samples, 3, 4), np.nan)
    v = np.full((n_samples, 3, 4), np.nan)
    r = np.full((n_samples, 3, 4), np.nan)
    density = np.full((n_samples, 4), np.nan)
    for sc in range(4):
        b[:, :, sc] = gap_aware_resample(*raw[sc]["B"], target)
        v[:, :, sc] = gap_aware_resample(*raw[sc]["V"], target)
        r[:, :, sc] = gap_aware_resample(*raw[sc]["R"], target)
        density[:, sc] = gap_aware_resample(*raw[sc]["N"], target)[:, 0]
    valid = np.ones(n_samples, dtype=bool)
    for sc in range(4):
        valid &= np.all(np.isfinite(b[:, :, sc]), axis=1)
        valid &= np.all(np.isfinite(v[:, :, sc]), axis=1)
        valid &= np.all(np.isfinite(r[:, :, sc]), axis=1)
        valid &= np.isfinite(density[:, sc])
    if np.count_nonzero(valid) < 20:
        raise RuntimeError("Fewer than 20 common gap-safe four-spacecraft samples")
    return {"tt": target[valid], "B": b[valid], "V": v[valid], "R": r[valid], "N": density[valid]}


def one_pass_median(field: np.ndarray, samples: int, tt2000: np.ndarray) -> np.ndarray:
    if samples <= 1:
        return field
    result = field.copy()
    split = break_rows(tt2000)
    starts = np.r_[0, split]
    stops = np.r_[split, field.shape[0]]
    for first, last in zip(starts, stops, strict=True):
        if last - first < 2:
            continue
        result[first:last] = median_filter(
            field[first:last], size=(samples, 1, 1), mode="nearest"
        )
    return result


def centered_time_mean(
    tt2000: np.ndarray, values: np.ndarray, window_seconds: float
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Average on the native timestamps in an exact centered time window.

    Each output at time ``t`` uses native samples in
    ``[t-window/2, t+window/2)``.  Burst gaps are never crossed, and samples
    within half a window of a gap edge are left as NaN so every retained value
    represents a complete physical-time window.
    """
    tt2000 = np.asarray(tt2000, dtype=np.int64)
    values = np.asarray(values, dtype=float)
    result = np.full_like(values, np.nan, dtype=float)
    if window_seconds <= 0 or tt2000.size < 2:
        return result, {
            "cadence_seconds": np.nan,
            "sample_rate_hz": np.nan,
            "nominal_window_samples": 0,
        }

    cadence = np.diff(tt2000).astype(float) / 1e9
    positive = cadence[np.isfinite(cadence) & (cadence > 0)]
    if positive.size == 0:
        return result, {
            "cadence_seconds": np.nan,
            "sample_rate_hz": np.nan,
            "nominal_window_samples": 0,
        }
    median_cadence = float(np.median(positive))
    half_window_ns = np.int64(round(0.5 * window_seconds * 1e9))
    nominal_samples = max(1, int(round(window_seconds / median_cadence)))

    gap_limit_ns = np.int64(round(max(10 * median_cadence, 1.0) * 1e9))
    split = np.flatnonzero(np.diff(tt2000) > gap_limit_ns) + 1
    starts = np.r_[0, split]
    stops = np.r_[split, tt2000.size]
    flat_values = values.reshape(values.shape[0], -1)
    flat_result = result.reshape(result.shape[0], -1)

    for first, last in zip(starts, stops, strict=True):
        if last - first < 2:
            continue
        segment_t = tt2000[first:last]
        complete = (
            (segment_t - segment_t[0] >= half_window_ns)
            & (segment_t[-1] - segment_t >= half_window_ns)
        )
        if not np.any(complete):
            continue
        left = np.searchsorted(segment_t, segment_t - half_window_ns, side="left")
        right = np.searchsorted(segment_t, segment_t + half_window_ns, side="left")
        for column in range(flat_values.shape[1]):
            segment_y = flat_values[first:last, column]
            valid = np.isfinite(segment_y)
            cumulative_sum = np.r_[0.0, np.cumsum(np.where(valid, segment_y, 0.0))]
            cumulative_count = np.r_[0, np.cumsum(valid.astype(np.int64))]
            sums = cumulative_sum[right] - cumulative_sum[left]
            counts = cumulative_count[right] - cumulative_count[left]
            # Require at least 80% of the nominal native samples to prevent a
            # sparse interval from masquerading as a complete 10 s average.
            usable = complete & (counts >= max(1, int(np.floor(0.8 * nominal_samples))))
            column_result = flat_result[first:last, column]
            column_result[usable] = sums[usable] / counts[usable]

    return result, {
        "cadence_seconds": median_cadence,
        "sample_rate_hz": 1.0 / median_cadence,
        "nominal_window_samples": nominal_samples,
    }


def irfu_solid_angle(vec_a: np.ndarray, vec_b: np.ndarray, vec_c: np.ndarray) -> np.ndarray:
    """Signed solid angle with the orientation used by ``irf.solidangle``.

    IRFU assigns the sign with ``dot(cross(VecC, VecB), VecA)``.  The atan2
    form below is mathematically equivalent to its spherical-cosine formula
    and remains stable when two mapped vectors are close together.
    """
    numerator = np.einsum("ij,ij->i", np.cross(vec_c, vec_b), vec_a)
    denominator = (
        1.0
        + np.einsum("ij,ij->i", vec_a, vec_b)
        + np.einsum("ij,ij->i", vec_b, vec_c)
        + np.einsum("ij,ij->i", vec_c, vec_a)
    )
    return 2.0 * np.arctan2(numerator, denominator)


def poincare_index_irfu(field: np.ndarray) -> np.ndarray:
    """Port of IRFU ``c_4_poincare_index`` for an ``(N, 3, 4)`` field.

    The four oriented faces follow the exact IRFU order: 123, 142, 134,
    and 243.  A value near +1 or -1 means the tetrahedron encloses an odd
    number of 3-D zeros; zero means zero or an even number.
    """
    field = np.asarray(field, dtype=float)
    if field.ndim != 3 or field.shape[1:] != (3, 4):
        raise ValueError(f"Expected field shape (N, 3, 4), got {field.shape}")
    norms = np.linalg.norm(field, axis=1)
    valid = np.all(np.isfinite(field), axis=(1, 2)) & np.all(norms > 0, axis=1)
    mapped = np.full_like(field, np.nan, dtype=float)
    mapped[valid] = field[valid] / norms[valid, None, :]
    areas = (
        irfu_solid_angle(mapped[:, :, 0], mapped[:, :, 1], mapped[:, :, 2])
        + irfu_solid_angle(mapped[:, :, 0], mapped[:, :, 3], mapped[:, :, 1])
        + irfu_solid_angle(mapped[:, :, 0], mapped[:, :, 2], mapped[:, :, 3])
        + irfu_solid_angle(mapped[:, :, 1], mapped[:, :, 3], mapped[:, :, 2])
    )
    index = areas / (4.0 * np.pi)
    index[~valid] = np.nan
    index[np.isclose(index, 0.0, atol=1e-10)] = 0.0
    near_integer = np.isfinite(index) & np.isclose(index, np.rint(index), atol=1e-8)
    index[near_integer] = np.rint(index[near_integer])
    return index


def classify_gradient(eigenvalues: np.ndarray) -> tuple[str, float]:
    scale = float(np.max(np.abs(eigenvalues)))
    if not np.isfinite(scale) or scale == 0:
        return "X", np.nan
    real = eigenvalues.real
    n_positive = np.count_nonzero(real > 1e-8 * scale)
    n_negative = np.count_nonzero(real < -1e-8 * scale)
    spiral = np.max(np.abs(eigenvalues.imag)) > 1e-6 * scale
    if spiral:
        if n_positive == 1 and n_negative == 2:
            return "As", -1.0
        if n_positive == 2 and n_negative == 1:
            return "Bs", 1.0
        if n_positive == 3:
            return "S+", 2.0
        if n_negative == 3:
            return "S-", -2.0
        return "O", np.nan
    if np.min(np.abs(real)) / scale < 1e-6:
        return "X", np.nan
    if n_positive == 1 and n_negative == 2:
        return "A", np.nan
    if n_positive == 2 and n_negative == 1:
        return "B", np.nan
    if n_positive == 3:
        return "S+", np.nan
    if n_negative == 3:
        return "S-", np.nan
    return "degenerate", np.nan


def analyze_affine_field(position: np.ndarray, field: np.ndarray) -> Analysis:
    samples = field.shape[0]
    gradient = np.full((samples, 3, 3), np.nan)
    eigenvalues = np.full((samples, 3), np.nan + 0j)
    divergence = np.full(samples, np.nan)
    curl_magnitude = np.full(samples, np.nan)
    xi = np.full(samples, np.nan)
    distance_km = np.full(samples, np.nan)
    distance_l = np.full(samples, np.nan)
    scale_l = np.full(samples, np.nan)
    gradient_condition = np.full(samples, np.nan)
    geometry_condition = np.full(samples, np.nan)
    topology = np.full(samples, "degenerate", dtype="U10")
    screw_sense = np.full(samples, np.nan)
    finite = np.zeros(samples, dtype=bool)
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

    for i in range(samples):
        r = position[i].T
        f = field[i].T
        if not np.all(np.isfinite(r)) or not np.all(np.isfinite(f)):
            continue
        r0, f0 = np.mean(r, axis=0), np.mean(f, axis=0)
        x, y = r - r0, f - f0
        if np.linalg.matrix_rank(x) < 3:
            continue
        g = np.linalg.lstsq(x, y, rcond=None)[0].T
        gradient[i] = g
        geometry_condition[i] = np.linalg.cond(x)
        gradient_condition[i] = np.linalg.cond(g)
        ev = np.linalg.eigvals(g)
        eigenvalues[i] = ev
        divergence[i] = np.trace(g)
        curl = np.array([g[2, 1] - g[1, 2], g[0, 2] - g[2, 0], g[1, 0] - g[0, 1]])
        curl_magnitude[i] = np.linalg.norm(curl)
        largest = np.max(np.abs(ev))
        if largest > 0:
            xi[i] = 100 * np.abs(np.sum(ev)) / largest
        separations = np.array([np.linalg.norm(r[a] - r[b]) for a, b in pairs])
        scale_l[i] = np.nanmean(separations)
        if np.linalg.cond(g) < 1e12:
            r_null = r0 - np.linalg.solve(g, f0)
            distance_km[i] = np.min(np.linalg.norm(r - r_null, axis=1))
            distance_l[i] = distance_km[i] / scale_l[i]
        topology[i], screw_sense[i] = classify_gradient(ev)
        finite[i] = np.isfinite(distance_l[i]) and np.isfinite(gradient_condition[i])

    return Analysis(
        gradient,
        eigenvalues,
        divergence,
        curl_magnitude,
        xi,
        distance_km,
        distance_l,
        scale_l,
        gradient_condition,
        geometry_condition,
        topology,
        screw_sense,
        finite,
    )


def safe_ratio_percent(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    result = np.full(numerator.shape, np.nan)
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
    result[valid] = 100 * numerator[valid] / denominator[valid]
    return result


def summary_dict(mag: Analysis, vel: Analysis, agreement: np.ndarray) -> dict[str, int]:
    return {
        "B_quality40": int(np.count_nonzero(mag.good)),
        "V_quality40": int(np.count_nonzero(vel.good)),
        "B_As_quality40": int(np.count_nonzero(mag.good & (mag.topology == "As"))),
        "B_Bs_quality40": int(np.count_nonzero(mag.good & (mag.topology == "Bs"))),
        "V_As_quality40": int(np.count_nonzero(vel.good & (vel.topology == "As"))),
        "V_Bs_quality40": int(np.count_nonzero(vel.good & (vel.topology == "Bs"))),
        "screw_same": int(np.count_nonzero(agreement == 1)),
        "screw_opposite": int(np.count_nonzero(agreement == -1)),
    }


def break_rows(tt2000: np.ndarray) -> np.ndarray:
    """Rows after burst discontinuities; plotting values there must be NaN."""
    if tt2000.size < 3:
        return np.array([], dtype=int)
    cadence = np.diff(tt2000).astype(float) / 1e9
    positive = cadence[np.isfinite(cadence) & (cadence > 0)]
    if positive.size == 0:
        return np.array([], dtype=int)
    limit = max(10 * np.median(positive), 1.0)
    return np.flatnonzero(cadence > limit) + 1


def broken(values: np.ndarray, rows: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=float).copy()
    result[rows, ...] = np.nan
    return result


def plot_type_markers(ax, x, y, topology, good, max_y):
    specifications = [
        ("A", "^", "#d92626", "none", "A"),
        ("B", ">", "#194ccc", "none", "B"),
        ("As", "^", "#d92626", "#d92626", r"$A_s$ screw-in"),
        ("Bs", ">", "#194ccc", "#194ccc", r"$B_s$ screw-out"),
        ("X", "x", "#222222", "none", "X"),
        ("O", "o", "#222222", "none", "O"),
        ("S+", "s", "#b826a8", "#b826a8", "source"),
        ("S-", "s", "#0099a6", "#0099a6", "sink"),
    ]
    for code, marker, edge, face, label in specifications:
        mask = good & (topology == code) & np.isfinite(y) & (y <= max_y)
        if not np.any(mask):
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=9,
            marker=marker,
            edgecolors=edge,
            facecolors=face,
            linewidths=0.55,
            alpha=0.65,
            label=label,
            zorder=3,
        )
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), fontsize=8, frameon=False)


def make_plot(
    times: list[datetime], aligned: dict[str, np.ndarray], mag: Analysis, vel: Analysis,
    agreement: np.ndarray, event_id: str, start: datetime, stop: datetime,
    quality_percent: float, reliable_l: float, max_l: float, smooth_samples: int,
    summary: dict[str, int],
    raw_times: list[datetime] | None = None,
    raw_aligned: dict[str, np.ndarray] | None = None,
    smooth_seconds: float = 0.0,
    poincare_b: np.ndarray | None = None,
    poincare_v: np.ndarray | None = None,
):
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})
    show_raw = raw_aligned is not None and raw_times is not None
    show_poincare = poincare_b is not None and poincare_v is not None
    panel_count = (7 if show_raw else 5) + (1 if show_poincare else 0)
    figure_height = (18.0 if show_raw else 14.0) if show_poincare else (16.0 if show_raw else 12.2)
    fig, axes = plt.subplots(
        panel_count, 1, figsize=(15.5, figure_height), sharex=True, constrained_layout=False
    )
    fig.subplots_adjust(
        left=0.09,
        right=0.80,
        top=0.955 if show_raw else 0.94,
        bottom=0.055 if show_raw else 0.07,
        hspace=0.10,
    )
    colors = ("#1959cf", "#198f42", "#d92626")

    def plot_b_panel(ax, panel_times, panel_data, ylabel):
        mean_field = np.nanmean(panel_data["B"], axis=2)
        panel_breaks = break_rows(panel_data["tt"])
        plot_field = broken(mean_field, panel_breaks)
        for component, color, label in zip(
            range(3), colors, (r"$B_x$", r"$B_y$", r"$B_z$"), strict=True
        ):
            ax.plot(panel_times, plot_field[:, component], color=color, lw=0.65, label=label)
        ax.plot(
            panel_times,
            np.linalg.norm(plot_field, axis=1),
            color="black",
            lw=0.85,
            label=r"$|B|$",
        )
        ax.set_ylabel(ylabel)
        ax.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), ncol=1, frameon=False)

    def plot_v_panel(ax, panel_times, panel_data, ylabel):
        mean_field = np.nanmean(panel_data["V"], axis=2)
        panel_breaks = break_rows(panel_data["tt"])
        plot_field = broken(mean_field, panel_breaks)
        for component, color, label in zip(
            range(3), colors, (r"$V_{ix}$", r"$V_{iy}$", r"$V_{iz}$"), strict=True
        ):
            ax.plot(panel_times, plot_field[:, component], color=color, lw=0.65, label=label)
        ax.set_ylabel(ylabel)
        ax.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), frameon=False)

    offset = 0
    if show_raw:
        plot_b_panel(axes[0], raw_times, raw_aligned, "Raw\n" + r"$B_{GSE}$ [nT]")
        plot_v_panel(axes[1], raw_times, raw_aligned, "Raw\n" + r"$V_{i,GSE}$ [km/s]")
        offset = 2

    mean_prefix = f"{smooth_seconds:g} s mean\n" if smooth_seconds > 0 else ""
    plot_b_panel(axes[offset], times, aligned, mean_prefix + r"$B_{GSE}$ [nT]")
    plot_v_panel(axes[offset + 1], times, aligned, mean_prefix + r"$V_{i,GSE}$ [km/s]")

    plot_breaks = break_rows(aligned["tt"])

    max_error_percent = 100.0
    b_error = np.maximum(mag.metric, mag.xi)
    b_error_axis = axes[offset + 2]
    v_error_axis = axes[offset + 3]
    poincare_axis = axes[offset + 4] if show_poincare else None
    type_axis = axes[offset + 5] if show_poincare else axes[offset + 4]
    b_error_axis.plot(
        times,
        broken(np.minimum(mag.metric, max_error_percent), plot_breaks),
        color="0.35",
        lw=0.55,
        label=r"$\eta_B$",
    )
    b_error_axis.plot(
        times,
        broken(np.minimum(mag.xi, max_error_percent), plot_breaks),
        color="#7b49a5",
        lw=0.55,
        label=r"$\xi_B$",
    )
    plot_type_markers(
        b_error_axis, np.asarray(times), b_error, mag.topology, mag.good, max_error_percent
    )
    b_error_axis.axhline(quality_percent, color="black", ls="--", lw=0.9)
    b_error_axis.set_ylim(0, max_error_percent)
    b_error_axis.set_ylabel("FOTE B error [%]")

    v_error_axis.plot(
        times,
        broken(np.minimum(vel.metric, max_error_percent), plot_breaks),
        color="0.35",
        lw=0.55,
        label=r"$\alpha_V$",
    )
    plot_type_markers(
        v_error_axis, np.asarray(times), vel.metric, vel.topology, vel.good, max_error_percent
    )
    v_error_axis.axhline(quality_percent, color="black", ls="--", lw=0.9)
    v_error_axis.set_ylim(0, max_error_percent)
    v_error_axis.set_ylabel(r"FOTE-V $V_i$ error [%]")

    if show_poincare:
        pi_b_plot = broken(poincare_b, plot_breaks)
        pi_v_plot = broken(poincare_v, plot_breaks)
        poincare_axis.axhline(0.0, color="0.55", ls="--", lw=0.7)
        poincare_axis.step(
            times, pi_b_plot, where="mid", color="#194ccc", lw=1.1, label=r"$PI_B$"
        )
        poincare_axis.step(
            times, pi_v_plot, where="mid", color="#d92626", lw=0.8,
            ls=(0, (4, 3)), label=r"$PI_V$"
        )
        b_detected = np.isfinite(poincare_b) & (np.abs(poincare_b) >= 0.5)
        v_detected = np.isfinite(poincare_v) & (np.abs(poincare_v) >= 0.5)
        poincare_axis.scatter(
            np.asarray(times)[b_detected], poincare_b[b_detected],
            s=12, marker="o", facecolors="none", edgecolors="#194ccc", linewidths=0.7,
        )
        poincare_axis.scatter(
            np.asarray(times)[v_detected], poincare_v[v_detected],
            s=13, marker="x", color="#d92626", linewidths=0.7,
        )
        poincare_axis.set_yticks([-1, 0, 1])
        poincare_axis.set_ylim(-1.25, 1.25)
        poincare_axis.set_ylabel("Poincare index")
        poincare_axis.legend(
            loc="center left", bbox_to_anchor=(1.005, 0.5), frameon=False
        )

    time_array = np.asarray(times)
    b_out = mag.good & (mag.topology == "Bs")
    b_in = mag.good & (mag.topology == "As")
    v_out = vel.good & (vel.topology == "Bs")
    v_in = vel.good & (vel.topology == "As")
    type_axis.scatter(time_array[b_out], np.full(np.count_nonzero(b_out), 3), s=12, marker=">", color="#194ccc", alpha=0.65)
    type_axis.scatter(time_array[b_in], np.full(np.count_nonzero(b_in), 2), s=12, marker="^", color="#d92626", alpha=0.65)
    type_axis.scatter(time_array[v_out], np.full(np.count_nonzero(v_out), 1), s=12, marker=">", color="#194ccc", alpha=0.65)
    type_axis.scatter(time_array[v_in], np.full(np.count_nonzero(v_in), 0), s=12, marker="^", color="#d92626", alpha=0.65)
    type_axis.set_yticks([0, 1, 2, 3], labels=[r"$V_i$ screw-in", r"$V_i$ screw-out", "B screw-in", "B screw-out"])
    type_axis.set_ylim(-0.6, 3.6)
    type_axis.set_ylabel("error < 40%")

    for ax in axes:
        ax.grid(True, which="both", color="0.88", lw=0.5)
        ax.set_xlim(start, stop)
    locator = mdates.AutoDateLocator(minticks=6, maxticks=10)
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S", tz=timezone.utc))
    axes[-1].set_xlabel(f"UTC on {start:%Y-%m-%d}")

    fig.suptitle(event_id, fontsize=14, fontweight="bold")
    return fig


def run_event(
    event_id: str,
    start_text: str,
    stop_text: str,
    data_root: Path,
    output_root: Path,
    smooth_samples: int = 7,
    quality_percent: float = 40.0,
    reliable_l: float = 2.0,
    max_l: float = 4.0,
    pdf_pages=None,
    figure_target: Path | None = None,
    save_detail: bool = True,
    smooth_seconds: float = 0.0,
) -> dict[str, object]:
    start, stop = parse_utc(start_text), parse_utc(stop_text)
    start_tt, stop_tt = to_tt2000(start), to_tt2000(stop)
    print(f"Loading {event_id}: {start.isoformat()} -- {stop.isoformat()}", flush=True)
    raw = []
    smoothing_diagnostics: list[dict[str, object]] = []
    for spacecraft in range(1, 5):
        data = load_spacecraft(data_root, spacecraft, start, stop, start_tt, stop_tt)
        raw.append(data)
        print(
            f"MMS{spacecraft}: B={data['B'][0].size}, Vi={data['V'][0].size}, Ni={data['N'][0].size}, R={data['R'][0].size}",
            flush=True,
        )
    raw_aligned = align_four_spacecraft(raw, start_tt, stop_tt) if smooth_seconds > 0 else None
    if smooth_seconds > 0:
        for spacecraft, data in enumerate(raw, start=1):
            for product in ("B", "V"):
                native_t, native_values = data[product]
                averaged, diagnostics = centered_time_mean(
                    native_t, native_values, smooth_seconds
                )
                data[product] = (native_t, averaged)
                smoothing_diagnostics.append(
                    {"spacecraft": spacecraft, "product": product, **diagnostics}
                )
    aligned = align_four_spacecraft(raw, start_tt, stop_tt)
    cadence = np.median(np.diff(aligned["tt"]).astype(float) / 1e9)
    print(f"Common gap-safe samples={aligned['tt'].size}, cadence={cadence:.3f} s", flush=True)
    if smooth_seconds > 0:
        for item in smoothing_diagnostics:
            print(
                f"MMS{item['spacecraft']} {item['product']}: "
                f"native cadence={item['cadence_seconds']:.9f} s, "
                f"{smooth_seconds:g} s nominal samples={item['nominal_window_samples']}",
                flush=True,
            )
    else:
        aligned["B"] = one_pass_median(aligned["B"], smooth_samples, aligned["tt"])
        aligned["V"] = one_pass_median(aligned["V"], smooth_samples, aligned["tt"])

    mag = analyze_affine_field(aligned["R"], aligned["B"])
    vel = analyze_affine_field(aligned["R"], aligned["V"])
    poincare_b = poincare_index_irfu(aligned["B"])
    poincare_v = poincare_index_irfu(aligned["V"])
    flux_field = aligned["V"] * aligned["N"][:, None, :]
    flux = analyze_affine_field(aligned["R"], flux_field)
    mag.metric = safe_ratio_percent(np.abs(mag.divergence), mag.curl_magnitude)
    vel.metric = safe_ratio_percent(np.abs(flux.divergence), flux.curl_magnitude)
    # Quality-only selection requested for this catalog run.  Null distance,
    # gradient condition number, and tetrahedron geometry are retained in the
    # exported diagnostics, but they do not determine which types are plotted.
    mag_eigen_finite = np.all(np.isfinite(mag.eigenvalues), axis=1)
    vel_eigen_finite = np.all(np.isfinite(vel.eigenvalues), axis=1)
    mag.good = (
        mag_eigen_finite
        & np.isfinite(mag.metric)
        & np.isfinite(mag.xi)
        & (mag.metric <= quality_percent)
        & (mag.xi <= quality_percent)
    )
    vel.good = (
        vel_eigen_finite
        & np.isfinite(vel.metric)
        & (vel.metric <= quality_percent)
    )
    mag.near = mag.good & (mag.distance_l <= reliable_l)
    vel.near = vel.good & (vel.distance_l <= reliable_l)
    agreement = np.full(aligned["tt"].shape, np.nan)
    spiral = (np.abs(mag.screw_sense) == 1) & (np.abs(vel.screw_sense) == 1)
    agreement[spiral & (mag.screw_sense == vel.screw_sense)] = 1
    agreement[spiral & (mag.screw_sense != vel.screw_sense)] = -1
    agreement[~(mag.good & vel.good)] = np.nan
    summary = summary_dict(mag, vel, agreement)
    b_pi_detected = np.isfinite(poincare_b) & (np.abs(poincare_b) >= 0.5)
    v_pi_detected = np.isfinite(poincare_v) & (np.abs(poincare_v) >= 0.5)
    summary.update(
        {
            "B_PI_nonzero": int(np.count_nonzero(b_pi_detected)),
            "V_PI_nonzero": int(np.count_nonzero(v_pi_detected)),
            "B_PI_and_quality40": int(np.count_nonzero(b_pi_detected & mag.good)),
            "V_PI_and_quality40": int(np.count_nonzero(v_pi_detected & vel.good)),
            "B_and_V_PI_nonzero": int(np.count_nonzero(b_pi_detected & v_pi_detected)),
        }
    )

    pdf_dir = output_root / "pdf"
    data_dir = output_root / "data"
    timeseries_dir = data_dir / "timeseries"
    mat_dir = data_dir / "mat"
    summary_dir = data_dir / "summaries"
    for folder in (pdf_dir, timeseries_dir, mat_dir, summary_dir):
        folder.mkdir(parents=True, exist_ok=True)
    compact_start = start.strftime("%Y%m%d%H%M%S")
    compact_stop = stop.strftime("%Y%m%d%H%M%S")
    smoothing_tag = (
        f"_{smooth_seconds:g}s_native_time_mean_rawpanels_poincare_quality40"
        if smooth_seconds > 0
        else ""
    )
    base = f"{event_id}_{compact_start}_{compact_stop}_FOTE_FOTEV_VI{smoothing_tag}"
    pdf_path = pdf_dir / f"{base}.pdf" if pdf_pages is None else figure_target
    csv_path = timeseries_dir / f"{base}_timeseries.csv"
    mat_path = mat_dir / f"{base}.mat"
    json_path = summary_dir / f"{base}_summary.json"
    times = [start + timedelta(seconds=float((tt - start_tt) / 1e9)) for tt in aligned["tt"]]
    raw_times = (
        [start + timedelta(seconds=float((tt - start_tt) / 1e9)) for tt in raw_aligned["tt"]]
        if raw_aligned is not None
        else None
    )
    fig = make_plot(
        times, aligned, mag, vel, agreement, event_id, start, stop,
        quality_percent, reliable_l, max_l, smooth_samples, summary,
        raw_times=raw_times, raw_aligned=raw_aligned, smooth_seconds=smooth_seconds,
        poincare_b=poincare_b, poincare_v=poincare_v,
    )
    if pdf_pages is None:
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    else:
        pdf_pages.savefig(fig, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    if save_detail:
        table = pd.DataFrame(
            {
                "TimeUTC": [time.isoformat().replace("+00:00", "Z") for time in times],
                "B_NullDistance_km": mag.distance_km,
                "B_NullDistance_L": mag.distance_l,
                "B_Type": mag.topology,
                "B_Eta_percent": mag.metric,
                "B_Xi_percent": mag.xi,
                "B_Good": mag.good,
                "B_Near": mag.near,
                "V_NullDistance_km": vel.distance_km,
                "V_NullDistance_L": vel.distance_l,
                "V_Type": vel.topology,
                "V_Alpha_percent": vel.metric,
                "V_Good": vel.good,
                "V_Near": vel.near,
                "B_PoincareIndex": poincare_b,
                "B_PI_Nonzero": b_pi_detected,
                "V_PoincareIndex": poincare_v,
                "V_PI_Nonzero": v_pi_detected,
                "ScrewAgreement": agreement,
            }
        )
        table.to_csv(csv_path, index=False)
        savemat(
            mat_path,
            {
                "tt2000": aligned["tt"], "B": aligned["B"], "Vi": aligned["V"],
                "R": aligned["R"], "Ni": aligned["N"], "B_gradient": mag.gradient,
                "Vi_gradient": vel.gradient, "B_distance_L": mag.distance_l,
                "Vi_distance_L": vel.distance_l, "B_type": mag.topology.astype(object),
                "Vi_type": vel.topology.astype(object), "B_eta_percent": mag.metric,
                "B_xi_percent": mag.xi, "Vi_alpha_percent": vel.metric,
                "B_poincare_index": poincare_b, "Vi_poincare_index": poincare_v,
                "agreement": agreement,
            },
            do_compression=True,
        )
    report = {
        "event_id": event_id,
        "start_utc": start_text,
        "stop_utc": stop_text,
        "velocity_field": "Vi",
        "smooth_samples": smooth_samples,
        "smooth_seconds": smooth_seconds,
        "smoothing_method": (
            "centered arithmetic mean on each product's native timestamps before common-time alignment"
            if smooth_seconds > 0
            else "one-pass median on the common timeline"
        ),
        "native_smoothing_diagnostics": smoothing_diagnostics,
        "raw_common_samples": int(raw_aligned["tt"].size) if raw_aligned is not None else 0,
        "common_samples": int(aligned["tt"].size),
        "median_cadence_seconds": float(cadence),
        "summary": summary,
        "figure_pdf": str(pdf_path),
        "timeseries": str(csv_path) if save_detail else "",
        "mat_file": str(mat_path) if save_detail else "",
    }
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("event_id")
    parser.add_argument("start_utc")
    parser.add_argument("stop_utc")
    parser.add_argument("--data-root", default=r"Z:\SPART-WORK\Data\MMS")
    parser.add_argument(
        "--output-root",
        default=str(Path(__file__).resolve().parents[2] / "outputs"),
    )
    parser.add_argument("--smooth-samples", type=int, default=7)
    parser.add_argument("--smooth-seconds", type=float, default=0.0)
    parser.add_argument("--quality-percent", type=float, default=40.0)
    parser.add_argument("--reliable-l", type=float, default=2.0)
    parser.add_argument("--max-l", type=float, default=4.0)
    args = parser.parse_args()
    run_event(
        args.event_id,
        args.start_utc,
        args.stop_utc,
        Path(args.data_root),
        Path(args.output_root),
        args.smooth_samples,
        args.quality_percent,
        args.reliable_l,
        args.max_l,
        smooth_seconds=args.smooth_seconds,
    )


if __name__ == "__main__":
    main()
