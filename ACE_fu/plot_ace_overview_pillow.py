#!/usr/bin/env python3
"""Headless ACE overview renderer for plot_ace_20260119_event.m."""

from __future__ import annotations

import math
import os
import struct
import sys
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFont


WIDTH = 980
HEIGHT = 1320
LEFT = 110
RIGHT = 895
TOP = 12
BOTTOM = 64
PANEL_COUNT = 7
PANEL_H = (HEIGHT - TOP - BOTTOM) // PANEL_COUNT

BLUE = (45, 78, 190)
GREEN = (89, 190, 0)
RED = (238, 48, 35)
BLACK = (0, 0, 0)
GRID = (218, 218, 218)
AXIS = (80, 80, 80)
TEXT = (45, 45, 45)
FLUX_COLORS = [
    (45, 78, 190),
    (89, 190, 0),
    (238, 48, 35),
    (0, 135, 160),
    (180, 60, 180),
    (165, 110, 0),
    (65, 65, 65),
    (235, 135, 25),
]


def load_font(size: int, bold: bool = False):
    candidates = []
    if os.name == "nt":
        candidates.extend([
            r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
            r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
        ])
    candidates.append("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except Exception:
            pass
    return ImageFont.load_default()


FONT_TICK = load_font(15)
FONT_LABEL = load_font(20)
FONT_PANEL = load_font(28, True)
FONT_LEGEND = load_font(18)
FONT_SMALL = load_font(13)


def read_mat_v4(path: str) -> dict[str, np.ndarray]:
    data: dict[str, np.ndarray] = {}
    with open(path, "rb") as handle:
        while True:
            header = handle.read(20)
            if not header:
                break
            if len(header) != 20:
                raise RuntimeError(f"Truncated MAT header in {path}")
            _mopt, mrows, ncols, imagf, namelen = struct.unpack("<5i", header)
            name = handle.read(namelen).rstrip(b"\0").decode("ascii", "replace")
            count = mrows * ncols
            raw = handle.read(count * 8)
            if len(raw) != count * 8:
                raise RuntimeError(f"Truncated MAT data in {path}:{name}")
            arr = np.frombuffer(raw, dtype="<f8").reshape((mrows, ncols), order="F").copy()
            data[name] = arr
            if imagf:
                handle.seek(count * 8, os.SEEK_CUR)
    return data


def as_vec(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr, dtype=float).reshape(-1)


def x_to_px(t: np.ndarray | float, plot0: float, plot1: float) -> np.ndarray | float:
    return LEFT + (np.asarray(t) - plot0) / (plot1 - plot0) * (RIGHT - LEFT)


def panel_box(index: int) -> tuple[int, int, int, int]:
    y0 = TOP + index * PANEL_H
    y1 = TOP + (index + 1) * PANEL_H
    if index == PANEL_COUNT - 1:
        y1 = HEIGHT - BOTTOM
    return LEFT, y0, RIGHT, y1


def nice_ticks(lo: float, hi: float, count: int = 5) -> list[float]:
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return [lo]
    span = abs(hi - lo)
    raw = span / max(count - 1, 1)
    base = 10 ** math.floor(math.log10(raw))
    step = base
    for mult in (1, 2, 5, 10):
        if raw <= mult * base:
            step = mult * base
            break
    start = math.ceil(lo / step) * step
    end = math.floor(hi / step) * step
    ticks = []
    value = start
    while value <= end + step * 0.1:
        ticks.append(value)
        value += step
    return ticks if len(ticks) >= 2 else [lo, hi]


def finite_limits(series: Iterable[np.ndarray], pad_frac: float = 0.06) -> tuple[float, float]:
    values = []
    for arr in series:
        a = np.asarray(arr, dtype=float)
        values.append(a[np.isfinite(a)])
    if not values or all(v.size == 0 for v in values):
        return 0.0, 1.0
    joined = np.concatenate([v for v in values if v.size])
    lo = float(np.nanmin(joined))
    hi = float(np.nanmax(joined))
    if lo == hi:
        delta = max(abs(lo) * 0.1, 1.0)
        return lo - delta, hi + delta
    pad = (hi - lo) * pad_frac
    return lo - pad, hi + pad


def format_tick(value: float) -> str:
    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.0f}"
    if abs(value) >= 1:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2g}"


def draw_rotated_label(img: Image.Image, text: str, x: int, y: int) -> None:
    bbox = ImageDraw.Draw(Image.new("RGB", (1, 1))).textbbox((0, 0), text, font=FONT_LABEL)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    label = Image.new("RGBA", (tw + 8, th + 8), (255, 255, 255, 0))
    draw = ImageDraw.Draw(label)
    draw.text((4, 4), text, fill=TEXT, font=FONT_LABEL)
    label = label.rotate(90, expand=True)
    img.paste(label, (x, y - label.height // 2), label)


def draw_base_axis(
    img: Image.Image,
    draw: ImageDraw.ImageDraw,
    index: int,
    ylabel: str,
    ylim: tuple[float, float],
    plot0: float,
    plot1: float,
    bottom_labels: bool = False,
) -> tuple[int, int, int, int]:
    left, top, right, bottom = panel_box(index)
    draw.rectangle([left, top, right, bottom], fill=(255, 255, 255), outline=AXIS)
    lo, hi = ylim
    for tick in nice_ticks(lo, hi):
        y = bottom - (tick - lo) / (hi - lo) * (bottom - top)
        if top <= y <= bottom:
            draw.line([left, y, right, y], fill=GRID, width=1)
            label = format_tick(tick)
            bbox = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((left - 8 - (bbox[2] - bbox[0]), y - 8), label, fill=TEXT, font=FONT_TICK)
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        x = float(x_to_px(tx, plot0, plot1))
        draw.line([x, top, x, bottom], fill=GRID, width=1)
        draw.line([x, bottom, x, bottom + 6], fill=AXIS, width=1)
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            bbox = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((x - (bbox[2] - bbox[0]) / 2, bottom + 10), label, fill=TEXT, font=FONT_TICK)
        else:
            draw.line([x, top, x, top + 6], fill=AXIS, width=1)
    draw_rotated_label(img, ylabel, 18, (top + bottom) // 2)
    return left, top, right, bottom


def y_to_px(values: np.ndarray, ylim: tuple[float, float], top: int, bottom: int) -> np.ndarray:
    lo, hi = ylim
    return bottom - (values - lo) / (hi - lo) * (bottom - top)


def draw_line_series(
    draw: ImageDraw.ImageDraw,
    t: np.ndarray,
    y: np.ndarray,
    color: tuple[int, int, int],
    ylim: tuple[float, float],
    top: int,
    bottom: int,
    plot0: float,
    plot1: float,
    width: int = 1,
) -> None:
    t = as_vec(t)
    y = as_vec(y)
    mask = np.isfinite(t) & np.isfinite(y) & (t >= plot0) & (t <= plot1)
    if not np.any(mask):
        return
    xpx = x_to_px(t[mask], plot0, plot1)
    ypx = y_to_px(y[mask], ylim, top, bottom)
    last = None
    for x, yy in np.column_stack([xpx, ypx]):
        xi = int(round(float(x)))
        yi = max(top, min(bottom, int(round(float(yy)))))
        if last is not None:
            draw.line([last, (xi, yi)], fill=color, width=width)
        last = (xi, yi)


def draw_inline_legend(draw: ImageDraw.ImageDraw, labels, colors, x: int, y: int) -> None:
    offset = 0
    for label, color in zip(labels, colors):
        draw.text((x + offset, y), label, fill=color, font=FONT_LEGEND)
        bbox = draw.textbbox((0, 0), label, font=FONT_LEGEND)
        offset += bbox[2] - bbox[0] + 12


def draw_panel_letter(draw: ImageDraw.ImageDraw, letter: str, top: int) -> None:
    draw.text((LEFT + 14, top + 4), letter, fill=BLACK, font=FONT_PANEL)


def draw_event_line(draw: ImageDraw.ImageDraw, tbmax: float, top: int, bottom: int, plot0: float, plot1: float) -> None:
    if not np.isfinite(tbmax) or tbmax < plot0 or tbmax > plot1:
        return
    x = int(round(float(x_to_px(tbmax, plot0, plot1))))
    y = top
    while y < bottom:
        draw.line([x, y, x, min(y + 9, bottom)], fill=(35, 35, 35), width=1)
        y += 16


def draw_zero_line(draw: ImageDraw.ImageDraw, ylim: tuple[float, float], top: int, bottom: int) -> None:
    lo, hi = ylim
    if lo <= 0 <= hi:
        yy = bottom - (0 - lo) / (hi - lo) * (bottom - top)
        x = LEFT
        while x < RIGHT:
            draw.line([x, yy, min(x + 10, RIGHT), yy], fill=(40, 40, 40), width=1)
            x += 18


def cmap(values: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    stops = np.array([
        [46, 45, 95],
        [30, 115, 190],
        [98, 190, 185],
        [240, 230, 70],
        [250, 125, 25],
        [185, 25, 25],
    ], dtype=float)
    u = np.clip((values - vmin) / max(vmax - vmin, 1e-9), 0, 1)
    p = u * (len(stops) - 1)
    i = np.floor(p).astype(int)
    i = np.clip(i, 0, len(stops) - 2)
    f = p - i
    return ((1 - f)[..., None] * stops[i] + f[..., None] * stops[i + 1]).astype(np.uint8)


def draw_colorbar(img, draw, x: int, top: int, bottom: int, vmin: float, vmax: float, label: str) -> None:
    h = bottom - top - 24
    y0 = top + 10
    grad = np.linspace(vmax, vmin, h)[:, None]
    colors = cmap(grad, vmin, vmax).reshape(h, 1, 3)
    bar = Image.fromarray(np.repeat(colors, 13, axis=1), "RGB")
    img.paste(bar, (x, y0))
    draw.rectangle([x, y0, x + 13, y0 + h], outline=(255, 255, 255), width=2)
    for val in [vmin, (vmin + vmax) / 2, vmax]:
        yy = y0 + h - (val - vmin) / max(vmax - vmin, 1e-9) * h
        draw.line([x + 13, yy, x + 18, yy], fill=TEXT, width=1)
        draw.text((x + 21, yy - 7), f"{val:.1f}", fill=(255, 255, 255), font=FONT_SMALL)
    draw.text((x - 3, y0 + h + 3), label, fill=(255, 255, 255), font=FONT_SMALL)


def draw_spectrogram(img, draw, index, t, flux, edges_ev, ylabel, plot0, plot1, tbmax, letter, cbar_label, bottom_labels):
    left, top, right, bottom = panel_box(index)
    draw.rectangle([left, top, right, bottom], fill=(38, 38, 82), outline=AXIS)
    log_edges = np.log10(np.asarray(edges_ev, dtype=float))
    ylo, yhi = float(log_edges[0]), float(log_edges[-1])
    valid = np.isfinite(flux) & (flux > 0)
    if np.any(valid):
        log_flux = np.log10(np.where(valid, flux, np.nan))
        vmin, vmax = np.nanpercentile(log_flux, [5, 95])
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = np.nanmin(log_flux), np.nanmax(log_flux)
            if vmin == vmax:
                vmax = vmin + 1
        x = np.asarray(x_to_px(t, plot0, plot1), dtype=float)
        for row in range(flux.shape[0] - 1):
            if x[row + 1] < left or x[row] > right:
                continue
            x0 = max(left, int(math.floor(x[row])))
            x1 = min(right, max(x0 + 1, int(math.ceil(x[row + 1]))))
            if x1 <= x0:
                continue
            for ch in range(flux.shape[1]):
                val = flux[row, ch]
                if not np.isfinite(val) or val <= 0:
                    continue
                yy0 = bottom - (log_edges[ch] - ylo) / (yhi - ylo) * (bottom - top)
                yy1 = bottom - (log_edges[ch + 1] - ylo) / (yhi - ylo) * (bottom - top)
                color = tuple(int(c) for c in cmap(np.array(np.log10(val)), vmin, vmax).reshape(3))
                draw.rectangle([x0, int(yy1), x1, int(yy0)], fill=color)
        draw_colorbar(img, draw, right - 76, top, bottom, float(vmin), float(vmax), cbar_label)
    else:
        draw.text((left + 20, top + 20), "No usable flux", fill=TEXT, font=FONT_LABEL)

    for p in range(math.floor(ylo), math.ceil(yhi) + 1):
        if ylo <= p <= yhi:
            yy = bottom - (p - ylo) / (yhi - ylo) * (bottom - top)
            draw.line([left, yy, right, yy], fill=(110, 110, 130), width=1)
            draw.text((left - 58, yy - 8), f"10^{p}", fill=TEXT, font=FONT_TICK)
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        xx = float(x_to_px(tx, plot0, plot1))
        draw.line([xx, top, xx, bottom], fill=(95, 95, 120), width=1)
        draw.line([xx, bottom, xx, bottom + 6], fill=AXIS, width=1)
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            bbox = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((xx - (bbox[2] - bbox[0]) / 2, bottom + 10), label, fill=TEXT, font=FONT_TICK)
    draw.rectangle([left, top, right, bottom], outline=AXIS)
    draw_rotated_label(img, ylabel, 18, (top + bottom) // 2)
    draw_panel_letter(draw, letter, top)
    draw_event_line(draw, tbmax, top, bottom, plot0, plot1)


def log_flux_limits(t: np.ndarray, flux: np.ndarray, plot0: float, plot1: float) -> tuple[float, float]:
    keep = (t >= plot0) & (t <= plot1)
    vals = flux[keep, :]
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return 0.0, 1.0
    logs = np.log10(vals)
    lo, hi = np.nanpercentile(logs, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo = float(np.nanmin(logs))
        hi = float(np.nanmax(logs))
    if lo == hi:
        return lo - 0.5, hi + 0.5
    pad = 0.08 * (hi - lo)
    return float(lo - pad), float(hi + pad)


def draw_log_flux_axis(img, draw, index, ylabel, ylim, plot0, plot1, bottom_labels=False):
    left, top, right, bottom = panel_box(index)
    draw.rectangle([left, top, right, bottom], fill=(255, 255, 255), outline=AXIS)
    lo, hi = ylim
    p0 = math.floor(lo)
    p1 = math.ceil(hi)
    for p in range(p0, p1 + 1):
        if lo <= p <= hi:
            yy = bottom - (p - lo) / (hi - lo) * (bottom - top)
            draw.line([left, yy, right, yy], fill=GRID, width=1)
            draw.text((left - 58, yy - 8), f"10^{p}", fill=TEXT, font=FONT_TICK)
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        xx = float(x_to_px(tx, plot0, plot1))
        draw.line([xx, top, xx, bottom], fill=GRID, width=1)
        draw.line([xx, bottom, xx, bottom + 6], fill=AXIS, width=1)
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            bbox = draw.textbbox((0, 0), label, font=FONT_TICK)
            draw.text((xx - (bbox[2] - bbox[0]) / 2, bottom + 10), label, fill=TEXT, font=FONT_TICK)
        else:
            draw.line([xx, top, xx, top + 6], fill=AXIS, width=1)
    draw_rotated_label(img, ylabel, 18, (top + bottom) // 2)
    return left, top, right, bottom


def draw_log_flux_lines(draw, t, flux, colors, ylim, top, bottom, plot0, plot1):
    for ch in range(flux.shape[1]):
        vals = flux[:, ch]
        yy = np.where(np.isfinite(vals) & (vals > 0), np.log10(vals), np.nan)
        draw_line_series(draw, t, yy, colors[ch % len(colors)], ylim, top, bottom, plot0, plot1, width=1)


def draw_flux_legend(draw, labels, colors, x, y, max_width=330):
    positions = []
    draw_x = x
    draw_y = y
    max_x = x
    for label in labels:
        bbox = draw.textbbox((0, 0), label, font=FONT_SMALL)
        w = bbox[2] - bbox[0]
        if draw_x + w > x + max_width:
            draw_x = x
            draw_y += 15
        positions.append((draw_x, draw_y, label))
        max_x = max(max_x, draw_x + w)
        draw_x += w + 10
    if positions:
        draw.rectangle([x - 4, y - 3, max_x + 4, positions[-1][1] + 15], fill=(255, 255, 255), outline=(210, 210, 210))
    for (draw_x, draw_y, label), color in zip(positions, colors):
        draw.text((draw_x, draw_y), label, fill=color, font=FONT_SMALL)


def draw_flux_line_panel(img, draw, index, t, flux, labels, ylabel, plot0, plot1, tbmax, letter, bottom_labels):
    ylim = log_flux_limits(t, flux, plot0, plot1)
    _, top, right, bottom = draw_log_flux_axis(img, draw, index, ylabel, ylim, plot0, plot1, bottom_labels)
    colors = FLUX_COLORS[:flux.shape[1]]
    draw_log_flux_lines(draw, t, flux, colors, ylim, top, bottom, plot0, plot1)
    draw_flux_legend(draw, labels, colors, right - 345, top + 7)
    draw_panel_letter(draw, letter, top)
    draw_event_line(draw, tbmax, top, bottom, plot0, plot1)


def render(out_png: str, plot0: float, plot1: float, tbmax: float, mfi_mat: str, swe_mat: str, epm_mat: str) -> None:
    mfi = read_mat_v4(mfi_mat)
    swe = read_mat_v4(swe_mat)
    epm = read_mat_v4(epm_mat)
    tm = as_vec(mfi["t"])
    ym = np.asarray(mfi["y"], dtype=float)
    ts = as_vec(swe["t"])
    ys = np.asarray(swe["y"], dtype=float)
    te = as_vec(epm["t"])
    ye = np.asarray(epm["y"], dtype=float)

    img = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(img)
    keep_m = (tm >= plot0) & (tm <= plot1)
    bx, by, bz, bmag = ym[:, 1], ym[:, 2], ym[:, 3], ym[:, 0]
    ylim_b = finite_limits([bx[keep_m], by[keep_m], bz[keep_m], bmag[keep_m]])
    _, top, right, bottom = draw_base_axis(img, draw, 0, "B [nT]", ylim_b, plot0, plot1)
    draw_line_series(draw, tm, bx, BLUE, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(draw, tm, by, GREEN, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(draw, tm, bz, RED, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(draw, tm, bmag, BLACK, ylim_b, top, bottom, plot0, plot1)
    draw_inline_legend(draw, ["B_x", "B_y", "B_z", "|B|"], [BLUE, GREEN, RED, BLACK], right - 195, top + 8)
    draw_panel_letter(draw, "a", top)
    draw_event_line(draw, tbmax, top, bottom, plot0, plot1)

    keep_s = (ts >= plot0) & (ts <= plot1)
    vp = ys[:, 0] if ys.shape[1] >= 1 else np.full_like(ts, np.nan)
    ylim_v = finite_limits([vp[keep_s]])
    _, top, right, bottom = draw_base_axis(img, draw, 1, "V_p [km/s]", ylim_v, plot0, plot1)
    draw_zero_line(draw, ylim_v, top, bottom)
    draw_line_series(draw, ts, vp, BLUE, ylim_v, top, bottom, plot0, plot1, width=2)
    draw_inline_legend(draw, ["V_p"], [BLUE], right - 50, top + 8)
    draw_panel_letter(draw, "b", top)
    draw_event_line(draw, tbmax, top, bottom, plot0, plot1)

    ti = ys[:, 1] / 11604.51812 if ys.shape[1] >= 2 else np.full_like(ts, np.nan)
    ylim_t = finite_limits([ti[keep_s]])
    _, top, right, bottom = draw_base_axis(img, draw, 2, "T_i [eV]", ylim_t, plot0, plot1)
    draw_line_series(draw, ts, ti, RED, ylim_t, top, bottom, plot0, plot1, width=2)
    draw_inline_legend(draw, ["T_i"], [RED], right - 45, top + 8)
    draw_panel_letter(draw, "c", top)
    draw_event_line(draw, tbmax, top, bottom, plot0, plot1)

    ions = ye[:, :8]
    electrons = ye[:, 8:12]
    ions[(ions <= 0) | (np.abs(ions) > 1e30)] = np.nan
    electrons[(electrons <= 0) | (np.abs(electrons) > 1e30)] = np.nan
    draw_spectrogram(img, draw, 3, te, electrons, [38e3, 53e3, 103e3, 175e3, 315e3], "E_e [eV]", plot0, plot1, tbmax, "d", "log10 J_e", False)
    draw_spectrogram(img, draw, 4, te, ions, [46e3, 68e3, 115e3, 195e3, 321e3, 580e3, 1060e3, 1900e3, 4800e3], "E_i [eV]", plot0, plot1, tbmax, "e", "log10 J_i", False)
    e_labels = ["e^- 38-53 keV", "e^- 53-103 keV", "e^- 103-175 keV", "e^- 175-315 keV"]
    i_labels = ["H^+ 46-68 keV", "H^+ 67-115 keV", "H^+ 115-195 keV", "H^+ 193-321 keV",
                "H^+ 315-580 keV", "H^+ 0.58-1.06 MeV", "H^+ 1.06-1.90 MeV", "H^+ 1.88-4.80 MeV"]
    draw_flux_line_panel(img, draw, 5, te, electrons, e_labels, "J_e", plot0, plot1, tbmax, "f", False)
    draw_flux_line_panel(img, draw, 6, te, ions, i_labels, "J_i", plot0, plot1, tbmax, "g", True)

    xlabel = "2026-01-19/20 UTC"
    bbox = draw.textbbox((0, 0), xlabel, font=FONT_LABEL)
    draw.text(((LEFT + RIGHT - (bbox[2] - bbox[0])) / 2, HEIGHT - 34), xlabel, fill=TEXT, font=FONT_LABEL)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    img.save(out_png)


def main(argv: list[str]) -> int:
    if len(argv) != 8:
        print("Usage: plot_ace_overview_pillow.py OUT_PNG PLOT0 PLOT1 TBMAX MFI_MAT SWE_MAT EPM_MAT", file=sys.stderr)
        return 2
    render(argv[1], float(argv[2]), float(argv[3]), float(argv[4]), argv[5], argv[6], argv[7])
    print(f"Saved PNG: {argv[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
