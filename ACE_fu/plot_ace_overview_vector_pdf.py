#!/usr/bin/env python3
"""Create an editable/vector PDF for the 2026-01-19 ACE overview."""

from __future__ import annotations

import math
import os
import struct
import sys
from typing import Iterable

import numpy as np
from reportlab.lib.colors import Color, black, white
from reportlab.pdfgen import canvas


WIDTH = 980.0
HEIGHT = 1320.0
LEFT = 110.0
RIGHT = 895.0
TOP = 12.0
BOTTOM = 64.0
PANEL_COUNT = 7
PANEL_H = (HEIGHT - TOP - BOTTOM) / PANEL_COUNT

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


def rgb(rgb255: tuple[int, int, int]) -> Color:
    return Color(rgb255[0] / 255.0, rgb255[1] / 255.0, rgb255[2] / 255.0)


def ypdf(y_top_origin: float) -> float:
    return HEIGHT - y_top_origin


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
            data[name] = np.frombuffer(raw, dtype="<f8").reshape((mrows, ncols), order="F").copy()
            if imagf:
                handle.seek(count * 8, os.SEEK_CUR)
    return data


def as_vec(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr, dtype=float).reshape(-1)


def x_to_px(t: np.ndarray | float, plot0: float, plot1: float):
    return LEFT + (np.asarray(t) - plot0) / (plot1 - plot0) * (RIGHT - LEFT)


def panel_box(index: int) -> tuple[float, float, float, float]:
    top = TOP + index * PANEL_H
    bottom = TOP + (index + 1) * PANEL_H
    if index == PANEL_COUNT - 1:
        bottom = HEIGHT - BOTTOM
    return LEFT, top, RIGHT, bottom


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


def draw_rotated_label(c: canvas.Canvas, text: str, x: float, y_top: float) -> None:
    c.saveState()
    c.translate(x, ypdf(y_top))
    c.rotate(90)
    c.setFont("Helvetica", 20)
    c.setFillColor(rgb(TEXT))
    c.drawCentredString(0, 0, text)
    c.restoreState()


def y_to_px(values: np.ndarray, ylim: tuple[float, float], top: float, bottom: float) -> np.ndarray:
    lo, hi = ylim
    return bottom - (values - lo) / (hi - lo) * (bottom - top)


def draw_base_axis(
    c: canvas.Canvas,
    index: int,
    ylabel: str,
    ylim: tuple[float, float],
    plot0: float,
    plot1: float,
    bottom_labels: bool = False,
) -> tuple[float, float, float, float]:
    left, top, right, bottom = panel_box(index)
    c.setFillColor(white)
    c.setStrokeColor(rgb(AXIS))
    c.rect(left, ypdf(bottom), right - left, bottom - top, stroke=1, fill=1)
    lo, hi = ylim
    c.setFont("Helvetica", 15)
    for tick in nice_ticks(lo, hi):
        y = bottom - (tick - lo) / (hi - lo) * (bottom - top)
        if top <= y <= bottom:
            c.setStrokeColor(rgb(GRID))
            c.line(left, ypdf(y), right, ypdf(y))
            label = format_tick(tick)
            c.setFillColor(rgb(TEXT))
            c.drawRightString(left - 8, ypdf(y) - 5, label)
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        x = float(x_to_px(tx, plot0, plot1))
        c.setStrokeColor(rgb(GRID))
        c.line(x, ypdf(top), x, ypdf(bottom))
        c.setStrokeColor(rgb(AXIS))
        c.line(x, ypdf(bottom), x, ypdf(bottom + 6))
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            c.setFillColor(rgb(TEXT))
            c.drawCentredString(x, ypdf(bottom + 24), label)
        else:
            c.line(x, ypdf(top), x, ypdf(top + 6))
    draw_rotated_label(c, ylabel, 28, (top + bottom) / 2)
    return left, top, right, bottom


def draw_line_series(
    c: canvas.Canvas,
    t: np.ndarray,
    y: np.ndarray,
    color: tuple[int, int, int],
    ylim: tuple[float, float],
    top: float,
    bottom: float,
    plot0: float,
    plot1: float,
    width: float = 0.75,
) -> None:
    t = as_vec(t)
    y = as_vec(y)
    mask = np.isfinite(t) & np.isfinite(y) & (t >= plot0) & (t <= plot1)
    if not np.any(mask):
        return
    xpx = np.asarray(x_to_px(t[mask], plot0, plot1), dtype=float)
    ypx = y_to_px(y[mask], ylim, top, bottom)
    c.setStrokeColor(rgb(color))
    c.setLineWidth(width)
    path = c.beginPath()
    started = False
    lastx = None
    for x, yy in zip(xpx, ypx):
        yy = max(top, min(bottom, float(yy)))
        if not started:
            path.moveTo(float(x), ypdf(yy))
            started = True
            lastx = x
            continue
        if lastx is not None and abs(float(x) - float(lastx)) > (RIGHT - LEFT) / 8:
            path.moveTo(float(x), ypdf(yy))
        else:
            path.lineTo(float(x), ypdf(yy))
        lastx = x
    c.drawPath(path, stroke=1, fill=0)


def draw_inline_legend(c: canvas.Canvas, labels, colors, x: float, y: float) -> None:
    c.setFont("Helvetica", 18)
    offset = 0.0
    for label, color in zip(labels, colors):
        c.setFillColor(rgb(color))
        c.drawString(x + offset, ypdf(y + 18), label)
        offset += c.stringWidth(label, "Helvetica", 18) + 12


def draw_panel_letter(c: canvas.Canvas, letter: str, top: float) -> None:
    c.setFont("Helvetica-Bold", 28)
    c.setFillColor(black)
    c.drawString(LEFT + 14, ypdf(top + 29), letter)


def draw_event_line(c: canvas.Canvas, tbmax: float, top: float, bottom: float, plot0: float, plot1: float) -> None:
    if not np.isfinite(tbmax) or tbmax < plot0 or tbmax > plot1:
        return
    x = float(x_to_px(tbmax, plot0, plot1))
    c.setStrokeColor(rgb((35, 35, 35)))
    c.setLineWidth(0.8)
    y = top
    while y < bottom:
        c.line(x, ypdf(y), x, ypdf(min(y + 9, bottom)))
        y += 16


def draw_zero_line(c: canvas.Canvas, ylim: tuple[float, float], top: float, bottom: float) -> None:
    lo, hi = ylim
    if lo <= 0 <= hi:
        yy = bottom - (0 - lo) / (hi - lo) * (bottom - top)
        c.setStrokeColor(rgb((40, 40, 40)))
        x = LEFT
        while x < RIGHT:
            c.line(x, ypdf(yy), min(x + 10, RIGHT), ypdf(yy))
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


def draw_colorbar(c: canvas.Canvas, x: float, top: float, bottom: float, vmin: float, vmax: float, label: str) -> None:
    h = bottom - top - 24
    y0 = top + 10
    steps = 80
    for i in range(steps):
        frac0 = i / steps
        frac1 = (i + 1) / steps
        val = vmax - (vmax - vmin) * (frac0 + frac1) / 2
        color = tuple(int(v) for v in cmap(np.array(val), vmin, vmax).reshape(3))
        c.setFillColor(rgb(color))
        yy0 = y0 + h * frac0
        yy1 = y0 + h * frac1
        c.rect(x, ypdf(yy1), 13, yy1 - yy0, stroke=0, fill=1)
    c.setStrokeColor(white)
    c.rect(x, ypdf(y0 + h), 13, h, stroke=1, fill=0)
    c.setFont("Helvetica", 11)
    c.setFillColor(white)
    for val in [vmin, (vmin + vmax) / 2, vmax]:
        yy = y0 + h - (val - vmin) / max(vmax - vmin, 1e-9) * h
        c.drawString(x + 20, ypdf(yy + 4), f"{val:.1f}")
    c.drawString(x - 3, ypdf(y0 + h + 15), label)


def draw_spectrogram(c, index, t, flux, edges_ev, ylabel, plot0, plot1, tbmax, letter, cbar_label, bottom_labels):
    left, top, right, bottom = panel_box(index)
    c.setFillColor(rgb((38, 38, 82)))
    c.setStrokeColor(rgb(AXIS))
    c.rect(left, ypdf(bottom), right - left, bottom - top, stroke=1, fill=1)
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
            x0 = max(left, float(x[row]))
            x1 = min(right, max(x0 + 0.2, float(x[row + 1])))
            if x1 <= x0:
                continue
            for ch in range(flux.shape[1]):
                val = flux[row, ch]
                if not np.isfinite(val) or val <= 0:
                    continue
                yy0 = bottom - (log_edges[ch] - ylo) / (yhi - ylo) * (bottom - top)
                yy1 = bottom - (log_edges[ch + 1] - ylo) / (yhi - ylo) * (bottom - top)
                color = tuple(int(v) for v in cmap(np.array(np.log10(val)), vmin, vmax).reshape(3))
                c.setFillColor(rgb(color))
                c.rect(x0, ypdf(yy0), x1 - x0, yy0 - yy1, stroke=0, fill=1)
        draw_colorbar(c, right - 76, top, bottom, float(vmin), float(vmax), cbar_label)
    c.setFont("Helvetica", 15)
    c.setFillColor(rgb(TEXT))
    c.setStrokeColor(rgb((110, 110, 130)))
    for p in range(math.floor(ylo), math.ceil(yhi) + 1):
        if ylo <= p <= yhi:
            yy = bottom - (p - ylo) / (yhi - ylo) * (bottom - top)
            c.line(left, ypdf(yy), right, ypdf(yy))
            c.drawRightString(left - 8, ypdf(yy) - 5, f"10^{p}")
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        xx = float(x_to_px(tx, plot0, plot1))
        c.setStrokeColor(rgb((95, 95, 120)))
        c.line(xx, ypdf(top), xx, ypdf(bottom))
        c.setStrokeColor(rgb(AXIS))
        c.line(xx, ypdf(bottom), xx, ypdf(bottom + 6))
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            c.setFillColor(rgb(TEXT))
            c.drawCentredString(xx, ypdf(bottom + 24), label)
    c.setStrokeColor(rgb(AXIS))
    c.rect(left, ypdf(bottom), right - left, bottom - top, stroke=1, fill=0)
    draw_rotated_label(c, ylabel, 28, (top + bottom) / 2)
    draw_panel_letter(c, letter, top)
    draw_event_line(c, tbmax, top, bottom, plot0, plot1)


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


def draw_log_flux_axis(c, index, ylabel, ylim, plot0, plot1, bottom_labels=False):
    left, top, right, bottom = panel_box(index)
    c.setFillColor(white)
    c.setStrokeColor(rgb(AXIS))
    c.rect(left, ypdf(bottom), right - left, bottom - top, stroke=1, fill=1)
    lo, hi = ylim
    c.setFont("Helvetica", 15)
    for p in range(math.floor(lo), math.ceil(hi) + 1):
        if lo <= p <= hi:
            yy = bottom - (p - lo) / (hi - lo) * (bottom - top)
            c.setStrokeColor(rgb(GRID))
            c.line(left, ypdf(yy), right, ypdf(yy))
            c.setFillColor(rgb(TEXT))
            c.drawRightString(left - 8, ypdf(yy) - 5, f"10^{p}")
    for hour in [0, 3, 6, 9, 12]:
        tx = plot0 + hour / 24.0
        x = float(x_to_px(tx, plot0, plot1))
        c.setStrokeColor(rgb(GRID))
        c.line(x, ypdf(top), x, ypdf(bottom))
        c.setStrokeColor(rgb(AXIS))
        c.line(x, ypdf(bottom), x, ypdf(bottom + 6))
        if bottom_labels:
            label = "00:00" if hour == 12 else f"{12 + hour:02d}:00"
            c.setFillColor(rgb(TEXT))
            c.drawCentredString(x, ypdf(bottom + 24), label)
        else:
            c.line(x, ypdf(top), x, ypdf(top + 6))
    draw_rotated_label(c, ylabel, 28, (top + bottom) / 2)
    return left, top, right, bottom


def draw_log_flux_lines(c, t, flux, colors, ylim, top, bottom, plot0, plot1):
    for ch in range(flux.shape[1]):
        vals = flux[:, ch]
        yy = np.where(np.isfinite(vals) & (vals > 0), np.log10(vals), np.nan)
        draw_line_series(c, t, yy, colors[ch % len(colors)], ylim, top, bottom, plot0, plot1, width=0.8)


def draw_flux_legend(c, labels, colors, x, y, max_width=330):
    c.setFont("Helvetica", 10)
    positions = []
    draw_x = x
    draw_y = y
    max_x = x
    for label in labels:
        w = c.stringWidth(label, "Helvetica", 10)
        if draw_x + w > x + max_width:
            draw_x = x
            draw_y += 14
        positions.append((draw_x, draw_y, label, w))
        max_x = max(max_x, draw_x + w)
        draw_x += w + 10
    if positions:
        c.setFillColor(white)
        c.setStrokeColor(rgb((210, 210, 210)))
        c.rect(x - 4, ypdf(positions[-1][1] + 14), max_x - x + 8, positions[-1][1] - y + 17, stroke=1, fill=1)
    for (draw_x, draw_y, label, _w), color in zip(positions, colors):
        c.setFillColor(rgb(color))
        c.drawString(draw_x, ypdf(draw_y + 10), label)


def draw_flux_line_panel(c, index, t, flux, labels, ylabel, plot0, plot1, tbmax, letter, bottom_labels):
    ylim = log_flux_limits(t, flux, plot0, plot1)
    _, top, right, bottom = draw_log_flux_axis(c, index, ylabel, ylim, plot0, plot1, bottom_labels)
    colors = FLUX_COLORS[:flux.shape[1]]
    draw_log_flux_lines(c, t, flux, colors, ylim, top, bottom, plot0, plot1)
    draw_flux_legend(c, labels, colors, right - 345, top + 7)
    draw_panel_letter(c, letter, top)
    draw_event_line(c, tbmax, top, bottom, plot0, plot1)


def render(out_pdf: str, plot0: float, plot1: float, tbmax: float, mfi_mat: str, swe_mat: str, epm_mat: str) -> None:
    mfi = read_mat_v4(mfi_mat)
    swe = read_mat_v4(swe_mat)
    epm = read_mat_v4(epm_mat)
    tm = as_vec(mfi["t"])
    ym = np.asarray(mfi["y"], dtype=float)
    ts = as_vec(swe["t"])
    ys = np.asarray(swe["y"], dtype=float)
    te = as_vec(epm["t"])
    ye = np.asarray(epm["y"], dtype=float)

    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    c = canvas.Canvas(out_pdf, pagesize=(WIDTH, HEIGHT), pageCompression=1)

    keep_m = (tm >= plot0) & (tm <= plot1)
    bx, by, bz, bmag = ym[:, 1], ym[:, 2], ym[:, 3], ym[:, 0]
    ylim_b = finite_limits([bx[keep_m], by[keep_m], bz[keep_m], bmag[keep_m]])
    _, top, right, bottom = draw_base_axis(c, 0, "B [nT]", ylim_b, plot0, plot1)
    draw_line_series(c, tm, bx, BLUE, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(c, tm, by, GREEN, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(c, tm, bz, RED, ylim_b, top, bottom, plot0, plot1)
    draw_line_series(c, tm, bmag, BLACK, ylim_b, top, bottom, plot0, plot1)
    draw_inline_legend(c, ["B_x", "B_y", "B_z", "|B|"], [BLUE, GREEN, RED, BLACK], right - 195, top + 8)
    draw_panel_letter(c, "a", top)
    draw_event_line(c, tbmax, top, bottom, plot0, plot1)

    keep_s = (ts >= plot0) & (ts <= plot1)
    vp = ys[:, 0] if ys.shape[1] >= 1 else np.full_like(ts, np.nan)
    ylim_v = finite_limits([vp[keep_s]])
    _, top, right, bottom = draw_base_axis(c, 1, "V_p [km/s]", ylim_v, plot0, plot1)
    draw_zero_line(c, ylim_v, top, bottom)
    draw_line_series(c, ts, vp, BLUE, ylim_v, top, bottom, plot0, plot1, width=1.2)
    draw_inline_legend(c, ["V_p"], [BLUE], right - 50, top + 8)
    draw_panel_letter(c, "b", top)
    draw_event_line(c, tbmax, top, bottom, plot0, plot1)

    ti = ys[:, 1] / 11604.51812 if ys.shape[1] >= 2 else np.full_like(ts, np.nan)
    ylim_t = finite_limits([ti[keep_s]])
    _, top, right, bottom = draw_base_axis(c, 2, "T_i [eV]", ylim_t, plot0, plot1)
    draw_line_series(c, ts, ti, RED, ylim_t, top, bottom, plot0, plot1, width=1.2)
    draw_inline_legend(c, ["T_i"], [RED], right - 45, top + 8)
    draw_panel_letter(c, "c", top)
    draw_event_line(c, tbmax, top, bottom, plot0, plot1)

    ions = ye[:, :8]
    electrons = ye[:, 8:12]
    ions[(ions <= 0) | (np.abs(ions) > 1e30)] = np.nan
    electrons[(electrons <= 0) | (np.abs(electrons) > 1e30)] = np.nan
    draw_spectrogram(c, 3, te, electrons, [38e3, 53e3, 103e3, 175e3, 315e3], "E_e [eV]", plot0, plot1, tbmax, "d", "log10 J_e", False)
    draw_spectrogram(c, 4, te, ions, [46e3, 68e3, 115e3, 195e3, 321e3, 580e3, 1060e3, 1900e3, 4800e3], "E_i [eV]", plot0, plot1, tbmax, "e", "log10 J_i", False)
    e_labels = ["e^- 38-53 keV", "e^- 53-103 keV", "e^- 103-175 keV", "e^- 175-315 keV"]
    i_labels = ["H^+ 46-68 keV", "H^+ 67-115 keV", "H^+ 115-195 keV", "H^+ 193-321 keV",
                "H^+ 315-580 keV", "H^+ 0.58-1.06 MeV", "H^+ 1.06-1.90 MeV", "H^+ 1.88-4.80 MeV"]
    draw_flux_line_panel(c, 5, te, electrons, e_labels, "J_e", plot0, plot1, tbmax, "f", False)
    draw_flux_line_panel(c, 6, te, ions, i_labels, "J_i", plot0, plot1, tbmax, "g", True)

    c.setFont("Helvetica", 20)
    c.setFillColor(rgb(TEXT))
    c.drawCentredString((LEFT + RIGHT) / 2, ypdf(HEIGHT - 34), "2026-01-19/20 UTC")
    c.showPage()
    c.save()


def main(argv: list[str]) -> int:
    if len(argv) != 8:
        print("Usage: plot_ace_overview_vector_pdf.py OUT_PDF PLOT0 PLOT1 TBMAX MFI_MAT SWE_MAT EPM_MAT", file=sys.stderr)
        return 2
    render(argv[1], float(argv[2]), float(argv[3]), float(argv[4]), argv[5], argv[6], argv[7])
    print(f"Saved PDF: {argv[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
