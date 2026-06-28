#!/usr/bin/env python3
"""Convert CDAWeb/HAPI CSV data to a small MATLAB v4 MAT cache.

The ACE plotting script uses this helper because the local MATLAB R2026a
batch process can crash while parsing large text files line by line.  MATLAB
Level-4 MAT files are simple and load quickly without optional Python
packages such as scipy.
"""

from __future__ import annotations

import csv
import math
import os
import struct
import sys
import urllib.request
from datetime import date


def matlab_datenum(value: str) -> float:
    value = value.strip()
    if value.endswith("Z"):
        value = value[:-1]
    day, clock = value.split("T", 1)
    year = int(day[0:4])
    month = int(day[5:7])
    day_num = int(day[8:10])
    hour = int(clock[0:2])
    minute = int(clock[3:5])
    second = float(clock[6:])
    return (
        date(year, month, day_num).toordinal()
        + 366
        + (hour * 3600.0 + minute * 60.0 + second) / 86400.0
    )


def parse_csv(csv_path: str) -> tuple[list[float], list[list[float]]]:
    times: list[float] = []
    values: list[list[float]] = []
    width = 0

    with open(csv_path, "r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row or not row[0] or row[0].startswith("#"):
                continue
            try:
                t = matlab_datenum(row[0])
            except Exception:
                continue

            vals: list[float] = []
            for item in row[1:]:
                try:
                    val = float(item)
                except Exception:
                    val = math.nan
                if abs(val) > 1.0e30:
                    val = math.nan
                vals.append(val)

            width = max(width, len(vals))
            times.append(t)
            values.append(vals)

    if not times or width == 0:
        raise RuntimeError(f"No numeric HAPI rows found in {csv_path}")

    for row in values:
        if len(row) < width:
            row.extend([math.nan] * (width - len(row)))

    return times, values


def write_doubles(handle, iterable) -> None:
    chunk = []
    for value in iterable:
        chunk.append(float(value))
        if len(chunk) >= 8192:
            handle.write(struct.pack(f"<{len(chunk)}d", *chunk))
            chunk.clear()
    if chunk:
        handle.write(struct.pack(f"<{len(chunk)}d", *chunk))


def write_matrix(handle, name: str, mrows: int, ncols: int, columns) -> None:
    encoded_name = name.encode("ascii") + b"\0"
    handle.write(struct.pack("<5i", 0, mrows, ncols, 0, len(encoded_name)))
    handle.write(encoded_name)
    write_doubles(handle, columns)


def write_mat_v4(mat_path: str, times: list[float], values: list[list[float]]) -> None:
    nrows = len(times)
    ncols = len(values[0])
    tmp_path = mat_path + ".tmp"
    with open(tmp_path, "wb") as handle:
        write_matrix(handle, "t", nrows, 1, times)
        write_matrix(
            handle,
            "y",
            nrows,
            ncols,
            (values[row][col] for col in range(ncols) for row in range(nrows)),
        )
    os.replace(tmp_path, mat_path)


def main(argv: list[str]) -> int:
    if len(argv) not in (3, 4):
        print("Usage: hapi_csv_to_mat_v4.py CSV_PATH MAT_PATH [URL]", file=sys.stderr)
        return 2

    csv_path = argv[1]
    mat_path = argv[2]
    url = argv[3] if len(argv) == 4 else ""

    if (not os.path.exists(csv_path) or os.path.getsize(csv_path) < 128) and url:
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        urllib.request.urlretrieve(url, csv_path)

    times, values = parse_csv(csv_path)
    os.makedirs(os.path.dirname(mat_path), exist_ok=True)
    write_mat_v4(mat_path, times, values)
    print(f"Wrote {mat_path} ({len(times)} rows, {len(values[0])} columns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
