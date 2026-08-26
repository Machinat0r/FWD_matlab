#!/usr/bin/env python3
"""Download only the COHO one-hour months needed by Case1 PPT events.

All transfers use the official NASA/SPDF archive directly.  The selected
month set is derived from each marked event day plus seven days on each
side.  Existing validated files are retained.
"""

from __future__ import annotations

import calendar
import datetime as dt
import json
import os
from pathlib import Path
import tempfile

# The configured Windows proxy truncates NASA/SPDF TLS handshakes on this
# computer.  Direct access to the same official host is available.
os.environ["NO_PROXY"] = "*"
os.environ["no_proxy"] = "*"

import download_voyager_files as voyager_download


PROGRAM_ROOT = Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
OUTPUT_ROOT = MONTHLY_ROOT / "Case1_Selected_Event_Data"
STAGE_ROOT = Path(tempfile.gettempdir()) / "Voyager_Case1_staging"
MANIFEST = OUTPUT_ROOT / "case1_ppt_coho_download_manifest.json"


EVENT_DATES = {
    1: [
        "2013-01-09", "2013-05-30", "2013-09-11", "2013-12-29",
        "2014-04-17", "2014-05-21", "2014-08-12", "2014-08-28",
        "2014-09-21", "2014-12-13",
        "2015-03-04", "2015-10-04", "2015-10-13", "2015-12-19",
        "2015-12-28",
        "2016-09-13", "2016-09-22",
        "2017-01-15", "2017-05-24", "2017-07-18", "2017-08-15",
        "2017-09-16", "2017-10-10", "2017-12-06",
        "2018-01-08", "2018-02-02", "2018-04-24", "2018-07-07",
        "2018-08-03", "2018-08-23", "2018-10-30", "2018-11-10",
        "2018-12-12",
        "2019-05-02", "2019-05-27", "2019-05-30", "2019-06-13",
        "2019-06-24", "2019-08-26",
        "2020-05-22", "2020-09-23", "2020-12-06",
        "2021-01-06", "2021-05-06", "2021-05-17",
    ],
    2: ["2019-02-24", "2019-09-20", "2019-11-13", "2020-11-21"],
}


def parse_dates(values: list[str]) -> list[dt.date]:
    return [dt.date.fromisoformat(value) for value in values]


def required_months() -> set[tuple[int, int, int]]:
    result: set[tuple[int, int, int]] = set()
    for spacecraft, strings in EVENT_DATES.items():
        for event_day in parse_dates(strings):
            start = event_day - dt.timedelta(days=7)
            stop = event_day + dt.timedelta(days=8)
            current = dt.date(start.year, start.month, 1)
            while current <= stop:
                result.add((spacecraft, current.year, current.month))
                if current.month == 12:
                    current = dt.date(current.year + 1, 1, 1)
                else:
                    current = dt.date(current.year, current.month + 1, 1)
    return result


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    STAGE_ROOT.mkdir(parents=True, exist_ok=True)
    months = required_months()
    discovered: list[dict[str, object]] = []

    for spacecraft, year in sorted({(sc, year) for sc, year, _ in months}):
        records = voyager_download.discover_coho_year(
            spacecraft,
            year,
            dt.date(year, 1, 1),
            dt.date(year, 12, 31),
        )
        for record in records:
            file_date = record["date"]
            if isinstance(file_date, str):
                file_date = dt.date.fromisoformat(file_date)
            if (spacecraft, file_date.year, file_date.month) in months:
                discovered.append(record)

    found_months = {
        (
            int(record["spacecraft"]),
            dt.date.fromisoformat(str(record["date"])).year,
            dt.date.fromisoformat(str(record["date"])).month,
        )
        for record in discovered
    }
    missing_months = sorted(months - found_months)
    if missing_months:
        raise RuntimeError(f"NASA listing lacks required COHO months: {missing_months}")

    completed: list[dict[str, object]] = []
    for record in discovered:
        completed.append(
            voyager_download.download_one(
                record,
                OUTPUT_ROOT,
                STAGE_ROOT,
                check_size=True,
                force=False,
            )
        )

    errors = [item for item in completed if item.get("status") == "error"]
    metadata = {
        "source": "NASA/SPDF Voyager COHO one-hour archive",
        "event_context_days_each_side": 7,
        "required_months": [
            {"spacecraft": sc, "year": year, "month": month}
            for sc, year, month in sorted(months)
        ],
        "records": completed,
        "errors": errors,
    }
    MANIFEST.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({
        "required_month_count": len(months),
        "file_count": len(completed),
        "downloaded": sum(item.get("status") == "downloaded" for item in completed),
        "skipped": sum(item.get("status") == "skipped" for item in completed),
        "errors": len(errors),
        "manifest": str(MANIFEST),
    }, indent=2))
    if errors:
        raise RuntimeError(f"{len(errors)} COHO downloads failed")


if __name__ == "__main__":
    main()
