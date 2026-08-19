"""Run the selected MMS KH windows into one multi-page FOTE/FOTE-V PDF."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
LOCAL_PACKAGES = SCRIPT_DIR.parent / "python_packages"
if LOCAL_PACKAGES.exists():
    sys.path.insert(0, str(LOCAL_PACKAGES))

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from kh_fote_event import run_event


def run_catalog(catalog: Path, data_root: Path, output_root: Path, smooth_samples: int = 7):
    rows = pd.read_csv(catalog, dtype={"EventID": str})
    pdf_dir = output_root / "pdf"
    data_dir = output_root / "data"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    tag = f"quality40_no_distance_{smooth_samples}pt"
    combined_pdf = pdf_dir / f"KH_all_FOTE_FOTEV_Vi_{tag}.pdf"
    status_csv = data_dir / f"KH_all_FOTE_FOTEV_Vi_{tag}_status.csv"
    status_json = data_dir / f"KH_all_FOTE_FOTEV_Vi_{tag}_status.json"
    statuses: list[dict[str, object]] = []

    with PdfPages(combined_pdf) as pdf:
        total = len(rows)
        for index, row in rows.iterrows():
            event_id = str(row["EventID"])
            mode = str(row.get("Mode", ""))
            print(f"BATCH {index + 1}/{total} {event_id} mode={mode}", flush=True)
            if mode == "no_local_file_index":
                statuses.append(
                    {
                        "EventID": event_id,
                        "PlotStartUTC": row.get("PlotStartUTC", ""),
                        "PlotEndUTC": row.get("PlotEndUTC", ""),
                        "Mode": mode,
                        "Status": "skipped_no_local_data",
                        "Error": "No local MMS file index for this event window",
                    }
                )
                pd.DataFrame(statuses).to_csv(status_csv, index=False)
                continue
            start = str(row["PlotStartUTC"]).replace(" ", "T") + "Z"
            stop = str(row["PlotEndUTC"]).replace(" ", "T") + "Z"
            try:
                report = run_event(
                    event_id,
                    start,
                    stop,
                    data_root,
                    output_root,
                    smooth_samples=smooth_samples,
                    pdf_pages=pdf,
                    figure_target=combined_pdf,
                    save_detail=True,
                )
                item = {
                    "EventID": event_id,
                    "PlotStartUTC": row["PlotStartUTC"],
                    "PlotEndUTC": row["PlotEndUTC"],
                    "Mode": mode,
                    "Status": "ok",
                    "Error": "",
                    "CommonSamples": report["common_samples"],
                    "CadenceSeconds": report["median_cadence_seconds"],
                }
                item.update(report["summary"])
                statuses.append(item)
            except Exception as exc:
                print(f"FAILED {event_id}: {exc}", flush=True)
                traceback.print_exc()
                statuses.append(
                    {
                        "EventID": event_id,
                        "PlotStartUTC": row["PlotStartUTC"],
                        "PlotEndUTC": row["PlotEndUTC"],
                        "Mode": mode,
                        "Status": "failed",
                        "Error": str(exc),
                    }
                )
            pd.DataFrame(statuses).to_csv(status_csv, index=False)

    report = {
        "catalog": str(catalog),
        "combined_pdf": str(combined_pdf),
        "smooth_samples": smooth_samples,
        "selection": "B: eta<=40% and xi<=40%; Vi: alpha<=40%; no distance filter",
        "total_events": int(len(rows)),
        "ok": sum(item["Status"] == "ok" for item in statuses),
        "skipped_no_local_data": sum(item["Status"] == "skipped_no_local_data" for item in statuses),
        "failed": sum(item["Status"] == "failed" for item in statuses),
        "status_csv": str(status_csv),
    }
    status_json.write_text(json.dumps({"run": report, "events": statuses}, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        default=r"C:\Users\Administrator\Documents\KH\selected_burst_windows_all.csv",
    )
    parser.add_argument("--data-root", default=r"Z:\SPART-WORK\Data\MMS")
    parser.add_argument(
        "--output-root",
        default=str(Path(__file__).resolve().parents[2] / "output"),
    )
    parser.add_argument("--smooth-samples", type=int, default=7)
    args = parser.parse_args()
    run_catalog(Path(args.catalog), Path(args.data_root), Path(args.output_root), args.smooth_samples)


if __name__ == "__main__":
    main()
