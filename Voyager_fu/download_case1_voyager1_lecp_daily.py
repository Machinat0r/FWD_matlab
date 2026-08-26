#!/usr/bin/env python3
"""Download the official V1 LECP daily sector data needed by Case1."""

from __future__ import annotations

import json
import os
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET


os.environ["NO_PROXY"] = "*"
os.environ["no_proxy"] = "*"

PROGRAM_ROOT = Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DATABASE_ROOT = MONTHLY_ROOT / "Case1_Selected_Event_Data" / "voyager1"
OUTPUT_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1d"
    / "l2"
    / "sectored_flux"
    / "2013-2021"
    / "voyager1_lecp_p1_sectored_daily_case1_20130102_20210525.cdf"
)
METADATA_FILE = OUTPUT_FILE.with_suffix(".json")

DATASET = "VOYAGER-1_LECP_LEV-2-DAILY-AVG"
START = "20130102T000000Z"
STOP = "20210526T000000Z"
VARIABLES = "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties"


def fetch(opener: urllib.request.OpenerDirector, url: str) -> bytes:
    with opener.open(url, timeout=240) as response:
        return response.read()


def main() -> None:
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    if OUTPUT_FILE.is_file() and OUTPUT_FILE.stat().st_size > 1000:
        print(json.dumps({
            "status": "existing",
            "file": str(OUTPUT_FILE),
            "bytes": OUTPUT_FILE.stat().st_size,
        }, indent=2))
        return

    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    request_url = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        f"{DATASET}/data/{START},{STOP}/{VARIABLES}?format=cdf"
    )
    xml_data = fetch(opener, request_url)
    root = ET.fromstring(xml_data)
    node = root.find(".//{http://cdaweb.gsfc.nasa.gov/schema}Name")
    if node is None or not node.text:
        raise RuntimeError(f"CDAWeb returned no generated CDF URL: {xml_data[:500]!r}")
    cdf_url = node.text.strip()
    payload = fetch(opener, cdf_url)
    if len(payload) < 1000:
        raise RuntimeError(f"Downloaded CDF is unexpectedly small: {len(payload)}")
    OUTPUT_FILE.write_bytes(payload)

    metadata = {
        "dataset": DATASET,
        "time_interval": f"{START}/{STOP}",
        "variables": VARIABLES.split(","),
        "request": request_url,
        "generated_cdf": cdf_url,
        "file": str(OUTPUT_FILE),
        "bytes": len(payload),
        "processing": (
            "Official NASA daily sector flux; no interpolation, fit, "
            "or added background correction"
        ),
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

