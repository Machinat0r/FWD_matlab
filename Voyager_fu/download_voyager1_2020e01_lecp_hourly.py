#!/usr/bin/env python3
"""Download Voyager 1 LECP hourly sectored P1 flux for event 2020-E01.

The request uses the NASA CDAWeb web service directly.  The generated CDF is
stored both beside the plotting programs and in the local Voyager event-data
database.  No browser session, resampling, interpolation, or time averaging is
performed by this downloader.
"""

from __future__ import annotations

import json
import pathlib
import urllib.request
import xml.etree.ElementTree as ET


PROGRAM_ROOT = pathlib.Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DOWNLOAD_ROOT = MONTHLY_ROOT / "Voyager1_LECP_Sectored_Hourly"
DATABASE_ROOT = MONTHLY_ROOT / "Voyager1_Selected_Event_Data" / "voyager1"

CDF_NAME = "voyager1_lecp_p1_sectored_hourly_20200622_20200908.cdf"
CONVENIENCE_FILE = DOWNLOAD_ROOT / CDF_NAME
DATABASE_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1hr"
    / "l2"
    / "sectored_flux"
    / "2020"
    / "06-09"
    / CDF_NAME
)
METADATA_FILE = DOWNLOAD_ROOT / "voyager1_2020e01_lecp_hourly_metadata.json"

DATASET = "VOYAGER-1_LECP_LEV-2-HOURLY-AVG"
START = "20200622T000000Z"
STOP = "20200909T000000Z"
VARIABLES = "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties"


def fetch(opener: urllib.request.OpenerDirector, url: str) -> bytes:
    with opener.open(url, timeout=180) as response:
        return response.read()


def generated_cdf_url(
    opener: urllib.request.OpenerDirector, request_url: str
) -> str:
    xml_data = fetch(opener, request_url)
    root = ET.fromstring(xml_data)
    node = root.find(".//{http://cdaweb.gsfc.nasa.gov/schema}Name")
    if node is None or not node.text:
        raise RuntimeError(f"CDAWeb returned no CDF URL: {xml_data[:500]!r}")
    return node.text.strip()


def main() -> None:
    DOWNLOAD_ROOT.mkdir(parents=True, exist_ok=True)
    DATABASE_FILE.parent.mkdir(parents=True, exist_ok=True)

    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    request_url = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        f"{DATASET}/data/{START},{STOP}/{VARIABLES}?format=cdf"
    )
    cdf_url = generated_cdf_url(opener, request_url)
    payload = fetch(opener, cdf_url)
    if len(payload) < 1000:
        raise RuntimeError(f"Downloaded CDF is unexpectedly small: {len(payload)} bytes")

    CONVENIENCE_FILE.write_bytes(payload)
    DATABASE_FILE.write_bytes(payload)

    metadata = {
        "dataset": DATASET,
        "time_interval": f"{START}/{STOP}",
        "variables": VARIABLES.split(","),
        "request": request_url,
        "generated_cdf": cdf_url,
        "convenience_file": str(CONVENIENCE_FILE),
        "database_file": str(DATABASE_FILE),
        "bytes": len(payload),
        "processing": "NASA hourly averages; no additional averaging or interpolation",
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
