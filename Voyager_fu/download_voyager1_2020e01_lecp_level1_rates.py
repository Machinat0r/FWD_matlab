#!/usr/bin/env python3
"""Download Voyager 1 LECP Level-1 sectored rates for event 2020-E01.

This diagnostic file preserves the public count rates.  Florinski et al.
(2008) used a separate LECP-team penetrating-particle background correction;
the paper does not publish that correction series or an automatic algorithm.
S8 can help diagnose the background, but its same-day value must not be
silently treated as the exact correction.  The CDF is obtained directly from
the official NASA CDAWeb API; no browser is used.
"""

from __future__ import annotations

import json
import pathlib
import urllib.request
import xml.etree.ElementTree as ET


PROGRAM_ROOT = pathlib.Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DOWNLOAD_ROOT = MONTHLY_ROOT / "Voyager1_LECP_Sectored_Level1_Rates"
DATABASE_ROOT = MONTHLY_ROOT / "Voyager1_Selected_Event_Data" / "voyager1"

CDF_NAME = "voyager1_lecp_p1_level1_sectored_rates_20200622_20200908.cdf"
CONVENIENCE_FILE = DOWNLOAD_ROOT / CDF_NAME
DATABASE_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "level1"
    / "sectored_rates"
    / "2020"
    / "06-09"
    / CDF_NAME
)
METADATA_FILE = DOWNLOAD_ROOT / "voyager1_2020e01_level1_rates_metadata.json"

DATASET = "VOYAGER-1_LECP_LEV-1-RATES"
START = "20200622T000000Z"
STOP = "20200909T000000Z"
VARIABLES = "FHDU_SectoredRates,FHDU_SectoredRateUncertainties"


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
        "paper_method": "Florinski, Decker, and le Roux (2008), JGR 113 A07103",
        "doi": "10.1029/2007JA012859",
        "dataset": DATASET,
        "time_interval": f"{START}/{STOP}",
        "variables": VARIABLES.split(","),
        "request": request_url,
        "generated_cdf": cdf_url,
        "convenience_file": str(CONVENIENCE_FILE),
        "database_file": str(DATABASE_FILE),
        "bytes": len(payload),
        "processing": "Original NASA Level-1 sectored rates; no averaging, interpolation, or background correction",
        "background_caveat": "S8 is a diagnostic; the paper's event-specific LECP-team correction is not public",
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
