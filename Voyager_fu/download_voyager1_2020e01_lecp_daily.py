#!/usr/bin/env python3
"""Download official daily Voyager 1 LECP P1 sectored flux for 2020-E01.

The file is requested directly from NASA CDAWeb, without browser automation.
CDAWeb's public Level-2 flux is preserved as supplied.  The LECP archive warns
that its automated public products do not contain the event-specific
penetrating-particle background corrections used in research analyses.
"""

from __future__ import annotations

import json
import pathlib
import urllib.request
import xml.etree.ElementTree as ET


PROGRAM_ROOT = pathlib.Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DOWNLOAD_ROOT = MONTHLY_ROOT / "Voyager1_LECP_Sectored_Daily"
DATABASE_ROOT = MONTHLY_ROOT / "Voyager1_Selected_Event_Data" / "voyager1"

CDF_NAME = "voyager1_lecp_p1_sectored_daily_20200622_20200908.cdf"
CONVENIENCE_FILE = DOWNLOAD_ROOT / CDF_NAME
DATABASE_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1d"
    / "l2"
    / "sectored_flux"
    / "2020"
    / "06-09"
    / CDF_NAME
)
METADATA_FILE = DOWNLOAD_ROOT / "voyager1_2020e01_lecp_daily_metadata.json"

DATASET = "VOYAGER-1_LECP_LEV-2-DAILY-AVG"
START = "20200622T000000Z"
STOP = "20200909T000000Z"
VARIABLES = "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties"


def fetch(opener: urllib.request.OpenerDirector, url: str) -> bytes:
    with opener.open(url, timeout=180) as response:
        return response.read()


def generated_cdf_url(
    opener: urllib.request.OpenerDirector, request_url: str
) -> str:
    root = ET.fromstring(fetch(opener, request_url))
    node = root.find(".//{http://cdaweb.gsfc.nasa.gov/schema}Name")
    if node is None or not node.text:
        raise RuntimeError("CDAWeb returned no generated CDF URL")
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
        raise RuntimeError(f"Downloaded CDF is unexpectedly small: {len(payload)}")

    CONVENIENCE_FILE.write_bytes(payload)
    DATABASE_FILE.write_bytes(payload)
    metadata = {
        "paper_method": (
            "Florinski, Decker, and le Roux (2008), "
            "doi:10.1029/2007JA012859"
        ),
        "dataset": DATASET,
        "time_interval": f"{START}/{STOP}",
        "variables": VARIABLES.split(","),
        "request": request_url,
        "generated_cdf": cdf_url,
        "convenience_file": str(CONVENIENCE_FILE),
        "database_file": str(DATABASE_FILE),
        "bytes": len(payload),
        "processing": (
            "Official NASA daily sector flux; no added averaging, "
            "interpolation, fit, or background correction"
        ),
        "background_caveat": (
            "Public automated LECP data are not event-specifically "
            "background corrected; S8 is retained only as a diagnostic"
        ),
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
