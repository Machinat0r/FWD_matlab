#!/usr/bin/env python3
"""Download Voyager 1 data for Florinski et al. (2008), Figure 4."""

from __future__ import annotations

import json
import pathlib
import urllib.request
import xml.etree.ElementTree as ET


PROGRAM_ROOT = pathlib.Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DATABASE_ROOT = MONTHLY_ROOT / "Voyager1_Selected_Event_Data" / "voyager1"

LECP_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1d"
    / "l2"
    / "sectored_flux"
    / "2004"
    / "08"
    / "voyager1_lecp_lev2_daily_sectored_20040801_20040831.cdf"
)
LECP_RATES_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "level1"
    / "sectored_rates"
    / "2004"
    / "08"
    / "voyager1_lecp_lev1_sectored_rates_20040801_20040831.cdf"
)
MAG_FILE = (
    DATABASE_ROOT
    / "coho"
    / "1hr"
    / "l2"
    / "merged_mag_plasma"
    / "2004"
    / "08"
    / "voyager1_coho1hr_merged_mag_plasma_20040801_v01.cdf"
)
METADATA_FILE = (
    MONTHLY_ROOT
    / "Voyager1_2004_Florinski_Figure4_PitchAngle_Spectrogram"
    / "voyager1_2004_figure4_download_metadata.json"
)


def fetch(opener: urllib.request.OpenerDirector, url: str) -> bytes:
    with opener.open(url, timeout=120) as response:
        return response.read()


def cdaweb_cdf_url(
    opener: urllib.request.OpenerDirector, request_url: str
) -> str:
    root = ET.fromstring(fetch(opener, request_url))
    node = root.find(".//{http://cdaweb.gsfc.nasa.gov/schema}Name")
    if node is None or not node.text:
        raise RuntimeError("CDAWeb did not return a generated CDF URL")
    return node.text.strip()


def main() -> None:
    LECP_FILE.parent.mkdir(parents=True, exist_ok=True)
    LECP_RATES_FILE.parent.mkdir(parents=True, exist_ok=True)
    MAG_FILE.parent.mkdir(parents=True, exist_ok=True)
    METADATA_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Direct API download; no browser session is used.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    lecp_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER-1_LECP_LEV-2-DAILY-AVG/data/"
        "20040801T000000Z,20040901T000000Z/"
        "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties?format=cdf"
    )
    lecp_url = cdaweb_cdf_url(opener, lecp_request)
    LECP_FILE.write_bytes(fetch(opener, lecp_url))

    # Level 2 has no P1 sector records on DOY 218-227. Level 1 contains
    # the original P1 sector count rates across the Figure 4 interval.
    lecp_rates_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER-1_LECP_LEV-1-RATES/data/"
        "20040801T000000Z,20040901T000000Z/"
        "FHDU_SectoredRates,FHDU_SectoredRateUncertainties?format=cdf"
    )
    lecp_rates_url = cdaweb_cdf_url(opener, lecp_rates_request)
    LECP_RATES_FILE.write_bytes(fetch(opener, lecp_rates_url))

    mag_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA/data/"
        "20040801T000000Z,20040901T000000Z/"
        "ABS_B,BR,BT,BN,V,protonDensity,protonTemp,"
        "protonFlux1_LECP,protonFlux2_LECP,protonFlux3_LECP?format=cdf"
    )
    mag_url = cdaweb_cdf_url(opener, mag_request)
    MAG_FILE.write_bytes(fetch(opener, mag_url))

    metadata = {
        "paper": "Florinski, Decker, and le Roux (2008), JGR, 113, A07103",
        "doi": "10.1029/2007JA012859",
        "figure": 4,
        "figure_days": "2004 DOY 221-224",
        "plot_interval": "2004-08-05/2004-08-15 (end exclusive)",
        "LECP_dataset": "VOYAGER-1_LECP_LEV-2-DAILY-AVG",
        "LECP_request": lecp_request,
        "LECP_generated_file": lecp_url,
        "LECP_database_file": str(LECP_FILE),
        "LECP_level1_rates_dataset": "VOYAGER-1_LECP_LEV-1-RATES",
        "LECP_level1_rates_request": lecp_rates_request,
        "LECP_level1_rates_generated_file": lecp_rates_url,
        "LECP_level1_rates_database_file": str(LECP_RATES_FILE),
        "MAG_dataset": "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA",
        "MAG_request": mag_request,
        "MAG_generated_file": mag_url,
        "MAG_database_file": str(MAG_FILE),
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print(f"LECP CDF bytes: {LECP_FILE.stat().st_size}")
    print(f"LECP Level-1 rates CDF bytes: {LECP_RATES_FILE.stat().st_size}")
    print(f"MAG CDF bytes: {MAG_FILE.stat().st_size}")


if __name__ == "__main__":
    main()
