#!/usr/bin/env python3
"""Download the Voyager 1 data needed for the 2002 DOY 221/225 PADs."""

from __future__ import annotations

import json
import pathlib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


PROGRAM_ROOT = pathlib.Path(__file__).resolve().parent
MONTHLY_ROOT = PROGRAM_ROOT / "Voyager_Interstellar_Monthly"
DOWNLOAD_ROOT = MONTHLY_ROOT / "Voyager1_LECP_Sectored_Daily"
DATABASE_ROOT = MONTHLY_ROOT / "Voyager1_Selected_Event_Data" / "voyager1"

LECP_CONVENIENCE_FILE = (
    DOWNLOAD_ROOT / "voyager1_lecp_p1_sectored_daily_20020808_20020815.cdf"
)
MAG_CSV_FILE = DOWNLOAD_ROOT / "voyager1_coho1hr_mag_20020808_20020815.csv"
METADATA_FILE = DOWNLOAD_ROOT / "voyager1_pad_2002_download_metadata.json"
LECP_DATABASE_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1d"
    / "l2"
    / "sectored_flux"
    / "2002"
    / "08"
    / "voyager1_lecp_lev2_daily_sectored_20020808_20020814.cdf"
)
MAG_DATABASE_FILE = (
    DATABASE_ROOT
    / "coho"
    / "1hr"
    / "l2"
    / "merged_mag_plasma"
    / "2002"
    / "08"
    / "voyager1_coho1hr_merged_mag_plasma_20020808_20020815.cdf"
)
LECP_SPECTROGRAM_DATABASE_FILE = (
    DATABASE_ROOT
    / "lecp"
    / "1d"
    / "l2"
    / "sectored_flux"
    / "2002"
    / "08"
    / "voyager1_lecp_lev2_daily_sectored_20020801_20020831.cdf"
)
MAG_MONTHLY_DATABASE_FILE = (
    DATABASE_ROOT
    / "coho"
    / "1hr"
    / "l2"
    / "merged_mag_plasma"
    / "2002"
    / "08"
    / "voyager1_coho1hr_merged_mag_plasma_20020801_v01.cdf"
)


def fetch(opener: urllib.request.OpenerDirector, url: str) -> bytes:
    with opener.open(url, timeout=120) as response:
        return response.read()


def cdaweb_cdf_url(
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
    LECP_DATABASE_FILE.parent.mkdir(parents=True, exist_ok=True)
    MAG_DATABASE_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Direct CDAWeb access avoids dependence on a browser download session.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    lecp_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER-1_LECP_LEV-2-DAILY-AVG/data/"
        "20020808T000000Z,20020815T000000Z/"
        "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties?format=cdf"
    )
    lecp_generated_url = cdaweb_cdf_url(opener, lecp_request)
    lecp_bytes = fetch(opener, lecp_generated_url)
    LECP_CONVENIENCE_FILE.write_bytes(lecp_bytes)
    LECP_DATABASE_FILE.write_bytes(lecp_bytes)

    mag_query = urllib.parse.urlencode(
        {
            "id": "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA",
            "parameters": "Time,ABS_B,BR,BT,BN",
            "time.min": "2002-08-08T00:00:00Z",
            "time.max": "2002-08-15T00:00:00Z",
            "format": "csv",
        }
    )
    mag_hapi_url = "https://cdaweb.gsfc.nasa.gov/hapi/data?" + mag_query
    MAG_CSV_FILE.write_bytes(fetch(opener, mag_hapi_url))

    mag_cdf_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA/data/"
        "20020808T000000Z,20020815T000000Z/ABS_B,BR,BT,BN?format=cdf"
    )
    mag_generated_url = cdaweb_cdf_url(opener, mag_cdf_request)
    MAG_DATABASE_FILE.write_bytes(fetch(opener, mag_generated_url))

    lecp_spectrogram_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER-1_LECP_LEV-2-DAILY-AVG/data/"
        "20020801T000000Z,20020901T000000Z/"
        "FHDU_SectoredFluxes,FHDU_SectoredFluxUncertainties?format=cdf"
    )
    lecp_spectrogram_url = cdaweb_cdf_url(opener, lecp_spectrogram_request)
    LECP_SPECTROGRAM_DATABASE_FILE.write_bytes(
        fetch(opener, lecp_spectrogram_url)
    )

    mag_monthly_request = (
        "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA/data/"
        "20020801T000000Z,20020901T000000Z/"
        "ABS_B,BR,BT,BN,V,protonDensity,protonTemp,"
        "protonFlux1_LECP,protonFlux2_LECP,protonFlux3_LECP?format=cdf"
    )
    mag_monthly_url = cdaweb_cdf_url(opener, mag_monthly_request)
    MAG_MONTHLY_DATABASE_FILE.write_bytes(fetch(opener, mag_monthly_url))

    metadata = {
        "time_interval": "2002-08-08T00:00:00Z/2002-08-15T00:00:00Z",
        "LECP_dataset": "VOYAGER-1_LECP_LEV-2-DAILY-AVG",
        "LECP_request": lecp_request,
        "LECP_generated_file": lecp_generated_url,
        "LECP_database_file": str(LECP_DATABASE_FILE),
        "MAG_dataset": "VOYAGER1_COHO1HR_MERGED_MAG_PLASMA",
        "MAG_HAPI_request": mag_hapi_url,
        "MAG_CDF_request": mag_cdf_request,
        "MAG_generated_CDF": mag_generated_url,
        "MAG_database_file": str(MAG_DATABASE_FILE),
        "LECP_spectrogram_request": lecp_spectrogram_request,
        "LECP_spectrogram_generated_file": lecp_spectrogram_url,
        "LECP_spectrogram_database_file": str(
            LECP_SPECTROGRAM_DATABASE_FILE
        ),
        "MAG_monthly_request": mag_monthly_request,
        "MAG_monthly_generated_file": mag_monthly_url,
        "MAG_monthly_database_file": str(MAG_MONTHLY_DATABASE_FILE),
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print(f"LECP CDF bytes: {LECP_DATABASE_FILE.stat().st_size}")
    print(f"MAG CDF bytes: {MAG_DATABASE_FILE.stat().st_size}")
    print(
        "LECP spectrogram CDF bytes: "
        f"{LECP_SPECTROGRAM_DATABASE_FILE.stat().st_size}"
    )
    print(f"MAG monthly CDF bytes: {MAG_MONTHLY_DATABASE_FILE.stat().st_size}")


if __name__ == "__main__":
    main()
