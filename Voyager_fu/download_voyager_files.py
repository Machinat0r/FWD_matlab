#!/usr/bin/env python3
"""Download Voyager field, plasma, and ephemeris files from NASA archives.

This backend intentionally uses only the Python standard library.  It is
normally called by ``VoyagerFilesDownload_NAS.m``, but it is also useful as a
standalone command-line downloader.
"""

from __future__ import annotations

import argparse
import calendar
import concurrent.futures
import datetime as dt
import hashlib
import html
from html.parser import HTMLParser
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


LOGGER = logging.getLogger("voyager_downloader")
SPDF_ROOT = "https://spdf.gsfc.nasa.gov/pub/data/voyager"
USER_AGENT = "Voyager-SPDF-downloader/1.0 (+https://spdf.gsfc.nasa.gov/)"
HTTP_TIMEOUT_SECONDS = 60
HTTP_RETRIES = 4
MIN_CDF_BYTES = 32
MIN_GENERIC_BYTES = 8
MAX_THREADS = 5
MANIFEST_NAME = "voyager_download_manifest.json"
DEFAULT_PRODUCTS = ("coho1hr", "position1day")
HIGHRES_PRODUCTS = (
    "mag2s",
    "mag48s_vim",
    "plasma_fine",
    "position1hr",
    "spice_spk",
)
OPTIONAL_PRODUCTS = ("mag48s",)
EXPERIMENTAL_PRODUCTS = (
    "mag2s_unreviewed_primary",
    "mag2s_unreviewed_secondary",
)
PRODUCT_ORDER = (
    DEFAULT_PRODUCTS + HIGHRES_PRODUCTS + OPTIONAL_PRODUCTS + EXPERIMENTAL_PRODUCTS
)
MISSION_START = {
    1: dt.date(1977, 9, 5),
    2: dt.date(1977, 8, 20),
}
LEGACY_MAG_LAST_YEAR = {1: 1991, 2: 1990}
PLS_SOLAR_WIND_LAST_YEAR = {1: 1980, 2: 2007}

PRODUCTS: Mapping[str, Mapping[str, object]] = {
    "coho1hr": {
        "description": (
            "One-hour COHO merged magnetic field and plasma CDFs "
            "(field, density, bulk speed, and temperature where available)"
        ),
        "cadence": "1 hour",
        "remote": "coho1hr_magplasma/YYYY/",
        "local": "coho/1hr/l2/merged_mag_plasma/YYYY/MM/",
    },
    "position1day": {
        "description": "One-day heliocentric Voyager position/ephemeris CDF",
        "cadence": "1 day",
        "remote": "helio1day/",
        "local": "ephemeris/1day/l2/helio_position/YYYY/MM/",
    },
    "mag2s": {
        "description": (
            "Highest-cadence legacy magnetic-field CDFs (1.92 s averages; "
            "calibrated but unreviewed)"
        ),
        "cadence": "1.92 seconds",
        "remote": "magnetic_fields_cdaweb/mag_2s/YYYY/",
        "local": "mag/1.92s/calibrated_unreviewed/YYYY/MM/",
    },
    "mag48s": {
        "description": (
            "Legacy 48-second magnetic-field averages "
            "(calibrated but unreviewed)"
        ),
        "cadence": "48 seconds",
        "remote": "magnetic_fields_cdaweb/mag_48s/YYYY/",
        "local": "mag/48s/calibrated_unreviewed/YYYY/MM/",
    },
    "mag48s_vim": {
        "description": "Reviewed Voyager Interstellar Mission 48-second MAG CDFs",
        "cadence": "48 seconds",
        "remote": "magnetic_fields_cdaweb/vim_48secmag/",
        "local": "mag/48s/reviewed_vim/YYYY/",
    },
    "mag2s_unreviewed_primary": {
        "description": (
            "Post-1991 1.92-second MAG CDFs from the primary/out-board "
            "sensor; NASA marks these data as generally not science quality"
        ),
        "cadence": "1.92 seconds",
        "remote": (
            "magnetic_fields_cdaweb/hires1991_2030/primary/mag_2s/YYYY/"
        ),
        "local": (
            "mag/1.92s/experimental_unreviewed_post1991/primary/YYYY/MM/"
        ),
    },
    "mag2s_unreviewed_secondary": {
        "description": (
            "Post-1991 1.92-second MAG CDFs from the secondary/in-board "
            "sensor; NASA marks these data as generally not science quality"
        ),
        "cadence": "1.92 seconds",
        "remote": (
            "magnetic_fields_cdaweb/hires1991_2030/secondary/mag_2s/YYYY/"
        ),
        "local": (
            "mag/1.92s/experimental_unreviewed_post1991/secondary/YYYY/MM/"
        ),
    },
    "plasma_fine": {
        "description": (
            "PLS high-resolution fitted plasma CDFs, including the "
            "Voyager 2 heliosheath extension"
        ),
        "cadence": "12 to 192 seconds (telemetry-mode dependent)",
        "remote": "plasma_cdaweb/hires_plasma/",
        "local": "pls/hires/l3/{solar_wind|heliosheath}/YYYY/",
    },
    "position1hr": {
        "description": "One-hour SPDF HelioWeb heliocentric position CDF",
        "cadence": "1 hour",
        "remote": "helio1hr/",
        "local": "ephemeris/1hr/l2/helio_position/YYYY/MM/",
    },
    "spice_spk": {
        "description": "JPL/NAIF SPICE spacecraft ephemeris kernels through 2100",
        "cadence": "continuous polynomial ephemeris",
        "remote": "https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/",
        "local": "ephemeris/spice/kernels/",
    },
}

COHO_RE_TEMPLATE = (
    r"^voyager{spacecraft}_coho1hr_merged_mag_plasma_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
POSITION_RE_TEMPLATE = (
    r"^voyager{spacecraft}_helio1day_position_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
POSITION_1HR_RE_TEMPLATE = (
    r"^voyager{spacecraft}_helio1hr_position_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
MAG_RE_TEMPLATE = (
    r"^voyager{spacecraft}_(?P<cadence>2s|48s)_mag_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
MAG_VIM_RE_TEMPLATE = (
    r"^voyager{spacecraft}_48s_mag-vim_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
MAG_POST1991_RE_TEMPLATE = (
    r"^voyager{spacecraft}_2s_mag-(?P<sensor>pri|sec)_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
PLS_RE_TEMPLATE = (
    r"^voyager{spacecraft}_pls_hires_plasma_data_"
    r"(?P<date>\d{{8}})_v(?P<version>\d+)\.cdf$"
)
PLS_HSH_RE = re.compile(
    r"^voyager2_pls_hires_plasma_data_hsh_"
    r"(?P<date>\d{8})_v(?P<version>\d+)\.cdf$",
    re.IGNORECASE,
)

NAIF_SPK_ROOT = "https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/"
NAIF_LSK_ROOT = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/"
NAIF_PLANET_SPK_ROOT = (
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
)
JPL_SSD_SPACECRAFT_ROOT = "https://ssd.jpl.nasa.gov/ftp/eph/spacecraft/"

SPICE_SHARED_FILES = (
    {
        "filename": "naif0012.tls",
        "url": urllib.parse.urljoin(NAIF_LSK_ROOT, "naif0012.tls"),
        "size": 5257,
        "sha256": (
            "678E32BDB5A744117A467CD9601CD6B373F0E9BC9BBDE1371D5EEE39600A039B"
        ),
        "format": "text",
        "quality": "NAIF leap-seconds kernel",
    },
    {
        "filename": "de440s.bsp",
        "url": urllib.parse.urljoin(NAIF_PLANET_SPK_ROOT, "de440s.bsp"),
        "size": 32726016,
        "sha256": (
            "C1C7FEEAB882263FC493A9D5A5B2DDD71B54826CDF65D8D17A76126B260A49F2"
        ),
        "format": "spice",
        "quality": "JPL DE440 short planetary ephemeris",
    },
)

SPICE_FILES = {
    1: (
        {
            "filename": "Voyager_1.a54206u_V0.2_merged.bsp",
            "url": urllib.parse.urljoin(
                NAIF_SPK_ROOT, "Voyager_1.a54206u_V0.2_merged.bsp"
            ),
            "size": 6374400,
            "sha256": (
                "47C6F2BE03668B50A1EFB5F96978A2B68B2B501DAE6F585841B1569BAA3F4311"
            ),
            "quality": "broad reconstructed Voyager 1 trajectory",
        },
        {
            "filename": "vgr1.x2100.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr1.x2100.bsp"),
            "size": 16286720,
            "sha256": (
                "AA8F37673B5EF60296D210566D01CC49F7B32711F37003393C3327616C9FC6D0"
            ),
            "quality": "2022 JPL fit and long-range extension through 2100",
        },
        {
            "filename": "vgr1_jup230.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr1_jup230.bsp"),
            "size": 325632,
            "sha256": (
                "E1EA3F72F19B15508BC45979771A36A97D02F33056B76867D444304CB82205C9"
            ),
            "quality": "Jupiter encounter reconstruction",
        },
        {
            "filename": "vgr1.sat427.bsp",
            "url": urllib.parse.urljoin(
                JPL_SSD_SPACECRAFT_ROOT, "vgr1.sat427.bsp"
            ),
            "size": 696320,
            "sha256": (
                "D45C923038BD4C0463A13BD0F9751BE5C6911E94C0C5D0FCB9D15B79B7D89659"
            ),
            "quality": "SAT427 Saturn encounter reconstruction",
        },
    ),
    2: (
        {
            "filename": "Voyager_2.m05016u.merged.bsp",
            "url": urllib.parse.urljoin(
                NAIF_SPK_ROOT, "Voyager_2.m05016u.merged.bsp"
            ),
            "size": 6447104,
            "sha256": (
                "CE66CBA12CF77BF3A1097F44142EF978F46656788ED08F9052238A102ED2E706"
            ),
            "quality": "broad reconstructed Voyager 2 trajectory",
        },
        {
            "filename": "vgr2.x2100.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr2.x2100.bsp"),
            "size": 15744000,
            "sha256": (
                "391F5361773322E2242FB1C81578AEB1FBF9134C3878B2F325E9B8B5F46AD7CB"
            ),
            "quality": "2022 JPL fit and long-range extension through 2100",
        },
        {
            "filename": "vgr2_jup230.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr2_jup230.bsp"),
            "size": 330752,
            "sha256": (
                "9C00BE3C83915F6C1FD8448D9266420E3C462E9D78E4C32B17145DAC81529D5A"
            ),
            "quality": "Jupiter encounter reconstruction",
        },
        {
            "filename": "vgr2.sat427.bsp",
            "url": urllib.parse.urljoin(
                JPL_SSD_SPACECRAFT_ROOT, "vgr2.sat427.bsp"
            ),
            "size": 834560,
            "sha256": (
                "E1977B3EC78704E5CE4892D2ED682BE04EEE469F57801430E6BB88E5A96FEF8C"
            ),
            "quality": "SAT427 Saturn encounter reconstruction",
        },
        {
            "filename": "vgr2.ura182.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr2.ura182.bsp"),
            "size": 1064960,
            "sha256": (
                "22D012AB426273FFB22710238CFE4F7D0C63850225CC9E46B1F5AA40308CD880"
            ),
            "quality": "DE442/URA182 Uranus encounter reconstruction",
        },
        {
            "filename": "vgr2_nep097.bsp",
            "url": urllib.parse.urljoin(NAIF_SPK_ROOT, "vgr2_nep097.bsp"),
            "size": 1469440,
            "sha256": (
                "1AEE50E5E4A5253C660D06F7359EB4DAC2FE4A5DAB5FDC821D2B868469BDC473"
            ),
            "quality": "DE440/NEP097 Neptune encounter reconstruction",
        },
    ),
}

# CDF 3.x and 2.x files begin with one of these magic values.  The next
# 32-bit word indicates an uncompressed or compressed CDF container.
CDF_FIRST_MAGICS = {b"\xcd\xf3\x00\x01", b"\xcd\xf2\x60\x02"}
CDF_SECOND_MAGICS = {b"\x00\x00\xff\xff", b"\xcc\xcc\x00\x01"}


class DiscoveryError(RuntimeError):
    """Raised when an SPDF directory cannot be discovered reliably."""


class DownloadError(RuntimeError):
    """Raised when a remote CDF cannot be downloaded or validated."""


class LinkCollector(HTMLParser):
    """Collect href values from a simple Apache-style directory listing."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.links: List[str] = []

    def handle_starttag(
        self, tag: str, attrs: List[Tuple[str, Optional[str]]]
    ) -> None:
        if tag.lower() != "a":
            return
        for name, value in attrs:
            if name.lower() == "href" and value:
                self.links.append(value)


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp suitable for a manifest."""

    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def log(message: str, *args: object, level: int = logging.INFO) -> None:
    """Write operational output to stderr through the logger."""

    LOGGER.log(level, message, *args)


def parse_bool(value: object) -> bool:
    """Parse MATLAB-friendly booleans as well as conventional CLI values."""

    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(
        "expected one of: 1, 0, true, false, yes, no, on, off"
    )


def parse_manifest_name(value: str) -> str:
    """Require a single safe JSON filename for the archive manifest."""

    name = value.strip()
    if (
        not name
        or name in {".", ".."}
        or re.search(r'[\\/:*?"<>|\x00]', name)
        or Path(name).name != name
    ):
        raise argparse.ArgumentTypeError(
            "--manifest-name must be a filename without directory components"
        )
    return name


def parse_date_value(value: str) -> dt.date:
    """Parse an ISO or compact calendar date."""

    stripped = value.strip()
    for date_format in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return dt.datetime.strptime(stripped, date_format).date()
        except ValueError:
            pass
    raise argparse.ArgumentTypeError(
        "dates must use YYYY-MM-DD (compact YYYYMMDD is also accepted)"
    )


def parse_date_range(value: str) -> Tuple[dt.date, dt.date]:
    """Parse an inclusive START/END interval."""

    parts = value.strip().split("/")
    if len(parts) != 2 or not all(part.strip() for part in parts):
        raise argparse.ArgumentTypeError("--date must be START/END")
    start = parse_date_value(parts[0])
    end = parse_date_value(parts[1])
    if start > end:
        raise argparse.ArgumentTypeError("--date START must not be after END")
    return start, end


def split_cli_values(values: Sequence[str]) -> List[str]:
    """Split repeated, comma-, semicolon-, or whitespace-separated values."""

    output: List[str] = []
    for value in values:
        output.extend(item for item in re.split(r"[,;\s]+", value.strip()) if item)
    return output


def parse_spacecraft(values: Sequence[str]) -> List[int]:
    """Normalize the spacecraft selection while retaining deterministic order."""

    raw = split_cli_values(values)
    if not raw:
        raise argparse.ArgumentTypeError("--spacecraft cannot be empty")
    result: List[int] = []
    for item in raw:
        try:
            number = int(item)
        except ValueError as exception:
            raise argparse.ArgumentTypeError(
                "--spacecraft accepts only 1 and/or 2"
            ) from exception
        if number not in (1, 2):
            raise argparse.ArgumentTypeError("--spacecraft accepts only 1 and/or 2")
        if number not in result:
            result.append(number)
    return result


def parse_products(values: Sequence[str]) -> List[str]:
    """Normalize requested products and aliases."""

    raw = [item.lower() for item in split_cli_values(values)]
    if not raw:
        raise argparse.ArgumentTypeError("--product cannot be empty")
    expanded: List[str] = []
    for item in raw:
        if item == "all":
            expanded.extend(PRODUCT_ORDER)
        elif item == "highres":
            expanded.extend(HIGHRES_PRODUCTS)
        elif item == "both":
            expanded.extend(DEFAULT_PRODUCTS)
        elif item in {"mag2s_unreviewed", "experimental_mag"}:
            expanded.extend(EXPERIMENTAL_PRODUCTS)
        else:
            expanded.append(item)
    unknown = [item for item in expanded if item not in PRODUCTS]
    if unknown:
        raise argparse.ArgumentTypeError(
            "unknown product(s): {0}; choices are {1}".format(
                ", ".join(unknown), ", ".join(PRODUCT_ORDER)
            )
        )
    return [product for product in PRODUCT_ORDER if product in expanded]


def make_request(url: str, method: str = "GET") -> urllib.request.Request:
    """Create an identity-encoded SPDF request."""

    return urllib.request.Request(
        url,
        method=method,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "*/*",
            "Accept-Encoding": "identity",
        },
    )


def force_powershell_http() -> bool:
    """Return whether the explicit Windows HTTPS fallback is enabled."""

    return str(
        os.environ.get("VOYAGER_FORCE_POWERSHELL_HTTP", "")
    ).strip().lower() in {"1", "true", "yes", "on"}


def fetch_listing_powershell(url: str) -> str:
    """Read an HTTPS directory index through Windows PowerShell."""

    powershell = (
        Path(os.environ.get("WINDIR", r"C:\Windows"))
        / "System32"
        / "WindowsPowerShell"
        / "v1.0"
        / "powershell.exe"
    )
    escaped_url = url.replace("'", "''")
    command = (
        "$ProgressPreference='SilentlyContinue'; "
        "$ErrorActionPreference='Stop'; "
        "[Console]::OutputEncoding="
        "[System.Text.UTF8Encoding]::new($false); "
        "$r=Invoke-WebRequest -UseBasicParsing -TimeoutSec 180 "
        "-Uri '{0}'; [Console]::Out.Write($r.Content)"
    ).format(escaped_url)
    last_error = ""
    for attempt in range(HTTP_RETRIES):
        try:
            completed = subprocess.run(
                [
                    str(powershell),
                    "-NoProfile",
                    "-NonInteractive",
                    "-Command",
                    command,
                ],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=240,
                check=False,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if completed.returncode == 0:
                return completed.stdout
            last_error = (
                completed.stderr or completed.stdout or
                "PowerShell Invoke-WebRequest failed"
            ).strip()
        except (OSError, subprocess.SubprocessError) as exception:
            last_error = str(exception)
        if attempt + 1 < HTTP_RETRIES:
            log(
                "PowerShell directory request failed (attempt %d/%d): %s",
                attempt + 1,
                HTTP_RETRIES,
                url,
                level=logging.WARNING,
            )
            time.sleep(retry_delay(attempt))
    raise DiscoveryError(
        "cannot read SPDF directory through PowerShell {0}: {1}".format(
            url, last_error
        )
    )


def retry_delay(attempt: int) -> float:
    """Small deterministic exponential delay for transient HTTP failures."""

    return min(8.0, 0.75 * (2**attempt))


def is_transient_http_error(exception: BaseException) -> bool:
    """Classify HTTP failures for retry."""

    if isinstance(exception, urllib.error.HTTPError):
        return exception.code in {408, 425, 429, 500, 502, 503, 504}
    return isinstance(exception, (urllib.error.URLError, TimeoutError, OSError))


def fetch_listing(url: str, allow_missing: bool = False) -> str:
    """Fetch and decode an SPDF directory index with retries."""

    hostname = (urllib.parse.urlsplit(url).hostname or "").lower()
    if (
        os.name == "nt"
        and force_powershell_http()
        and hostname == "spdf.gsfc.nasa.gov"
    ):
        return fetch_listing_powershell(url)

    last_error: Optional[BaseException] = None
    for attempt in range(HTTP_RETRIES):
        try:
            with urllib.request.urlopen(
                make_request(url), timeout=HTTP_TIMEOUT_SECONDS
            ) as response:
                encoding = response.headers.get_content_charset() or "utf-8"
                return response.read().decode(encoding, errors="replace")
        except urllib.error.HTTPError as exception:
            if allow_missing and exception.code == 404:
                return ""
            last_error = exception
            if not is_transient_http_error(exception):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exception:
            last_error = exception
        if attempt + 1 < HTTP_RETRIES:
            log(
                "Directory request failed (attempt %d/%d): %s",
                attempt + 1,
                HTTP_RETRIES,
                url,
                level=logging.WARNING,
            )
            time.sleep(retry_delay(attempt))
    raise DiscoveryError("cannot read SPDF directory {0}: {1}".format(url, last_error))


def basenames_from_listing(page: str) -> List[str]:
    """Extract safe file basenames from an Apache directory listing."""

    collector = LinkCollector()
    try:
        collector.feed(page)
    except Exception:
        # Apache indexes are regular HTML, but the regex fallback below also
        # handles malformed mirrors or proxy-generated pages.
        pass

    candidates: List[str] = []
    for link in collector.links:
        parsed = urllib.parse.urlsplit(html.unescape(link))
        filename = urllib.parse.unquote(parsed.path.rsplit("/", 1)[-1])
        if filename and "/" not in filename and "\\" not in filename:
            candidates.append(filename)

    return sorted(set(candidates))


def filenames_from_listing(page: str) -> List[str]:
    """Extract CDF basenames, with a fallback for malformed HTML."""

    candidates = [
        name for name in basenames_from_listing(page) if name.lower().endswith(".cdf")
    ]
    # Preserve CDF files even if an intermediary stripped the anchor href.
    candidates.extend(
        urllib.parse.unquote(html.unescape(match))
        for match in re.findall(
            r"""(?i)([a-z0-9][a-z0-9_.-]*\.cdf)(?=["'<\s])""", page
        )
    )
    return sorted(set(candidates))


def years_from_listing(page: str) -> List[int]:
    """Extract four-digit year subdirectories from an archive index."""

    collector = LinkCollector()
    try:
        collector.feed(page)
        collector.close()
    except Exception:
        pass

    years = set()
    for link in collector.links:
        parsed = urllib.parse.urlsplit(html.unescape(link))
        path = urllib.parse.unquote(parsed.path).rstrip("/")
        name = path.rsplit("/", 1)[-1]
        if re.fullmatch(r"(?:19|20)\d{2}", name):
            years.add(int(name))

    years.update(
        int(match)
        for match in re.findall(
            r"""(?i)href\s*=\s*["']((?:19|20)\d{2})/["']""", page
        )
    )
    return sorted(years)


def month_intersects(year: int, month: int, start: dt.date, end: dt.date) -> bool:
    """Return whether a calendar month intersects the inclusive date range."""

    month_start = dt.date(year, month, 1)
    month_end = dt.date(year, month, calendar.monthrange(year, month)[1])
    return month_start <= end and month_end >= start


def make_file_record(
    spacecraft: int,
    product: str,
    filename: str,
    url: str,
    file_date: dt.date,
    version: int,
    file_format: str = "cdf",
    relative_path: Optional[Path] = None,
    quality: Optional[str] = None,
    sensor: Optional[str] = None,
    expected_size: Optional[int] = None,
    expected_sha256: Optional[str] = None,
) -> Dict[str, object]:
    """Build the stable public record used for list output and manifests."""

    year = file_date.year
    month = file_date.month
    if relative_path is not None:
        relative = relative_path
    elif product == "coho1hr":
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "coho"
            / "1hr"
            / "l2"
            / "merged_mag_plasma"
            / "{0:04d}".format(year)
            / "{0:02d}".format(month)
            / filename
        )
    elif product == "position1day":
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "ephemeris"
            / "1day"
            / "l2"
            / "helio_position"
            / "{0:04d}".format(year)
            / "{0:02d}".format(month)
            / filename
        )
    elif product == "position1hr":
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "ephemeris"
            / "1hr"
            / "l2"
            / "helio_position"
            / "{0:04d}".format(year)
            / "{0:02d}".format(month)
            / filename
        )
    elif product in {"mag2s", "mag48s"}:
        cadence = "1.92s" if product == "mag2s" else "48s"
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "mag"
            / cadence
            / "calibrated_unreviewed"
            / "{0:04d}".format(year)
            / "{0:02d}".format(month)
            / filename
        )
    elif product == "mag48s_vim":
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "mag"
            / "48s"
            / "reviewed_vim"
            / "{0:04d}".format(year)
            / filename
        )
    elif product in EXPERIMENTAL_PRODUCTS:
        sensor = (
            "primary"
            if product == "mag2s_unreviewed_primary"
            else "secondary"
        )
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "mag"
            / "1.92s"
            / "experimental_unreviewed_post1991"
            / sensor
            / "{0:04d}".format(year)
            / "{0:02d}".format(month)
            / filename
        )
    elif product == "plasma_fine":
        relative = (
            Path("voyager{0}".format(spacecraft))
            / "pls"
            / "hires"
            / "l3"
            / "solar_wind"
            / "{0:04d}".format(year)
            / filename
        )
    else:
        raise ValueError("relative path is required for product {0}".format(product))
    record: Dict[str, object] = {
        "spacecraft": spacecraft,
        "product": product,
        "date": file_date.isoformat(),
        "year": year,
        "month": month,
        "version": version,
        "filename": filename,
        "url": url,
        "relative_path": relative.as_posix(),
        "format": file_format,
        "size": None,
        "expected_size": expected_size,
    }
    if quality:
        record["quality"] = quality
    if sensor:
        record["sensor"] = sensor
    if product in EXPERIMENTAL_PRODUCTS:
        record.update(
            {
                "sensor_location": (
                    "out-board" if sensor == "primary" else "in-board"
                ),
                "cadence_seconds": 1.92,
                "source_collection": "hires1991_2030",
                "science_quality": False,
                "quality_flag_semantics": "0=good, 1=bad",
            }
        )
    if expected_sha256:
        record["expected_sha256"] = expected_sha256.upper()
    return record


def discover_coho_year(
    spacecraft: int, year: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover highest-version COHO files for one spacecraft/year."""

    start = max(start, MISSION_START[spacecraft])
    if end < start:
        return []

    directory_url = (
        "{0}/voyager{1}/coho1hr_magplasma/{2:04d}/".format(
            SPDF_ROOT, spacecraft, year
        )
    )
    page = fetch_listing(directory_url, allow_missing=True)
    if not page:
        log("No COHO directory for Voyager %d in %04d", spacecraft, year)
        return []

    expression = re.compile(
        COHO_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE
    )
    best_by_month: Dict[
        Tuple[int, int], Tuple[Tuple[int, int, int, str], int, str, dt.date]
    ] = {}
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date.year != year:
            continue
        if not month_intersects(file_date.year, file_date.month, start, end):
            continue
        version = int(match.group("version"))
        key = (file_date.year, file_date.month)
        # A few launch-month directories contain both a nominal day-01
        # monthly file and a second file beginning on the actual launch/data
        # start day.  Keep exactly one file per month: day 01 wins; otherwise
        # use the earliest start day; finally choose the highest version.
        priority = (
            1 if file_date.day == 1 else 0,
            -file_date.day,
            version,
            filename,
        )
        current = best_by_month.get(key)
        if current is None or priority > current[0]:
            best_by_month[key] = (priority, version, filename, file_date)

    records: List[Dict[str, object]] = []
    for _, (_, version, filename, file_date) in sorted(best_by_month.items()):
        records.append(
            make_file_record(
                spacecraft=spacecraft,
                product="coho1hr",
                filename=filename,
                url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
                file_date=file_date,
                version=version,
            )
        )
    return records


def discover_position(
    spacecraft: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover the highest-version full-mission heliocentric position CDF."""

    if end < max(start, MISSION_START[spacecraft]):
        return []
    directory_url = "{0}/voyager{1}/helio1day/".format(SPDF_ROOT, spacecraft)
    page = fetch_listing(directory_url)
    expression = re.compile(
        POSITION_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE
    )
    candidates: List[Tuple[int, str, dt.date]] = []
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date > end:
            continue
        candidates.append((int(match.group("version")), filename, file_date))
    if not candidates:
        return []
    version, filename, file_date = max(candidates, key=lambda item: (item[0], item[1]))
    return [
        make_file_record(
            spacecraft=spacecraft,
            product="position1day",
            filename=filename,
            url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
            file_date=file_date,
            version=version,
        )
    ]


def discover_position_1hr(
    spacecraft: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover the full-mission one-hour HelioWeb position CDF."""

    if end < max(start, MISSION_START[spacecraft]):
        return []
    directory_url = "{0}/voyager{1}/helio1hr/".format(SPDF_ROOT, spacecraft)
    page = fetch_listing(directory_url)
    expression = re.compile(
        POSITION_1HR_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE
    )
    candidates: List[Tuple[int, str, dt.date]] = []
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date <= end:
            candidates.append((int(match.group("version")), filename, file_date))
    if not candidates:
        return []
    version, filename, file_date = max(candidates, key=lambda item: (item[0], item[1]))
    return [
        make_file_record(
            spacecraft=spacecraft,
            product="position1hr",
            filename=filename,
            url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
            file_date=file_date,
            version=version,
            quality="SPDF HelioWeb hourly ephemeris",
        )
    ]


def discover_legacy_mag_year(
    spacecraft: int,
    product: str,
    year: int,
    start: dt.date,
    end: dt.date,
) -> List[Dict[str, object]]:
    """Discover highest-version legacy 1.92 s or 48 s MAG CDFs."""

    directory_name = "mag_2s" if product == "mag2s" else "mag_48s"
    cadence_name = "2s" if product == "mag2s" else "48s"
    directory_url = (
        "{0}/voyager{1}/magnetic_fields_cdaweb/{2}/{3:04d}/".format(
            SPDF_ROOT, spacecraft, directory_name, year
        )
    )
    page = fetch_listing(directory_url, allow_missing=True)
    if not page:
        return []
    expression = re.compile(MAG_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE)
    best_by_date: Dict[dt.date, Tuple[int, str]] = {}
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match or match.group("cadence").lower() != cadence_name:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date.year != year:
            continue
        version = int(match.group("version"))
        current = best_by_date.get(file_date)
        if current is None or (version, filename) > current:
            best_by_date[file_date] = (version, filename)

    selected_dates = [
        file_date for file_date in best_by_date if start <= file_date <= end
    ]
    preceding_dates = [file_date for file_date in best_by_date if file_date < start]
    if preceding_dates:
        predecessor = max(preceding_dates)
        # Most granules span five days; the latest legacy blocks are longer.
        if predecessor + dt.timedelta(days=31) >= start:
            selected_dates.append(predecessor)

    quality = "calibrated legacy MAG; not individually reviewed"
    return [
        make_file_record(
            spacecraft=spacecraft,
            product=product,
            filename=filename,
            url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
            file_date=file_date,
            version=version,
            quality=quality,
        )
        for file_date in sorted(set(selected_dates))
        for version, filename in [best_by_date[file_date]]
    ]


def discover_mag_vim(
    spacecraft: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover reviewed annual/partial-year 48-second VIM MAG CDFs."""

    directory_url = (
        "{0}/voyager{1}/magnetic_fields_cdaweb/vim_48secmag/".format(
            SPDF_ROOT, spacecraft
        )
    )
    page = fetch_listing(directory_url)
    expression = re.compile(
        MAG_VIM_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE
    )
    best_by_date: Dict[dt.date, Tuple[int, str]] = {}
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date.year < start.year or file_date.year > end.year:
            continue
        version = int(match.group("version"))
        current = best_by_date.get(file_date)
        if current is None or (version, filename) > current:
            best_by_date[file_date] = (version, filename)
    return [
        make_file_record(
            spacecraft=spacecraft,
            product="mag48s_vim",
            filename=filename,
            url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
            file_date=file_date,
            version=version,
            quality="calibrated and reviewed VIM magnetic field",
        )
        for file_date, (version, filename) in sorted(best_by_date.items())
    ]


def discover_post1991_mag(
    spacecraft: int,
    product: str,
    start: dt.date,
    end: dt.date,
) -> List[Dict[str, object]]:
    """Discover unreviewed post-1991 1.92-second MAG CDFs for one sensor."""

    if product not in EXPERIMENTAL_PRODUCTS:
        raise ValueError("unsupported post-1991 MAG product {0}".format(product))

    sensor = (
        "primary"
        if product == "mag2s_unreviewed_primary"
        else "secondary"
    )
    sensor_code = "pri" if sensor == "primary" else "sec"
    root_url = (
        "{0}/voyager{1}/magnetic_fields_cdaweb/"
        "hires1991_2030/{2}/mag_2s/"
    ).format(SPDF_ROOT, spacecraft, sensor)
    root_page = fetch_listing(root_url)
    available_years = [
        year
        for year in years_from_listing(root_page)
        if start.year - 1 <= year <= end.year
    ]

    expression = re.compile(
        MAG_POST1991_RE_TEMPLATE.format(spacecraft=spacecraft),
        re.IGNORECASE,
    )
    best_by_date: Dict[dt.date, Tuple[int, str, str]] = {}
    for year in available_years:
        directory_url = urllib.parse.urljoin(root_url, "{0:04d}/".format(year))
        page = fetch_listing(directory_url)
        for filename in filenames_from_listing(page):
            match = expression.fullmatch(filename)
            if not match or match.group("sensor").lower() != sensor_code:
                continue
            try:
                file_date = dt.datetime.strptime(
                    match.group("date"), "%Y%m%d"
                ).date()
            except ValueError:
                continue
            if file_date.year != year:
                continue
            version = int(match.group("version"))
            current = best_by_date.get(file_date)
            candidate = (version, filename, directory_url)
            if current is None or candidate[:2] > current[:2]:
                best_by_date[file_date] = candidate

    selected_dates = [
        file_date for file_date in best_by_date if start <= file_date <= end
    ]
    preceding_dates = [file_date for file_date in best_by_date if file_date < start]
    if preceding_dates:
        predecessor = max(preceding_dates)
        # These files normally contain consecutive 30-day data blocks.
        if predecessor + dt.timedelta(days=31) >= start:
            selected_dates.append(predecessor)

    quality = (
        "NASA calibrated but unreviewed post-1991 {0} MAG; generally not "
        "science quality; retain and apply the CDF point-quality flags "
        "(0=good, 1=bad)"
    ).format(
        "primary/out-board" if sensor == "primary" else "secondary/in-board"
    )
    return [
        make_file_record(
            spacecraft=spacecraft,
            product=product,
            filename=filename,
            url=urllib.parse.urljoin(
                directory_url, urllib.parse.quote(filename)
            ),
            file_date=file_date,
            version=version,
            quality=quality,
            sensor=sensor,
        )
        for file_date in sorted(set(selected_dates))
        for version, filename, directory_url in [best_by_date[file_date]]
    ]


def discover_plasma_fine_year(
    spacecraft: int, year: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover one year of PLS high-resolution solar-wind CDFs."""

    if year < start.year or year > end.year:
        return []
    directory_url = (
        "{0}/voyager{1}/plasma_cdaweb/hires_plasma/{2:04d}/".format(
            SPDF_ROOT, spacecraft, year
        )
    )
    page = fetch_listing(directory_url, allow_missing=True)
    if not page:
        return []
    expression = re.compile(
        PLS_RE_TEMPLATE.format(spacecraft=spacecraft), re.IGNORECASE
    )
    best_by_date: Dict[dt.date, Tuple[int, str]] = {}
    for filename in filenames_from_listing(page):
        match = expression.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date.year != year:
            continue
        version = int(match.group("version"))
        current = best_by_date.get(file_date)
        if current is None or (version, filename) > current:
            best_by_date[file_date] = (version, filename)

    return [
        make_file_record(
            spacecraft=spacecraft,
            product="plasma_fine",
            filename=filename,
            url=urllib.parse.urljoin(directory_url, urllib.parse.quote(filename)),
            file_date=file_date,
            version=version,
            quality=(
                "PLS Level-3 fitted solar-wind parameters; nominal sample "
                "spacing 12 to 192 seconds"
            ),
        )
        for file_date, (version, filename) in sorted(best_by_date.items())
    ]


def discover_plasma_hsh(
    start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Discover the Voyager 2 high-resolution heliosheath PLS CDFs."""

    sheath_start = dt.date(2007, 8, 27)
    sheath_end = dt.date(2018, 11, 5)
    if start > sheath_end or end < sheath_start:
        return []
    sheath_url = (
        "{0}/voyager2/plasma_cdaweb/hires_plasma/heliosheath/".format(
            SPDF_ROOT
        )
    )
    sheath_page = fetch_listing(sheath_url)
    best_by_date: Dict[dt.date, Tuple[int, str]] = {}
    for filename in filenames_from_listing(sheath_page):
        match = PLS_HSH_RE.fullmatch(filename)
        if not match:
            continue
        try:
            file_date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
        except ValueError:
            continue
        if file_date.year < start.year or file_date.year > end.year:
            continue
        version = int(match.group("version"))
        current = best_by_date.get(file_date)
        if current is None or (version, filename) > current:
            best_by_date[file_date] = (version, filename)

    records: List[Dict[str, object]] = []
    for file_date, (version, filename) in sorted(best_by_date.items()):
        relative = (
            Path("voyager2")
            / "pls"
            / "hires"
            / "l3"
            / "heliosheath"
            / "{0:04d}".format(file_date.year)
            / filename
        )
        records.append(
            make_file_record(
                spacecraft=2,
                product="plasma_fine",
                filename=filename,
                url=urllib.parse.urljoin(
                    sheath_url, urllib.parse.quote(filename)
                ),
                file_date=file_date,
                version=version,
                relative_path=relative,
                quality=(
                    "PLS Level-3 heliosheath fits with per-record uncertainties"
                ),
            )
        )
    return records


def discover_spice(
    spacecraft: int, start: dt.date, end: dt.date
) -> List[Dict[str, object]]:
    """Return the verified JPL/NAIF Voyager SPICE kernel stack."""

    if end < max(start, MISSION_START[spacecraft]):
        return []
    base = Path("voyager{0}".format(spacecraft)) / "ephemeris" / "spice"
    records: List[Dict[str, object]] = []
    for entry in SPICE_SHARED_FILES + SPICE_FILES[spacecraft]:
        filename = str(entry["filename"])
        records.append(
            make_file_record(
                spacecraft=spacecraft,
                product="spice_spk",
                filename=filename,
                url=str(entry["url"]),
                file_date=MISSION_START[spacecraft],
                version=1,
                file_format=str(entry.get("format", "spice")),
                relative_path=base / "kernels" / filename,
                quality=str(entry["quality"]),
                expected_size=int(entry["size"]),
                expected_sha256=str(entry["sha256"]),
            )
        )
    return records


def discover_files(
    spacecraft: Sequence[int],
    products: Sequence[str],
    start: dt.date,
    end: dt.date,
    threads: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Discover all requested remote files, parallelizing independent indexes."""

    jobs: List[Tuple[str, int, Optional[int]]] = []
    for number in spacecraft:
        if "coho1hr" in products:
            first_year = max(start.year, 1977)
            if end.year >= first_year:
                jobs.extend(
                    ("coho1hr", number, year)
                    for year in range(first_year, end.year + 1)
                )
        if "position1day" in products:
            jobs.append(("position1day", number, None))
        if "position1hr" in products:
            jobs.append(("position1hr", number, None))
        for magnetic_product in ("mag2s", "mag48s"):
            if magnetic_product not in products:
                continue
            first_year = max(start.year - 1, MISSION_START[number].year)
            last_year = min(end.year, LEGACY_MAG_LAST_YEAR[number])
            if first_year <= last_year:
                jobs.extend(
                    (magnetic_product, number, year)
                    for year in range(first_year, last_year + 1)
                )
        if "mag48s_vim" in products:
            jobs.append(("mag48s_vim", number, None))
        for experimental_product in EXPERIMENTAL_PRODUCTS:
            if experimental_product in products:
                jobs.append((experimental_product, number, None))
        if "plasma_fine" in products:
            first_year = max(start.year, MISSION_START[number].year)
            last_year = min(end.year, PLS_SOLAR_WIND_LAST_YEAR[number])
            if first_year <= last_year:
                jobs.extend(
                    ("plasma_fine", number, year)
                    for year in range(first_year, last_year + 1)
                )
            if number == 2:
                jobs.append(("plasma_fine_hsh", number, None))
        if "spice_spk" in products:
            jobs.append(("spice_spk", number, None))

    if not jobs:
        return [], []

    records: List[Dict[str, object]] = []
    errors: List[Dict[str, object]] = []
    worker_count = min(max(1, threads), len(jobs))

    def run_job(job: Tuple[str, int, Optional[int]]) -> List[Dict[str, object]]:
        product, number, year = job
        if product == "coho1hr":
            assert year is not None
            return discover_coho_year(number, year, start, end)
        if product in {"mag2s", "mag48s"}:
            assert year is not None
            return discover_legacy_mag_year(
                number, product, year, start, end
            )
        if product == "position1day":
            return discover_position(number, start, end)
        if product == "position1hr":
            return discover_position_1hr(number, start, end)
        if product == "mag48s_vim":
            return discover_mag_vim(number, start, end)
        if product in EXPERIMENTAL_PRODUCTS:
            return discover_post1991_mag(
                number, product, start, end
            )
        if product == "plasma_fine":
            assert year is not None
            return discover_plasma_fine_year(
                number, year, start, end
            )
        if product == "plasma_fine_hsh":
            return discover_plasma_hsh(start, end)
        if product == "spice_spk":
            return discover_spice(number, start, end)
        raise ValueError("unsupported discovery product {0}".format(product))

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=worker_count, thread_name_prefix="voyager-discovery"
    ) as executor:
        future_jobs = {executor.submit(run_job, job): job for job in jobs}
        for future in concurrent.futures.as_completed(future_jobs):
            job = future_jobs[future]
            try:
                records.extend(future.result())
            except Exception as exception:
                product, number, year = job
                label = "Voyager {0} {1}".format(number, product)
                if year is not None:
                    label += " {0:04d}".format(year)
                error = {
                    "stage": "discovery",
                    "spacecraft": number,
                    "product": product,
                    "year": year,
                    "label": label,
                    "error": str(exception),
                }
                errors.append(error)
                log("Discovery failed for %s: %s", label, exception, level=logging.ERROR)
    records.sort(
        key=lambda item: (
            int(item["spacecraft"]),
            PRODUCT_ORDER.index(str(item["product"])),
            str(item["date"]),
            str(item["filename"]),
        )
    )
    return records, errors


def remote_content_length(url: str) -> Optional[int]:
    """Read Content-Length using HEAD, with a range fallback if unsupported."""

    last_error: Optional[BaseException] = None
    for attempt in range(HTTP_RETRIES):
        try:
            with urllib.request.urlopen(
                make_request(url, method="HEAD"), timeout=HTTP_TIMEOUT_SECONDS
            ) as response:
                value = response.headers.get("Content-Length")
                return int(value) if value and value.isdigit() else None
        except urllib.error.HTTPError as exception:
            if exception.code in {405, 501}:
                request = make_request(url)
                request.add_header("Range", "bytes=0-0")
                with urllib.request.urlopen(
                    request, timeout=HTTP_TIMEOUT_SECONDS
                ) as response:
                    content_range = response.headers.get("Content-Range", "")
                    match = re.search(r"/(\d+)$", content_range)
                    if match:
                        return int(match.group(1))
                    value = response.headers.get("Content-Length")
                    return int(value) if value and value.isdigit() else None
            last_error = exception
            if not is_transient_http_error(exception):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exception:
            last_error = exception
        if attempt + 1 < HTTP_RETRIES:
            log(
                "HEAD failed (attempt %d/%d): %s",
                attempt + 1,
                HTTP_RETRIES,
                url,
                level=logging.WARNING,
            )
            time.sleep(retry_delay(attempt))
    raise DownloadError("cannot read remote size for {0}: {1}".format(url, last_error))


def validate_cdf(path: Path, expected_size: Optional[int] = None) -> Tuple[bool, str]:
    """Perform inexpensive size and CDF magic checks."""

    try:
        size = path.stat().st_size
    except OSError as exception:
        return False, "cannot stat file: {0}".format(exception)
    if size < MIN_CDF_BYTES:
        return False, "file is only {0} bytes".format(size)
    if expected_size is not None and size != expected_size:
        return False, "size {0} does not match remote {1}".format(size, expected_size)
    try:
        with path.open("rb") as stream:
            header = stream.read(8)
    except OSError as exception:
        return False, "cannot read file: {0}".format(exception)
    if len(header) < 8:
        return False, "CDF header is truncated"
    if header[:4] not in CDF_FIRST_MAGICS or header[4:8] not in CDF_SECOND_MAGICS:
        return False, "CDF magic is invalid ({0})".format(header.hex())
    return True, ""


def validate_file(
    path: Path,
    expected_size: Optional[int] = None,
    file_format: str = "cdf",
    expected_sha256: Optional[str] = None,
) -> Tuple[bool, str]:
    """Validate a downloaded CDF, SPICE kernel, or accompanying text file."""

    normalized_format = file_format.strip().lower()
    if normalized_format == "cdf":
        valid, reason = validate_cdf(path, expected_size)
        if not valid:
            return valid, reason
        if expected_sha256:
            return validate_sha256(path, expected_sha256)
        return True, ""

    try:
        size = path.stat().st_size
    except OSError as exception:
        return False, "cannot stat file: {0}".format(exception)
    if size < MIN_GENERIC_BYTES:
        return False, "file is only {0} bytes".format(size)
    if expected_size is not None and size != expected_size:
        return False, "size {0} does not match remote {1}".format(size, expected_size)
    try:
        with path.open("rb") as stream:
            header = stream.read(512)
    except OSError as exception:
        return False, "cannot read file: {0}".format(exception)

    if normalized_format == "spice":
        if not header.startswith(b"DAF/SPK "):
            return False, "SPICE SPK signature is invalid ({0})".format(
                header[:8].hex()
            )
    elif normalized_format in {"text", "ascii"}:
        if b"\x00" in header:
            return False, "text file contains NUL bytes"
        try:
            header.decode("utf-8")
        except UnicodeDecodeError:
            try:
                header.decode("ascii")
            except UnicodeDecodeError:
                return False, "text file header is not UTF-8/ASCII"
    else:
        return False, "unsupported file format {0}".format(file_format)
    if expected_sha256:
        return validate_sha256(path, expected_sha256)
    return True, ""


def validate_sha256(path: Path, expected_sha256: str) -> Tuple[bool, str]:
    """Compare one file with a trusted SHA-256 digest."""

    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            while True:
                block = stream.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
    except OSError as exception:
        return False, "cannot hash file: {0}".format(exception)
    actual = digest.hexdigest().upper()
    expected = expected_sha256.strip().upper()
    if actual != expected:
        return False, "SHA-256 {0} does not match expected {1}".format(
            actual, expected
        )
    return True, ""


def unlink_quietly(path: Path) -> None:
    """Remove a known temporary file and suppress cleanup-only failures."""

    try:
        path.unlink()
    except FileNotFoundError:
        pass
    except OSError as exception:
        log("Could not remove temporary file %s: %s", path, exception, level=logging.WARNING)


def requires_powershell_transfer(url: str) -> bool:
    """Select the Windows HTTPS stack for JPL or an explicit SPDF fallback."""

    if os.name != "nt":
        return False
    hostname = (urllib.parse.urlsplit(url).hostname or "").lower()
    if hostname in {"naif.jpl.nasa.gov", "ssd.jpl.nasa.gov"}:
        return True
    return force_powershell_http() and hostname == "spdf.gsfc.nasa.gov"


def download_to_stage_powershell(
    url: str,
    stage_file: Path,
    expected_size: Optional[int],
    file_format: str,
    expected_sha256: Optional[str],
) -> int:
    """Download through Windows PowerShell and apply the normal validators."""

    powershell = (
        Path(os.environ.get("WINDIR", r"C:\Windows"))
        / "System32"
        / "WindowsPowerShell"
        / "v1.0"
        / "powershell.exe"
    )
    last_error = ""
    for attempt in range(HTTP_RETRIES):
        partial = stage_file.with_name(
            ".part-ps-{0}-{1}".format(os.getpid(), uuid.uuid4().hex[:12])
        )
        escaped_url = url.replace("'", "''")
        escaped_partial = str(partial).replace("'", "''")
        command = (
            "$ProgressPreference='SilentlyContinue'; "
            "$ErrorActionPreference='Stop'; "
            "Invoke-WebRequest -UseBasicParsing -TimeoutSec 180 "
            "-Uri '{0}' -OutFile '{1}'"
        ).format(escaped_url, escaped_partial)
        try:
            completed = subprocess.run(
                [
                    str(powershell),
                    "-NoProfile",
                    "-NonInteractive",
                    "-Command",
                    command,
                ],
                capture_output=True,
                text=True,
                timeout=240,
                check=False,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if completed.returncode != 0:
                details = (completed.stderr or completed.stdout or "").strip()
                raise DownloadError(
                    details or "PowerShell Invoke-WebRequest failed"
                )
            valid, reason = validate_file(
                partial, expected_size, file_format, expected_sha256
            )
            if not valid:
                raise DownloadError(reason)
            os.replace(str(partial), str(stage_file))
            return stage_file.stat().st_size
        except (OSError, subprocess.SubprocessError, DownloadError) as exception:
            last_error = str(exception)
            unlink_quietly(partial)
            if attempt + 1 < HTTP_RETRIES:
                log(
                    "PowerShell download failed (attempt %d/%d): %s (%s)",
                    attempt + 1,
                    HTTP_RETRIES,
                    url,
                    exception,
                    level=logging.WARNING,
                )
                time.sleep(retry_delay(attempt))
    raise DownloadError(
        "cannot download {0} through PowerShell: {1}".format(url, last_error)
    )


def download_to_stage(
    url: str,
    stage_file: Path,
    expected_size: Optional[int],
    file_format: str = "cdf",
    expected_sha256: Optional[str] = None,
) -> int:
    """Download one file through a unique partial file and atomically stage it."""

    stage_file.parent.mkdir(parents=True, exist_ok=True)
    if requires_powershell_transfer(url):
        return download_to_stage_powershell(
            url,
            stage_file,
            expected_size,
            file_format,
            expected_sha256,
        )
    last_error: Optional[BaseException] = None
    for attempt in range(HTTP_RETRIES):
        partial = stage_file.with_name(
            ".part-{0}-{1}".format(os.getpid(), uuid.uuid4().hex[:12])
        )
        try:
            request = make_request(url)
            with urllib.request.urlopen(
                request, timeout=HTTP_TIMEOUT_SECONDS
            ) as response, partial.open("wb") as output:
                response_size_text = response.headers.get("Content-Length")
                response_size = (
                    int(response_size_text)
                    if response_size_text and response_size_text.isdigit()
                    else None
                )
                if expected_size is None:
                    expected_size = response_size
                shutil.copyfileobj(response, output, length=1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            valid, reason = validate_file(
                partial, expected_size, file_format, expected_sha256
            )
            if not valid:
                raise DownloadError(reason)
            os.replace(str(partial), str(stage_file))
            return stage_file.stat().st_size
        except (
            urllib.error.HTTPError,
            urllib.error.URLError,
            TimeoutError,
            OSError,
            DownloadError,
        ) as exception:
            last_error = exception
            unlink_quietly(partial)
            if isinstance(exception, urllib.error.HTTPError) and not is_transient_http_error(
                exception
            ):
                break
            if attempt + 1 < HTTP_RETRIES:
                log(
                    "Download failed (attempt %d/%d): %s (%s)",
                    attempt + 1,
                    HTTP_RETRIES,
                    url,
                    exception,
                    level=logging.WARNING,
                )
                time.sleep(retry_delay(attempt))
    raise DownloadError("cannot download {0}: {1}".format(url, last_error))


def ensure_archive_parent(parent: Path) -> None:
    """Create an archive directory, including a Windows WebDAV/RaiDrive fallback."""

    try:
        parent.mkdir(parents=True, exist_ok=True)
        return
    except OSError as direct_error:
        if os.name != "nt" or getattr(direct_error, "winerror", None) != 5:
            raise

        powershell = (
            Path(os.environ.get("WINDIR", r"C:\Windows"))
            / "System32"
            / "WindowsPowerShell"
            / "v1.0"
            / "powershell.exe"
        )
        escaped_parent = str(parent).replace("'", "''")
        command = (
            "$ErrorActionPreference='Stop'; "
            "New-Item -ItemType Directory -Path '{0}' -Force | Out-Null"
        ).format(escaped_parent)
        completed = subprocess.run(
            [
                str(powershell),
                "-NoProfile",
                "-NonInteractive",
                "-Command",
                command,
            ],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if completed.returncode != 0 or not parent.is_dir():
            details = (completed.stderr or completed.stdout or "").strip()
            raise PermissionError(
                "cannot create archive directory {0}: {1}; PowerShell fallback: {2}".format(
                    parent, direct_error, details or "directory was not created"
                )
            ) from direct_error
        log("Created archive directory through PowerShell fallback: %s", parent)


def publish_from_stage(
    stage_file: Path,
    target: Path,
    expected_size: Optional[int],
    file_format: str = "cdf",
    expected_sha256: Optional[str] = None,
) -> int:
    """Copy a validated staged file and atomically publish it in the archive."""

    ensure_archive_parent(target.parent)
    temporary = target.with_name(
        ".tmp-{0}-{1}".format(os.getpid(), uuid.uuid4().hex[:12])
    )
    try:
        with stage_file.open("rb") as source, temporary.open("wb") as destination:
            shutil.copyfileobj(source, destination, length=1024 * 1024)
            destination.flush()
            os.fsync(destination.fileno())
        valid, reason = validate_file(
            temporary, expected_size, file_format, expected_sha256
        )
        if not valid:
            raise DownloadError("publish validation failed: {0}".format(reason))
        os.replace(str(temporary), str(target))
        return target.stat().st_size
    finally:
        unlink_quietly(temporary)


def download_one(
    source_record: Mapping[str, object],
    output_root: Path,
    stage_root: Path,
    check_size: bool,
    force: bool,
) -> Dict[str, object]:
    """Download, validate, and atomically publish one discovered file."""

    record: Dict[str, object] = dict(source_record)
    relative_path = Path(str(record["relative_path"]))
    target = output_root / relative_path
    stage_file = (
        stage_root
        / "voyager{0}".format(record["spacecraft"])
        / str(record["product"])
        / str(record["filename"])
    )
    record["target"] = str(target)
    record["stage"] = str(stage_file)
    file_format = str(record.get("format", "cdf"))
    expected_sha256 = (
        str(record["expected_sha256"])
        if record.get("expected_sha256")
        else None
    )

    try:
        size_hint = record.get("expected_size")
        expected_size = (
            int(size_hint) if isinstance(size_hint, int) and size_hint > 0 else None
        )
        if check_size and expected_size is None:
            try:
                expected_size = remote_content_length(str(record["url"]))
            except DownloadError as exception:
                record["size_check_warning"] = str(exception)
                log(
                    "Remote size unavailable; continuing with response-size "
                    "and format validation: %s (%s)",
                    record["url"],
                    exception,
                    level=logging.WARNING,
                )
        record["expected_size"] = expected_size
        if not force and target.is_file():
            valid, reason = validate_file(
                target, expected_size, file_format, expected_sha256
            )
            if valid:
                size = target.stat().st_size
                record.update(
                    {
                        "status": "skipped",
                        "size": size,
                        "detail": "existing complete file",
                    }
                )
                if expected_sha256:
                    record["sha256"] = expected_sha256.upper()
                log("Skip complete file: %s", target)
                return record
            log(
                "Existing file will be replaced because validation failed: %s (%s)",
                target,
                reason,
                level=logging.WARNING,
            )

        used_stage = False
        if not force and stage_file.is_file():
            staged_valid, _ = validate_file(
                stage_file, expected_size, file_format, expected_sha256
            )
            used_stage = staged_valid
        if not used_stage:
            download_to_stage(
                str(record["url"]),
                stage_file,
                expected_size,
                file_format,
                expected_sha256,
            )

        size = publish_from_stage(
            stage_file,
            target,
            expected_size,
            file_format,
            expected_sha256,
        )
        record.update(
            {
                "status": "downloaded",
                "size": size,
                "detail": (
                    "published from an existing complete staging file"
                    if used_stage
                    else "downloaded, validated, and published"
                ),
            }
        )
        if expected_sha256:
            record["sha256"] = expected_sha256.upper()
        log("Published: %s", target)
        return record
    except Exception as exception:
        record.update(
            {
                "status": "failed",
                "size": target.stat().st_size if target.is_file() else None,
                "error": str(exception),
            }
        )
        log("Failed: %s (%s)", record["url"], exception, level=logging.ERROR)
        return record


def atomic_write_json(path: Path, value: object) -> None:
    """Atomically write UTF-8 JSON next to its final destination."""

    ensure_archive_parent(path.parent)
    temporary = path.with_name(
        ".json-{0}-{1}.tmp".format(os.getpid(), uuid.uuid4().hex[:12])
    )
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2, sort_keys=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(str(temporary), str(path))
    finally:
        unlink_quietly(temporary)


def build_manifest(
    args: argparse.Namespace,
    start: dt.date,
    end: dt.date,
    spacecraft: Sequence[int],
    products: Sequence[str],
    started_at: str,
    files: Sequence[Mapping[str, object]],
    discovery_errors: Sequence[Mapping[str, object]],
) -> Dict[str, object]:
    """Create the complete download manifest."""

    counts = {"discovered": len(files), "downloaded": 0, "skipped": 0, "failed": 0}
    bytes_published = 0
    errors: List[Dict[str, object]] = [dict(item) for item in discovery_errors]
    for record in files:
        status = str(record.get("status", ""))
        if status in counts:
            counts[status] += 1
        size = record.get("size")
        if status in {"downloaded", "skipped"} and isinstance(size, int):
            bytes_published += size
        if status == "failed":
            errors.append(
                {
                    "url": record.get("url"),
                    "target": record.get("target"),
                    "error": record.get("error", "unknown failure"),
                }
            )
    counts["discovery_failed"] = len(discovery_errors)
    status = "completed" if not errors else "completed_with_errors"
    return {
        "schema_version": 2,
        "status": status,
        "started_at": started_at,
        "finished_at": utc_now(),
        "source": (
            "NASA Space Physics Data Facility (SPDF) and JPL/NAIF"
            if "spice_spk" in products
            else "NASA Space Physics Data Facility (SPDF)"
        ),
        "options": {
            "date": "{0}/{1}".format(start.isoformat(), end.isoformat()),
            "spacecraft": list(spacecraft),
            "products": list(products),
            "out": str(Path(args.out)),
            "stage": str(Path(args.stage)),
            "manifest_name": args.manifest_name,
            "threads": args.threads,
            "check_size": args.check_size,
            "force": args.force,
            "selection_policy": "highest_version_per_granule_date",
        },
        "summary": {
            **counts,
            "successful": counts["downloaded"] + counts["skipped"],
            "bytes_selected": bytes_published,
        },
        "files": list(files),
        "errors": errors,
    }


def make_parser() -> argparse.ArgumentParser:
    """Build the public command-line interface."""

    today = dt.datetime.now(dt.timezone.utc).date()
    parser = argparse.ArgumentParser(
        description=(
            "Download Voyager 1/2 magnetic field, plasma, and ephemeris "
            "products from NASA/SPDF and JPL/NAIF."
        )
    )
    parser.add_argument(
        "--date",
        default="1977-08-20/{0}".format(today.isoformat()),
        metavar="START/END",
        help="inclusive date interval (default: mission start through today)",
    )
    parser.add_argument(
        "--spacecraft",
        action="append",
        default=None,
        metavar="LIST",
        help="1, 2, or a comma-separated list (default: 1,2)",
    )
    parser.add_argument(
        "--product",
        action="append",
        default=None,
        metavar="LIST",
        help=(
            "comma-separated products, or aliases both/highres/"
            "mag2s_unreviewed/experimental_mag/all "
            "(default: coho1hr,position1day)"
        ),
    )
    parser.add_argument(
        "--out",
        default=r"Z:\SPART-WORK\Data\Voyager",
        help=r"archive root (default: Z:\SPART-WORK\Data\Voyager)",
    )
    parser.add_argument(
        "--stage",
        default=str(Path(tempfile.gettempdir()) / "Voyager_staging"),
        help="local staging root (default: system temporary directory)",
    )
    parser.add_argument(
        "--manifest-name",
        default=MANIFEST_NAME,
        type=parse_manifest_name,
        metavar="FILE",
        help="manifest filename written under --out",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=MAX_THREADS,
        help="simultaneous discovery/download workers, maximum 5 (default: 5)",
    )
    parser.add_argument(
        "--check-size",
        nargs="?",
        const=True,
        default=True,
        type=parse_bool,
        metavar="BOOL",
        help="compare files with HTTP Content-Length (default: true)",
    )
    parser.add_argument(
        "--force",
        nargs="?",
        const=True,
        default=False,
        type=parse_bool,
        metavar="BOOL",
        help="redownload complete files (default: false)",
    )
    parser.add_argument(
        "--list-only",
        nargs="?",
        const=True,
        default=False,
        type=parse_bool,
        metavar="BOOL",
        help="discover and print files without writing anything",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit machine-readable JSON to stdout",
    )
    parser.add_argument(
        "--show-products",
        action="store_true",
        help="show supported products and exit",
    )
    return parser


def show_products(as_json: bool) -> None:
    """Print the supported product catalogue."""

    catalogue = [
        {"name": name, **dict(PRODUCTS[name])}
        for name in PRODUCT_ORDER
    ]
    if as_json:
        json.dump(catalogue, sys.stdout, ensure_ascii=False, indent=2)
        sys.stdout.write("\n")
        return
    for item in catalogue:
        print(
            "{name}: {description}\n"
            "  cadence: {cadence}\n"
            "  remote:  {remote}\n"
            "  local:   {local}".format(**item)
        )


def list_record_for_output(
    record: Mapping[str, object], output_root: Path
) -> Dict[str, object]:
    """Return the list-only schema, including the prospective target."""

    output = dict(record)
    output["target"] = str(output_root / Path(str(record["relative_path"])))
    return output


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    parser = make_parser()
    args = parser.parse_args(argv)

    if args.show_products:
        show_products(args.json)
        return 0
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.threads > MAX_THREADS:
        log(
            "Requested %d workers; limiting NASA/SPDF concurrency to %d",
            args.threads,
            MAX_THREADS,
            level=logging.WARNING,
        )
        args.threads = MAX_THREADS

    try:
        start, end = parse_date_range(args.date)
        spacecraft = parse_spacecraft(args.spacecraft or ["1,2"])
        products = parse_products(args.product or ["coho1hr,position1day"])
    except argparse.ArgumentTypeError as exception:
        parser.error(str(exception))

    started_at = utc_now()
    try:
        discovered, discovery_errors = discover_files(
            spacecraft=spacecraft,
            products=products,
            start=start,
            end=end,
            threads=args.threads,
        )
    except DiscoveryError as exception:
        log("Discovery failed: %s", exception, level=logging.ERROR)
        if args.json:
            json.dump(
                {"status": "discovery_failed", "error": str(exception)},
                sys.stdout,
                ensure_ascii=False,
            )
            sys.stdout.write("\n")
        return 2

    output_root = Path(args.out)
    if args.list_only:
        listed = [
            list_record_for_output(record, output_root) for record in discovered
        ]
        if args.json:
            json.dump(listed, sys.stdout, ensure_ascii=False, indent=2)
            sys.stdout.write("\n")
        else:
            for record in listed:
                print(
                    "{spacecraft}\t{product}\t{date}\t{url}\t{relative_path}".format(
                        **record
                    )
                )
        return 2 if discovery_errors else 0

    stage_root = Path(args.stage)
    ensure_archive_parent(output_root)
    stage_root.mkdir(parents=True, exist_ok=True)
    log(
        "Discovered %d file(s); downloading with %d worker(s)",
        len(discovered),
        min(args.threads, max(1, len(discovered))),
    )
    completed: List[Dict[str, object]] = []
    if discovered:
        worker_count = min(args.threads, len(discovered))
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=worker_count, thread_name_prefix="voyager-download"
        ) as executor:
            futures = [
                executor.submit(
                    download_one,
                    record,
                    output_root,
                    stage_root,
                    args.check_size,
                    args.force,
                )
                for record in discovered
            ]
            for future in concurrent.futures.as_completed(futures):
                completed.append(future.result())
    completed.sort(
        key=lambda item: (
            int(item["spacecraft"]),
            PRODUCT_ORDER.index(str(item["product"])),
            str(item["date"]),
            str(item["filename"]),
        )
    )
    manifest = build_manifest(
        args=args,
        start=start,
        end=end,
        spacecraft=spacecraft,
        products=products,
        started_at=started_at,
        files=completed,
        discovery_errors=discovery_errors,
    )
    manifest_path = output_root / args.manifest_name
    atomic_write_json(manifest_path, manifest)
    log("Manifest: %s", manifest_path)
    if args.json:
        json.dump(manifest, sys.stdout, ensure_ascii=False, indent=2)
        sys.stdout.write("\n")
    return 0 if manifest["status"] == "completed" else 3


if __name__ == "__main__":
    raise SystemExit(main())
