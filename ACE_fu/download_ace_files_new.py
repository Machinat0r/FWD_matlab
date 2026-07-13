#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACE CDAWeb crawler and downloader.

This script follows the same split used by the MMS downloader:
1. crawl CDAWeb directory listings and choose the newest CDF version per day;
2. skip complete local files;
3. download missing/incomplete files with optional multi-thread Range requests.
"""

from __future__ import annotations

import argparse
import datetime as dt
import html
import http.client
import json
import os
import random
import re
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.error import HTTPError as UrlHTTPError
from urllib.error import URLError
from urllib.parse import urljoin
from urllib.request import Request, build_opener, ProxyHandler


# Prefer dependencies vendored next to this script.  This keeps the downloader
# self-contained when MATLAB launches a Python interpreter with no site-wide
# packages installed.
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR_DIR = SCRIPT_DIR / "_vendor"
if VENDOR_DIR.is_dir():
    sys.path.insert(0, str(VENDOR_DIR))


try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    HAVE_REQUESTS = True
except Exception:
    requests = None
    HTTPAdapter = None
    Retry = None
    HAVE_REQUESTS = False

try:
    from cdflib import CDF as CdflibCDF

    HAVE_CDFLIB = True
except Exception:
    CdflibCDF = None
    HAVE_CDFLIB = False

try:
    import warnings
    from urllib3.exceptions import NotOpenSSLWarning

    warnings.filterwarnings("ignore", category=NotOpenSSLWarning)
except Exception:
    pass


CDAWEB_ACE_ROOT = "https://cdaweb.gsfc.nasa.gov/pub/data/ace/"


@dataclass(frozen=True)
class Product:
    key: str
    path: str
    note: str = ""

    @property
    def base_url(self) -> str:
        return urljoin(CDAWEB_ACE_ROOT, self.path.strip("/") + "/")


PRODUCTS: Dict[str, Product] = {
    "mfi_h0": Product("mfi_h0", "mag/level_2_cdaweb/mfi_h0", "MAG/MFI H0"),
    "mfi_h1": Product("mfi_h1", "mag/level_2_cdaweb/mfi_h1", "MAG/MFI H1"),
    "mfi_h2": Product("mfi_h2", "mag/level_2_cdaweb/mfi_h2", "MAG/MFI H2"),
    "mfi_h3": Product("mfi_h3", "mag/level_2_cdaweb/mfi_h3", "MAG/MFI H3"),
    "mfi_k0": Product("mfi_k0", "mag/level_2_cdaweb/mfi_k0", "MAG/MFI key parameter"),
    "mfi_k1": Product("mfi_k1", "mag/level_2_cdaweb/mfi_k1", "MAG/MFI key parameter"),
    "mfi_k2": Product("mfi_k2", "mag/level_2_cdaweb/mfi_k2", "MAG/MFI key parameter"),
    "swe_h0": Product("swe_h0", "swepam/level_2_cdaweb/swe_h0", "SWEPAM H0"),
    "swe_h2": Product("swe_h2", "swepam/level_2_cdaweb/swe_h2", "SWEPAM H2"),
    "swe_k0": Product("swe_k0", "swepam/level_2_cdaweb/swe_k0", "SWEPAM key parameter"),
    "swe_k1": Product("swe_k1", "swepam/level_2_cdaweb/swe_k1", "SWEPAM key parameter"),
    "epm_h1": Product("epm_h1", "epam/level_2_cdaweb/epm_h1", "EPAM H1"),
    "epm_h2": Product("epm_h2", "epam/level_2_cdaweb/epm_h2", "EPAM H2"),
    "epm_h3": Product("epm_h3", "epam/level_2_cdaweb/epm_h3", "EPAM H3"),
    "epm_k0": Product("epm_k0", "epam/level_2_cdaweb/epm_k0", "EPAM key parameter"),
    "epm_k1": Product("epm_k1", "epam/level_2_cdaweb/epm_k1", "EPAM key parameter"),
    "sis_h1": Product("sis_h1", "sis/level_2_cdaweb/sis_h1", "SIS H1"),
    "sis_h2": Product("sis_h2", "sis/level_2_cdaweb/sis_h2", "SIS H2"),
    "sis_k0": Product("sis_k0", "sis/level_2_cdaweb/sis_k0", "SIS key parameter"),
}

ALIASES = {
    "mag": "mfi_h0",
    "mfi": "mfi_h0",
    "mfi_h": "mfi_h0",
    "fgm": "mfi_h0",
    "swe": "swe_h0",
    "swepam": "swe_h0",
    "plasma": "swe_h0",
    "epam": "epm_h1",
    "epm": "epm_h1",
    "sis": "sis_h1",
}

_progress_lock = threading.Lock()
_cdflib_warning_lock = threading.Lock()
_cdflib_warning_emitted = False


class ProgressReporter:
    """Thread-safe, log-friendly progress reporting with strong throttling."""

    def __init__(self, total_bytes: int, min_interval: float = 2.0, min_percent: float = 5.0):
        self.total_bytes = max(1, int(total_bytes))
        self.min_interval = float(min_interval)
        self.min_percent = float(min_percent)
        self.current_bytes = 0
        self.started = time.time()
        self.last_report_time = self.started
        self.last_report_percent = 0.0

    def update(self, nbytes: int) -> None:
        with _progress_lock:
            self.current_bytes += int(nbytes)
            now = time.time()
            percent = min(100.0, self.current_bytes / self.total_bytes * 100.0)
            complete = self.current_bytes >= self.total_bytes
            enough_time = now - self.last_report_time >= self.min_interval
            enough_progress = percent - self.last_report_percent >= self.min_percent
            if not complete and not (enough_time and enough_progress):
                return

            elapsed = max(now - self.started, 0.001)
            speed = self.current_bytes / elapsed / 1024 / 1024
            print(
                f"Downloaded: {self.current_bytes / 1024 / 1024:.2f} MB "
                f"[{percent:6.2f}%] {speed:.2f} MB/s",
                flush=True,
            )
            self.last_report_time = now
            self.last_report_percent = percent


class SimpleHTTPError(RuntimeError):
    pass


class UrllibResponse:
    def __init__(self, response):
        self._response = response
        self.status_code = response.getcode()
        self.url = response.geturl()
        self.headers = dict(response.headers.items())

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def close(self) -> None:
        self._response.close()

    @property
    def text(self) -> str:
        data = self._response.read()
        charset = self._response.headers.get_content_charset() or "utf-8"
        return data.decode(charset, errors="replace")

    def iter_content(self, chunk_size: int = 1024 * 256):
        while True:
            chunk = self._response.read(chunk_size)
            if not chunk:
                break
            yield chunk

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise SimpleHTTPError(f"HTTP {self.status_code} for {self.url}")


class UrllibSession:
    def __init__(self, trust_env: bool = True):
        proxy_handler = ProxyHandler() if trust_env else ProxyHandler({})
        self._opener = build_opener(proxy_handler)
        self.headers = {}

    def get(self, url, headers=None, stream=False, allow_redirects=True, timeout=None):
        return self._open("GET", url, headers, timeout)

    def head(self, url, headers=None, allow_redirects=True, timeout=None):
        return self._open("HEAD", url, headers, timeout)

    def _open(self, method: str, url: str, headers, timeout):
        merged_headers = dict(self.headers)
        if headers:
            merged_headers.update(headers)
        req = Request(url, headers=merged_headers, method=method)
        try:
            response = self._opener.open(req, timeout=_timeout_seconds(timeout))
        except UrlHTTPError as exc:
            response = exc
        except URLError as exc:
            raise OSError(str(exc)) from exc
        return UrllibResponse(response)


def _timeout_seconds(timeout) -> float:
    if isinstance(timeout, tuple):
        return float(timeout[-1])
    if timeout is None:
        return 45.0
    return float(timeout)


if not HAVE_REQUESTS:
    class _RequestsExceptions:
        ChunkedEncodingError = OSError
        ConnectionError = OSError
        HTTPError = SimpleHTTPError
        ReadTimeout = TimeoutError
        SSLError = OSError

    class _RequestsShim:
        exceptions = _RequestsExceptions

    requests = _RequestsShim()


def _build_session(trust_env: bool = True) -> requests.Session:
    if not HAVE_REQUESTS:
        s = UrllibSession(trust_env)
        s.headers.update(
            {
                "User-Agent": (
                    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0 Safari/537.36"
                ),
                "Accept": "*/*",
            }
        )
        return s

    s = requests.Session()
    s.trust_env = trust_env

    retry = Retry(
        total=6,
        connect=6,
        read=6,
        status=6,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["HEAD", "GET"]),
        raise_on_status=False,
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=64, pool_maxsize=64)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/120.0 Safari/537.36"
            ),
            "Accept": "*/*",
        }
    )
    return s


_SESSION_ENV = _build_session(True)
_SESSION_NOENV = _build_session(False)


def _sleep_backoff(attempt: int, cap: float = 8.0) -> None:
    base = min(cap, 0.5 * (2 ** (attempt - 1)))
    time.sleep(base + random.uniform(0, 0.25))


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_date(value: str) -> dt.date:
    value = str(value).strip()
    if "T" in value:
        value = value.split("T", 1)[0]
    value = value.replace(".", "-")
    if re.fullmatch(r"\d{8}", value):
        return dt.datetime.strptime(value, "%Y%m%d").date()
    return dt.datetime.strptime(value, "%Y-%m-%d").date()


def parse_date_range(value: str) -> Tuple[dt.date, dt.date]:
    parts = re.split(r"/|,", str(value).strip())
    parts = [p for p in parts if p]
    if not parts:
        raise ValueError("date is empty")
    start = parse_date(parts[0])
    end = parse_date(parts[-1]) if len(parts) > 1 else start
    if end < start:
        start, end = end, start
    return start, end


def iter_years(start: dt.date, end: dt.date) -> Iterable[int]:
    for year in range(start.year, end.year + 1):
        yield year


def normalize_product(value: str) -> str:
    key = str(value).strip().lower().replace("-", "_")
    key = ALIASES.get(key, key)
    if key not in PRODUCTS:
        raise KeyError(
            f"Unknown ACE product '{value}'. Use --show-products to list supported products."
        )
    return key


def parse_product_list(value: str) -> List[str]:
    products: List[str] = []
    for item in re.split(r"[,;\s]+", str(value).strip()):
        if not item:
            continue
        products.append(normalize_product(item))
    if not products:
        products.append("mfi_h0")
    return list(dict.fromkeys(products))


def fetch_html(url: str, timeout: Tuple[float, float] = (8, 45)) -> str:
    last_err: Optional[BaseException] = None
    for sess in (_SESSION_ENV, _SESSION_NOENV):
        for attempt in range(1, 7):
            try:
                with sess.get(url, timeout=timeout, allow_redirects=True) as r:
                    if r.status_code == 404:
                        return ""
                    r.raise_for_status()
                    return r.text
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                requests.exceptions.SSLError,
                requests.exceptions.HTTPError,
                http.client.IncompleteRead,
                http.client.RemoteDisconnected,
                OSError,
                TimeoutError,
            ) as exc:
                last_err = exc
                if attempt == 6:
                    break
                _sleep_backoff(attempt)
    raise RuntimeError(f"Failed to read directory: {url}. Last error: {last_err}")


def extract_cdf_links(index_html: str) -> List[str]:
    links = re.findall(r'href=["\']([^"\']+\.cdf)["\']', index_html, flags=re.IGNORECASE)
    names = [os.path.basename(html.unescape(link)) for link in links]
    return sorted(dict.fromkeys(names))


def filename_date_version(filename: str) -> Optional[Tuple[dt.date, int]]:
    m = re.search(r"_(\d{8})_v(\d+)\.cdf$", filename, flags=re.IGNORECASE)
    if not m:
        return None
    return dt.datetime.strptime(m.group(1), "%Y%m%d").date(), int(m.group(2))


def crawl_product(product_key: str, start: dt.date, end: dt.date) -> List[dict]:
    product = PRODUCTS[product_key]
    by_date: Dict[dt.date, dict] = {}

    for year in iter_years(start, end):
        year_url = urljoin(product.base_url, f"{year}/")
        index_html = fetch_html(year_url)
        if not index_html:
            print(f"Warning: remote year directory not found: {year_url}", file=sys.stderr)
            continue

        for filename in extract_cdf_links(index_html):
            parsed = filename_date_version(filename)
            if parsed is None:
                continue
            day, version = parsed
            if day < start or day > end:
                continue
            old = by_date.get(day)
            if old is None or version > old["version"]:
                by_date[day] = {
                    "product": product.key,
                    "date": day.isoformat(),
                    "year": year,
                    "version": version,
                    "filename": filename,
                    "url": urljoin(year_url, filename),
                }

    return [by_date[k] for k in sorted(by_date)]


def crawl_products(product_keys: Sequence[str], start: dt.date, end: dt.date) -> List[dict]:
    files: List[dict] = []
    for product_key in product_keys:
        files.extend(crawl_product(product_key, start, end))
    return sorted(files, key=lambda item: (item["product"], item["date"], item["filename"]))


def _head_with_retries(
    session: requests.Session, url: str, timeout: Tuple[float, float], attempts: int = 6
) -> Tuple[str, dict]:
    last_err: Optional[BaseException] = None
    for i in range(1, attempts + 1):
        try:
            with session.head(url, allow_redirects=True, timeout=timeout) as r:
                r.raise_for_status()
                return r.url, dict(r.headers)
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            requests.exceptions.HTTPError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            OSError,
            TimeoutError,
        ) as exc:
            last_err = exc
            if i == attempts:
                break
            _sleep_backoff(i)
    raise RuntimeError(f"HEAD failed for {url}. Last error: {last_err}")


def _probe_size_by_range(
    session: requests.Session, url: str, timeout: Tuple[float, float], attempts: int = 6
) -> Tuple[str, int]:
    last_err: Optional[BaseException] = None
    headers = {"Range": "bytes=0-0", "Accept-Encoding": "identity"}
    for i in range(1, attempts + 1):
        try:
            with session.get(
                url, headers=headers, stream=True, allow_redirects=True, timeout=timeout
            ) as r:
                if r.status_code != 206:
                    raise RuntimeError(f"Range probe expected 206, got {r.status_code}")
                cr = r.headers.get("Content-Range", "")
                match = re.fullmatch(r"bytes\s+0-0/(\d+)", cr.strip(), flags=re.IGNORECASE)
                if match is None:
                    raise RuntimeError(f"Missing/invalid Content-Range: {cr}")
                total = int(match.group(1))
                probe_bytes = 0
                for chunk in r.iter_content(chunk_size=16):
                    if not chunk:
                        continue
                    probe_bytes += len(chunk)
                    if probe_bytes > 1:
                        raise RuntimeError(
                            f"Range probe returned {probe_bytes} bytes; expected exactly 1"
                        )
                if probe_bytes != 1:
                    raise RuntimeError(
                        f"Range probe returned {probe_bytes} bytes; expected exactly 1"
                    )
                return r.url, total
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            OSError,
            TimeoutError,
            RuntimeError,
            ValueError,
        ) as exc:
            last_err = exc
            if i == attempts:
                break
            _sleep_backoff(i)
    raise RuntimeError(f"Range probe failed for {url}. Last error: {last_err}")


def resolve_url_and_size(url: str, timeout: Tuple[float, float] = (8, 45)) -> Tuple[str, int, requests.Session, bool]:
    for sess in (_SESSION_ENV, _SESSION_NOENV):
        try:
            final_url, headers = _head_with_retries(sess, url, timeout=timeout, attempts=6)
            cl = headers.get("Content-Length") or headers.get("content-length")
            if cl is not None:
                content_size = int(cl)
                if content_size > 0:
                    accept_ranges = (
                        headers.get("Accept-Ranges")
                        or headers.get("accept-ranges")
                        or ""
                    ).lower()
                    return final_url, content_size, sess, accept_ranges == "bytes"
        except Exception:
            continue

    for sess in (_SESSION_ENV, _SESSION_NOENV):
        try:
            final_url, total = _probe_size_by_range(sess, url, timeout=timeout, attempts=6)
            if total > 0:
                return final_url, total, sess, True
        except Exception:
            continue

    raise RuntimeError(f"Failed to resolve Content-Length for {url}")


def _info_value(info, name: str, default=None):
    if isinstance(info, dict):
        if name in info:
            return info[name]
        lowered = {str(key).lower(): value for key, value in info.items()}
        return lowered.get(name.lower(), default)
    return getattr(info, name, default)


def _cdf_variable_names(cdf) -> List[str]:
    info = cdf.cdf_info()
    names: List[str] = []
    for field in ("rVariables", "zVariables"):
        values = _info_value(info, field, []) or []
        for value in values:
            if isinstance(value, bytes):
                value = value.decode("utf-8", errors="replace")
            names.append(str(value))
    return names


def _value_is_nonempty(value) -> bool:
    if value is None:
        return False
    size = getattr(value, "size", None)
    if size is not None:
        try:
            return int(size) > 0
        except (TypeError, ValueError):
            pass
    try:
        return len(value) > 0
    except TypeError:
        return True


def _cdf_variable_has_records(cdf, var_name: str) -> bool:
    try:
        inquiry = cdf.varinq(var_name)
        last_record = _info_value(inquiry, "Last_Rec", None)
        if last_record is not None and int(last_record) < 0:
            return False
    except Exception:
        # Some older cdflib releases do not expose Last_Rec consistently;
        # reading one record below is the authoritative fallback.
        pass

    try:
        value = cdf.varget(var_name, startrec=0, endrec=0)
    except TypeError:
        value = cdf.varget(var_name)
    return _value_is_nonempty(value)


def _warn_cdflib_unavailable() -> None:
    global _cdflib_warning_emitted
    with _cdflib_warning_lock:
        if _cdflib_warning_emitted:
            return
        print(
            "Warning: cdflib is unavailable; only file size and a fast non-zero "
            "CDF-header check will be performed. MFI Epoch/Magnitude records are "
            "not being content-validated.",
            file=sys.stderr,
        )
        _cdflib_warning_emitted = True


def validate_local_cdf(
    path: Path,
    product_key: str,
    expected_size: Optional[int] = None,
    use_cdflib: bool = True,
) -> Tuple[bool, str]:
    """Validate size, a non-zero CDF header, and optionally CDF contents."""
    try:
        if not path.is_file():
            return False, "file does not exist"
        actual_size = path.stat().st_size
        if expected_size is not None and actual_size != expected_size:
            return False, f"size mismatch: {actual_size} != {expected_size}"
        if actual_size < 16:
            return False, f"file is too short to be a CDF: {actual_size} bytes"
        with path.open("rb") as handle:
            header = handle.read(16)
        if len(header) < 16:
            return False, "could not read the first 16 bytes"
        if not any(header):
            return False, "the first 16 bytes are all zero"
    except OSError as exc:
        return False, f"basic file validation failed: {exc}"

    if not use_cdflib:
        return True, "basic CDF validation passed (cdflib validation disabled)"
    if not HAVE_CDFLIB:
        _warn_cdflib_unavailable()
        return True, "basic CDF validation passed (cdflib unavailable)"

    cdf = None
    try:
        cdf = CdflibCDF(str(path))
        names = _cdf_variable_names(cdf)
        if not names:
            return False, "cdflib found no CDF variables"

        if str(product_key).lower().startswith("mfi_"):
            by_lower = {name.lower(): name for name in names}
            missing = [name for name in ("Epoch", "Magnitude") if name.lower() not in by_lower]
            if missing:
                return False, f"MFI CDF is missing variable(s): {', '.join(missing)}"
            for required in ("Epoch", "Magnitude"):
                actual_name = by_lower[required.lower()]
                if not _cdf_variable_has_records(cdf, actual_name):
                    return False, f"MFI variable {actual_name} has no records"
        return True, "cdflib content validation passed"
    except Exception as exc:
        return False, f"cdflib could not validate the CDF: {exc}"
    finally:
        if cdf is not None:
            try:
                cdf.close()
            except Exception:
                pass


def _new_part_path(final_path: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=f".{final_path.name}.", suffix=".part", dir=str(final_path.parent)
    )
    os.close(descriptor)
    return Path(name)


def download_chunk(args) -> Tuple[int, int, int]:
    (
        session,
        url,
        start_byte,
        end_byte,
        content_size,
        local_file_path,
        progress_callback,
    ) = args
    timeout = (8, 45)
    current_pos = start_byte
    total_written = 0
    max_retries = 10
    attempt = 0
    base_headers = {"Accept-Encoding": "identity"}

    while current_pos <= end_byte:
        headers = dict(base_headers)
        headers["Range"] = f"bytes={current_pos}-{end_byte}"
        try:
            with session.get(
                url, headers=headers, stream=True, timeout=timeout, allow_redirects=True
            ) as response:
                if response.status_code != 206:
                    raise RuntimeError(
                        f"Expected 206 for range request, got {response.status_code}"
                    )
                content_range = response.headers.get("Content-Range", "")
                match = re.fullmatch(
                    r"bytes\s+(\d+)-(\d+)/(\d+)", content_range.strip(), flags=re.IGNORECASE
                )
                if match is None:
                    raise RuntimeError(f"Missing/invalid Content-Range: {content_range}")
                cr_start, cr_end, cr_total = (int(value) for value in match.groups())
                if cr_start != current_pos:
                    raise RuntimeError(f"Content-Range start mismatch: {cr_start} != {current_pos}")
                if cr_end != end_byte:
                    raise RuntimeError(f"Content-Range end mismatch: {cr_end} != {end_byte}")
                if cr_total != content_size:
                    raise RuntimeError(f"Content-Range total mismatch: {cr_total} != {content_size}")
                if cr_end < cr_start:
                    raise RuntimeError(f"Invalid Content-Range endpoints: {content_range}")

                response_written = 0
                with open(local_file_path, "r+b") as file:
                    file.seek(current_pos)
                    for chunk in response.iter_content(chunk_size=1024 * 256):
                        if not chunk:
                            continue
                        next_pos = current_pos + len(chunk)
                        if next_pos > cr_end + 1 or next_pos > end_byte + 1:
                            raise RuntimeError(
                                f"Range body exceeds declared endpoint {cr_end}: "
                                f"next byte would be {next_pos - 1}"
                            )
                        file.write(chunk)
                        current_pos = next_pos
                        response_written += len(chunk)
                        total_written += len(chunk)
                        progress_callback(len(chunk))
                    file.flush()

                declared_bytes = cr_end - cr_start + 1
                if response_written != declared_bytes:
                    raise RuntimeError(
                        f"Range body length mismatch for {cr_start}-{cr_end}: "
                        f"wrote {response_written}, expected {declared_bytes}"
                    )
                if current_pos != cr_end + 1:
                    raise RuntimeError(
                        f"Range write stopped at byte {current_pos - 1}, expected {cr_end}"
                    )
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            OSError,
            TimeoutError,
            RuntimeError,
        ) as exc:
            attempt += 1
            if attempt >= max_retries:
                raise RuntimeError(
                    f"Failed downloading range {start_byte}-{end_byte} after "
                    f"{max_retries} retries. Last error: {exc}"
                ) from exc
            _sleep_backoff(attempt, cap=6.0)

    expected_written = end_byte - start_byte + 1
    if total_written != expected_written:
        raise RuntimeError(
            f"Range {start_byte}-{end_byte} wrote {total_written} bytes; "
            f"expected {expected_written}"
        )
    return start_byte, end_byte, total_written


def stream_download(
    session: requests.Session, final_url: str, local_file_path: Path, expected_size: Optional[int]
) -> None:
    headers = {"Accept-Encoding": "identity"}
    max_retries = 6
    last_error: Optional[BaseException] = None

    for attempt in range(1, max_retries + 1):
        current_size = 0
        reporter = ProgressReporter(expected_size) if expected_size is not None else None
        try:
            with session.get(
                final_url, headers=headers, stream=True, allow_redirects=True, timeout=(8, 45)
            ) as response:
                response.raise_for_status()
                with open(local_file_path, "wb") as file:
                    for chunk in response.iter_content(chunk_size=1024 * 256):
                        if not chunk:
                            continue
                        next_size = current_size + len(chunk)
                        if expected_size is not None and next_size > expected_size:
                            raise RuntimeError(
                                f"Stream exceeded expected size: {next_size} > {expected_size}"
                            )
                        file.write(chunk)
                        current_size = next_size
                        if reporter is not None:
                            reporter.update(len(chunk))
                    file.flush()

            if expected_size is not None and current_size != expected_size:
                raise RuntimeError(
                    f"Stream size mismatch: wrote {current_size}, expected {expected_size}"
                )
            return
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ConnectionError,
            requests.exceptions.HTTPError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            OSError,
            TimeoutError,
            RuntimeError,
        ) as exc:
            last_error = exc
            if attempt >= max_retries:
                break
            print(
                f"\nStream download attempt {attempt}/{max_retries} failed: {exc}; retrying.",
                file=sys.stderr,
            )
            _sleep_backoff(attempt, cap=6.0)

    raise RuntimeError(
        f"Stream download failed after {max_retries} attempts. Last error: {last_error}"
    ) from last_error


def range_download(
    session: requests.Session,
    final_url: str,
    local_file_path: Path,
    content_size: int,
    num_threads: int,
) -> None:
    num_threads = max(1, int(num_threads))
    chunk_size = content_size // num_threads
    if chunk_size <= 0:
        num_threads = 1
        chunk_size = content_size

    local_file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(local_file_path, "wb") as file:
        file.truncate(content_size)

    reporter = ProgressReporter(content_size)

    def progress_callback(nbytes: int) -> None:
        reporter.update(nbytes)

    futures = []
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_byte = i * chunk_size
            end_byte = start_byte + chunk_size - 1 if i < num_threads - 1 else content_size - 1
            args = (
                session,
                final_url,
                start_byte,
                end_byte,
                content_size,
                str(local_file_path),
                progress_callback,
            )
            futures.append((executor.submit(download_chunk, args), start_byte, end_byte))

        total_written = 0
        for future, expected_start, expected_end in futures:
            actual_start, actual_end, written = future.result()
            if (actual_start, actual_end) != (expected_start, expected_end):
                raise RuntimeError(
                    f"Range result mismatch: got {actual_start}-{actual_end}, "
                    f"expected {expected_start}-{expected_end}"
                )
            expected_written = expected_end - expected_start + 1
            if written != expected_written:
                raise RuntimeError(
                    f"Range {expected_start}-{expected_end} wrote {written} bytes; "
                    f"expected {expected_written}"
                )
            total_written += written

    if total_written != content_size:
        raise RuntimeError(
            f"Combined ranges wrote {total_written} bytes; expected {content_size}"
        )
    actual_size = local_file_path.stat().st_size
    if actual_size != content_size:
        raise RuntimeError(
            f"Range output size mismatch: {actual_size} != {content_size}"
        )


def target_path(out_dir: Path, item: dict, keep_tree: bool) -> Path:
    if keep_tree:
        return out_dir / item["product"] / str(item["year"]) / item["filename"]
    return out_dir / item["filename"]


def is_complete(
    path: Path,
    size: Optional[int],
    check_size: bool,
    product_key: str = "",
    validate_cdf: bool = True,
) -> bool:
    expected_size = size if check_size else None
    valid, _message = validate_local_cdf(
        path, product_key, expected_size=expected_size, use_cdflib=validate_cdf
    )
    return valid


def download_one(
    item: dict,
    out_dir: Path,
    threads: int,
    check_size: bool,
    keep_tree: bool,
    validate_cdf: bool = True,
) -> str:
    local_path = target_path(out_dir, item, keep_tree)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    final_url, content_size, session, _range_hint = resolve_url_and_size(item["url"])
    item["remote_bytes"] = content_size
    item["local_path"] = str(local_path)

    existing_expected_size = content_size if check_size else None
    existing_valid, existing_message = validate_local_cdf(
        local_path,
        item.get("product", ""),
        expected_size=existing_expected_size,
        use_cdflib=validate_cdf,
    )
    if existing_valid:
        return "skip"

    if local_path.exists():
        print(
            f"Existing file failed validation, re-downloading: {local_path.name} "
            f"({existing_message})"
        )

    range_ok = False
    if threads > 1:
        try:
            _final2, total2 = _probe_size_by_range(session, final_url, timeout=(8, 45), attempts=3)
            if total2 > 0:
                final_url = _final2
                content_size = total2
                item["remote_bytes"] = content_size
                range_ok = True
        except Exception:
            range_ok = False

    part_path = _new_part_path(local_path)
    item["part_path"] = str(part_path)
    try:
        if threads > 1 and range_ok:
            range_download(session, final_url, part_path, content_size, threads)
        else:
            stream_download(session, final_url, part_path, content_size)

        actual_size = part_path.stat().st_size
        if actual_size != content_size:
            raise RuntimeError(
                f"Strict size check failed for {local_path.name}: "
                f"{actual_size} != {content_size}"
            )

        valid, validation_message = validate_local_cdf(
            part_path,
            item.get("product", ""),
            expected_size=content_size,
            use_cdflib=validate_cdf,
        )
        if not valid:
            raise RuntimeError(
                f"Downloaded CDF failed validation for {local_path.name}: "
                f"{validation_message}"
            )

        os.replace(part_path, local_path)
        return "download"
    finally:
        if part_path.exists():
            try:
                part_path.unlink()
            except OSError as exc:
                print(f"Warning: could not remove temporary file {part_path}: {exc}", file=sys.stderr)


def show_products() -> None:
    for key in sorted(PRODUCTS):
        product = PRODUCTS[key]
        print(f"{key:8s}  {product.path:42s}  {product.note}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download ACE CDF data from NASA CDAWeb public directories."
    )
    parser.add_argument("--date", required=False, help="YYYY-MM-DD/YYYY-MM-DD or YYYYMMDD/YYYYMMDD")
    parser.add_argument(
        "--product",
        default="mfi_h0",
        help="Comma-separated products, e.g. mfi_h0,swe_h0,epm_h1. Default: mfi_h0",
    )
    parser.add_argument("--out", default=".", help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Download threads per file")
    parser.add_argument("--check-size", default="1", help="1/0, verify local size with Content-Length")
    parser.add_argument(
        "--validate-cdf",
        default="1",
        help="1/0, validate CDF contents with cdflib when available (default: 1)",
    )
    parser.add_argument("--keep-tree", default="0", help="1/0, save as product/year/filename")
    parser.add_argument("--list-only", action="store_true", help="Only crawl and print matched files")
    parser.add_argument("--json", action="store_true", help="Print list-only result as JSON")
    parser.add_argument("--show-products", action="store_true", help="Show supported products and exit")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.show_products:
        show_products()
        return 0

    if not args.date:
        parser.error("--date is required unless --show-products is used")

    start, end = parse_date_range(args.date)
    product_keys = parse_product_list(args.product)
    out_dir = Path(args.out).expanduser()
    keep_tree = parse_bool(args.keep_tree)
    check_size = parse_bool(args.check_size)
    validate_cdf = parse_bool(args.validate_cdf)

    files = crawl_products(product_keys, start, end)

    if args.list_only:
        if args.json:
            sys.stdout.write(json.dumps(files, ensure_ascii=False, indent=2) + "\n")
        else:
            for item in files:
                print(f"{item['product']} {item['date']} {item['filename']} {item['url']}")
            print(f"Matched {len(files)} file(s).")
        return 0

    if not files:
        print("No ACE files matched this request.")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    print("========== ACE download started ==========")
    print(f"Date range: {start.isoformat()} / {end.isoformat()}")
    print(f"Products: {', '.join(product_keys)}")
    print(f"Output: {out_dir}")
    print(f"Files: {len(files)}")

    downloaded = 0
    skipped = 0
    failures: List[dict] = []
    for idx, item in enumerate(files, start=1):
        print(f"[{idx}/{len(files)}] {item['filename']}")
        t0 = time.time()
        try:
            status = download_one(
                item,
                out_dir,
                args.threads,
                check_size,
                keep_tree,
                validate_cdf=validate_cdf,
            )
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            error_text = f"{type(exc).__name__}: {exc}"
            item["error"] = error_text
            failures.append(
                {
                    "product": item.get("product", ""),
                    "date": item.get("date", ""),
                    "filename": item.get("filename", ""),
                    "error": error_text,
                }
            )
            print(
                f"FAILED [{idx}/{len(files)}] {item['filename']}: {error_text}",
                file=sys.stderr,
                flush=True,
            )
            continue
        dt_s = max(time.time() - t0, 0.001)
        local = Path(item.get("local_path", target_path(out_dir, item, keep_tree)))
        mb = local.stat().st_size / 1024 / 1024 if local.exists() else 0.0
        if status == "skip":
            skipped += 1
            print(f"Already exists: {local}")
        else:
            downloaded += 1
            print(f"Saved: {local} ({mb:.2f} MB, {mb / dt_s:.2f} MB/s)")

    print("========== ACE download finished ==========")
    print(
        f"Downloaded: {downloaded}; skipped: {skipped}; "
        f"failed: {len(failures)}; total: {len(files)}"
    )
    if failures:
        print(
            "One or more files failed; all remaining files were attempted. "
            "See the FAILED records above.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
