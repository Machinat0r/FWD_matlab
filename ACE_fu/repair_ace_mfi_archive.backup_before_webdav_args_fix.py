#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Repair an ACE/MFI CDF archive with the best product available per day.

The existing ``download_ace_files_new.py`` remains the source of truth for
CDAWeb directory crawling, URL/size resolution and HTTP transfer.  This
orchestrator adds the archive-safety behavior that a long repair run needs:

* evaluate daily candidates in measured cadence order H3 > H0 > H1 > H2;
* validate existing and newly downloaded CDF files with ``cdflib``;
* download ``*.part`` only on local NTFS, never on Z:/WebDAV;
* serially copy to a sibling ``*.uploading``, validate, then replace;
* retain an existing final file until a replacement has passed validation;
* union valid hourly coverage and fall back until all 24 UTC hours are covered;
* append one manifest row per candidate and continue after file errors.

Expected destination layout below ``--out-root``::

    ace/mfi/h3/l2/YYYY/MM/ac_h3_mfi_YYYYMMDD_vNN.cdf
    ace/mfi/h2/l2/YYYY/MM/ac_h2_mfi_YYYYMMDD_vNN.cdf
    ace/mfi/h1/l2/YYYY/MM/ac_h1_mfi_YYYYMMDD_vNN.cdf
    ace/mfi/h0/l2/YYYY/MM/ac_h0_mfi_YYYYMMDD_vNN.cdf

This file intentionally does not run anything when imported.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import logging
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DOWNLOADER = SCRIPT_DIR / "download_ace_files_new.py"
DEFAULT_OUT_ROOT = Path(r"Z:\SPART-WORK\Data\ACE")
DEFAULT_VENDOR_DIR = SCRIPT_DIR / "_vendor"
DEFAULT_STAGING_DIR = SCRIPT_DIR / "ace_repair_staging"
DEFAULT_STATE_DIR = SCRIPT_DIR / "ace_repair_state"

# Measured cadences, not numeric product-name order:
# H3=1 s, H0=16 s, H1=240 s, H2=3600 s.
PRODUCT_PRIORITY: Tuple[str, ...] = ("mfi_h3", "mfi_h0", "mfi_h1", "mfi_h2")
PRODUCT_RANK = {product: rank for rank, product in enumerate(PRODUCT_PRIORITY, 1)}
PRODUCT_TO_LEVEL = {
    "mfi_h3": "h3",
    "mfi_h2": "h2",
    "mfi_h1": "h1",
    "mfi_h0": "h0",
}
EXPECTED_SAMPLES_PER_HOUR = {
    "mfi_h3": 3600,
    "mfi_h0": 225,
    "mfi_h1": 15,
    "mfi_h2": 1,
}

MANIFEST_FIELDS: Tuple[str, ...] = (
    "run_id",
    "timestamp_utc",
    "date",
    "candidate_rank",
    "product",
    "filename",
    "url",
    "remote_bytes",
    "local_path",
    "status",
    "validation",
    "epoch_records",
    "magnitude_records",
    "q_flag_records",
    "finite_sample_count",
    "expected_samples_per_hour",
    "coverage_hours",
    "coverage_mask",
    "hourly_sample_counts",
    "missing_hours",
    "attempts",
    "elapsed_seconds",
    "error",
)

EMPTY_COVERAGE_MASK: Tuple[bool, ...] = (False,) * 24
EMPTY_HOURLY_COUNTS: Tuple[int, ...] = (0,) * 24


@dataclass(frozen=True)
class ValidationResult:
    ok: bool
    reason: str
    epoch_records: int = 0
    magnitude_records: int = 0
    q_flag_records: int = 0
    finite_sample_count: int = 0
    expected_samples_per_hour: int = 0
    coverage_mask: Tuple[bool, ...] = EMPTY_COVERAGE_MASK
    hourly_sample_counts: Tuple[int, ...] = EMPTY_HOURLY_COUNTS

    @property
    def coverage_hours(self) -> int:
        return sum(bool(value) for value in self.coverage_mask)

    @property
    def missing_hours(self) -> Tuple[int, ...]:
        return tuple(hour for hour, covered in enumerate(self.coverage_mask) if not covered)


@dataclass(frozen=True)
class CandidateResult:
    success: bool
    status: str
    remote_bytes: int
    local_path: str
    validation: ValidationResult
    attempts: int
    elapsed_seconds: float
    error: str = ""


class ManifestWriter:
    """Append CSV rows and flush each row so a stopped run remains auditable."""

    def __init__(self, path: Path) -> None:
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        existed = path.is_file() and path.stat().st_size > 0
        self._handle = path.open("a", encoding="utf-8-sig", newline="")
        self._writer = csv.DictWriter(
            self._handle, fieldnames=list(MANIFEST_FIELDS), extrasaction="ignore"
        )
        if not existed:
            self._writer.writeheader()
            self._handle.flush()

    def write(self, row: Mapping[str, Any]) -> None:
        output = {field: row.get(field, "") for field in MANIFEST_FIELDS}
        self._writer.writerow(output)
        self._handle.flush()

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "ManifestWriter":
        return self

    def __exit__(self, exc_type, exc, traceback) -> bool:
        self.close()
        return False


def parse_date(value: str) -> dt.date:
    text = str(value).strip().replace(".", "-")
    if re.fullmatch(r"\d{8}", text):
        return dt.datetime.strptime(text, "%Y%m%d").date()
    return dt.date.fromisoformat(text)


def iter_days(start: dt.date, end: dt.date) -> Iterable[dt.date]:
    day = start
    while day <= end:
        yield day
        day += dt.timedelta(days=1)


def _field(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def load_downloader(path: Path) -> Any:
    """Import the existing downloader without copying its crawler logic."""

    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"ACE downloader not found: {path}")
    module_name = "ace_existing_downloader_backend"
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import ACE downloader: {path}")
    module = importlib.util.module_from_spec(spec)
    # dataclasses and some import-time helpers expect the module to be present.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    required = (
        "crawl_products",
        "resolve_url_and_size",
        "_probe_size_by_range",
        "range_download",
        "stream_download",
    )
    missing = [name for name in required if not callable(getattr(module, name, None))]
    if missing:
        raise AttributeError(
            f"Downloader is missing required callable(s): {', '.join(missing)}"
        )
    return module


def load_cdf_dependencies(vendor_dir: Path) -> Tuple[Any, Any, Any]:
    """Load numpy and cdflib, preferring the explicitly vendored directory."""

    vendor_dir = vendor_dir.expanduser().resolve()
    if vendor_dir.is_dir():
        sys.path.insert(0, str(vendor_dir))
    try:
        import numpy as np
        from cdflib import CDF, cdfepoch
    except Exception as exc:  # pragma: no cover - depends on deployment
        raise RuntimeError(
            "cdflib/numpy are required for integrity validation. "
            f"Vendor directory checked: {vendor_dir}. Import error: {exc}"
        ) from exc
    return np, CDF, cdfepoch


def destination_path(out_root: Path, item: Mapping[str, Any]) -> Path:
    product = str(item["product"])
    if product not in PRODUCT_TO_LEVEL:
        raise ValueError(f"Unsupported MFI product: {product}")
    day = parse_date(str(item["date"]))
    return (
        out_root
        / "ace"
        / "mfi"
        / PRODUCT_TO_LEVEL[product]
        / "l2"
        / f"{day.year:04d}"
        / f"{day.month:02d}"
        / str(item["filename"])
    )


def staging_path(staging_dir: Path, item: Mapping[str, Any]) -> Path:
    """Return a local-only .part path, namespaced to avoid candidate collisions."""

    product = str(item["product"])
    day = parse_date(str(item["date"]))
    filename = Path(str(item["filename"])).name
    root = staging_dir.expanduser().resolve()
    candidate = (
        root
        / product
        / f"{day.year:04d}"
        / f"{day.month:02d}"
        / f"{filename}.part"
    ).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Unsafe staging path outside staging root: {candidate}") from exc
    return candidate


def require_local_ntfs_staging(staging_dir: Path) -> Path:
    """Refuse staging on mapped/network/non-NTFS drives."""

    staging_dir = staging_dir.expanduser().resolve()
    staging_dir.mkdir(parents=True, exist_ok=True)
    if os.name != "nt":
        raise RuntimeError("Safe staging requires Windows drive-type and NTFS checks")

    import ctypes

    root = staging_dir.anchor
    if not root:
        raise RuntimeError(f"Cannot determine staging drive: {staging_dir}")
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    drive_type = int(kernel32.GetDriveTypeW(str(root)))
    drive_fixed = 3
    if drive_type != drive_fixed:
        raise RuntimeError(
            f"--staging-dir must be on a local fixed drive, not WebDAV/network: "
            f"{staging_dir} (drive type {drive_type})"
        )

    serial = ctypes.c_ulong()
    max_component = ctypes.c_ulong()
    flags = ctypes.c_ulong()
    filesystem = ctypes.create_unicode_buffer(64)
    ok = kernel32.GetVolumeInformationW(
        str(root),
        None,
        0,
        ctypes.byref(serial),
        ctypes.byref(max_component),
        ctypes.byref(flags),
        filesystem,
        len(filesystem),
    )
    if not ok:
        raise OSError(ctypes.get_last_error(), f"Cannot inspect staging filesystem: {root}")
    if filesystem.value.upper() != "NTFS":
        raise RuntimeError(
            f"--staging-dir must be NTFS: {staging_dir} "
            f"(found {filesystem.value or 'unknown'})"
        )
    return staging_dir


def _variable_names(cdf: Any) -> List[str]:
    info = cdf.cdf_info()
    zvars = list(_field(info, "zVariables", []) or [])
    rvars = list(_field(info, "rVariables", []) or [])
    return [str(name) for name in zvars + rvars]


def _record_count(cdf: Any, variable: str) -> int:
    info = cdf.varinq(variable)
    last_record = int(_field(info, "Last_Rec", -1))
    return max(last_record + 1, 0)


def _full_epoch_values(
    cdf: Any, variable: str, records: int, np: Any, cdfepoch: Any
) -> Any:
    """Read and decode every Epoch record."""

    raw = cdf.varget(variable)
    decoded = np.asarray(cdfepoch.to_datetime(raw)).astype("datetime64[ns]").reshape(-1)
    if decoded.size != records:
        raise ValueError(
            f"Epoch record mismatch: metadata={records}, decoded={decoded.size}"
        )
    return decoded


def _full_magnitude_validity(cdf: Any, variable: str, records: int, np: Any) -> Any:
    """Read every Magnitude record and return one validity flag per record."""

    values = np.asarray(cdf.varget(variable), dtype=float)
    if values.ndim == 0:
        values = values.reshape(1, 1)
    elif values.shape[0] == records:
        values = values.reshape(records, -1)
    elif records == 1:
        values = values.reshape(1, -1)
    else:
        raise ValueError(
            f"Magnitude record mismatch: metadata={records}, array_shape={values.shape}"
        )

    if values.shape[0] != records:
        raise ValueError(
            f"Magnitude record mismatch: metadata={records}, decoded={values.shape[0]}"
        )
    valid_elements = np.isfinite(values)

    try:
        attributes = cdf.varattsget(variable)
    except Exception:
        attributes = {}

    # Fill values such as -1e31 are finite, so isfinite alone is insufficient.
    fill = _field(attributes, "FILLVAL", None)
    if fill is not None:
        try:
            fill_value = float(np.asarray(fill).reshape(-1)[0])
            valid_elements &= values != fill_value
        except (TypeError, ValueError, IndexError):
            pass

    # The archive's usable physical range is explicitly constrained to
    # 0..500 nT.  This rejects finite sentinels and implausible spikes.
    valid_elements &= (values >= 0.0) & (values <= 500.0)
    for attribute_name, operator in (("VALIDMIN", "min"), ("VALIDMAX", "max")):
        limit = _field(attributes, attribute_name, None)
        if limit is None:
            continue
        try:
            scalar_limit = float(np.asarray(limit).reshape(-1)[0])
            if operator == "min":
                valid_elements &= values >= scalar_limit
            else:
                valid_elements &= values <= scalar_limit
        except (TypeError, ValueError, IndexError):
            pass
    return np.any(valid_elements, axis=1)


def _full_q_flag_validity(cdf: Any, variable: str, records: int, np: Any) -> Any:
    """Read every Q_FLAG record; a record is usable only when every flag is 0."""

    values = np.asarray(cdf.varget(variable))
    if values.ndim == 0:
        values = values.reshape(1, 1)
    elif values.shape[0] == records:
        values = values.reshape(records, -1)
    elif records == 1:
        values = values.reshape(1, -1)
    else:
        raise ValueError(
            f"Q_FLAG record mismatch: metadata={records}, array_shape={values.shape}"
        )
    if values.shape[0] != records:
        raise ValueError(
            f"Q_FLAG record mismatch: metadata={records}, decoded={values.shape[0]}"
        )
    try:
        numeric = np.asarray(values, dtype=float)
        return np.all(np.isfinite(numeric) & (numeric == 0), axis=1)
    except (TypeError, ValueError):
        return np.all(values == 0, axis=1)


def _hourly_coverage(
    epochs: Any, valid_magnitude: Any, expected_day: dt.date, np: Any
) -> Tuple[Tuple[bool, ...], Tuple[int, ...]]:
    """Count valid target-day records in each UTC hour."""

    if epochs.size != valid_magnitude.size:
        raise ValueError(
            f"Epoch/Magnitude length mismatch: {epochs.size} != {valid_magnitude.size}"
        )
    day_start = np.datetime64(expected_day.isoformat(), "ns")
    day_end = day_start + np.timedelta64(1, "D")
    valid_epoch = ~np.isnat(epochs)
    in_day = valid_epoch & (epochs >= day_start) & (epochs < day_end)
    usable = in_day & valid_magnitude
    if np.any(usable):
        hours = ((epochs[usable] - day_start) / np.timedelta64(1, "h")).astype(int)
        counts_array = np.bincount(hours, minlength=24)[:24]
    else:
        counts_array = np.zeros(24, dtype=int)
    counts = tuple(int(value) for value in counts_array.tolist())
    mask = tuple(value > 0 for value in counts)
    return mask, counts


def validate_cdf(
    path: Path,
    expected_size: Optional[int],
    expected_day: dt.date,
    product: Optional[str],
    np: Any,
    CDF: Any,
    cdfepoch: Any,
) -> ValidationResult:
    """Fully read science variables and derive target-day hourly coverage."""

    expected_per_hour = EXPECTED_SAMPLES_PER_HOUR.get(str(product or ""), 0)

    if not path.is_file():
        return ValidationResult(
            False, "file does not exist", expected_samples_per_hour=expected_per_hour
        )
    actual_size = path.stat().st_size
    if expected_size is not None:
        if expected_size <= 0:
            return ValidationResult(
                False,
                f"invalid expected size: {expected_size}",
                expected_samples_per_hour=expected_per_hour,
            )
        if actual_size != expected_size:
            return ValidationResult(
                False,
                f"size mismatch: local={actual_size}, expected={expected_size}",
                expected_samples_per_hour=expected_per_hour,
            )

    cdf = None
    try:
        cdf = CDF(str(path))
        names = _variable_names(cdf)
        name_lookup = {name.casefold(): name for name in names}
        missing = [name for name in ("Epoch", "Magnitude") if name.casefold() not in name_lookup]
        if missing:
            return ValidationResult(
                False,
                f"missing variable(s): {', '.join(missing)}",
                expected_samples_per_hour=expected_per_hour,
            )

        epoch_name = name_lookup["epoch"]
        magnitude_name = name_lookup["magnitude"]
        q_flag_name = name_lookup.get("q_flag")
        epoch_records = _record_count(cdf, epoch_name)
        magnitude_records = _record_count(cdf, magnitude_name)
        q_flag_records = _record_count(cdf, q_flag_name) if q_flag_name else 0
        if epoch_records <= 0 or magnitude_records <= 0:
            return ValidationResult(
                ok=False,
                reason="Epoch or Magnitude has no records",
                epoch_records=epoch_records,
                magnitude_records=magnitude_records,
                q_flag_records=q_flag_records,
                expected_samples_per_hour=expected_per_hour,
            )

        if epoch_records != magnitude_records:
            return ValidationResult(
                ok=False,
                reason=f"Epoch/Magnitude record mismatch: {epoch_records} != {magnitude_records}",
                epoch_records=epoch_records,
                magnitude_records=magnitude_records,
                q_flag_records=q_flag_records,
                expected_samples_per_hour=expected_per_hour,
            )
        if q_flag_name and q_flag_records != epoch_records:
            return ValidationResult(
                ok=False,
                reason=f"Epoch/Magnitude/Q_FLAG record mismatch: "
                f"{epoch_records}/{magnitude_records}/{q_flag_records}",
                epoch_records=epoch_records,
                magnitude_records=magnitude_records,
                q_flag_records=q_flag_records,
                expected_samples_per_hour=expected_per_hour,
            )

        epochs = _full_epoch_values(
            cdf, epoch_name, epoch_records, np=np, cdfepoch=cdfepoch
        )
        valid_magnitude = _full_magnitude_validity(
            cdf, magnitude_name, magnitude_records, np=np
        )
        if q_flag_name:
            valid_magnitude &= _full_q_flag_validity(
                cdf, q_flag_name, q_flag_records, np=np
            )
        finite_count = int(np.count_nonzero(valid_magnitude))
        coverage_mask, hourly_counts = _hourly_coverage(
            epochs, valid_magnitude, expected_day, np=np
        )
        target_count = sum(hourly_counts)
        if finite_count <= 0 or target_count <= 0:
            return ValidationResult(
                ok=False,
                reason=f"no valid Magnitude records for target date {expected_day.isoformat()}",
                epoch_records=epoch_records,
                magnitude_records=magnitude_records,
                q_flag_records=q_flag_records,
                finite_sample_count=finite_count,
                expected_samples_per_hour=expected_per_hour,
                coverage_mask=coverage_mask,
                hourly_sample_counts=hourly_counts,
            )
        coverage_hours = sum(coverage_mask)
        return ValidationResult(
            ok=True,
            reason=f"ok; hourly_coverage={coverage_hours}/24; "
            f"expected_samples_per_hour={expected_per_hour}",
            epoch_records=epoch_records,
            magnitude_records=magnitude_records,
            q_flag_records=q_flag_records,
            finite_sample_count=finite_count,
            expected_samples_per_hour=expected_per_hour,
            coverage_mask=coverage_mask,
            hourly_sample_counts=hourly_counts,
        )
    except Exception as exc:
        return ValidationResult(
            False,
            f"CDF validation error: {type(exc).__name__}: {exc}",
            expected_samples_per_hour=expected_per_hour,
        )
    finally:
        close = getattr(cdf, "close", None)
        if callable(close):
            try:
                close()
            except Exception:
                pass


def _download_to_staging(
    backend: Any,
    item: Mapping[str, Any],
    part_path: Path,
    threads: int,
    resolved: Tuple[str, int, Any, bool],
) -> Tuple[int, str]:
    """Transfer one already-resolved candidate into a local NTFS .part file."""

    final_url, remote_bytes, session, range_hint = resolved
    remote_bytes = int(remote_bytes)
    if remote_bytes <= 0:
        raise RuntimeError(f"Remote Content-Length is invalid: {remote_bytes}")

    if part_path.exists():
        part_path.unlink()
    part_path.parent.mkdir(parents=True, exist_ok=True)

    # HEAD already reported Accept-Ranges in the common CDAWeb path.  Avoid a
    # redundant bytes=0-0 request for every one of ~10,000 daily files; probe
    # only when HEAD did not advertise range support.
    range_ok = threads > 1 and bool(range_hint)
    if threads > 1 and not range_ok:
        try:
            probed_url, probed_bytes = backend._probe_size_by_range(
                session, final_url, timeout=(8, 45), attempts=3
            )
            probed_bytes = int(probed_bytes)
            if probed_bytes != remote_bytes:
                raise RuntimeError(
                    f"HEAD/range size disagreement: {remote_bytes} != {probed_bytes}"
                )
            final_url = probed_url
            range_ok = True
        except Exception as exc:
            logging.getLogger("ace_mfi_repair").warning(
                "Range probe unavailable for %s; using stream download (%s)",
                item["filename"],
                exc,
            )

    if threads > 1 and range_ok and range_hint:
        backend.range_download(session, final_url, part_path, remote_bytes, threads)
    elif threads > 1 and range_ok:
        # Some servers omit Accept-Ranges on HEAD but prove range support above.
        backend.range_download(session, final_url, part_path, remote_bytes, threads)
    else:
        backend.stream_download(session, final_url, part_path, remote_bytes)

    actual_bytes = part_path.stat().st_size if part_path.is_file() else -1
    if actual_bytes != remote_bytes:
        raise RuntimeError(
            f"Downloaded size mismatch: local={actual_bytes}, remote={remote_bytes}"
        )
    return remote_bytes, final_url


def _copy_staging_to_archive(
    staging_file: Path,
    final_path: Path,
    expected_day: dt.date,
    product: str,
    np: Any,
    CDF: Any,
    cdfepoch: Any,
) -> ValidationResult:
    """Serially copy local staging to sibling .uploading, validate, then replace."""

    uploading_path = final_path.with_name(final_path.name + ".uploading")
    expected_size = staging_file.stat().st_size
    ensure_archive_parent(final_path.parent)
    upload_validated = False

    try:
        uploaded = validate_cdf(
            uploading_path, expected_size, expected_day, product,
            np=np, CDF=CDF, cdfepoch=cdfepoch,
        )
        if not uploaded.ok:
            if uploading_path.exists():
                uploading_path.unlink()
            # Deliberately single-threaded.  No downloader worker ever receives
            # Z:/WebDAV, whose only temporary object is this .uploading file.
            with staging_file.open("rb") as source, uploading_path.open("wb") as target:
                while True:
                    chunk = source.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    target.write(chunk)
                target.flush()
            uploaded = validate_cdf(
                uploading_path, expected_size, expected_day, product,
                np=np, CDF=CDF, cdfepoch=cdfepoch,
            )
            if not uploaded.ok:
                raise RuntimeError(f".uploading validation failed: {uploaded.reason}")
        upload_validated = True

        # Same-directory replace is the only operation that changes the final
        # archive name.  WebDAV may reject overwrite semantics, so if direct
        # replace fails, revalidate the old final and delete it only when it is
        # still conclusively bad before retrying the rename.
        try:
            os.replace(uploading_path, final_path)
        except OSError as direct_replace_error:
            current_final = validate_cdf(
                final_path, None, expected_day, product,
                np=np, CDF=CDF, cdfepoch=cdfepoch,
            )
            if current_final.ok:
                raise RuntimeError(
                    "Direct replace failed and current final is valid; refusing to delete it"
                ) from direct_replace_error
            if final_path.exists():
                final_path.unlink()
            try:
                os.replace(uploading_path, final_path)
            except OSError as fallback_replace_error:
                raise RuntimeError(
                    f"WebDAV replace failed after invalid-final removal: {fallback_replace_error}"
                ) from direct_replace_error

        # ``uploaded`` was fully parsed from the Z:/WebDAV .uploading object.
        # A same-directory rename changes only the name, not its bytes.  Check
        # existence, exact size and non-zero CDF header after the rename, then
        # retain the already-computed full validation result.  This removes a
        # redundant second full WebDAV read while preserving the content gate.
        if not final_path.is_file():
            raise RuntimeError("final file is absent after WebDAV rename")
        final_size = final_path.stat().st_size
        if final_size != expected_size:
            raise RuntimeError(
                f"final size after WebDAV rename is {final_size}, expected {expected_size}"
            )
        with final_path.open("rb") as handle:
            final_header = handle.read(16)
        if len(final_header) != 16 or not any(final_header):
            raise RuntimeError("final CDF header is missing or all zero after WebDAV rename")
        return uploaded
    except Exception:
        # Retain an already-validated .uploading when rename failed.  A retry
        # can install it without another Z:/WebDAV copy.
        if uploading_path.exists() and not upload_validated:
            try:
                uploading_path.unlink()
            except OSError:
                pass
        raise


def ensure_archive_parent(parent: Path) -> None:
    """Create an archive directory, with a RaiDrive-compatible fallback.

    RaiDrive/WebDAV can reject Python's recursive ``os.mkdir`` even though the
    same signed-in user can create the directory through PowerShell.  Existing
    directories take the normal fast path; the helper is invoked only after a
    genuine access-denied result.
    """

    try:
        parent.mkdir(parents=True, exist_ok=True)
        return
    except OSError as direct_error:
        if os.name != "nt" or getattr(direct_error, "winerror", None) != 5:
            raise

        powershell = Path(os.environ.get(
            "WINDIR", r"C:\Windows"
        )) / "System32" / "WindowsPowerShell" / "v1.0" / "powershell.exe"
        command = (
            "$ErrorActionPreference='Stop'; "
            "New-Item -ItemType Directory -LiteralPath $args[0] -Force | Out-Null"
        )
        completed = subprocess.run(
            [str(powershell), "-NoProfile", "-NonInteractive", "-Command", command, str(parent)],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if completed.returncode != 0 or not parent.is_dir():
            detail = (completed.stderr or completed.stdout).strip()
            raise PermissionError(
                f"Python and PowerShell could not create archive directory {parent}: "
                f"{detail or direct_error}"
            ) from direct_error
        logging.getLogger("ace_mfi_repair").info(
            "Created WebDAV archive directory through PowerShell fallback: %s",
            parent,
        )


def process_candidate(
    backend: Any,
    item: Mapping[str, Any],
    out_root: Path,
    staging_dir: Path,
    threads: int,
    max_attempts: int,
    np: Any,
    CDF: Any,
    cdfepoch: Any,
) -> CandidateResult:
    """Validate/repair one candidate while preserving final-file atomicity."""

    started = time.monotonic()
    day = parse_date(str(item["date"]))
    product = str(item["product"])
    final_path = destination_path(out_root, item)
    part_path = staging_path(staging_dir, item)
    uploading_path = final_path.with_name(final_path.name + ".uploading")
    remote_bytes = 0
    last_error = ""
    last_validation = ValidationResult(False, "not validated")

    # Content-first validation avoids a HEAD request for every already-good
    # archive file.  Full Epoch/Magnitude/Q_FLAG reads make a size-only skip
    # unnecessary and unsafe.
    existing = validate_cdf(
        final_path,
        None,
        day,
        product,
        np=np,
        CDF=CDF,
        cdfepoch=cdfepoch,
    )
    last_validation = existing
    if existing.ok:
        for stale in (part_path, uploading_path):
            if stale.exists():
                try:
                    stale.unlink()
                except OSError as cleanup_exc:
                    logging.getLogger("ace_mfi_repair").warning(
                        "Valid final retained, but stale temporary file could not be removed: %s",
                        cleanup_exc,
                    )
        return CandidateResult(
            True,
            "valid_existing",
            0,
            str(final_path),
            existing,
            0,
            time.monotonic() - started,
        )

    for attempt in range(1, max_attempts + 1):
        try:
            staged = validate_cdf(
                part_path,
                None,
                day,
                product,
                np=np,
                CDF=CDF,
                cdfepoch=cdfepoch,
            )
            reused_staging = staged.ok
            if not staged.ok:
                if part_path.exists():
                    part_path.unlink()
                # HEAD is reached only after both final and local staging fail
                # full content validation.
                resolved = backend.resolve_url_and_size(item["url"])
                _, remote_bytes, _, _ = resolved
                remote_bytes = int(remote_bytes)
                remote_bytes, _ = _download_to_staging(
                    backend, item, part_path, threads=threads, resolved=resolved
                )
                staged = validate_cdf(
                    part_path,
                    remote_bytes,
                    day,
                    product,
                    np=np,
                    CDF=CDF,
                    cdfepoch=cdfepoch,
                )
                if not staged.ok:
                    raise RuntimeError(f"local staging validation failed: {staged.reason}")
            else:
                remote_bytes = part_path.stat().st_size

            last_validation = staged
            installed = _copy_staging_to_archive(
                part_path,
                final_path,
                day,
                product,
                np=np,
                CDF=CDF,
                cdfepoch=cdfepoch,
            )
            try:
                part_path.unlink()
            except OSError as cleanup_exc:
                logging.getLogger("ace_mfi_repair").warning(
                    "Installed final file, but local staging cleanup failed: %s",
                    cleanup_exc,
                )
            return CandidateResult(
                True,
                "installed_from_staging" if reused_staging else "downloaded",
                remote_bytes,
                str(final_path),
                installed,
                attempt,
                time.monotonic() - started,
            )
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            # Keep a fully valid local staging file after an upload/rename
            # failure so the next attempt does not redownload it.
            if part_path.exists():
                retained = validate_cdf(
                    part_path,
                    None,
                    day,
                    product,
                    np=np,
                    CDF=CDF,
                    cdfepoch=cdfepoch,
                )
                if retained.ok:
                    last_validation = retained
                else:
                    try:
                        part_path.unlink()
                    except OSError as cleanup_exc:
                        last_error += f"; invalid local .part cleanup failed: {cleanup_exc}"
            if uploading_path.exists():
                upload_retained = validate_cdf(
                    uploading_path,
                    part_path.stat().st_size if part_path.exists() else None,
                    day,
                    product,
                    np=np,
                    CDF=CDF,
                    cdfepoch=cdfepoch,
                )
                if upload_retained.ok:
                    logging.getLogger("ace_mfi_repair").warning(
                        "Retaining validated .uploading for rename retry: %s",
                        uploading_path,
                    )
                else:
                    try:
                        uploading_path.unlink()
                    except OSError as cleanup_exc:
                        last_error += f"; stale .uploading cleanup failed: {cleanup_exc}"
            if attempt < max_attempts:
                delay = min(2 ** (attempt - 1), 30)
                logging.getLogger("ace_mfi_repair").warning(
                    "Candidate %s attempt %d/%d failed; retrying in %ds: %s",
                    item["filename"],
                    attempt,
                    max_attempts,
                    delay,
                    last_error,
                )
                time.sleep(delay)

    return CandidateResult(
        False,
        "candidate_failed",
        remote_bytes,
        str(final_path),
        last_validation,
        max_attempts,
        time.monotonic() - started,
        last_error,
    )


def group_candidates(items: Sequence[Mapping[str, Any]]) -> Dict[dt.date, List[dict]]:
    grouped: Dict[dt.date, List[dict]] = {}
    for raw in items:
        item = dict(raw)
        product = str(item.get("product", ""))
        if product not in PRODUCT_RANK:
            continue
        day = parse_date(str(item["date"]))
        grouped.setdefault(day, []).append(item)
    for candidates in grouped.values():
        candidates.sort(
            key=lambda item: (
                PRODUCT_RANK[str(item["product"])],
                -int(item.get("version", 0)),
                str(item.get("filename", "")),
            )
        )
    return grouped


def _utc_timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def manifest_row(
    run_id: str,
    day: dt.date,
    rank: Any,
    item: Optional[Mapping[str, Any]],
    result: Optional[CandidateResult],
    status: Optional[str] = None,
    error: str = "",
    validation_override: Optional[ValidationResult] = None,
) -> Dict[str, Any]:
    item = item or {}
    validation = validation_override or (
        result.validation if result is not None else ValidationResult(False, "")
    )
    return {
        "run_id": run_id,
        "timestamp_utc": _utc_timestamp(),
        "date": day.isoformat(),
        "candidate_rank": rank,
        "product": item.get("product", ""),
        "filename": item.get("filename", ""),
        "url": item.get("url", ""),
        "remote_bytes": result.remote_bytes if result is not None else "",
        "local_path": result.local_path if result is not None else "",
        "status": status or (result.status if result is not None else ""),
        "validation": validation.reason,
        "epoch_records": validation.epoch_records,
        "magnitude_records": validation.magnitude_records,
        "q_flag_records": validation.q_flag_records,
        "finite_sample_count": validation.finite_sample_count,
        "expected_samples_per_hour": validation.expected_samples_per_hour,
        "coverage_hours": validation.coverage_hours,
        "coverage_mask": "".join("1" if value else "0" for value in validation.coverage_mask),
        "hourly_sample_counts": ",".join(
            str(value) for value in validation.hourly_sample_counts
        ),
        "missing_hours": ",".join(f"{hour:02d}" for hour in validation.missing_hours),
        "attempts": result.attempts if result is not None else "",
        "elapsed_seconds": (
            f"{result.elapsed_seconds:.3f}" if result is not None else ""
        ),
        "error": error or (result.error if result is not None else ""),
    }


def configure_logging(log_path: Path, verbose: bool) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("ace_mfi_repair")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    stream = logging.StreamHandler(sys.stdout)
    stream.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Repair ACE/MFI CDF files with H3>H0>H1>H2 hourly fallback."
    )
    parser.add_argument(
        "--downloader-path", type=Path, default=DEFAULT_DOWNLOADER,
        help="Existing download_ace_files_new.py to import."
    )
    parser.add_argument(
        "--out-root", type=Path, default=DEFAULT_OUT_ROOT,
        help="Root containing ace/mfi/hX/l2/YYYY/MM."
    )
    parser.add_argument("--start", type=parse_date, default=dt.date(1997, 9, 2))
    parser.add_argument("--end", type=parse_date, default=dt.date.today())
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--crawl-attempts", type=int, default=2)
    parser.add_argument("--candidate-attempts", type=int, default=3)
    parser.add_argument(
        "--staging-dir", type=Path, default=DEFAULT_STAGING_DIR,
        help="Local fixed NTFS directory for all .part downloads."
    )
    parser.add_argument(
        "--state-dir", type=Path, default=DEFAULT_STATE_DIR,
        help="Local directory for manifests and text logs."
    )
    parser.add_argument(
        "--vendor-dir", type=Path, default=DEFAULT_VENDOR_DIR,
        help="Directory containing vendored cdflib/numpy packages."
    )
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--log", type=Path, default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.end < args.start:
        raise SystemExit("--end must be on or after --start")
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    if args.crawl_attempts < 1:
        raise SystemExit("--crawl-attempts must be >= 1")
    if args.candidate_attempts < 1:
        raise SystemExit("--candidate-attempts must be >= 1")

    run_id = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    staging_dir = require_local_ntfs_staging(args.staging_dir)
    state_dir = args.state_dir.expanduser().resolve()
    manifest_path = args.manifest or state_dir / f"ACE_MFI_repair_manifest_{run_id}.csv"
    log_path = args.log or state_dir / f"ACE_MFI_repair_{run_id}.log"
    logger = configure_logging(log_path, args.verbose)

    logger.info("Loading downloader: %s", args.downloader_path)
    backend = load_downloader(args.downloader_path)
    np, CDF, cdfepoch = load_cdf_dependencies(args.vendor_dir)
    logger.info(
        "Crawling %s to %s for %s",
        args.start,
        args.end,
        ">".join(PRODUCT_PRIORITY),
    )
    # Crawl each product independently for clear diagnostics.  If any product
    # remains unavailable after its retries, abort before archive writes: a
    # crawl outage must never be mistaken for a legitimate low-cadence fallback.
    remote_items: List[dict] = []
    crawl_failures: List[str] = []
    for product in PRODUCT_PRIORITY:
        product_crawled = False
        for attempt in range(1, args.crawl_attempts + 1):
            try:
                product_items = backend.crawl_products((product,), args.start, args.end)
                remote_items.extend(product_items)
                logger.info("Crawl %s: %d daily file(s)", product, len(product_items))
                product_crawled = True
                break
            except Exception as exc:
                logger.error(
                    "Crawl %s attempt %d/%d failed: %s: %s",
                    product,
                    attempt,
                    args.crawl_attempts,
                    type(exc).__name__,
                    exc,
                )
                if attempt < args.crawl_attempts:
                    time.sleep(min(2 ** (attempt - 1), 30))
        if not product_crawled:
            crawl_failures.append(product)
    if crawl_failures:
        logger.critical(
            "Aborting before archive writes because product crawl failed: %s",
            ",".join(crawl_failures),
        )
        return 3
    if not remote_items:
        logger.critical(
            "Crawl succeeded but returned no remote files; cannot derive mission/publication bounds"
        )
        return 3

    candidates_by_day = group_candidates(remote_items)
    remote_dates = sorted(candidates_by_day)
    remote_start = remote_dates[0]
    latest_available = remote_dates[-1]
    logger.info(
        "Remote availability inferred from crawl manifest: %s through %s",
        remote_start,
        latest_available,
    )

    counters = {
        "requested_days": 0,
        "remote_interval_days": 0,
        "day_complete": 0,
        "candidate_valid": 0,
        "downloaded": 0,
        "installed_from_staging": 0,
        "valid_existing": 0,
        "no_remote_candidate": 0,
        "day_failed": 0,
        "source_data_gap": 0,
        "fallback_completed": 0,
        "pre_mission": 0,
        "publication_lag": 0,
    }
    with ManifestWriter(manifest_path) as manifest:
        for day in iter_days(args.start, args.end):
            counters["requested_days"] += 1
            if day < remote_start:
                counters["pre_mission"] += 1
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="pre_mission",
                        error=f"Earlier than first remotely listed date {remote_start}",
                    )
                )
                continue
            if day > latest_available:
                counters["publication_lag"] += 1
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="publication_lag",
                        error=f"Later than latest remotely listed date {latest_available}",
                    )
                )
                continue

            counters["remote_interval_days"] += 1
            candidates = candidates_by_day.get(day, [])
            if not candidates:
                counters["no_remote_candidate"] += 1
                counters["source_data_gap"] += 1
                empty_summary = ValidationResult(
                    ok=False,
                    reason="source_data_gap; no remote H3/H0/H1/H2 candidate",
                )
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="source_data_gap",
                        error="No H3/H0/H1/H2 file was listed by CDAWeb",
                        validation_override=empty_summary,
                    )
                )
                logger.warning("%s: no remote MFI candidate", day)
                continue

            completed = False
            any_valid_candidate = False
            union_mask = [False] * 24
            union_counts = [0] * 24
            for rank, item in enumerate(candidates, 1):
                try:
                    result = process_candidate(
                        backend,
                        item,
                        args.out_root,
                        staging_dir,
                        args.threads,
                        args.candidate_attempts,
                        np=np,
                        CDF=CDF,
                        cdfepoch=cdfepoch,
                    )
                except Exception as exc:
                    # A programming or unforeseen per-file exception must not
                    # abort a decades-long repair run.
                    result = CandidateResult(
                        False,
                        "candidate_failed",
                        0,
                        str(destination_path(args.out_root, item)),
                        ValidationResult(False, "unhandled candidate exception"),
                        0,
                        0.0,
                        f"{type(exc).__name__}: {exc}",
                    )
                manifest.write(manifest_row(run_id, day, rank, item, result))
                if result.success:
                    any_valid_candidate = True
                    counters["candidate_valid"] += 1
                    counters[result.status] += 1
                    for hour in range(24):
                        union_mask[hour] = (
                            union_mask[hour] or result.validation.coverage_mask[hour]
                        )
                        union_counts[hour] += result.validation.hourly_sample_counts[hour]
                    logger.info(
                        "%s: %s via %s (rank %d, coverage %d/24): %s",
                        day,
                        result.status,
                        item["product"],
                        rank,
                        result.validation.coverage_hours,
                        item["filename"],
                    )
                    if all(union_mask):
                        completed = True
                        if rank > 1:
                            counters["fallback_completed"] += 1
                        break
                    logger.info(
                        "%s: continuing to lower cadence; union missing UTC hour(s): %s",
                        day,
                        ",".join(f"{hour:02d}" for hour, value in enumerate(union_mask) if not value),
                    )
                    continue
                logger.error(
                    "%s: candidate %s failed: %s",
                    day,
                    item["product"],
                    result.error or result.validation.reason,
                )

            union_validation = ValidationResult(
                ok=completed,
                reason=(
                    "day_complete; coverage_union=24/24"
                    if completed
                    else (
                        f"source_data_gap; coverage_union={sum(union_mask)}/24"
                        if any_valid_candidate
                        else "day_failed; no candidate completed validation/transfer"
                    )
                ),
                finite_sample_count=sum(union_counts),
                coverage_mask=tuple(union_mask),
                hourly_sample_counts=tuple(union_counts),
            )
            if completed:
                counters["day_complete"] += 1
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="day_complete",
                        validation_override=union_validation,
                    )
                )
            elif any_valid_candidate:
                counters["source_data_gap"] += 1
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="source_data_gap",
                        error="All H3/H0/H1/H2 candidates exhausted; hourly union is incomplete",
                        validation_override=union_validation,
                    )
                )
                logger.warning(
                    "%s: source_data_gap; missing UTC hour(s): %s",
                    day,
                    ",".join(f"{hour:02d}" for hour, value in enumerate(union_mask) if not value),
                )
            else:
                counters["day_failed"] += 1
                manifest.write(
                    manifest_row(
                        run_id,
                        day,
                        "",
                        None,
                        None,
                        status="day_failed",
                        error="Every listed H3/H0/H1/H2 candidate failed validation or transfer",
                        validation_override=union_validation,
                    )
                )

    logger.info("Repair finished: %s", ", ".join(f"{k}={v}" for k, v in counters.items()))
    logger.info("Manifest: %s", manifest_path)
    logger.info("Log: %s", log_path)
    return 0 if counters["day_complete"] == counters["remote_interval_days"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
