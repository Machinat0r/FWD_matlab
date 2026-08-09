#!/usr/bin/env python3
"""Convert selected NASA Voyager CDF variables to a MATLAB-safe binary cache."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROGRAM_ROOT = Path(__file__).resolve().parent
VENDORED_PACKAGES = PROGRAM_ROOT / "python_packages"
if VENDORED_PACKAGES.is_dir():
    sys.path.insert(0, str(VENDORED_PACKAGES))

try:
    import numpy as np
    from cdflib import CDF, cdfepoch
except ImportError as exc:  # pragma: no cover - exercised by MATLAB setup checks
    raise SystemExit(
        "Missing Python dependency. Install cdflib and numpy, or retain the "
        f"bundled python_packages directory. Original error: {exc}"
    ) from exc


BRIDGE_VERSION = "1.0.0"

PROFILE_VARIABLES: dict[str, tuple[str, ...]] = {
    "mag48s": (
        "Epoch",
        "Epoch_ephem",
        "spacecraftID",
        "F1",
        "BR",
        "BT",
        "BN",
        "dF",
        "dBR",
        "dBT",
        "dBN",
        "Radius",
        "hg_lat",
        "hg_lon",
        "hgi_lon",
    ),
    "coho": (
        "Epoch",
        "heliocentricDistance",
        "heliographicLatitude",
        "heliographicLongitude",
        "ABS_B",
        "F",
        "BR",
        "BT",
        "BN",
        "V",
        "elevAngle",
        "azimuthAngle",
        "protonDensity",
        "protonTemp",
    ),
    "pls": (
        "Epoch",
        "chi2",
        "V_rtn",
        "xV_rtn",
        "V",
        "xV",
        "dens",
        "xdens",
        "w",
        "xw",
        "T",
        "xT",
    ),
    "mag2s": (
        "Epoch",
        "Epoch2",
        "F1",
        "F2",
        "B1",
        "B2",
        "B3",
        "RMSB1",
        "RMSB2",
        "RMSB3",
        "quality_flag",
        "dataConfidence",
        "numDetailPoints",
        "scDistance",
    ),
}

PROFILE_PATTERNS: dict[str, tuple[re.Pattern[str], ...]] = {
    "coho": (
        re.compile(r"^protonFlux\d+_LECP$"),
        re.compile(r"^protonFlux\d+_CRS$"),
    )
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--cache-root", required=True)
    parser.add_argument("--profile", required=True, choices=sorted(PROFILE_VARIABLES))
    return parser.parse_args()


def safe_json(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [safe_json(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return safe_json(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    return str(value)


def number_attribute(attributes: dict[str, Any], name: str) -> float | None:
    raw = attributes.get(name)
    if raw is None:
        return None
    try:
        return float(np.asarray(raw).reshape(-1)[0])
    except (TypeError, ValueError, OverflowError):
        return None


def source_identity(source: Path, profile: str) -> tuple[os.stat_result, str]:
    stat = source.stat()
    identity = "|".join(
        (
            BRIDGE_VERSION,
            str(source.resolve()),
            str(stat.st_size),
            str(stat.st_mtime_ns),
            profile,
        )
    )
    return stat, hashlib.sha256(identity.encode("utf-8")).hexdigest()


def stage_source(source: Path, cache_root: Path, stat: os.stat_result, key: str) -> Path:
    stage_dir = cache_root / "cdf"
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged = stage_dir / f"{key}{source.suffix.lower()}"
    if staged.is_file() and staged.stat().st_size == stat.st_size:
        return staged

    temporary = staged.with_suffix(staged.suffix + f".{os.getpid()}.partial")
    if temporary.exists():
        temporary.unlink()
    shutil.copyfile(source, temporary)
    if temporary.stat().st_size != stat.st_size:
        temporary.unlink(missing_ok=True)
        raise OSError(f"Incomplete staging copy for {source}")
    os.replace(temporary, staged)
    return staged


def selected_variables(cdf: CDF, profile: str) -> list[str]:
    info = cdf.cdf_info()
    available = list(info.zVariables) + list(info.rVariables)
    requested = set(PROFILE_VARIABLES[profile])
    patterns = PROFILE_PATTERNS.get(profile, ())
    selected: list[str] = []
    for name in available:
        if name in requested or any(pattern.match(name) for pattern in patterns):
            selected.append(name)
    return selected


def is_epoch(description: str) -> bool:
    upper = description.upper()
    return upper.startswith("CDF_EPOCH") or upper == "CDF_TIME_TT2000"


def epoch_to_unix_seconds(raw: np.ndarray) -> np.ndarray:
    converted = np.asarray(cdfepoch.to_datetime(raw)).astype("datetime64[ns]")
    epoch = np.datetime64("1970-01-01T00:00:00", "ns")
    return (converted - epoch).astype("timedelta64[ns]").astype(np.float64) / 1.0e9


def numeric_array(raw: np.ndarray, attributes: dict[str, Any]) -> np.ndarray:
    output = np.asarray(raw, dtype=np.float64)
    fill = number_attribute(attributes, "FILLVAL")
    valid_min = number_attribute(attributes, "VALIDMIN")
    valid_max = number_attribute(attributes, "VALIDMAX")
    if fill is not None:
        output[output == fill] = np.nan
    if valid_min is not None:
        output[output < valid_min] = np.nan
    if valid_max is not None:
        output[output > valid_max] = np.nan
    output[~np.isfinite(output) | (np.abs(output) >= 1.0e30)] = np.nan
    return output


def write_cache(
    source: Path,
    staged: Path,
    cache_root: Path,
    profile: str,
    stat: os.stat_result,
    key: str,
) -> Path:
    converted_dir = cache_root / "converted"
    converted_dir.mkdir(parents=True, exist_ok=True)
    metadata_file = converted_dir / f"{key}.json"
    binary_file = converted_dir / f"{key}.bin"
    if metadata_file.is_file() and binary_file.is_file():
        try:
            metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
            if (
                metadata.get("bridge_version") == BRIDGE_VERSION
                and metadata.get("source_size") == stat.st_size
                and metadata.get("profile") == profile
                and binary_file.stat().st_size == metadata.get("binary_size")
            ):
                return metadata_file
        except (OSError, ValueError, TypeError):
            pass

    cdf = CDF(str(staged))
    variables = selected_variables(cdf, profile)
    temp_binary = binary_file.with_suffix(f".bin.{os.getpid()}.partial")
    temp_metadata = metadata_file.with_suffix(f".json.{os.getpid()}.partial")
    variable_metadata: list[dict[str, Any]] = []
    offset = 0

    with temp_binary.open("wb") as stream:
        for name in variables:
            inquiry = cdf.varinq(name)
            attributes = cdf.varattsget(name)
            raw = np.asarray(cdf.varget(name))
            if is_epoch(inquiry.Data_Type_Description):
                array = epoch_to_unix_seconds(raw)
                time_variable = True
            elif np.issubdtype(raw.dtype, np.number):
                array = numeric_array(raw, attributes)
                time_variable = False
            else:
                continue
            if array.ndim == 0:
                array = array.reshape(1)
            array = np.asarray(array, dtype="<f8", order="F")
            flat = array.ravel(order="F")
            flat.tofile(stream)
            nbytes = int(flat.size * 8)
            variable_metadata.append(
                {
                    "name": name,
                    "shape": list(array.shape),
                    "offset_bytes": offset,
                    "count": int(flat.size),
                    "is_time": time_variable,
                    "cdf_type": inquiry.Data_Type_Description,
                    "attributes": safe_json(attributes),
                }
            )
            offset += nbytes

    metadata = {
        "schema_version": 1,
        "bridge_version": BRIDGE_VERSION,
        "source_file": str(source),
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
        "profile": profile,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "binary_file": binary_file.name,
        "binary_size": offset,
        "variables": variable_metadata,
    }
    temp_metadata.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(temp_binary, binary_file)
    os.replace(temp_metadata, metadata_file)
    return metadata_file


def main() -> int:
    args = parse_args()
    source = Path(args.source).expanduser()
    cache_root = Path(args.cache_root).expanduser()
    if not source.is_file():
        raise FileNotFoundError(source)
    cache_root.mkdir(parents=True, exist_ok=True)
    stat, key = source_identity(source, args.profile)
    staged = stage_source(source, cache_root, stat, key)
    metadata = write_cache(source, staged, cache_root, args.profile, stat, key)
    print(f"CACHE_JSON|{metadata}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # concise diagnostics for MATLAB system()
        print(f"ERROR|{type(exc).__name__}|{exc}", file=sys.stderr)
        raise
