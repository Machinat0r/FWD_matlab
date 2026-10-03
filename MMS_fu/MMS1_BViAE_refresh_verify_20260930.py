"""核对临时重画图件的记录、尺寸与原 B/AE 区域；不修改任何图件。"""
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
from PIL import Image

root = Path(r"Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815")
audit = root / "Vi_completeness_20260930"
redraw = audit / "redraw"
source = Path(r"C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\B_Vi_AE")
staging = Path(os.environ["TEMP"]) / "MMS1_BViAE_refresh_20260930" / "plots" / "B_Vi_AE"
expected = json.loads((audit / "local" / "after_download" / "windows_audit.json").read_text(encoding="utf-8"))
gap_records = json.loads((audit / "local" / "after_download" / "window_gaps.json").read_text(encoding="utf-8"))
assert len(expected) == 26
assert len(list(staging.glob("*.png"))) == 26
rows = []
for entry in expected:
    ident = entry["id"]
    old_record = json.loads((root / (ident + "_overview.json")).read_text(encoding="utf-8"))
    new_record = json.loads((redraw / (ident + "_overview.json")).read_text(encoding="utf-8"))
    assert new_record["complete"] and not new_record["error"], ident
    assert new_record["nativeCounts"][1] == [entry["burstGsmAnyFinite"], entry["fastGsmAnyFinite"]], ident
    assert new_record["nativeCounts"][0] == old_record["nativeCounts"][0], ident
    unchanged_keys = ["id", "startUTC", "endUTC", "spacecraft", "coordinateSystem", "midpointXYZ_RE",
                      "nativeCountNames", "modeOrder", "aeRows", "aeValid", "aeSources",
                      "aeCadenceSeconds", "timeLabelRotation", "eventMarkers"]
    for key in unchanged_keys:
        assert new_record[key] == old_record[key], (ident, key)
    assert new_record["aeCadenceSeconds"] == 60
    assert new_record["panelAvailable"] == [True, True, True]
    assert len(new_record["png"]) == 1
    new_path = staging / (ident + "_MMS1_B_Vi_AE.png")
    old_path = source / new_path.name
    assert Path(new_record["png"][0]).resolve() == new_path.resolve()
    with Image.open(new_path) as im:
        im.verify()
    new_image = np.array(Image.open(new_path).convert("RGB"))
    old_image = np.array(Image.open(old_path).convert("RGB"))
    same_size = new_image.shape == old_image.shape
    row = {"id": ident, "newPNG": str(new_path), "sourcePNG": str(old_path),
           "newSize": [new_image.shape[1], new_image.shape[0]], "oldSize": [old_image.shape[1], old_image.shape[0]],
           "sameSize": same_size, "nativeVi": new_record["nativeCounts"][1],
           "allRecordChecksPassed": True,
           "newSHA256": hashlib.sha256(new_path.read_bytes()).hexdigest(),
           "oldSHA256": hashlib.sha256(old_path.read_bytes()).hexdigest()}
    if same_size:
        # 原图高度约919：这些区域包含完整标题/B曲线及AE曲线/时间标签。
        # 避开Vi panel及相邻边框，以便检查重画后未改变的区域。
        h = new_image.shape[0]
        b_end = round(h * 0.326)
        ae_begin = round(h * 0.643)
        row["B_region_pixel_equal"] = bool(np.array_equal(new_image[:b_end], old_image[:b_end]))
        row["AE_region_pixel_equal"] = bool(np.array_equal(new_image[ae_begin:], old_image[ae_begin:]))
    # 图中真实空档的独立像素检查：寻找Vi左右坐标轴，避开顶部图例。
    h, w = new_image.shape[:2]
    ylo, yhi = round(h * 0.402), round(h * 0.610)
    region = new_image[ylo:yhi]
    gray = (region[:, :, 0] == region[:, :, 1]) & (region[:, :, 1] == region[:, :, 2])
    gray &= (region[:, :, 0] >= 100) & (region[:, :, 0] <= 200)
    columns = np.where(gray.sum(axis=0) >= (yhi - ylo) * 0.9)[0]
    left_columns = columns[(columns > w * 0.05) & (columns < w * 0.15)]
    right_columns = columns[columns > w * 0.95]
    assert len(left_columns) > 0 and len(right_columns) > 0, (ident, columns.tolist())
    left, right = int(left_columns[-1]), int(right_columns[0])
    row["detectedViAxesPixels"] = [left, right, ylo, yhi]
    row["rasterGapChecks"] = []
    t0 = datetime.fromisoformat(new_record["startUTC"].replace("Z", "+00:00")).timestamp()
    t1 = datetime.fromisoformat(new_record["endUTC"].replace("Z", "+00:00")).timestamp()
    for gap in gap_records:
        if gap["id"] != ident:
            continue
        gs = datetime.fromisoformat(gap["startUTC"].replace("Z", "+00:00")).timestamp()
        ge = datetime.fromisoformat(gap["endUTC"].replace("Z", "+00:00")).timestamp()
        x0 = int(np.ceil(left + (gs - t0) / (t1 - t0) * (right - left))) + 3
        x1 = int(np.floor(left + (ge - t0) / (t1 - t0) * (right - left))) - 3
        check = {"startUTC": gap["startUTC"], "endUTC": gap["endUTC"], "interiorWidthPixels": x1 - x0}
        if x1 > x0:
            block = region[:, x0:x1].astype(np.int16)
            colored = block.max(axis=2) - block.min(axis=2) > 80
            check["coloredPixels"] = int(colored.sum())
            assert check["coloredPixels"] == 0, (ident, check)
        else:
            check["status"] = "gap too narrow for a raster check; verified in native samples"
        row["rasterGapChecks"].append(check)
    rows.append(row)

result = {"complete": True, "count": len(rows), "recordChecksPassed": True,
          "sameSizeCount": sum(r["sameSize"] for r in rows),
          "sameSizeBPixelEqualCount": sum(r.get("B_region_pixel_equal", False) for r in rows),
          "sameSizeAEPixelEqualCount": sum(r.get("AE_region_pixel_equal", False) for r in rows),
          "sourceRecordsAndPNGUntouchedByThisVerifier": True, "images": rows}
(redraw / "redraw_verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({k: v for k, v in result.items() if k != "images"}, ensure_ascii=False))
for row in rows:
    if not row["sameSize"] or not row.get("B_region_pixel_equal", False) or not row.get("AE_region_pixel_equal", False):
        print(json.dumps(row, ensure_ascii=False))
