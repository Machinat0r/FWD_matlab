"""核验输入未被另行修改后，仅覆盖指定26张PNG和用户指定PPT。
原文件备份置于TEMP，交付和完整性记录置于Z盘derived。
"""
import datetime
import hashlib
import json
import os
import shutil
from pathlib import Path

work = Path(os.environ['TEMP']) / 'MMS1_BViAE_refresh_20260930'
build = work / 'build'
manifest = json.loads((build / 'replacement_manifest.json').read_text('utf-8'))
proof = json.loads((build / 'package_preservation.json').read_text('utf-8'))
result = Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815')
audit = Path(r'Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930\redraw')
target = result / 'MMS1_B_Vi_AE_20260720_0815_75percent.pptx'
candidate = work / 'output' / target.name
backup = work / 'backup'
backup.mkdir(exist_ok=True)
png_backup = backup / 'B_Vi_AE'
png_backup.mkdir(exist_ok=True)

def sha(file):
    return hashlib.sha256(Path(file).read_bytes()).hexdigest()

assert sha(target) == manifest['source_sha256'], '原PPT已被用户另行修改，停止覆盖'
assert sha(candidate) == proof['final_sha256'], '候选稿验证后发生变化'
assert proof['changed_media_count'] == 26 and proof['other_parts_identical'] == 801
redraw = json.loads((audit / 'redraw_verification.json').read_text('utf-8'))
assert redraw['complete'] and redraw['count'] == 26 and redraw['recordChecksPassed']
assert redraw['sameSizeBPixelEqualCount'] == 26 and redraw['sameSizeAEPixelEqualCount'] == 26
original_hashes = {p.name: sha(p) for p in (result / 'B_Vi_AE').glob('*.png')}
assert len(original_hashes) == 162
for item in manifest['matches']:
    old = Path(item['old_png'])
    assert old.resolve().parent == (result / 'B_Vi_AE').resolve()
    assert sha(old) == item['old_sha256'], '原图已被另行修改，停止覆盖'
    shutil.copy2(old, png_backup / old.name)
shutil.copy2(target, backup / target.name)

changed = []
try:
    # PPT只保存已通过PowerPoint渲染、结构和内容差异核验的原样字节。
    shutil.copyfile(candidate, target)
    changed.append(target)
    for item in manifest['matches']:
        shutil.copyfile(item['new_png'], item['old_png'])
        changed.append(Path(item['old_png']))
    after_hashes = {p.name: sha(p) for p in (result / 'B_Vi_AE').glob('*.png')}
    expected = {Path(m['old_png']).name for m in manifest['matches']}
    differences = {name for name in original_hashes if original_hashes[name] != after_hashes[name]}
    assert differences == expected
    assert sha(target) == proof['final_sha256']
    for item in manifest['matches']:
        assert sha(item['old_png']) == sha(item['new_png'])
except Exception:
    for path in changed:
        origin = backup / path.name if path == target else png_backup / path.name
        shutil.copyfile(origin, path)
    raise

receipt = {
    'complete': True,
    'published_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'ppt': str(target), 'ppt_sha256': sha(target),
    'original_ppt_sha256': manifest['source_sha256'],
    'slides': 162, 'updated_png_count': 26, 'unchanged_png_count': 136,
    'changed_slide_numbers': [m['slide'] for m in manifest['matches']],
    'other_ppt_parts_identical': 801,
    'images': [{'window_id': m['window_id'], 'slide': m['slide'], 'path': m['old_png'],
                'old_sha256': m['old_sha256'], 'new_sha256': sha(m['old_png'])} for m in manifest['matches']],
    'backup_path': str(backup),
    'unchanged_png_sha256': {k: v for k, v in original_hashes.items() if k not in expected},
}
for name in ['validation.json', 'package_preservation.json', 'replacement_manifest.json']:
    shutil.copy2(build / name, audit / name)
(audit / 'delivery_receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in receipt.items() if k not in ['images','unchanged_png_sha256']}, ensure_ascii=False))
