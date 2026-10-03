"""汇总此次CDF核查；用户检查表输出Recovery，来源及验证记录保存在Z盘。"""
import collections
import csv
import datetime as dt
import json
from pathlib import Path

audit = Path(r'Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930')
result = Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815')
baseline = json.loads((audit / 'local' / 'windows_audit.json').read_text('utf-8'))
new = json.loads((audit / 'local' / 'after_download' / 'windows_audit.json').read_text('utf-8'))
before_summary = json.loads((audit / 'local' / 'summary.json').read_text('utf-8'))
new_summary = json.loads((audit / 'local' / 'after_download' / 'summary.json').read_text('utf-8'))
download = json.loads((audit / 'cdaweb_download_verification.json').read_text('utf-8'))
cda = json.loads((audit / 'cdaweb' / 'summary.json').read_text('utf-8'))
sdc_path = json.loads((audit / 'sdc' / 'latest_audit.json').read_text('utf-8'))['report']
sdc = json.loads(Path(sdc_path).read_text('utf-8'))
assert len(baseline) == 162 and len(new) == 26 and download['complete'] and len(download['files']) == 49
assert before_summary['originalCountMismatches'] == 0
for section in (before_summary, new_summary):
    for key in ('fileReadErrors', 'windowReadErrors', 'databaseDirectVelocityMismatches', 'databaseDirectCountMismatches', 'coordinateValidCountLoss'):
        assert section[key] == 0, (key, section[key])
with (audit / 'cdaweb' / 'file_inventory.csv').open(encoding='utf-8-sig', newline='') as stream:
    inventory = [f for f in csv.DictReader(stream) if f['intersects_requested_interval'] == 'True']
assert len(inventory) == 204
for f in inventory:
    assert Path(f['local_path']).is_file() and Path(f['local_path']).stat().st_size == int(f['website_length']), f['filename']

updates = {row['id']: row for row in new}
categories = {'covered': '整段覆盖', 'partial_gap': '局部空档', 'entire_window_blank': '全段无当前公开L2速度数据'}
rows = []
current = []
for old in baseline:
    latest = updates.get(old['id'], old)
    current.append(latest)
    needs_update = old['id'] in updates
    if needs_update:
        assert old['fastGsmAnyFinite'] == 0 and latest['fastGsmAnyFinite'] > 0
    rows.append({
        '图号/PPT页码': int(old['id'][1:4]),
        '图片文件名': old['id'] + '_MMS1_B_Vi_AE.png',
        '窗口开始UTC': old['startUTC'], '窗口结束UTC': old['endUTC'],
        '原图Vi有效点数': old['oldBurstAnyFinite'] + old['oldFastAnyFinite'],
        '补齐后可读取Vi点数': latest['burstGsmAnyFinite'] + latest['fastGsmAnyFinite'],
        '补齐后数据覆盖': categories[latest['gapCategory']],
        '首个有效速度UTC': latest['firstValidUTC'], '最后有效速度UTC': latest['lastValidUTC'],
        '需重画图片并更新PPT': '是' if needs_update else '否',
        '当前图片状态': '仍为核查前版本',
    })
output = result / 'MMS1_Vi_completeness_check_20260930.csv'
with output.open('w', newline='', encoding='utf-8-sig') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
counts = collections.Counter(row['gapCategory'] for row in current)
assert counts == {'covered': 71, 'partial_gap': 42, 'entire_window_blank': 49}
summary = {
    'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
    'sdc_fast_files': sdc['products'][0]['latest_count'], 'sdc_burst_files': sdc['products'][1]['latest_count'],
    'cdaweb_files_requested_interval': cda['requested_latest_counts'],
    'initial_local_files': 155, 'missing_files_confirmed_and_downloaded': 49,
    'current_local_files_checked': len(inventory), 'download_bytes': sum(f['bytes_verified'] for f in download['files']),
    'new_valid_velocity_samples': new_summary['fastValidSamples'],
    'affected_windows': sorted(int(key[1:4]) for key in updates),
    'coverage_categories_after_download': dict(counts),
    'new_coverage_hours': new_summary['coverageHours'],
    'current_total_valid_velocity_samples': sum(w['fastGsmAnyFinite'] + w['burstGsmAnyFinite'] for w in current),
    'old_file_read_errors': before_summary['fileReadErrors'], 'new_file_read_errors': new_summary['fileReadErrors'],
    'coordinate_valid_sample_loss': before_summary['coordinateValidCountLoss'] + new_summary['coordinateValidCountLoss'],
    'images_redrawn': False, 'ppt_updated': False, 'user_check_table': str(output),
    'conclusion': 'CDAWeb contains 49 files absent from the original local collection. They provide finite velocity samples for 26 previously blank plots. All 49 CDFs are now archived and verified. Existing PNG/PPT files still show the old data coverage.',
    'availability_caveat': 'Query results establish present public availability. File LastModified does not establish when a file first became available.',
}
(audit / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
