"""复用现有SDC接口，对照当前FPI离子矩文件清单；不下载或改变原CDF。"""
import collections
import datetime as dt
import importlib.util
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
from pathlib import Path

code = Path(__file__).parent
spec = importlib.util.spec_from_file_location('existing_mms_download', code / 'MMS_event_overview_20260923_download.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
old_root = d.ROOT / 'derived' / 'MMS1_tail_20260720_0815'
audit_root = old_root / 'Vi_completeness_20260930' / 'sdc'
audit_root.mkdir(parents=True, exist_ok=True)
# 每次运行建立新的查询目录，避免沿用历史缓存。
if len(sys.argv) > 1:
    d.AUDIT = Path(sys.argv[1])  # 仅用于续跑本次已保存、仍需补完的查询。
    assert d.AUDIT.is_dir() and d.AUDIT.parent == audit_root
else:
    d.AUDIT = audit_root / dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    d.AUDIT.mkdir(parents=True, exist_ok=False)
windows = json.loads((old_root / 'windows.json').read_text('utf-8'))
old_manifest = json.loads((old_root / 'science_manifest.json').read_text('utf-8'))
old_selected = {f['file_name'] for f in old_manifest['files'] if '_fpi_' in f['file_name']}
report = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'query_directory': str(d.AUDIT), 'products': []}
all_latest = []
for mode in ('fast', 'brst'):
    files = d.query(('fpi', mode, 'dis-moms'), '2026-07-19', '2026-08-17', (1,))
    latest = {}
    for f in files:
        stem, version = f['file_name'].rsplit('_v', 1)
        version = tuple(map(int, version.removesuffix('.cdf').split('.')))
        if stem not in latest or version > latest[stem][0]:
            latest[stem] = (version, f)
    rows = sorted((f for _, f in latest.values()), key=lambda f: f['file_name'])
    params = dict(sc_id='mms1', instrument_id='fpi', data_rate_mode=mode, data_level='l2', descriptor='dis-moms', start_date='2026-07-19', end_date='2026-08-17')
    url = d.BASE + 'file_names/science?' + urllib.parse.urlencode(params)
    raw_names = None
    names_error = None
    for attempt in range(3):
        try:
            with d.open_url(url) as response:
                raw_names = response.read().decode('utf-8')
            break
        except (http.client.RemoteDisconnected, urllib.error.URLError, TimeoutError) as exc:
            names_error = repr(exc)
            if attempt < 2:
                time.sleep(attempt + 1)
    names = sorted(set(re.findall(r'mms1_fpi_[^,/\\\s]+\.cdf', raw_names))) if raw_names is not None else None
    if raw_names is not None:
        (d.AUDIT / f'{mode}_file_names.txt').write_text(raw_names, 'utf-8')
    metadata_names = {f['file_name'] for f in files}
    if names is not None:
        assert set(names) == metadata_names, (mode, 'file_info/file_names disagree', len(names), len(metadata_names))
    inventory = []
    for f in rows:
        local = d.local_path(f['file_name'])
        start = d.timestamp(f['timetag']).timestamp()
        # 仅供定位候选；最终覆盖由MATLAB读取CDF样本判定。
        estimated_end = start + (7200 if mode == 'fast' else 3600)
        candidate_windows = [w['id'] for w in windows if start < w['endEpoch'] and estimated_end > w['startEpoch']]
        exists = local.is_file()
        actual_size = local.stat().st_size if exists else None
        entry = {**f, 'path': str(local), 'exists': exists, 'local_size': actual_size,
                 'size_matches': actual_size == int(f['file_size']), 'in_old_manifest': f['file_name'] in old_selected,
                 'candidate_windows_by_filename_only': candidate_windows}
        inventory.append(entry)
        all_latest.append(entry)
    old_query_path = old_root / 'queries' / f'fpi_{mode}_dis-moms_1_2026-07-19_2026-08-17.json'
    old_query = json.loads(old_query_path.read_text('utf-8')) if old_query_path.is_file() else {}
    old_names = {f['file_name'] for f in old_query.get('files', [])}
    summary = {'mode': mode, 'all_versions_count': len(files), 'latest_count': len(rows),
               'file_names_match_file_info': set(names) == metadata_names if names is not None else None,
               'file_names_error': names_error if names is None else None, 'names_url': url,
               'by_day': dict(sorted(collections.Counter(f['timetag'][:10] for f in rows).items())),
               'first_file': rows[0]['file_name'] if rows else None,
               'last_file': rows[-1]['file_name'] if rows else None,
               'missing_or_wrong_size': [f['file_name'] for f in inventory if not f['size_matches']],
               'missing_or_wrong_size_candidate': [f['file_name'] for f in inventory if not f['size_matches'] and f['candidate_windows_by_filename_only']],
               'new_since_old_query': sorted(metadata_names - old_names),
               'old_query_utc': old_query.get('_queried_utc'), 'inventory': inventory}
    report['products'].append(summary)
    d.save(d.AUDIT / 'sdc_audit.json', report)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('inventory', 'by_day')}, ensure_ascii=False), flush=True)
report['completed_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
d.save(d.AUDIT / 'sdc_audit.json', report)
d.save(audit_root / 'latest_audit.json', {'report': str(d.AUDIT / 'sdc_audit.json'), 'completed_utc': report['completed_utc']})
print('AUDIT_SAVED', d.AUDIT, flush=True)
