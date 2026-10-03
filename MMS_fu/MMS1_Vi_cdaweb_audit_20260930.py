"""仅核查 CDAWeb 原始文件清单及本地文件；不下载或处理 CDF。"""
import csv
import datetime as dt
import json
from pathlib import Path
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(r'Z:\SPART-WORK\Data\MMS')
OUT = ROOT / 'derived/MMS1_tail_20260720_0815/Vi_completeness_20260930/cdaweb'
OUT.mkdir(parents=True, exist_ok=True)
BASE = 'https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/'
QUERIES = [
    ('fast_full', 'FAST', '20260719T000000Z', '20260817T000000Z'),
    ('brst_full', 'BRST', '20260719T000000Z', '20260817T000000Z'),
    ('fast_late_check', 'FAST', '20260809T000000Z', '20260816T000000Z'),
    ('fast_positive', 'FAST', '20260725T000000Z', '20260726T000000Z'),
    ('brst_positive', 'BRST', '20250725T000000Z', '20250726T000000Z'),
]

def query(spec):
    name, mode, start, end = spec
    url = BASE + f'MMS1_FPI_{mode}_L2_DIS-MOMS/orig_data/{start},{end}'
    stamp = dt.datetime.now(dt.timezone.utc).isoformat()
    request = urllib.request.Request(url, headers={'Accept': 'application/json'})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            raw = response.read()
            meta = {'query_utc': stamp, 'url': url, 'status': response.status,
                    'headers': dict(response.headers), 'bytes': len(raw)}
        (OUT / f'{name}.response.json').write_bytes(raw)
        obj = json.loads(raw)
        meta['response_keys'] = list(obj)
        rows = obj.get('FileDescription', [])
        meta['file_count'] = len(rows)
        meta['non_file_fields'] = {k: v for k, v in obj.items() if k != 'FileDescription'}
        meta['pagination_headers'] = {k:v for k,v in meta['headers'].items() if k.lower() in ('link','content-range','x-total-count')}
    except Exception as exc:
        meta = {'query_utc': stamp, 'url': url, 'error': repr(exc)}
        rows = []
    (OUT / f'{name}.metadata.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    return name, meta, rows

if '--map-only' in sys.argv:
    audit = json.loads((OUT/'summary.json').read_text())
    windows = json.loads((OUT.parent.parent/'windows.json').read_text())
    missing = audit['missing_requested_files']
    rows = []
    for window in windows:
        hits = []
        for file in missing:
            a = dt.datetime.fromisoformat(file['start_utc'].replace('Z','+00:00')).timestamp()
            b = dt.datetime.fromisoformat(file['end_utc'].replace('Z','+00:00')).timestamp()
            if b >= window['startEpoch'] and a < window['endEpoch']:
                hits.append(file)
        if hits:
            plotted = json.loads((OUT.parent.parent/f"{window['id']}_overview.json").read_text())
            rows.append({'window_number':window['number'],'window_id':window['id'],
                'window_start_utc':window['startUTC'],'window_end_utc':window['endUTC'],
                'existing_Vi_native_counts':plotted['nativeCounts'][plotted['nativeCountNames'].index('Vi')],
                'missing_file_count':len(hits),'missing_files':[x['filename'] for x in hits],
                'file_metadata_first_start':min(x['start_utc'] for x in hits),
                'file_metadata_last_end':max(x['end_utc'] for x in hits)})
    result = {'missing_file_count':len(missing),'window_count':len(rows),
        'window_numbers':[x['window_number'] for x in rows],
        'note':'The overlaps use CDAWeb FileDescription StartTime/EndTime. They identify windows to recheck; actual finite bulk-velocity records require reading the CDF.',
        'last_modified_note':'CDAWeb FileDescription LastModified records a modification timestamp. It does not establish the first public availability time.',
        'windows':rows}
    (OUT/'missing_files_window_map.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    with (OUT/'missing_files_window_map.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result,indent=2))
    sys.exit(0)

if '--checks-only' in sys.argv:
    checks = [
        ('fast_check_week1', 'FAST', '20260719T000000Z', '20260726T000000Z'),
        ('fast_check_week2', 'FAST', '20260726T000000Z', '20260802T000000Z'),
        ('fast_check_week3', 'FAST', '20260802T000000Z', '20260809T000000Z'),
        ('fast_check_week4', 'FAST', '20260809T000000Z', '20260817T000000Z'),
        ('brst_positive_20151016', 'BRST', '20151016T000000Z', '20151017T000000Z'),
    ]
    with ThreadPoolExecutor(max_workers=3) as pool:
        checked = dict((name, (meta, rows)) for name, meta, rows in pool.map(query, checks))
    full_names = {x['Name'] for x in json.loads((OUT/'fast_full.response.json').read_text())['FileDescription']}
    split_names = {x['Name'] for name, (meta, rows) in checked.items() if name.startswith('fast_check') for x in rows}
    result = {'checks':{name:meta for name,(meta,rows) in checked.items()},
        'full_query_unique_files':len(full_names),'weekly_query_unique_files':len(split_names),
        'weekly_missing_from_full_query':sorted(split_names-full_names),
        'full_missing_from_weekly_queries':sorted(full_names-split_names),
        'conclusion':'The full-range fast response agrees exactly with the union of four weekly responses.' if full_names==split_names else 'File lists differ.'}
    (OUT/'independent_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    sys.exit(0)

with ThreadPoolExecutor(max_workers=3) as pool:
    results = dict((name, (meta, rows)) for name, meta, rows in pool.map(query, QUERIES))

inventory = []
daily = []
for mode in ('fast', 'brst'):
    rows = results[f'{mode}_full'][1]
    latest = {}
    for item in rows:
        filename = item['Name'].rsplit('/', 1)[-1]
        match = re.search(r'_(\d{14})_v(\d+)\.(\d+)\.(\d+)\.cdf$', filename)
        if not match:
            raise ValueError(filename)
        start, *version = match.groups()
        version = tuple(map(int, version))
        if start not in latest or version > latest[start][0]:
            latest[start] = (version, item)
    for start, (version, item) in sorted(latest.items()):
        filename = item['Name'].rsplit('/', 1)[-1]
        sub = ROOT / f'mms1/fpi/{mode}/l2/dis-moms/{start[:4]}/{start[4:6]}'
        if mode == 'brst':
            sub = sub / start[6:8]
        local = sub / filename
        in_requested = item['EndTime'] >= '2026-07-20T00:00:00.000Z' and item['StartTime'] < '2026-08-16T00:00:00.000Z'
        record = {'mode':mode,'date':f'{start[:4]}-{start[4:6]}-{start[6:8]}',
            'filename': filename,'source_url':item['Name'],'start_utc':item['StartTime'],
            'end_utc':item['EndTime'],'website_length':item.get('Length'),
            'website_last_modified':item.get('LastModified'), 'local_path':str(local),
            'local_exists':local.is_file(), 'local_length':local.stat().st_size if local.is_file() else None,
            'intersects_requested_interval':in_requested}
        record['same_length'] = record['local_length'] == record['website_length']
        inventory.append(record)
    for offset in range(27):
        day = dt.date(2026,7,20) + dt.timedelta(days=offset)
        current = [x for x in inventory if x['mode']==mode and x['date']==str(day)]
        daily.append({'date':str(day),'mode':mode,'website_latest_file_count':len(current),
            'local_missing_count':sum(not x['local_exists'] for x in current),
            'length_mismatch_count':sum(x['local_exists'] and not x['same_length'] for x in current),
            'first_start_utc':min((x['start_utc'] for x in current),default=''),
            'last_end_utc':max((x['end_utc'] for x in current),default='')})

for filename, rows in [('file_inventory.csv',inventory),('daily_counts.csv',daily)]:
    with (OUT / filename).open('w',newline='',encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else ['empty'])
        writer.writeheader()
        writer.writerows(rows)
summary = {'query_results':{name:meta for name,(meta,rows) in results.items()},
    'requested_interval_utc':['2026-07-20T00:00:00Z','2026-08-16T00:00:00Z'],
    'requested_latest_counts':{mode:sum(x['mode']==mode and x['intersects_requested_interval'] for x in inventory) for mode in ('fast','brst')},
    'missing_requested_files':[x for x in inventory if x['intersects_requested_interval'] and not x['local_exists']],
    'wrong_length_requested_files':[x for x in inventory if x['intersects_requested_interval'] and x['local_exists'] and not x['same_length']],
    'daily_counts':daily}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({'query_results':{k: {'status':v[0].get('status'),'file_count':v[0].get('file_count'),'error':v[0].get('error')} for k,v in results.items()},
    'requested_latest_counts':summary['requested_latest_counts'],
    'missing_requested_files':summary['missing_requested_files'],
    'wrong_length_requested_files':summary['wrong_length_requested_files'],
    'daily_counts':daily,'output':str(OUT)},indent=2))
