"""五事件四星电流图：复用既有SDC查询/下载程序，原CDF归档到Z盘。"""
import concurrent.futures as cf
import datetime as dt
import importlib.util
import json
from pathlib import Path

code = Path(__file__).parent
spec = importlib.util.spec_from_file_location('existing_mms_download', code / 'MMS_event_overview_20260923_download.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
d.AUDIT = d.ROOT / 'derived' / 'MMS_current_5events_20261001'
d.AUDIT.mkdir(parents=True, exist_ok=True)
windows = json.loads((d.ROOT / 'derived' / 'MMS1_tail_20260720_0815' / 'windows.json').read_text('utf-8'))
chosen = [windows[i - 1] for i in (3, 25, 50, 51, 82)]
d.save(d.AUDIT / 'windows.json', chosen)
days = ('2026-07-20', '2026-07-23', '2026-07-24', '2026-07-28', '2026-08-02')
products = (('fgm', 'srvy', ''), ('fgm', 'brst', ''), ('mec', 'srvy', 'epht89d'))
latest = {}
queries = []
for day in days:
    end = (dt.date.fromisoformat(day) + dt.timedelta(days=1)).isoformat()
    for product in products:
        files = d.query(product, day, end, (1, 2, 3, 4))
        queries.append({'day': day, 'product': product, 'count': len(files)})
        for f in files:
            stem, ver = f['file_name'].rsplit('_v', 1)
            version = tuple(map(int, ver.removesuffix('.cdf').split('.')))
            if stem not in latest or version > latest[stem][0]:
                latest[stem] = (version, f)
        d.save(d.AUDIT / 'inventory_progress.json', {'queries': queries})
        print('QUERY', day, product, len(files), flush=True)
selected = []
for _, f in latest.values():
    # 日文件覆盖当日；burst以官网文件起止时间为候选依据。
    if '_brst_' in f['file_name']:
        a = d.timestamp(f['timetag']).timestamp()
        stop = f.get('end_time')
        b = d.timestamp(stop).timestamp() if stop else a + 3600
        if not any(a <= w['endEpoch'] and b >= w['startEpoch'] for w in chosen):
            continue
    selected.append({**f, 'file_size': int(f['file_size']), 'path': str(d.local_path(f['file_name']))})
selected.sort(key=lambda f: f['file_name'])
d.save(d.AUDIT / 'manifest.json', {'windows': chosen, 'queries': queries, 'files': selected})
for folder in sorted({str(Path(f['path']).parent) for f in selected}):
    Path(folder).mkdir(parents=True, exist_ok=True)
print('FILES', len(selected), 'MISSING', sum(not Path(f['path']).is_file() for f in selected), flush=True)
results = []
with cf.ThreadPoolExecutor(max_workers=4) as pool:
    jobs = {pool.submit(d.download, f): f for f in selected}
    for future in cf.as_completed(jobs):
        r = future.result()
        results.append(r)
        d.save(d.AUDIT / 'download_progress.json', {'done': len(results), 'total': len(selected), 'results': results})
        print('CDF', len(results), '/', len(selected), r['status'], r['name'], flush=True)
d.save(d.AUDIT / 'download_complete.json', {'complete': all(r['status'] != 'failed' for r in results), 'results': results})
if any(r['status'] == 'failed' for r in results):
    raise RuntimeError('下载失败，请查看download_complete.json')
print('DOWNLOAD_FINISHED', flush=True)
