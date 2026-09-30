"""复用现有下载函数，仅补充20个事件的原始MEC星历CDF。"""
import concurrent.futures as cf
import datetime as dt
import json
from pathlib import Path
import MMS_event_overview_20260923_download as d
d.AUDIT=d.ROOT/'derived'/'event_orbits_20260924'
d.AUDIT.mkdir(parents=True,exist_ok=True)
events=json.loads((d.CODE/'MMS_event_overview_20260924_events.json').read_text('utf-8'))
days=sorted({e['start'][:10] for e in events})
selected={}
coverage=[]
for day in days:
    end=(dt.date.fromisoformat(day)+dt.timedelta(days=1)).isoformat()
    candidates=d.query(('mec','srvy','epht89d'),day,end)
    found={f['file_name'].split('_')[0] for f in candidates}
    if len(found)<4:
        missing=tuple(sc for sc in range(1,5) if 'mms'+str(sc) not in found)
        candidates+=d.query(('mec','srvy','epht89q'),day,end,missing)
    latest={}
    for f in candidates:
        name=f['file_name'];stem,ver=name.rsplit('_v',1)
        if name.split('_')[-2][:8]!=day.replace('-',''):continue
        version=tuple(map(int,ver.removesuffix('.cdf').split('.')))
        if stem not in latest or version>latest[stem][0]:latest[stem]=(version,f)
    for _,f in latest.values():
        f=dict(f);f['path']=str(d.local_path(f['file_name']))
        selected[f['file_name']]=f
    coverage.append(dict(day=day,files=[f['file_name'] for _,f in latest.values()]))
    print('MEC',day,len(latest),flush=True)
d.save(d.AUDIT/'mec_manifest.json',dict(coverage=coverage,files=list(selected.values())))
for f in selected.values():
    for folder in reversed(Path(f['path']).parents):
        if not folder.exists():folder.mkdir()
with cf.ThreadPoolExecutor(max_workers=3) as pool:
    results=list(pool.map(d.download,selected.values()))
d.save(d.AUDIT/'mec_download_results.json',results)
print('DOWNLOAD_COMPLETE',len(results),'failed',sum(r['status']=='failed' for r in results),'bytes',sum(r.get('bytes',0) for r in results),flush=True)
