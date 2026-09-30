"""下载MEC暂缺日期的官方DEFEPH原文件，供IRFU读取。"""
import json,urllib.parse,concurrent.futures as cf
from pathlib import Path
import MMS_event_overview_20260923_download as d
a=d.ROOT/'derived'/'event_orbits_20260924'
query=json.loads((a/'ancillary_query.json').read_text('utf-8'))
days={'2026-08-27','2026-08-28','2026-08-29','2026-09-15','2026-09-16'}
chosen={}
for f in query['files']:
    if f['start_date'][:10] not in days:continue
    stem,version=f['file_name'].rsplit('.V',1)
    if stem not in chosen or int(version)>int(chosen[stem]['file_name'].rsplit('.V',1)[1]):chosen[stem]=f
files=[]
for f in chosen.values():
    f=dict(f);sc=f['file_name'].split('_')[0].lower()
    p=d.ROOT/sc/'ancillary'/'defeph'/'2026'/f['file_name']
    for q in reversed(p.parents):
        if not q.exists():q.mkdir()
    f['path']=str(p);files.append(f)
def fetch(f):
    p=Path(f['path'])
    if not p.exists() or p.stat().st_size!=f['file_size']:
        b=d.open_url(d.BASE+'download/ancillary?'+urllib.parse.urlencode({'file':f['file_name']})).read()
        assert len(b)==f['file_size'],f['file_name']
        assert b.startswith(b'Definitive Orbit Ephemeris')
        assert b'ReferenceFrame = Mean of J2000' in b[:2000]
        p.write_bytes(b)
    return dict(name=f['file_name'],bytes=p.stat().st_size)
with cf.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(fetch,files))
d.save(a/'defeph_manifest.json',dict(files=files,results=results,source='MMS SDC ancillary DEFEPH',referenceFrame='Mean of J2000'))
print('DEFEPH_COMPLETE',len(files),sum(r['bytes'] for r in results),flush=True)
