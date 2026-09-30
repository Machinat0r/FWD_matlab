"""保存NASA SSC返回的原始GSE星历响应；科学读取与绘图由MATLAB完成。"""
import json,datetime as dt,urllib.request
from pathlib import Path
import MMS_event_overview_20260923_download as d
root=d.ROOT/'ancillary'/'sscweb'/'2026'
for p in reversed([root,*root.parents]):
    if not p.exists():p.mkdir()
events=json.loads((d.CODE/'MMS_event_overview_20260924_events.json').read_text('utf-8'))
records=[]
for e in events[14:]:
    t=d.timestamp(e['start'])+(d.timestamp(e['end'])-d.timestamp(e['start']))/2
    span=','.join((t+dt.timedelta(minutes=k)).strftime('%Y%m%dT%H%M%SZ') for k in (-4,4))
    url='https://sscweb.gsfc.nasa.gov/WS/sscr/2/locations/mms1,mms2,mms3,mms4/'+span+'/gse/'
    dest=root/(e['id']+'_gse.json')
    if not dest.exists():
        with urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json'}),timeout=60) as r:b=r.read()
        assert b'"SUCCESS"' in b and b'"mms4"' in b
        dest.write_bytes(b)
    records.append(dict(event=e['id'],time=t.isoformat(),path=str(dest),url=url,source='NASA SSCWeb',coordinateSystem='GSE',units='km'))
    print('SSC',e['id'],dest.stat().st_size,flush=True)
d.save(d.ROOT/'derived'/'event_orbits_20260924'/'ssc_manifest.json',records)
