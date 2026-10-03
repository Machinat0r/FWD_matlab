"""2019-08-05 16:24-16:25电流：复用本地FGM，仅补缺失MEC及OMNI原CDF。"""
from pathlib import Path
import importlib.util,json,datetime as dt,os
code=Path(__file__).parent
spec=importlib.util.spec_from_file_location('existing_mms_download',code/'MMS_event_overview_20260923_download.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
d.AUDIT=d.ROOT/'derived/MMS_current_20190805_1624_1625_native_20261002'
d.AUDIT.mkdir(parents=True,exist_ok=True)
window=[{'id':'MMS_20190805_1624_1625','startUTC':'2019-08-05T16:24:00Z','endUTC':'2019-08-05T16:25:00Z'}]
d.save(d.AUDIT/'windows.json',window)
missing=[];reused=[]
for ic in range(1,5):
    folder=d.ROOT/f'mms{ic}/mec/srvy/l2/epht89d/2019/08'
    files=[p for p in folder.glob(f'mms{ic}_mec_srvy_l2_epht89d_20190805_v*.cdf') if p.stat().st_size>10000]
    if files: reused.extend(str(p) for p in files)
    else: missing.append(ic)
results=[];selected=[]
if missing:
    files=d.query(('mec','srvy','epht89d'),'2019-08-05','2019-08-06',tuple(missing))
    for ic in missing:
        candidates=[f for f in files if f['file_name'].startswith(f'mms{ic}_mec_srvy_l2_epht89d_20190805_')]
        assert candidates,f'Missing official MEC for MMS{ic}'
        f=max(candidates,key=lambda f:tuple(map(int,f['file_name'].rsplit('_v',1)[1].removesuffix('.cdf').split('.'))))
        item={**f,'file_size':int(f['file_size']),'path':str(d.local_path(f['file_name']))}
        Path(item['path']).parent.mkdir(parents=True,exist_ok=True)
        selected.append(item);r=d.download(item);results.append(r)
        print('MEC',ic,r['status'],r['name'],flush=True)
        assert r['status']!='failed',r
omni=d.ROOT/'ancillary/omni/hro_1min/2019/omni_hro_1min_20190801_v01.cdf'
omni.parent.mkdir(parents=True,exist_ok=True)
url='https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/2019/'+omni.name
if not omni.exists():
    tmp=omni.with_suffix('.cdf.part')
    with d.open_url(url) as response,tmp.open('wb') as out:
        expected=response.headers.get('Content-Length')
        while chunk:=response.read(1024*1024):out.write(chunk)
    if expected:assert tmp.stat().st_size==int(expected)
    assert tmp.read_bytes()[:4] in (bytes.fromhex('cdf30001'),bytes.fromhex('cdf26002'),bytes.fromhex('0000ffff'))
    os.replace(tmp,omni);omni_status='downloaded'
else:omni_status='existing'
d.save(d.AUDIT/'download_complete.json',{'complete':True,'MECResults':results,'reusedMEC':reused,'OMNI':{'path':str(omni),'url':url,'status':omni_status}})
d.save(d.AUDIT/'manifest.json',{'windows':window,'MECFiles':selected,'reusedMEC':reused,'OMNI':str(omni),'FGM':'reuse existing four-spacecraft local burst CDF; actual time coverage verified by MATLAB'})
print('PREPARED',len(results),'MEC files; OMNI',omni_status,flush=True)
