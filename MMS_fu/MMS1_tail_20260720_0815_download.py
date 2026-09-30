"""复用现有下载器，原始CDF按产品层级归档；不做科学数据转换。"""
import importlib.util,json,sys,datetime as dt
from pathlib import Path
import concurrent.futures as cf
code=Path(r'C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
spec=importlib.util.spec_from_file_location('existing_mms_download',code/'MMS_event_overview_20260923_download.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
d.AUDIT=d.ROOT/'derived'/'MMS1_tail_20260720_0815'
d.AUDIT.mkdir(parents=True,exist_ok=True)
phase=sys.argv[1] if len(sys.argv)>1 else 'mec'
products=[('mec','srvy','epht89d')] if phase=='mec' else [('fgm','srvy',''),('fgm','brst',''),('fpi','fast','dis-moms'),('fpi','brst','dis-moms')]
if phase!='mec':
 windows=json.loads((d.AUDIT/'windows.json').read_text('utf-8'))
selected={};coverage=[]
for product in products:
 files=d.query(product,'2026-07-19','2026-08-17',(1,))
 latest={}
 for f in files:
  stem,version=f['file_name'].rsplit('_v',1);ver=tuple(map(int,version.removesuffix('.cdf').split('.')))
  if stem not in latest or ver>latest[stem][0]:latest[stem]=(ver,f)
 rows=sorted((f for _,f in latest.values()),key=lambda f:f['timetag'])
 chosen=[]
 for i,f in enumerate(rows):
  t=d.timestamp(f['timetag']).timestamp()
  duration=3600 if product[1]=='brst' else (7200 if product[0]=='fpi' else 86400)
  end=min(t+duration,d.timestamp(rows[i+1]['timetag']).timestamp()) if i+1<len(rows) else t+duration
  if phase=='mec' or any(t<w['endEpoch'] and end>w['startEpoch'] for w in windows):
   name=f['file_name'];selected[name]={**f,'path':str(d.local_path(name))};chosen.append(name)
 coverage.append({'product':product,'available':len(files),'selected':chosen})
 print(product,'available',len(files),'selected',len(chosen),flush=True)
d.save(d.AUDIT/(phase+'_manifest.json'),{'files':list(selected.values()),'coverage':coverage,'note':'Candidate overlap by file timetags. Actual samples checked with IRFU in MATLAB.'})
print('CANDIDATES',len(selected),'GB',sum(int(f['file_size']) for f in selected.values())/1e9,flush=True)
for directory in sorted({str(Path(f['path']).parent) for f in selected.values()}):Path(directory).mkdir(parents=True,exist_ok=True)
results=[]
with cf.ThreadPoolExecutor(max_workers=4) as pool:
 pending={pool.submit(d.download,f):f for f in selected.values()}
 for future in cf.as_completed(pending):
  result=future.result();results.append(result)
  d.save(d.AUDIT/(phase+'_progress.json'),{'done':len(results),'total':len(selected),'results':results})
  print('CDF',len(results),'/',len(selected),result['status'],result['name'],flush=True)
d.save(d.AUDIT/(phase+'_complete.json'),{'complete':all(r['status']!='failed' for r in results),'results':results})
assert all(r['status']!='failed' for r in results)
print('DOWNLOAD_FINISHED',flush=True)
