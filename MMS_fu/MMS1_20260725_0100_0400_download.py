"""调用已有下载器查询本次MMS1原始CDF；科学读取和绘图继续使用MATLAB。"""
import importlib.util,datetime as dt,json
from pathlib import Path
import concurrent.futures as cf
code=Path(r'C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
spec=importlib.util.spec_from_file_location('existing_mms_download',code/'MMS_event_overview_20260923_download.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
d.AUDIT=d.ROOT/'derived'/'MMS1_20260725_0100_0400'
d.AUDIT.mkdir(parents=True,exist_ok=True)
a=d.timestamp('2026-07-25T01:00:00Z');b=d.timestamp('2026-07-25T04:00:00Z')
selected={};coverage=[]
for product in d.PRODUCTS:
 files=d.query(product,'2026-07-25','2026-07-26',(1,))
 latest={}
 for f in files:
  stem,version=f['file_name'].rsplit('_v',1)
  ver=tuple(map(int,version.removesuffix('.cdf').split('.')))
  if stem not in latest or ver>latest[stem][0]:latest[stem]=(ver,f)
 rows=sorted((f for _,f in latest.values()),key=lambda f:f['timetag'])
 duration=dt.timedelta(hours=1 if product[1]=='brst' else (2 if product[0]=='fpi' else 24))
 chosen=[]
 for i,f in enumerate(rows):
  t=d.timestamp(f['timetag'])
  end=min(t+duration,d.timestamp(rows[i+1]['timetag'])) if i+1<len(rows) else t+duration
  if t<b and end>a:
   name=f['file_name'];selected[name]={**f,'path':str(d.local_path(name))};chosen.append(name)
 coverage.append({'product':product,'filesInDay':len(files),'candidates':chosen})
 print(product,'available',len(files),'selected',len(chosen),flush=True)
manifest={'startUTC':'2026-07-25T01:00:00Z','endUTC':'2026-07-25T04:00:00Z','spacecraft':1,'files':list(selected.values()),'coverage':coverage,'note':'Candidate file selection by timetag; actual observation intervals checked in MATLAB.'}
d.save(d.AUDIT/'manifest.json',manifest)
print('CANDIDATES',len(selected),'BYTES',sum(int(f['file_size']) for f in selected.values()),flush=True)
for directory in sorted({str(Path(f['path']).parent) for f in selected.values()}):Path(directory).mkdir(parents=True,exist_ok=True)
results=[]
with cf.ThreadPoolExecutor(max_workers=4) as pool:
 pending={pool.submit(d.download,f):f for f in selected.values()}
 for future in cf.as_completed(pending):
  result=future.result();results.append(result)
  d.save(d.AUDIT/'download_progress.json',{'done':len(results),'total':len(selected),'results':results})
  print('CDF',len(results),'/',len(selected),result['status'],result['name'],flush=True)
d.save(d.AUDIT/'download_complete.json',{'complete':all(r['status']!='failed' for r in results),'results':results})
assert all(r['status']!='failed' for r in results)
print('DOWNLOAD_FINISHED',flush=True)
