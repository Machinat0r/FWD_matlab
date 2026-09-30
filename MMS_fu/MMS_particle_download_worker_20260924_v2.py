"""按卫星复用已有download函数。四颗卫星路径互不重叠；每个进程4个下载线程。"""
import importlib.util,json,sys,concurrent.futures as cf,datetime as dt
from pathlib import Path
code=Path(__file__).parent
spec=importlib.util.spec_from_file_location('existing_mms_download',code/'MMS_event_overview_20260923_download.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
audit=Path(r'Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2\particles')
sc=int(sys.argv[1]);items=[f for f in json.loads((audit/'manifest.json').read_text(encoding='utf8'))['files'] if f['file_name'].startswith(f'mms{sc}_')]
results=[]
with cf.ThreadPoolExecutor(max_workers=4) as pool:
 pending={pool.submit(m.download,item):item for item in items}
 for future in cf.as_completed(pending):
  r=future.result();r['completedUTC']=dt.datetime.now(dt.timezone.utc).isoformat();results.append(r)
  m.save(audit/f'download_mms{sc}.json',{'done':len(results),'total':len(items),'failed':sum(r['status']=='failed' for r in results),'results':results})
  print('CDF',sc,len(results),'/',len(items),r['status'],r['name'],flush=True)
m.save(audit/f'download_mms{sc}_complete.json',{'complete':all(r['status']!='failed' for r in results),'results':results})