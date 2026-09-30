"""等候每个事件的所有候选CDF完成，再调用原样MATLAB脚本；最多两个MATLAB进程。"""
import json,os,time,subprocess,datetime as dt,hashlib
from pathlib import Path
code=Path(r'C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
audit=Path(r'Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2')
work=Path(os.environ['TEMP'])/'MMS_original_style_20260924_v2'
manifest=json.loads((audit/'particles'/'manifest.json').read_text(encoding='utf8'))
files={f['file_name']:f for f in manifest['files']}
source=code/'Overview_events_original_style_20260924_v2.m';source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
needed={}
for sc in range(1,5):
 for ie in range(1,21):
  needed[(sc,ie)]={n for c in manifest['candidate_coverage'] if c['event']==f'EV{ie:02d}' and c['product'].startswith(f'mms{sc}_') for n in c['candidate_files']}
pending=set(needed);running=[];finished=[];errors=[]
group=max([0]+[int(p.name.split('group')[1][:2]) for p in work.glob('overview_final_group??_mms?.log')])
for sc,ie in list(pending):
 q=audit/f'EV{ie:02d}_MMS{sc}_overview.json'
 if not q.is_file():continue
 r=json.loads(q.read_text(encoding='utf8'))
 if not r.get('complete') or not r.get('particlePanelsChecked') or r.get('error'):continue
 names=needed[(sc,ie)];ready=all(Path(files[n]['path']).is_file() and Path(files[n]['path']).stat().st_size==int(files[n]['file_size']) for n in names)
 if not ready:continue
 fresh=max([source.stat().st_mtime]+[Path(files[n]['path']).stat().st_mtime for n in names])
 if (r.get('status')=='no L2 data' and not names) or (r.get('png') and Path(r['png']).is_file() and Path(r['png']).stat().st_mtime>=fresh):
  pending.remove((sc,ie));finished.append([sc,ie])
print('REUSE_VERIFIED',len(finished),flush=True)
while pending or running:
 assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash,'Source changed during final run'
 for job in list(running):
  rc=job['process'].poll()
  if rc is None:continue
  job['handle'].close();running.remove(job)
  for ie in job['events']:
   q=audit/f'EV{ie:02d}_MMS{job["sc"]}_overview.json'
   if not q.is_file():errors.append([job['sc'],ie,'No output']);continue
   logtext=Path(job['log']).read_text(encoding='utf8',errors='replace')
   if f'DONE EV{ie:02d} MMS{job["sc"]} complete=1' not in logtext:errors.append([job['sc'],ie,'No completion in MATLAB log']);continue
   r=json.loads(q.read_text(encoding='utf8'))
   if not r.get('complete') or r.get('error') or not r.get('particlePanelsChecked'):errors.append([job['sc'],ie,r.get('error')]);continue
   finished.append([job['sc'],ie])
  print('FINISHED',job['sc'],job['events'],'rc',rc,flush=True)
 if errors:
  (work/'final_driver_errors.json').write_text(json.dumps(errors),encoding='utf8')
  break
 if len(running)<2 and pending:
  have={n:Path(f['path']).is_file() and Path(f['path']).stat().st_size==int(f['file_size']) for n,f in files.items()}
  for sc in range(1,5):
   if len(running)>=2:break
   if any(j['sc']==sc for j in running):continue
   ready=sorted(ie for s,ie in pending if s==sc and all(have[n] for n in needed[(s,ie)]))
   if not ready:continue
   group+=1;log=work/f'overview_final_group{group:02d}_mms{sc}.log';err=work/f'overview_final_group{group:02d}_mms{sc}.stdout'
   batch=f"addpath('{code.as_posix()}');EventList=[{' '.join(map(str,ready))}];Spacecraft={sc};Overview_events_original_style_20260924_v2"
   handle=err.open('w',encoding='utf8');start=time.time()
   proc=subprocess.Popen([r'C:\Matlab\bin\win64\MATLAB.exe','-batch',batch,'-logfile',str(log)],cwd=code,stdout=handle,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
   running.append(dict(sc=sc,events=ready,process=proc,handle=handle,start=start,log=str(log)))
   pending-={(sc,ie) for ie in ready}
   print('STARTED',sc,ready,'pid',proc.pid,flush=True)
   break  # stagger MATLAB initialization to avoid concurrent IRFU datastore writes
 state={'complete':not pending and not running and not errors,'finished':finished,'pending':sorted(pending),'running':[{'sc':j['sc'],'events':j['events'],'pid':j['process'].pid,'log':j['log']} for j in running],'errors':errors,'sourceHash':source_hash,'updatedUTC':dt.datetime.now(dt.timezone.utc).isoformat()}
 (work/'final_driver_state.json').write_text(json.dumps(state,indent=2),encoding='utf8')
 if pending or running:time.sleep(30)
print('FINAL_DRIVER_FINISHED',len(finished),'errors',errors,flush=True)