"""仅查询/下载官方原始 CDF；科学数据读取和绘图由 MATLAB 完成。"""
import concurrent.futures as cf
import datetime as dt
import json
import os
from pathlib import Path
import re
import threading
import time
import urllib.request
import urllib.parse
import urllib.error

CODE = Path(__file__).parent
ROOT = Path(r'Z:\SPART-WORK\Data\MMS')
AUDIT = ROOT / 'derived' / 'event_overview_20260923'
AUDIT.mkdir(parents=True, exist_ok=True)
BASE = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/'
EVENTS = json.loads((CODE / 'MMS_event_overview_20260923_events.json').read_text('utf-8'))
PRODUCTS = [('fgm','srvy',''),('fgm','brst',''),('fpi','fast','dis-moms'),
            ('fpi','fast','des-moms'),('fpi','brst','dis-moms'),('fpi','brst','des-moms'),
            ('edp','fast','dce'),('edp','brst','dce')]
lock = threading.Lock()
next_request = 0.0

def timestamp(s):
    return dt.datetime.fromisoformat(s.replace('Z','+00:00')).replace(tzinfo=dt.timezone.utc)

def save(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), 'utf-8')
    os.replace(tmp, path)

def open_url(url):
    global next_request
    for attempt in range(8):
        with lock:
            delay = max(0, next_request - time.monotonic())
            next_request = max(time.monotonic(), next_request) + 1.0
        if delay:
            time.sleep(delay)
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent':'MMS-event-overview/1.0'}), timeout=180)
        except urllib.error.HTTPError as exc:
            if exc.code not in (429,500,502,503,504):
                raise
            retry = exc.headers.get('Retry-After','')
            pause = float(retry) if retry.isdigit() else min(120, 8 * 2**attempt)
            with lock:
                next_request = max(next_request, time.monotonic()+pause)
            print('HTTP retry', exc.code, pause, flush=True)
        except (TimeoutError, urllib.error.URLError):
            if attempt == 7:
                raise
            time.sleep(min(60, 3*2**attempt))
    raise RuntimeError('Request retry limit: '+url)

def query(product, begin, end, spacecraft=(1,2,3,4)):
    inst, mode, desc = product
    params = dict(sc_id=','.join('mms'+str(n) for n in spacecraft), instrument_id=inst,
                  data_rate_mode=mode, data_level='l2', start_date=begin, end_date=end)
    if desc:
        params['descriptor'] = desc
    key = '_'.join([inst, mode, desc or 'none', ''.join(map(str,spacecraft)), begin, end]).replace(':','')
    cache = AUDIT/'queries'/(key+'.json')
    cache.parent.mkdir(exist_ok=True)
    if cache.exists():
        data = json.loads(cache.read_text('utf-8'))
    else:
        url = BASE+'file_info/science?'+urllib.parse.urlencode(params)
        with open_url(url) as response:
            data = json.load(response)
        data['_url'] = url
        data['_queried_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        save(cache,data)
    if data.get('files_were_truncated'):
        if len(spacecraft)>1:
            return sum((query(product,begin,end,(sc,)) for sc in spacecraft),[])
        names_url=BASE+'file_names/science?'+urllib.parse.urlencode(params)
        with open_url(names_url) as response:
            names=sorted(set(re.findall(r'mms[1-4]_[^,/\\\s]+\.cdf',response.read().decode())))
        known={f['file_name']:f for f in data.get('files',[])}
        missing=[n for n in names if n not in known]
        for offset in range(0,len(missing),10):
            batch=missing[offset:offset+10]
            url=BASE+'file_info/science?'+urllib.parse.urlencode({'files':','.join(batch)})
            with open_url(url) as response:
                chunk=json.load(response)
            if chunk.get('files_were_truncated'):
                raise RuntimeError('Metadata chunk still truncated '+key)
            known.update({f['file_name']:f for f in chunk.get('files',[])})
        if set(names)-set(known):
            raise RuntimeError('Missing named-file metadata '+key)
        data['files']=list(known.values())
        data['files_were_truncated']=False
        data['_resolved_via_file_names']=True
        data['_names_url']=names_url
        save(cache,data)
    files=data.get('files',[])
    pattern = re.compile(r'^mms[1-4]_'+re.escape(inst+'_'+mode+'_l2'+('_'+desc if desc else ''))+r'_\d{8,14}_v')
    if any(not pattern.match(f['file_name']) for f in files):
        raise RuntimeError('Unexpected product in query '+key)
    return files

def local_path(name):
    parts=name.split('_')
    tag=parts[-2]
    path=ROOT.joinpath(*parts[:-2],tag[:4],tag[4:6])
    if parts[2]=='brst':
        path=path/tag[6:8]
    return path/name

def download(item):
    name=item['file_name']; dest=Path(item['path']); expected=item['file_size']
    if dest.exists() and dest.stat().st_size==expected:
        return dict(name=name,status='existing',bytes=expected)
    # Parent directories are prepared sequentially to avoid SMB mkdir races.
    tmp=dest.with_suffix('.cdf.part')
    for attempt in range(4):
        try:
            url=BASE+'download/science?'+urllib.parse.urlencode({'file':name})
            with open_url(url) as response, tmp.open('wb') as output:
                while chunk:=response.read(1024*1024):
                    output.write(chunk)
            if tmp.stat().st_size!=expected:
                raise RuntimeError(f'Size mismatch {name}: {tmp.stat().st_size}/{expected}')
            with tmp.open('rb') as f:
                magic=f.read(4)
            if magic not in (bytes.fromhex('cdf30001'),bytes.fromhex('cdf26002'),bytes.fromhex('0000ffff')):
                raise RuntimeError('Unexpected CDF header '+name+' '+magic.hex())
            os.replace(tmp,dest)
            return dict(name=name,status='downloaded',bytes=expected)
        except Exception as exc:
            if attempt==3:
                return dict(name=name,status='failed',error=str(exc))
            time.sleep(4*(attempt+1))

def main():
    windows=[];days=set()
    for event in EVENTS:
        a=timestamp(event['start'])-dt.timedelta(minutes=10)
        b=timestamp(event['end'])+dt.timedelta(minutes=10)
        windows.append((a,b))
        day=a.date()
        while day<=b.date():
            days.add(day);day+=dt.timedelta(days=1)
    all_files={}; queries=[]
    jobs=[(p,day.isoformat(),(day+dt.timedelta(days=1)).isoformat()) for day in sorted(days) for p in PRODUCTS]
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        pending={pool.submit(query,*job):job for job in jobs}
        for future in cf.as_completed(pending):
            job=pending[future]
            files=future.result()
            for f in files:
                all_files[f['file_name']]=f
            queries.append(dict(product=job[0],day=job[1],files=len(files)))
            save(AUDIT/'inventory_progress.json',dict(done=len(queries),total=len(jobs),queries=queries))
            print('Inventory',len(queries),'/',len(jobs),job,len(files),flush=True)
    # 同一产品/时间只取数值版本号最高者；目录查询证据保留全部版本。
    latest={}
    for name,f in all_files.items():
        stem,version=name.rsplit('_v',1)
        ver=tuple(map(int,version.removesuffix('.cdf').split('.')))
        if stem not in latest or ver>latest[stem][0]:
            latest[stem]=(ver,f)
    groups={}
    for _,f in latest.values():
        key=f['file_name'].rsplit('_',2)[0]
        groups.setdefault(key,[]).append(f)
    selected={};coverage=[]
    for key,files in groups.items():
        files.sort(key=lambda f:f['timetag'])
        for ev,(a,b) in zip(EVENTS,windows):
            # 文件起始时间用于选候选；真正覆盖必须在 MATLAB 读取 Epoch 后验证。
            duration=dt.timedelta(hours=1 if '_brst_' in key else (2 if '_fpi_fast_' in key else 24))
            chosen=[]
            for i,f in enumerate(files):
                t=timestamp(f['timetag'])
                estimated_end=min(t+duration,timestamp(files[i+1]['timetag'])) if i+1<len(files) else t+duration
                if t<=b and estimated_end>=a:
                    name=f['file_name'];selected[name]={**f,'path':str(local_path(name))}
                    chosen.append(name)
            coverage.append(dict(event=ev['id'],product=key,candidate_files=chosen))
    selected_list=sorted(selected.values(),key=lambda f: ('_brst_' in f['file_name'], '_edp_' in f['file_name'], '_fgm_' in f['file_name'], f['file_name']))
    save(AUDIT/'manifest.json',dict(events=EVENTS,products=PRODUCTS,files=selected_list,candidate_coverage=coverage,
        note='File time tags identify candidates only. Read CDF record times to confirm coverage.'))
    print('Download',len(selected_list),'CDF, bytes',sum(f['file_size'] for f in selected_list),flush=True)
    for directory in sorted({str(Path(f['path']).parent) for f in selected_list}):
        Path(directory).mkdir(parents=True,exist_ok=True)
    results=[]
    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        pending={pool.submit(download,item):item for item in selected_list}
        for future in cf.as_completed(pending):
            result=future.result()
            result['completedUTC']=dt.datetime.now(dt.timezone.utc).isoformat()
            results.append(result)
            save(AUDIT/'download_progress.json',dict(done=len(results),total=len(selected_list),failed=sum(r['status']=='failed' for r in results),results=results))
            print('CDF',len(results),'/',len(selected_list),result['status'],result['name'],flush=True)
    save(AUDIT/'download_complete.json',dict(complete=all(r['status']!='failed' for r in results),results=results))

if __name__=='__main__':
    main()
