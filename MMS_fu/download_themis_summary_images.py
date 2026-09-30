"""下载 THEMIS Summary 页的全部卫星图种，按仪器或图种分类，日期保留在官方文件名中。
仅复制图片；不读取或改写科学数据。清单和运行日志放在临时目录。
"""
from __future__ import annotations
import argparse, collections, concurrent.futures, datetime as dt, hashlib
import http.client, io, json, os, re, ssl, threading, time, faulthandler
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin, urlsplit, unquote
from PIL import Image

HOST = 'themis.ssl.berkeley.edu'
ORIGIN = 'https://' + HOST
ROOT = ORIGIN + '/themisdata/overplots/'
ORBIT = ORIGIN + '/themisdata/thg/l0/asi/'
LOCAL = threading.local()
WORK = Path(os.environ['TEMP']) / 'codex_themis_download_20260919'
WORK.mkdir(parents=True, exist_ok=True)
STACK_LOG = (WORK / 'thread_stacks.log').open('w', encoding='utf-8')
faulthandler.dump_traceback_later(60, repeat=True, file=STACK_LOG)

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            value = dict(attrs).get('href', '')
            if value: self.links.append(value)

def request(url):
    parts = urlsplit(url)
    if parts.scheme != 'https' or parts.hostname != HOST:
        raise ValueError('Unapproved download origin: ' + url)
    last = None
    for attempt in range(5):
        try:
            if getattr(LOCAL, 'connection', None) is None:
                LOCAL.connection = http.client.HTTPSConnection(HOST, timeout=12, context=ssl.create_default_context())
            LOCAL.connection.request('GET', parts.path + ('?' + parts.query if parts.query else ''), headers={'User-Agent':'THEMIS-summary-image-archive/1.0', 'Accept-Encoding':'identity'})
            resp = LOCAL.connection.getresponse()
            data = resp.read()
            headers = dict(resp.getheaders())
            if resp.status != 200: raise IOError('HTTP %s: %s' % (resp.status, url))
            length = resp.getheader('Content-Length')
            if length is not None and len(data) != int(length): raise IOError('Incomplete response')
            LOCAL.reuse_count = getattr(LOCAL, 'reuse_count', 0) + 1
            if LOCAL.reuse_count >= 24:
                LOCAL.connection.close(); LOCAL.connection = None; LOCAL.reuse_count = 0
            return data, headers
        except Exception as exc:
            last = exc
            with (WORK / 'network_retries.jsonl').open('a', encoding='utf-8') as retry_log:
                retry_log.write(json.dumps({'url':url,'attempt':attempt+1,'error':str(exc)})+'\n')
            try: LOCAL.connection.close()
            except Exception: pass
            LOCAL.connection = None
            if attempt < 4: time.sleep(min(2 ** attempt, 8))
    raise last

def listing(url):
    data, _ = request(url)
    parser = Links(); parser.feed(data.decode('utf-8', errors='replace'))
    return sorted(set(parser.links))

def directories(url, pattern):
    return [x for x in listing(url) if re.fullmatch(pattern, x)]

def parse_date(name):
    match = re.search(r'(20\d{2})[-.]?(\d{2})[-.]?(\d{2})', name)
    if not match: return None
    try: return dt.date(*map(int, match.groups()))
    except ValueError: return None

def category(name):
    if re.match(r'th[a-e]_l2_', name):
        return name[:3] + '/' + name.split('_')[2]
    if name.startswith('thm_tohban_'):
        return 'multi_probe/' + name.split('_')[2]
    if '_wave_survey_' in name:
        return 'fgm_wave_survey/' + name.split('_wave_survey_')[1].rsplit('.',1)[0]
    if '-1days-memory-' in name: return 'memory/' + name.split('-1days-memory-')[1].split('.')[0]
    if name.startswith('goes_'): return 'goes/' + name.split('_')[1]
    if name.startswith(('metop', 'noaa')): return 'poes_metop/' + name.split('_')[0]
    if name.startswith('kompsat_'): return 'kompsat'
    return 'other_satellite'

def instrument_path(cat, name):
    parts=cat.split('/')
    group=parts[0]
    if re.fullmatch(r'th[a-e]',group):
        satellite='THEMIS_'+group[-1].upper()
        if parts[1]=='overview': folder='Overview/'+satellite
        else: folder='ESA/'+satellite+'/'+{'moms':'Space_Moments','gmoms':'Ground_Moments'}[parts[1]]
    elif group=='multi_probe':
        product=parts[1]
        if product.startswith('esa'): folder='ESA/Multi_Probe/'+product
        elif product.startswith('sst'): folder='SST/Multi_Probe/'+product
        elif product=='fgm': folder='FGM/Multi_Probe'
        elif product.startswith('fftfbk'): folder='DFB/'+product
        else: folder='Multi_Probe/'+product
    elif group=='fgm_wave_survey': folder='FGM/Wave_Survey/'+parts[1]
    elif group=='memory': folder='Engineering/Memory/'+parts[1]
    elif group in ('goes','poes_metop'): folder='Overview/'+parts[1].upper()
    elif group=='kompsat': folder='Overview/KOMPSAT'
    elif group=='orbits': folder='Orbits/'+parts[1]
    else: folder='Other_Satellite_Images/'+cat
    return folder+'/'+name

def discover(start, destination, workers):
    all_rows=[]; dirs=[]
    years=directories(ROOT, r'\d{4}/')
    for ys in years:
        year=int(ys[:-1])
        if year < start.year: continue
        for ms in directories(urljoin(ROOT, ys), r'\d{2}/'):
            month=int(ms[:-1])
            if (year,month) < (start.year,start.month): continue
            month_url=urljoin(ROOT, ys+ms)
            for ds in directories(month_url, r'\d{2}/'):
                day=dt.date(year,month,int(ds[:-1]))
                if day >= start: dirs.append((day,urljoin(month_url,ds)))
    print(json.dumps({'phase':'directory_inventory','day_directories':len(dirs),'first':str(min(x[0] for x in dirs)),'last':str(max(x[0] for x in dirs))}),flush=True)
    def one_day(item):
        day,url=item
        rows=[]
        for name in listing(url):
            if '/' in name or not re.search(r'\.(png|gif|jpe?g)$', name, re.I): continue
            if parse_date(name) != day: raise ValueError('Filename date mismatch: '+name)
            cat=category(name)
            rows.append({'date':str(day),'category':cat,'url':urljoin(url,name),'relative_path':instrument_path(cat, unquote(name))})
        return rows
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for rows in pool.map(one_day,dirs): all_rows.extend(rows)
    for ys in directories(ORBIT, r'\d{4}/'):
        year=int(ys[:-1])
        if year < start.year: continue
        for ms in directories(urljoin(ORBIT,ys), r'\d{2}/'):
            month=int(ms[:-1])
            if (year,month) < (start.year,start.month): continue
            url=urljoin(ORBIT,ys+ms)
            for name in listing(url):
                if '/' in name or not re.match(r'(orbit_|themis_foot_|themis_south_).*\.gif$',name): continue
                day=parse_date(name)
                if day is None or day < start: continue
                cat='orbits/' + re.split(r'_20\d{2}-',name)[0]
                all_rows.append({'date':str(day),'category':cat,'url':urljoin(url,name),'relative_path':instrument_path(cat, unquote(name))})
    all_rows.sort(key=lambda r:(r['date'],r['category'],r['relative_path']))
    assert len(all_rows)==len({r['url'] for r in all_rows})==len({r['relative_path'] for r in all_rows})
    inventory={'source':ORIGIN+'/summary.php','created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'start':str(start),'destination':str(destination),'day_directories':[str(d) for d,_ in dirs],'scope':'All satellite observation, multi-probe, FGM wave survey, memory, GOES, POES/METOP, KOMPSAT, orbit and footprint images; all listed time intervals. Ground magnetometer and all-sky ground camera products are outside satellite scope.','rows':all_rows}
    (WORK/'inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2),encoding='utf-8')
    summary={'phase':'inventory_complete','files':len(all_rows),'first':min(r['date'] for r in all_rows),'last':max(r['date'] for r in all_rows),'categories':dict(collections.Counter(r['category'] for r in all_rows))}
    (WORK/'inventory_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary),flush=True)
    return inventory

def validate(data):
    with Image.open(io.BytesIO(data)) as im:
        fmt=im.format; size=im.size
        im.verify()
    with Image.open(io.BytesIO(data)) as im:
        im.load()
    return fmt,size

def download(row,destination):
    target=destination/row['relative_path']
    try:
        if not target.resolve().is_relative_to(destination.resolve()): raise ValueError('Unsafe target path')
        if target.exists():
            data=target.read_bytes(); fmt,size=validate(data)
            return dict(row,status='existing_verified',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),format=fmt,width=size[0],height=size[1])
        data,headers=request(row['url'])
        fmt,size=validate(data)
        target.parent.mkdir(parents=True,exist_ok=True)
        partial=target.with_name(target.name+'.part')
        partial.write_bytes(data)
        partial.replace(target)
        return dict(row,status='downloaded',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),format=fmt,width=size[0],height=size[1],last_modified=headers.get('Last-Modified',''))
    except Exception as exc:
        return dict(row,status='failed',error=str(exc))

def run_download(inventory,destination,workers):
    rows=inventory['rows']; destination.mkdir(parents=True,exist_ok=True)
    total=len(rows); done=0; nbytes=0; failed=[]; counts=collections.Counter(); beginning=time.monotonic(); last=beginning
    status_path=WORK/'progress.json'
    def progress(final=False):
        elapsed=time.monotonic()-beginning
        state={'phase':'complete' if final else 'downloading','total':total,'completed':done,'verified':done-len(failed),'failed':len(failed),'bytes':nbytes,'elapsed_seconds':round(elapsed,1),'files_per_second':round(done/max(elapsed,1),2),'eta_seconds':round((total-done)*elapsed/max(done,1)),'by_group':dict(counts),'destination':str(destination)}
        try:
            status_tmp = status_path.with_suffix('.tmp')
            status_tmp.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
            status_tmp.replace(status_path)
            print(json.dumps(state,ensure_ascii=False),flush=True)
        except OSError as exc:
            with (WORK/'progress_errors.log').open('a',encoding='utf-8') as error_log:
                error_log.write(repr(exc)+'\n')
    progress()
    with (WORK/'download_results.jsonl').open('w',encoding='utf-8') as log, concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(download,r,destination):r for r in rows}
        for future in concurrent.futures.as_completed(futures):
            try:
                result=future.result()
            except Exception as exc:
                result=dict(futures[future],status='failed',error=repr(exc))
            done+=1
            log.write(json.dumps(result,ensure_ascii=False)+'\n'); log.flush()
            if result['status']=='failed': failed.append(result)
            else:
                nbytes+=result['bytes']; counts[result['category'].split('/')[0]]+=1
            now=time.monotonic()
            if now-last >= 15: progress(); last=now
    (WORK/'failed.json').write_text(json.dumps(failed,ensure_ascii=False,indent=2),encoding='utf-8')
    progress(True)
    return 1 if failed else 0

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--start',default='2026-07-20')
    parser.add_argument('--destination',default='C:/Users/Administrator/Documents/KH/SMILE_THEMIS')
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--inventory-only',action='store_true')
    parser.add_argument('--reuse-inventory',action='store_true')
    args=parser.parse_args(); destination=Path(args.destination)
    if args.reuse_inventory: inventory=json.loads((WORK/'inventory.json').read_text(encoding='utf-8'))
    else: inventory=discover(dt.date.fromisoformat(args.start),destination,args.workers)
    if not args.inventory_only: raise SystemExit(run_download(inventory,destination,args.workers))


