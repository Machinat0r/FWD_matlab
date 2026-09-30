"""交付前检查：逐事件时间窗、文件、坐标系和所有本地链接。"""
import json,datetime as dt,os,argparse
from pathlib import Path
from collections import Counter
from html.parser import HTMLParser
from urllib.parse import unquote
from pypdf import PdfReader
CODE=Path(__file__).parent
AUDIT=Path(r'Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923')
OUT=Path(r'C:\Users\Administrator\Documents\KH\MMS_event_overviews_20260721_20260916')
EVENTS=json.loads((CODE/'MMS_event_overview_20260923_events.json').read_text('utf-8'))
REF=json.loads((AUDIT/'official_reference_index.json').read_text('utf-8'))
def parse(t):return dt.datetime.fromisoformat(t.replace('Z','+00:00')).replace(tzinfo=None)
def window(e):return parse(e['start'])-dt.timedelta(minutes=10),parse(e['end'])+dt.timedelta(minutes=10)
def main(final=False):
    errors=[];results=[]
    for p in AUDIT.glob('EV*_overview.json'):
        r=json.loads(p.read_text('utf-8'));results.append(r);e=EVENTS[int(r['event'][2:])-1];a,b=window(e)
        if (parse(r['startUTC']),parse(r['endUTC']))!=(a,b):errors.append(p.name+': time window')
        if r.get('renderVersion')!=3:errors.append(p.name+': original-program redraw pending')
        if r['readErrors']:errors.append(p.name+': read errors')
        for key in ['png','pdfPage']:
            if not Path(r[key]).is_file():errors.append(p.name+': missing '+key)
        if int(r['event'][2:])<=12:
            expected=6 if r['spacecraft']==4 else 8
            if r['statistics']['panelsAvailable']!=expected:errors.append(p.name+': panel count')
            if r.get('renderVersion')!=3:errors.append(p.name+': old render')
            if r['spacecraft']==4 and 'DSL XY' not in r['coordinateSystem']:errors.append(p.name+': E coordinates')
        elif r['statistics']['panelsAvailable']!=1:errors.append(p.name+': unexpected late L2 availability')
    for rr in REF:
        e=EVENTS[int(rr['event'][2:])-1];a,b=window(e)
        spans=sorted(set((parse(q['startUTC']),parse(q['endUTC'])) for q in rr['quicklook']))
        cursor=a
        for s,t in spans:
            if s<=cursor:cursor=max(cursor,t)
        if cursor<b:errors.append(f"{rr['event']}/MMS{rr['spacecraft']}: reference does not cover requested context")
        for q in rr['quicklook']+rr['burst']:
            if not Path(q['path']).is_file():errors.append('Missing reference: '+q['path'])
    summary=dict(nativeFigures=len(results),panels=dict(Counter(r['statistics']['panelsAvailable'] for r in results)),referenceSets=len(REF),errors=errors)
    if final:
        if len(results)!=60:errors.append('Expected 60 native figures')
        download=json.loads((AUDIT/'download_complete.json').read_text('utf-8'))
        manifest=json.loads((AUDIT/'manifest.json').read_text('utf-8'))
        if not download['complete'] or len(download['results'])!=len(manifest['files']):errors.append('Incomplete CDF download')
        summary.update(cdfFiles=len(download['results']),cdfBytes=sum(f['file_size'] for f in manifest['files']))
        d=json.loads((AUDIT/'delivery_manifest.json').read_text('utf-8'))
        pdf=PdfReader(d['pdf'])
        if len(pdf.pages)!=d['pages']:errors.append('PDF page count')
        if len(d['pageIndex'])!=80:errors.append('Expected 80 event/spacecraft sets')
        if any(x['lastPage']<x['firstPage'] for x in d['pageIndex']):errors.append('Empty event/spacecraft PDF section')
        class Links(HTMLParser):
            def handle_starttag(self,tag,attrs):
                for key,value in attrs:
                    if key not in ('href','src') or value.startswith(('#','http:','https:','data:')):continue
                    if not (OUT/unquote(value)).is_file():errors.append('Broken local link: '+value)
        Links().feed((OUT/'index.html').read_text('utf-8'))
        summary.update(pdfPages=len(pdf.pages),bookmarkedEvents=sum(isinstance(x,dict) for x in pdf.outline),eventSpacecraftSets=len(d['pageIndex']))
        (AUDIT/'delivery_validation.json').write_text(json.dumps(summary,indent=2),'utf-8')
    print(json.dumps(summary,ensure_ascii=False))
    if errors:raise SystemExit(1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--final',action='store_true');main(p.parse_args().final)