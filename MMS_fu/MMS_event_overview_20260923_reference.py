"""已下载官网图片的事件索引。原PNG原样复制，保留官网坐标系、时间窗和像素。"""
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil

CODE=Path(__file__).parent
OUT=Path(r'C:\Users\Administrator\Documents\KH\MMS_event_overviews_20260721_20260916')
WORK=Path(os.environ['TEMP'])/'MMS_overview_20260923'
AUDIT=Path(r'Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923')
events=json.loads((CODE/'MMS_event_overview_20260923_events.json').read_text('utf-8'))
images=json.loads((Path(os.environ['TEMP'])/'MMS_plots_20260919'/'expected_manifest.json').read_text('utf-8-sig'))

def parse(s):
    return dt.datetime.fromisoformat(s.replace('Z','+00:00')).replace(tzinfo=None)

def decode(item):
    if item['kind']=='quicklook':
        hit=re.search(r'_(\d{8})_(\d{4})_(\d{4})\.png$',item['name'])
        if not hit:return None
        a=dt.datetime.strptime(hit[1]+hit[2],'%Y%m%d%H%M')
        return a,a+dt.timedelta(minutes=int(hit[3])),int(hit[3])
    hit=re.search(r'_(\d{8})_(\d{6})\.png$',item['name'])
    if not hit:return None
    a=dt.datetime.strptime(hit[1]+hit[2],'%Y%m%d%H%M%S')
    return a,None,None

def copy_image(item,folder):
    target=folder/item['name']
    target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists():
        shutil.copy2(item['path'],target)
    return str(target)

index=[]
for event in events:
    a=parse(event['start'])-dt.timedelta(minutes=10)
    b=parse(event['end'])+dt.timedelta(minutes=10)
    for sc in range(1,5):
        folder=OUT/(event['id']+'_'+event['start'][:10])/'official_reference'/('MMS'+str(sc))
        record=dict(event=event['id'],spacecraft=sc,requestedStartUTC=a.isoformat(),requestedEndUTC=b.isoformat(),quicklook=[],burst=[])
        for plot in ('all_mms'+str(sc)+'_summ','fpi_mms'+str(sc)+'_summ'):
            eligible=[(item,decode(item)) for item in images if item['plot']==plot and item['kind']=='quicklook' and Path(item['path']).exists()]
            covering=[(item,span) for item,span in eligible if span[0]<=a and span[1]>=b]
            if covering and min(x[1][2] for x in covering)==120:
                selected=[min(covering,key=lambda x:x[1][2])]
            else:
                selected=sorted([(item,span) for item,span in eligible if span[2]==120 and span[0]<b and span[1]>a],key=lambda x:x[1][0])
            for item,span in selected:
                record['quicklook'].append(dict(path=copy_image(item,folder/'quicklook'),plot=plot,startUTC=span[0].isoformat(),endUTC=span[1].isoformat(),
                    sourcePath=item['path'],note='Unmodified official quicklook PNG. Coordinates and data levels are as labeled in the source image.'))
        # burst原图仅按起始时刻选入窗口内，不把图名推定为完整持续区间。
        for item in images:
            if item['kind']!='burst' or item['spacecraft']!='MMS'+str(sc) or item['instrument']!='综合图':continue
            span=decode(item)
            if span and a<=span[0]<=b and Path(item['path']).exists():
                record['burst'].append(dict(path=copy_image(item,folder/'burst'),startUTC=span[0].isoformat(),sourcePath=item['path'],
                    note='Start time inside requested window; exact end and gaps remain as shown in the original image.'))
        index.append(record)
WORK.mkdir(exist_ok=True)
(AUDIT/'official_reference_index.json').write_text(json.dumps(index,indent=2,ensure_ascii=False),'utf-8')
print(json.dumps(dict(eventSpacecraftSets=len(index),quicklook=sum(len(x['quicklook']) for x in index),burst=sum(len(x['burst']) for x in index))))
