"""核对这次重绘与已交付版本的时间、CDF来源及数据覆盖，检查轴标签和panel间距。"""
import json,hashlib,sys
from pathlib import Path
code=Path(__file__).parent
old=Path(r'Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923')
new=Path(r'Z:\SPART-WORK\Data\MMS\derived\event_overview_20260924')
expected=['B [nT]','Vi [km/s]','Ve [km/s]','E [mV/m]','N [cm^{-3}]','Ti [eV]','Te [eV]','Ei(ev)','Ee(ev)']
errors=[];records=[]
for p in sorted(new.glob('EV*_overview.json')):
 r=json.loads(p.read_text('utf-8'));o=json.loads((old/p.name).read_text('utf-8'))
 for key in ['event','spacecraft','startUTC','endUTC','coordinateSystem','sources']:
  if r[key]!=o[key]:errors.append([p.name,key,'differs from baseline'])
 if r['readErrors'] or not r['complete']:errors.append([p.name,'read failed'])
 if r['renderVersion']!=4:errors.append([p.name,'wrong render version'])
 l=r['layout']
 if l['axisCount']!=9 or l['modeTimeline'] is not False or l['ylabels']!=expected:errors.append([p.name,'axis/labels/timeline'])
 pos=l['positions'];gaps=[pos[k][1]-(pos[k+1][1]+pos[k+1][3]) for k in range(8)]
 if any(abs(g-.002)>1e-10 for g in gaps):errors.append([p.name,'spacing',gaps])
 if r['statistics']['coverage']!=o['statistics']['coverage']:errors.append([p.name,'coverage changed'])
 a=r['statistics']['panelRecords'];b=o['statistics']['panelRecords']
 if a[:5]!=b[:5] or a[7:]!=b[6:] or a[5]+a[6]!=b[5]:errors.append([p.name,'record counts changed',a,b])
 png=Path(r['png'])
 if not png.exists() or png.stat().st_size<20000:errors.append([p.name,'missing PNG'])
 records.append({'event':r['event'],'spacecraft':r['spacecraft'],'panels':r['statistics']['panelsAvailable'],'png':r['png']})
report={'figures':len(records),'expected':60,'labels':expected,'errors':errors,'records':records}
(new/'redraw_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),'utf-8')
print(json.dumps({'figures':len(records),'expected':60,'errors':errors},ensure_ascii=False))
if errors or ('--complete' in sys.argv and len(records)!=60):sys.exit(1)
