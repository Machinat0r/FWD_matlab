"""编制带书签的PDF及本地图册；官网PNG保持原样。"""
import argparse,datetime as dt,html,json,os
from pathlib import Path
from collections import defaultdict
from urllib.parse import quote
from PIL import Image
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib.colors import HexColor
CODE=Path(__file__).parent
OUT=Path(r'C:\Users\Administrator\Documents\KH\MMS_event_overviews_20260721_20260916')
WORK=Path(os.environ['TEMP'])/'MMS_overview_20260923'
AUDIT=Path(r'Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923')
EVENTS=json.loads((CODE/'MMS_event_overview_20260923_events.json').read_text('utf-8'))
REF=json.loads((AUDIT/'official_reference_index.json').read_text('utf-8'))
REFMAP={(r['event'],r['spacecraft']):r for r in REF}
PDFNAME='MMS_20_events_MMS1-4_overviews.pdf'
INK=HexColor('#183047');MUTED=HexColor('#526778')
def rel(p):return quote(Path(p).relative_to(OUT).as_posix())
def stamp(t):return t.replace('T',' ').replace('Z','')
def expanded(e):
    a=dt.datetime.fromisoformat(e['start'].replace('Z','+00:00'))-dt.timedelta(minutes=10)
    b=dt.datetime.fromisoformat(e['end'].replace('Z','+00:00'))+dt.timedelta(minutes=10)
    return a.strftime('%Y-%m-%d %H:%M'),b.strftime('%Y-%m-%d %H:%M')
def basename(e,sc):return f"{e['id']}_{e['start'][:10].replace('-','')}_MMS{sc}_overview"
def native_result(e,sc):
    p=AUDIT/(basename(e,sc)+'.json')
    return json.loads(p.read_text('utf-8')) if p.exists() else None
def txt(c,x,y,text,size=12,color=INK,bold=False):
    c.setFillColor(color);c.setFont('Helvetica-Bold' if bold else 'Helvetica',size);c.drawString(x,y,text)
def fit(c,path,x,y,w,h):
    im=Image.open(path);iw,ih=im.size;s=min(w/iw,h/ih);ww,hh=iw*s,ih*s
    c.drawImage(ImageReader(im),x+(w-ww)/2,y+h-hh,width=ww,height=hh,mask='auto');im.close()
def reference_pages():
    dest=WORK/'reference_pages';dest.mkdir(exist_ok=True);output={}
    for e in EVENTS:
        for sc in range(1,5):
            groups=defaultdict(list)
            for q in REFMAP[e['id'],sc]['quicklook']:groups[q['startUTC'],q['endUTC']].append(q)
            pages=[]
            for num,((start,end),items) in enumerate(sorted(groups.items()),1):
                path=dest/f"{e['id']}_MMS{sc}_reference_{num:02d}.pdf"
                c=canvas.Canvas(str(path),pagesize=(1584,1260))
                txt(c,46,1216,f"{e['id']} | MMS{sc} | Official quicklook reference",23,bold=True)
                a,b=expanded(e);txt(c,46,1185,f"Requested context: {a} to {b} UTC",13)
                txt(c,46,1162,f"Original image interval: {stamp(start)} to {stamp(end)} UTC",12,MUTED)
                txt(c,46,1139,"Unchanged source images. Coordinates remain as labeled, including DMPA / DBCS / DSL.",12,MUTED)
                if len(items)==1:fit(c,items[0]['path'],55,70,1470,1040)
                else:
                    for j,q in enumerate(sorted(items,key=lambda q:q['plot'])):fit(c,q['path'],44+j*766,70,742,1040)
                txt(c,46,38,"MMS Science Data Center | Original quicklook PNGs archived locally; these pages are not CDF redraws.",11,MUTED)
                c.linkURL('https://lasp.colorado.edu/mms/sdc/public/quicklook/',(40,23,990,50),relative=0)
                c.showPage();c.save();pages.append(str(path))
            output[f"{e['id']}_MMS{sc}"]=pages
    (AUDIT/'reference_pdf_pages.json').write_text(json.dumps(output,indent=2),'utf-8')
    return output
def preface():
    path=WORK/'atlas_preface.pdf';c=canvas.Canvas(str(path),pagesize=(842,595))
    txt(c,48,526,'MMS Event Overviews',31,bold=True)
    txt(c,48,486,'20 user-selected intervals | MMS1-4 | July-September 2026',17)
    txt(c,48,446,'Data availability checked on 23 September 2026',12,MUTED)
    lines=[
        'Each requested interval is extended by 10 minutes at both ends (UTC).',
        'A single supplied time is the center of a 20-minute context window.',
        'The 11 August 15:00 +/-5 min interval becomes 14:45-15:15 UTC.','',
        'MATLAB redraws adapt the user script Overview_download.m (L2 CDF via IRFU).',
        'B/V and MMS1-3 E are in GSM; MMS4 E retains native DSL XY.',
        'Burst takes priority in its observed coverage; survey/fast supplies the background.',
        'Panels: B, Vi, Ve, E, ni/ne, Ti/Te, ion/electron differential energy-flux spectra.',
        'No smoothing or gap interpolation. Dense curves retain original display extrema.',
        'Public release limitations:',
        '- No requested FPI burst-moment L2 CDF was listed at retrieval.',
        '- MMS4 electron L2 moments were unavailable; its E panel uses native DSL XY only.',
        '- On 11 August only magnetic-field L2 data were available for these panels.',
        '- On 28-29 August and 16 September the required public L2 CDF were unavailable.','',
        'Official quicklook pages supplement incomplete later dates. Their coordinates',
        'and original time windows remain visible; original PNGs are also included.',
        'Available official burst overview PNGs are linked in index.html.','',
        'Event descriptions are user-supplied candidates, not classifications from this run.']
    y=405
    for line in lines:txt(c,48,y,line,11);y-=17
    txt(c,48,34,'Source: MMS Science Data Center public API and locally archived official overview images.',9,MUTED)
    c.linkURL('https://lasp.colorado.edu/mms/sdc/public/about/how-to/',(42,24,780,47),relative=0);c.showPage()
    for block in range(2):
        txt(c,40,550,f'Event index | {block*10+1}-{block*10+10}',22,bold=True)
        for x,label in [(40,'ID / date'),(175,'User interval (UTC)'),(315,'Plotted context (UTC)'),(542,'Available L2 redraw')]:txt(c,x,516,label,10,MUTED,True)
        y=484
        for e in EVENTS[block*10:block*10+10]:
            a,b=expanded(e);idx=int(e['id'][2:]);core=e['start'][11:16]
            if e['start']!=e['end']:core+='-'+e['end'][11:16]
            avail='B / particles / E; MMS4 electron gap' if idx<=12 else ('B only; official reference included' if idx<=15 else 'Official reference; CDF unavailable')
            context=a[11:]+'-'+b[11:] if a[:10]==b[:10] else a[5:]+' / '+b[5:]
            txt(c,40,y,e['id']+' '+e['start'][:10],10,bold=True);txt(c,175,y,core,10)
            txt(c,315,y,context,10);txt(c,542,y,avail,9)
            txt(c,175,y-15,e['label'][:96],9,MUTED);c.setStrokeColor(HexColor('#dde5ec'));c.line(40,y-24,802,y-24);y-=44
        txt(c,40,30,'Use PDF bookmarks or index.html to open an event and spacecraft.',10,MUTED);c.showPage()
    c.save();return path
def make_html(results):
    sections=[]
    for e in EVENTS:
        a,b=expanded(e);cards=[]
        for sc in range(1,5):
            r=results.get((e['id'],sc));rr=REFMAP[e['id'],sc]
            if r:
                status='CDF 重绘' if r['statistics']['panelsAvailable']>1 else '仅磁场 L2 可用'
                native=f'<h4>{status}</h4><a href="{rel(r["png"])}"><img loading="lazy" src="{rel(r["png"])}" alt="{e["id"]} MMS{sc}"></a>'
            else:native='<p class="notice">所需公开 CDF 暂未查到；下面为官网原始参考图。</p>'
            ql=''.join(f'<a class="reference" href="{rel(q["path"])}"><img loading="lazy" src="{rel(q["path"])}" alt="{html.escape(q["plot"])}"><span>{html.escape(q["plot"])}<br>{stamp(q["startUTC"])} — {stamp(q["endUTC"])} UTC</span></a>' for q in rr['quicklook'])
            burst=''.join(f'<a href="{rel(q["path"])}">{q["startUTC"][11:]} UTC</a>' for q in rr['burst'])
            cards.append(f'<div class="sc sc{sc}"><h3>MMS{sc}</h3>{native}<details><summary>官网 quicklook 原图（{len(rr["quicklook"])} 张）</summary><div class="refs">{ql}</div></details><details><summary>官网 burst 原图（{len(rr["burst"])} 张）</summary><div class="burstlinks">{burst or "本次官网清单未找到窗口内起始图。"}</div></details></div>')
        sections.append(f'<section id="{e["id"]}"><h2>{e["id"]} · {e["start"][:10]}</h2><p>{html.escape(e["label"])}</p><p class="muted">绘图窗口：{a} — {b} UTC</p>{"".join(cards)}</section>')
    nav=' '.join(f'<a href="#{e["id"]}">{e["id"]}</a>' for e in EVENTS)
    doc=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MMS 20 个事件 overview</title>
<style>body{{font:16px/1.55 system-ui,"Microsoft YaHei",sans-serif;margin:0;background:#f3f6f9;color:#193044}}header,main{{max-width:1250px;margin:auto;padding:28px}}header{{padding-bottom:0}}h1{{font-size:30px}}h2{{margin-bottom:4px}}a{{color:#125d9c}}.muted{{color:#5d6e7e}}.notice{{border-left:4px solid #d68c23;background:#fff6e8;padding:14px}}nav{{display:flex;flex-wrap:wrap;gap:10px;margin:18px 0}}nav a{{background:white;padding:6px 12px;border-radius:5px}}.toolbar{{position:sticky;top:0;background:#193044;color:white;padding:12px 28px;z-index:2}}select{{font-size:16px;padding:5px 20px;margin-left:12px}}section{{background:white;border:1px solid #d8e2ea;border-radius:8px;padding:24px;margin-bottom:28px;scroll-margin-top:74px}}.sc>a img{{width:100%;max-width:1100px}}.sc{{display:none}}.sc1{{display:block}}details{{margin:18px 0;border-top:1px solid #d8e2ea;padding-top:14px}}summary{{cursor:pointer;font-weight:600}}.refs{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px;margin-top:18px}}.reference img{{width:100%}}.reference span{{display:block;font-size:13px}}.burstlinks{{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}}.burstlinks a{{padding:5px 10px;background:#fff2dd;border-radius:4px}}@media(max-width:750px){{.refs{{grid-template-columns:1fr}}header,main{{padding:15px}}section{{padding:14px}}}}</style>
<header><h1>MMS 20 个事件 overview</h1><p>每个事件覆盖指定时段前后各 10 分钟；MMS1–4 分别处理，时间为 UTC。</p><p><a href="{PDFNAME}">打开汇总 PDF</a> · <a href="README_说明.md">数据覆盖与说明</a></p><p class="notice">原始 CDF 的公开发布尚不齐全。MMS4 电场使用原生 DSL XY 分量。8 月 11 日的重绘仅有磁场；8 月 28–29 日和 9 月 16 日提供官网参考图。MMS4 部分电子 L2 数据缺失。官网图保留其原有坐标系及时间范围。</p><nav>{nav}</nav></header>
<div class="toolbar">选择卫星 <select onchange="document.querySelectorAll('.sc').forEach(x=>x.style.display=x.classList.contains('sc'+this.value)?'block':'none')"><option value="1">MMS1</option><option value="2">MMS2</option><option value="3">MMS3</option><option value="4">MMS4</option></select></div><main>{"".join(sections)}</main></html>'''
    (OUT/'index.html').write_text(doc,'utf-8')
def make_readme(results,page_count):
    all8=sum(r['statistics']['panelsAvailable']==8 for r in results.values())
    bonly=sum(r['statistics']['panelsAvailable']==1 for r in results.values())
    lines=['# MMS 20 个候选事件 overview','',
        f'20 个事件，MMS1–4 分别处理。汇总 PDF 共 {page_count} 页，附事件/卫星书签；index.html 可切换卫星并打开原图。','',
        f'生成 {len(results)} 张 MATLAB L2 重绘图，其中 {all8} 张八个面板均有记录、{bonly} 张仅磁场有 L2 记录。缺测或未发布的量明确留空。','',
        '## 时间与物理量','',
        '- 时间按 UTC。单个时刻前后各 10 分钟；区间两端各加 10 分钟。',
        '- 8 月 11 日 15:00 ±5 分钟按 14:55–15:05 处理，最终窗口 14:45–15:15。',
        '- MATLAB 图的 B、Vi、Ve 及 MMS1–3 的 E 使用 GSM；MMS4 电场使用 CDF 原生 DSL XY 两个分量并单独标注，未推算缺失分量。其余为密度 ni/ne；标量温度 Ti/Te=(Tparallel+2Tperpendicular)/3；FPI 离子/电子全向差分能量通量谱。',
        '- 不平滑、不跨缺口插值。burst 按实际观测覆盖优先，其余采用 survey/fast。底部灰条为 survey/fast，橙条为 burst。',
        '- 密集曲线为显示保留各小块内原始极值点；CDF 保留原始分辨率。电场分段读取，原生记录数量与覆盖单独核对。',
        '- 对数能谱不显示非正通量；CDF 填充值留空。本次没有新增物理筛选判据，没有进行 FOTE 或事件分类。','',
        '## 原程序复用','',
        '- 绘图程序 Overview_download_events_20260923.m 从原 Overview_download.m 复制，保留 B、Vi/Ve、N、能谱绘图段，电场段复用 Overview_download_mms4.m。',
        '- 保留 c_eval、irf_subplot、irf_plot、irf_legend 和 irf_spectrogram 写法及原 b/g/r、jet 配色；温度采用原有标量温度公式。',
        '- 必要修改为事件/卫星参数化、burst 与 survey/fast 选择、缺测/二维电场适配、逐记录能量表、去掉旧事件专用坐标轴范围、输出位置。原程序保持不动。','',
        '## 当前公开数据限制','',
        '- 本次未找到所需 FPI burst L2 矩产品。可用粒子 L2 数据使用 fast；磁场/电场在有 burst 时优先使用 burst。',
        '- 早期日期 MMS4 电子 L2 矩缺失，Ve/ne/Te/电子能谱保留对应缺口。MMS4 电场 CDF 仅提供 DSL XY 二维变量，图中显示这两个原生分量。',
        '- 8 月 28–29 日的 MMS3 FPI 官网参考图也未在现有官网清单中找到，不能补齐这些事件的电子速度。',
        '- 8 月 11 日三个事件：所需 L2 产品仅磁场可用；另附官网综合和 FPI quicklook 原图。',
        '- 8 月 28–29 日、9 月 16 日五个事件：所需公开 CDF 未查到，以官网原图供检查。尚未完成这些日期的全物理量 CDF 重绘。',
        '- 官网 PNG 未改写。坐标系沿用原图 DMPA、DBCS、DSL 等标注，不能直接按 GSM 分量解读。',
        '- 官网固定窗口可能更宽；所选 quicklook 图覆盖事件前后至少 10 分钟。跨日事件使用相邻窗口。',
        '- burst 参考 PNG 按起始时刻位于目标窗口内选入；结束时刻及缺口以原图为准。该列表不宣称穷尽跨边界 burst。','',
        '## 事件索引','','| ID | 用户日期/时间 UTC | 本次窗口 UTC | L2 重绘范围 |','|---|---|---|---|']
    for e in EVENTS:
        a,b=expanded(e);idx=int(e['id'][2:])
        status='磁场/粒子/电场；MMS4 电子缺口' if idx<=12 else ('仅磁场，另附官网图' if idx<=15 else '官网参考图；CDF 未公开齐全')
        core=stamp(e['start']) if e['start']==e['end'] else stamp(e['start'])+' — '+e['end'][11:19]
        lines.append(f'| {e["id"]} | {core} | {a} — {b} | {status} |')
    lines+=['','## 保存位置','',
        '- 新图及说明：当前文件夹。官网参考原图保留原名，在各事件的 official_reference/MMS*/ 中。',
        r'- 原 CDF：Z:\SPART-WORK\Data\MMS 标准产品层级。',
        r'- 查询清单、CDF 来源、真实覆盖统计：Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923。',
        r'- 主程序：C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\MMS_event_overview_20260923.m。',
        '- 下载、索引和打包程序位于同一代码目录，以 MMS_event_overview_20260923 开头。','',
        '来源：[MMS SDC API](https://lasp.colorado.edu/mms/sdc/public/about/how-to/)；[MMS Quicklook](https://lasp.colorado.edu/mms/sdc/public/quicklook/)。可用性按 2026-09-23 本次查询。']
    (OUT/'README_说明.md').write_text('\n'.join(lines)+'\n','utf-8')
def finalize(refpages):
    results={(e['id'],sc):r for e in EVENTS for sc in range(1,5) if (r:=native_result(e,sc))}
    if len(results)!=60:raise RuntimeError(f'Wait for 60 L2 redraws; currently {len(results)}.')
    if any(r['readErrors'] for r in results.values()):raise RuntimeError('CDF read errors remain.')
    writer=PdfWriter();writer.append(str(preface()),import_outline=False);report=[]
    for e in EVENTS:
        mark=writer.add_outline_item(e['id']+' '+e['start'][:10]+' '+e['start'][11:16],len(writer.pages))
        for sc in range(1,5):
            p=len(writer.pages);writer.add_outline_item('MMS'+str(sc),p,parent=mark);r=results.get((e['id'],sc))
            if r:writer.append(r['pdfPage'],import_outline=False)
            if int(e['id'][2:])>=13:
                for path in refpages[f"{e['id']}_MMS{sc}"]:writer.append(path,import_outline=False)
            report.append(dict(event=e['id'],spacecraft=sc,firstPage=p+1,lastPage=len(writer.pages),native=bool(r)))
    writer.add_metadata({'/Title':'MMS 20 events: MMS1-4 overviews','/Subject':'L2 CDF redraws and explicitly labeled official references'})
    target=OUT/PDFNAME
    with target.open('wb') as f:writer.write(f)
    check=PdfReader(str(target))
    if len(check.pages)!=len(writer.pages):raise RuntimeError('PDF page count changed.')
    make_html(results);make_readme(results,len(check.pages))
    receipt=dict(pdf=str(target),pages=len(check.pages),nativeFigures=len(results),eventSpacecraftSets=len(report),pageIndex=report)
    (AUDIT/'delivery_manifest.json').write_text(json.dumps(receipt,indent=2),'utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k!='pageIndex'}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');args=ap.parse_args()
    refs=reference_pages();print('Prepared reference PDF pages:',sum(map(len,refs.values())),flush=True)
    if not args.prepare:finalize(refs)
