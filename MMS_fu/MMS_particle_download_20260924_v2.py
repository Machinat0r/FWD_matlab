"""复用已有官方CDF下载程序，仅更换仪器列表与本次记录目录；科学读取/绘图仍由MATLAB完成。"""
import importlib.util
from pathlib import Path
import json
code=Path(__file__).parent
spec=importlib.util.spec_from_file_location('existing_mms_download',code/'MMS_event_overview_20260923_download.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.AUDIT=Path(r'Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2\particles')
m.AUDIT.mkdir(exist_ok=True)
m.EVENTS=json.loads((code/'MMS_event_overview_20260924_events.json').read_text(encoding='utf8'))
m.PRODUCTS=[(inst,mode,desc) for inst,descs in [('feeps',['electron','ion']),('hpca',['ion','moments'])] for mode in ['brst','srvy'] for desc in descs]
m.main()