"""合并MATLAB矢量PDF并检查五事件交付；不改变或计算科学数据。"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
from pypdf import PdfReader, PdfWriter
from PIL import Image

record = Path(r'Z:\SPART-WORK\Data\MMS\derived\MMS_current_5events_20261001')
original = record.parent / 'MMS1_tail_20260720_0815'
output = Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS_current_5events_20261001')
qa = Path(os.environ['TEMP']) / 'MMS_current_5events_20261001' / 'qa'
qa.mkdir(parents=True, exist_ok=True)
windows = json.loads((record / 'windows.json').read_text('utf-8'))
writer = PdfWriter()
reports = []
for w in windows:
    r = json.loads((record / (w['id'] + '_current.json')).read_text('utf-8'))
    old = json.loads((original / (w['id'] + '_overview.json')).read_text('utf-8'))
    assert r['complete'] and not r['error'], r
    assert r['originalBViCountsMatch'] and all(r['panelAvailable'])
    assert r['startUTC'] == old['startUTC'] == w['startUTC']
    assert r['endUTC'] == old['endUTC'] == w['endUTC']
    assert r['aeRows'] == old['aeRows'] and r['aeValid'] == old['aeValid']
    assert r['coordinateSystem'] == 'GSM' and r['currentUnit'] == 'nA/m^2'
    assert r['timeAveraging'] and r['currentAveragingSeconds']==60
    assert not r['originalPanelsTimeAveraging'] and not r['smoothing'] and not r['qualityFilter']
    assert not r['fixedBackgroundField'] and r['minuteAverage']['decompositionBeforeAveraging']
    assert all(v for v,m in zip(r['instantaneousReferenceVerified'],r['current']) if m['finiteCurrentRows']>0)
    assert len(r['current']) == 2
    assert r['parallelSigned']
    assert r['parallelOpacity'] == 1 and r['perpendicularOpacity'] == 0.5
    assert r['panelOrder'] == ['B', 'Vi', '|J|', 'Jparallel,Jperp', 'AE']
    pdf, png = Path(r['pdf']), Path(r['png'])
    reader = PdfReader(pdf)
    assert len(reader.pages) == 1 and pdf.stat().st_size > 10000
    text = reader.pages[0].extract_text()
    assert 'AE' in text and 'Vi' in text and 'J' in text
    writer.add_page(reader.pages[0])
    with Image.open(png) as im:
        size = im.size
        assert size[0] > 1000 and size[1] > 1000
    reports.append({'id': w['id'], 'startUTC': r['startUTC'], 'endUTC': r['endUTC'],
                    'pngSize': size, 'png': str(png), 'pdf': str(pdf),
                    'pngSHA256': hashlib.sha256(png.read_bytes()).hexdigest(),
                    'pdfSHA256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
                    'baselineCountsMatch': True, 'aeCountsMatch': True,
                    'currentValidRowsByMode': [m['finiteCurrentRows'] for m in r['current']],
                    'parallelNegativeRowsByMode': [m.get('parallelNegativeRows',0) for m in r['current']],
                    'panelOrder': r['panelOrder'], 'parallelSigned': r['parallelSigned'],
                    'perpendicularDefinition': r['perpendicularDefinition'],
                    'parallelOpacity': r['parallelOpacity'],
                    'perpendicularOpacity': r['perpendicularOpacity'],
                    'currentAveragingSeconds':r['currentAveragingSeconds'],
                    'minuteAverageRows':r['minuteAverage']['rows'],
                    'minuteAverageFiniteRows':r['minuteAverage']['finiteRows'],
                    'instantaneousReferenceVerified':r['instantaneousReferenceVerified']})
combined = output / 'MMS_5events_B_Vi_J_AE_20261001.pdf'
writer.add_metadata({'/Title': 'MMS four-spacecraft current: five events, 1-minute means', '/Author': 'MATLAB / IRFU'})
with combined.open('wb') as f:
    writer.write(f)
assert len(PdfReader(combined).pages) == 5
poppler = Path(r'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin')
subprocess.run([str(poppler / 'pdftoppm.exe'), '-r', '120', '-png', str(combined), str(qa / 'page')], check=True)
rendered = sorted(qa.glob('page-*.png'))
assert len(rendered) == 5
for p in rendered:
    with Image.open(p) as im:
        assert im.size[0] > 1000 and im.size[1] > 800
algorithm = json.loads((record / 'algorithm_verification.json').read_text('utf-8'))
assert algorithm['passed']
report = {'complete': True, 'combinedPDF': str(combined), 'pages': 5,
          'combinedSHA256': hashlib.sha256(combined.read_bytes()).hexdigest(),
          'events': reports, 'algorithmVerification': algorithm,
          'qaRenderedPages': [str(p) for p in rendered], 'visualReview': 'pending'}
(record / 'delivery_verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), 'utf-8')
print(json.dumps({'combined': str(combined), 'events': reports, 'rendered': [str(p) for p in rendered]}, ensure_ascii=False))
