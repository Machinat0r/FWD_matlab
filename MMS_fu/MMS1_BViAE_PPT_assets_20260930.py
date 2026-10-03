"""检查PPT图片引用，并将经Artifact Tool替换验证的新图片移回原始包。
只允许26个指定PNG媒体部件变化，其余PPT内部文件逐字节保留。
"""
import hashlib
import json
import os
import posixpath
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

work = Path(os.environ['TEMP']) / 'MMS1_BViAE_refresh_20260930'
build = work / 'build'
source = build / 'source.pptx'
result = Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815')
audit = Path(r'Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930')
ns = {'p': 'http://schemas.openxmlformats.org/presentationml/2006/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
phase = sys.argv[1]

if phase == 'inspect':
    changed = json.loads((audit / 'local' / 'after_download' / 'windows_audit.json').read_text('utf-8'))
    targets = {}
    for row in changed:
        png = result / 'B_Vi_AE' / (row['id'] + '_MMS1_B_Vi_AE.png')
        digest = hashlib.sha256(png.read_bytes()).hexdigest()
        assert digest not in targets
        targets[digest] = {'window_id': row['id'], 'old_png': str(png), 'old_sha256': digest,
                           'new_png': str(work / 'plots' / 'B_Vi_AE' / png.name)}
    assert len(targets) == 26
    matches = []
    slides = []
    with zipfile.ZipFile(source) as z:
        pres = ET.fromstring(z.read('ppt/presentation.xml'))
        size = pres.find('p:sldSz', ns).attrib
        relations = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('ppt/_rels/presentation.xml.rels'))}
        for number, item in enumerate(pres.find('p:sldIdLst', ns), 1):
            target = relations[item.get('{' + ns['r'] + '}id')]
            part = target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join('ppt', target))
            root = ET.fromstring(z.read(part))
            rel_part = posixpath.join(posixpath.dirname(part), '_rels', posixpath.basename(part) + '.rels')
            rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read(rel_part))}
            pics = root.findall('.//p:pic', ns)
            slides.append({'slide': number, 'part': part, 'pictures': len(pics),
                           'shapes': len(root.findall('.//p:sp', ns)), 'text_runs': len(root.findall('.//a:t', ns))})
            for pic in pics:
                rid = pic.find('.//a:blip', ns).get('{' + ns['r'] + '}embed')
                if rid is None:
                    continue
                target = rels[rid]
                media = target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join(posixpath.dirname(part), target))
                digest = hashlib.sha256(z.read(media)).hexdigest()
                if digest not in targets:
                    continue
                nv = pic.find('p:nvPicPr/p:cNvPr', ns)
                matches.append({**targets[digest], 'slide': number, 'slide_part': part,
                                'media_part': media, 'picture_name': nv.get('name'), 'picture_id': nv.get('id')})
        assert {r['old_sha256'] for r in matches} == set(targets), 'PPT中未能精确匹配全部待更新PNG'
        assert len({r['media_part'] for r in matches}) == 26
        state = {'source_path': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                 'slide_count': len(slides), 'slide_size_emu': size,
                 'matches': matches, 'slides': slides, 'original_parts': len(z.namelist())}
    (build / 'replacement_manifest.json').write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'slides': len(slides), 'matched_pictures': len(matches), 'changed_media': 26,
                      'affected_slides': [m['slide'] for m in matches], 'slide_size_emu': size,
                      'extra_text_runs': sum(s['text_runs'] for s in slides)}, ensure_ascii=False))

elif phase == 'assemble':
    state = json.loads((build / 'replacement_manifest.json').read_text('utf-8'))
    replacements = {}
    with zipfile.ZipFile(build / 'artifact_candidate.pptx') as authored:
        exported_images = {hashlib.sha256(authored.read(n)).hexdigest() for n in authored.namelist() if n.startswith('ppt/media/')}
    for item in state['matches']:
        blob = Path(item['new_png']).read_bytes()
        digest = hashlib.sha256(blob).hexdigest()
        assert digest in exported_images, '新图必须已经过Artifact Tool图片替换与导出'
        assert digest != item['old_sha256'], '目标图未发生变化'
        replacements[item['media_part']] = blob
    assert len(replacements) == 26
    with zipfile.ZipFile(source) as original, zipfile.ZipFile(build / 'candidate.pptx', 'w') as output:
        output.comment = original.comment
        for info in original.infolist():
            output.writestr(info, replacements.get(info.filename, original.read(info.filename)))
    print('ASSEMBLED_ONLY_26_MEDIA_PARTS')

elif phase == 'verify':
    state = json.loads((build / 'replacement_manifest.json').read_text('utf-8'))
    final = Path(sys.argv[2])
    allowed = {m['media_part'] for m in state['matches']}
    with zipfile.ZipFile(source) as old, zipfile.ZipFile(final) as new:
        assert old.namelist() == new.namelist(), 'PPT内部文件列表或顺序改变'
        assert old.comment == new.comment
        differences = [name for name in old.namelist() if old.read(name) != new.read(name)]
        assert set(differences) == allowed, ('除指定图片外有其他改动', differences)
        for m in state['matches']:
            assert new.read(m['media_part']) == Path(m['new_png']).read_bytes()
        proof = {'source_sha256': state['source_sha256'], 'final_sha256': hashlib.sha256(final.read_bytes()).hexdigest(),
                 'slide_count': state['slide_count'], 'changed_parts': differences,
                 'changed_media_count': len(differences), 'other_parts_identical': len(old.namelist()) - len(differences),
                 'all_slide_xml_and_relationships_identical': True,
                 'all_notes_layouts_masters_and_properties_identical': True,
                 'affected_slides': [m['slide'] for m in state['matches']]}
    (build / 'package_preservation.json').write_text(json.dumps(proof, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in proof.items() if k != 'changed_parts'}, ensure_ascii=False))
else:
    raise ValueError(phase)
