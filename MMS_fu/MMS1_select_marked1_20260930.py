"""按原PPT中的独立文字标记1提取页面，保留全部原页面内容。"""
import hashlib
import json
import os
import posixpath
import shutil
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

root = Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815')
work = Path(os.environ['TEMP']) / 'MMS1_marked1_20260930'
build = work / 'build'
build.mkdir(parents=True, exist_ok=True)
source = root / 'MMS1.pptx'
snapshot = build / 'source.pptx'
ns = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def resolved(base, target):
    return target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join(posixpath.dirname(base), target))

if sys.argv[1] == 'inspect':
    assert not snapshot.exists(), '已有快照，避免覆盖'
    shutil.copy2(source, snapshot)
    with zipfile.ZipFile(snapshot) as z:
        pres = ET.fromstring(z.read('ppt/presentation.xml'))
        rels = {r.get('Id'):r for r in ET.fromstring(z.read('ppt/_rels/presentation.xml.rels'))}
        slides = []
        for number, elem in enumerate(pres.find('p:sldIdLst',ns),1):
            rid = elem.get('{'+ns['r']+'}id')
            part = resolved('ppt/presentation.xml',rels[rid].get('Target'))
            xml = ET.fromstring(z.read(part))
            texts = [''.join(t.itertext()).strip() for t in xml.findall('.//a:t',ns)]
            shape_texts = [''.join(t.text or '' for t in sp.findall('.//a:t',ns)).strip() for sp in xml.findall('.//p:sp',ns)]
            selected = '1' in shape_texts
            slides.append({'original_slide':number,'part':part,'rid':rid,'slide_id':elem.get('id'),
                           'texts':texts,'shape_texts':shape_texts,'selected':selected})
        manifest = {'source':str(source),'snapshot':str(snapshot),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                    'source_slide_count':len(slides),'slide_size':pres.find('p:sldSz',ns).attrib,
                    'selected_pages':[s['original_slide'] for s in slides if s['selected']],
                    'slides':slides}
    (build/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in manifest.items() if k!='slides'},ensure_ascii=False))
    print(json.dumps([s for s in slides if s['texts']],ensure_ascii=False))
elif sys.argv[1] == 'assemble':
    # Artifact Tool完成页选择后，复用原始页面和依赖部件，避免转换改变排版。
    from lxml import etree as LX
    manifest = json.loads((build/'manifest.json').read_text('utf-8'))
    chosen = [s for s in manifest['slides'] if s['selected']]
    with zipfile.ZipFile(build/'artifact_candidate.pptx') as authored:
        authored_pres = LX.fromstring(authored.read('ppt/presentation.xml'))
        assert len(authored_pres.find('p:sldIdLst', ns)) == len(chosen)
    with zipfile.ZipFile(snapshot) as z:
        blobs = {n:z.read(n) for n in z.namelist()}
        pres = LX.fromstring(blobs['ppt/presentation.xml'])
        sldlist = pres.find('p:sldIdLst',ns)
        keep_ids = {s['rid'] for s in chosen}
        for item in list(sldlist):
            if item.get('{'+ns['r']+'}id') not in keep_ids:
                sldlist.remove(item)
        # 本源稿没有自定义放映和分节；如未来出现则停止并显式处理。
        assert not pres.xpath('//*[local-name()="custShowLst" or local-name()="sectionLst"]')
        blobs['ppt/presentation.xml'] = LX.tostring(pres,xml_declaration=True,encoding='UTF-8',standalone=True)
        rels = LX.fromstring(blobs['ppt/_rels/presentation.xml.rels'])
        for relation in list(rels):
            if relation.get('Type').endswith('/slide') and relation.get('Id') not in keep_ids:
                rels.remove(relation)
        blobs['ppt/_rels/presentation.xml.rels'] = LX.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)
        roots = LX.fromstring(blobs['_rels/.rels'])
        for relation in list(roots):
            if relation.get('Type').endswith('/thumbnail'):
                roots.remove(relation)
        blobs['_rels/.rels'] = LX.tostring(roots,xml_declaration=True,encoding='UTF-8',standalone=True)
        # 只保留被选页面和共享母版所需的依赖，排除未选页面的原图。
        keep = {'_rels/.rels'}
        queue = [resolved('',r.get('Target')) for r in roots if r.get('TargetMode') != 'External']
        while queue:
            part = queue.pop()
            if part in keep:
                continue
            assert part in blobs, part
            keep.add(part)
            relpart = posixpath.join(posixpath.dirname(part),'_rels',posixpath.basename(part)+'.rels')
            if relpart in blobs:
                keep.add(relpart)
                for r in LX.fromstring(blobs[relpart]):
                    if r.get('TargetMode') != 'External':
                        queue.append(resolved(part,r.get('Target')))
        expected_slides = {s['part'] for s in chosen}
        assert {p for p in keep if p.startswith('ppt/slides/') and p.endswith('.xml')} == expected_slides
        if 'docProps/app.xml' in keep:
            app = LX.fromstring(blobs['docProps/app.xml'])
            for el in list(app):
                local = LX.QName(el).localname
                if local in ['HeadingPairs','TitlesOfParts']:
                    app.remove(el)
                elif local in ['Slides','Notes']:
                    el.text = str(len(chosen))
                elif local == 'HiddenSlides':
                    el.text = '0'
            blobs['docProps/app.xml'] = LX.tostring(app,xml_declaration=True,encoding='UTF-8',standalone=True)
        types = LX.fromstring(blobs['[Content_Types].xml'])
        for el in list(types):
            if LX.QName(el).localname == 'Override' and el.get('PartName').lstrip('/') not in keep:
                types.remove(el)
        keep.add('[Content_Types].xml')
        blobs['[Content_Types].xml'] = LX.tostring(types,xml_declaration=True,encoding='UTF-8',standalone=True)
        for part in expected_slides:
            assert blobs[part] == z.read(part)
        with zipfile.ZipFile(build/'candidate.pptx','w') as output:
            for info in z.infolist():
                if info.filename in keep:
                    output.writestr(info,blobs[info.filename])
        print(json.dumps({'slides':len(chosen),'parts':len(keep),'original_slide_xml_preserved':True}))
elif sys.argv[1] == 'publish':
    from PIL import Image, ImageChops
    manifest = json.loads((build/'manifest.json').read_text('utf-8'))
    candidate = work/'output'/'MMS1_marked1.pptx'
    final = root/'MMS1_标注1.pptx'
    assert not final.exists(), '已有同名输出，不覆盖'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest['source_sha256']
    for i in range(1,len(manifest['selected_pages'])+1):
        before = Image.open(work/'before'/f'slide_{i}.png').convert('RGB')
        after = Image.open(work/'after'/f'slide_{i}.png').convert('RGB')
        assert before.size == after.size and ImageChops.difference(before,after).getbbox() is None, f'第{i}页显示有变化'
    with zipfile.ZipFile(snapshot) as original, zipfile.ZipFile(candidate) as new:
        for s in manifest['slides']:
            if s['selected']:
                assert original.read(s['part']) == new.read(s['part'])
    with final.open('xb') as output:
        output.write(candidate.read_bytes())
    assert final.read_bytes() == candidate.read_bytes()
    receipt = {'output':str(final),'source':str(source),'selected_pages':manifest['selected_pages'],
               'page_count':len(manifest['selected_pages']),'native_render_pixel_equal':True,
               'source_unchanged':True,'sha256':hashlib.sha256(final.read_bytes()).hexdigest()}
    (build/'delivery.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False))
