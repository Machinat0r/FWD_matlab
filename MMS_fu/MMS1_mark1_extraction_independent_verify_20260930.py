"""只读独立验证9页标注1提取稿：页序、原部件字节及PowerPoint原生像素。"""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import posixpath
import re
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image

WORK=Path(os.environ['TEMP'])/'MMS1_marked1_20260930'
SOURCE=WORK/'build/source.pptx'
CANDIDATE=WORK/'output/MMS1_marked1.pptx'
OUT=WORK/'build/independent_verification.json'
EXPECTED=[37,50,67,87,96,110,130,134,135]
NS={'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

def resolve(part,target):
    return target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join(posixpath.dirname(part),target))

def relationships(archive,part):
    name=posixpath.join(posixpath.dirname(part),'_rels',posixpath.basename(part)+'.rels')
    return {r.attrib['Id']:r.attrib for r in ET.fromstring(archive.read(name))}

def ordered_parts(archive):
    part='ppt/presentation.xml'
    root=ET.fromstring(archive.read(part))
    rels=relationships(archive,part)
    return [resolve(part,rels[n.attrib['{'+NS['r']+'}id']]['Target']) for n in root.find('p:sldIdLst',NS)]

def images(archive,part):
    root=ET.fromstring(archive.read(part))
    rels=relationships(archive,part)
    return [resolve(part,rels[n.attrib['{'+NS['r']+'}embed']]['Target']) for n in root.findall('.//a:blip',NS) if '{'+NS['r']+'}embed' in n.attrib]

manifest=json.loads((WORK/'build/manifest.json').read_text(encoding='utf-8'))
checks={}
source_hash=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
checks['source_matches_manifest_hash']=source_hash==manifest['source_sha256']
checks['manifest_selection_matches_independent_selection']=manifest['selected_pages']==EXPECTED
results=[]
with zipfile.ZipFile(SOURCE) as src,zipfile.ZipFile(CANDIDATE) as dest:
    original_parts,new_parts=ordered_parts(src),ordered_parts(dest)
    checks['source_162_pages_candidate_9_pages']=len(original_parts)==162 and len(new_parts)==9
    checks['selected_order_preserved']=new_parts==[original_parts[i-1] for i in EXPECTED]
    selected_media=set()
    for index,(original_number,new_part) in enumerate(zip(EXPECTED,new_parts),start=1):
        source_part=original_parts[original_number-1]
        old_xml,new_xml=src.read(source_part),dest.read(new_part)
        old_images,new_images=images(src,source_part),images(dest,new_part)
        selected_media.update(new_images)
        texts=[n.text or '' for n in ET.fromstring(new_xml).findall('.//a:t',NS)]
        all_image_bytes_equal=len(old_images)==len(new_images) and all(src.read(a)==dest.read(b) for a,b in zip(old_images,new_images))
        with Image.open(WORK/f'before/slide_{index}.png') as old_img,Image.open(WORK/f'after/slide_{index}.png') as new_img:
            old_pixels,new_pixels=np.asarray(old_img.convert('RGBA')),np.asarray(new_img.convert('RGBA'))
            size_ok=old_img.size==new_img.size==(1700,919)
            changed_pixels=int(np.any(old_pixels!=new_pixels,axis=2).sum()) if old_pixels.shape==new_pixels.shape else None
        results.append({'output_page':index,'source_page':original_number,'source_part':source_part,'candidate_part':new_part,
            'slide_xml_byte_identical':old_xml==new_xml,'image_count':len(new_images),'image_parts':new_images,
            'all_image_bytes_identical':all_image_bytes_equal,'texts':texts,'exact_marker_1_preserved':texts==['1'],
            'native_render_size_1700x919':size_ok,'changed_native_render_pixels':changed_pixels,
            'native_render_pixel_identical':changed_pixels==0})
    all_slide_parts={p for p in dest.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml',p)}
    all_media_parts={p for p in dest.namelist() if p.startswith('ppt/media/') and not p.endswith('/')}
    all_other_image_parts=[p for p in dest.namelist() if not p.startswith('ppt/media/') and re.search(r'\.(png|jpe?g|gif|bmp|wmf|emf|tiff?|webp)$',p,re.I)]
    checks['exactly_selected_9_slide_parts_no_unselected_slides']=all_slide_parts==set(new_parts) and len(all_slide_parts)==9
    checks['exactly_selected_image_parts_no_unselected_media']=all_media_parts==selected_media
    checks['no_other_image_parts_or_thumbnail']=not all_other_image_parts
    checks['all_9_slide_xml_byte_identical']=len(results)==9 and all(r['slide_xml_byte_identical'] for r in results)
    checks['all_9_image_bytes_unchanged']=len(results)==9 and all(r['all_image_bytes_identical'] for r in results)
    checks['all_9_markers_1_preserved']=len(results)==9 and all(r['exact_marker_1_preserved'] for r in results)
    checks['all_9_native_render_pixels_identical']=len(results)==9 and all(r['native_render_pixel_identical'] and r['native_render_size_1700x919'] for r in results)
report={'verified_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source':str(SOURCE),'candidate':str(CANDIDATE),
    'source_sha256':source_hash,'candidate_sha256':hashlib.sha256(CANDIDATE.read_bytes()).hexdigest(),
    'selected_source_pages':EXPECTED,'checks':checks,'passed':all(checks.values()),
    'candidate_slide_part_count':len(all_slide_parts),'candidate_media_part_count':len(all_media_parts),
    'other_image_parts':all_other_image_parts,'pages':results,'record':str(OUT)}
OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
assert report['passed'],'Independent extraction verification failed; see JSON.'
