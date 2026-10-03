"""独立只读核查PPT局部换图：ZIP部件、原生渲染、图片几何及用户文字。"""
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
from PIL import Image

WORK = Path(r'C:\Users\Administrator\AppData\Local\Temp\MMS1_BViAE_refresh_20260930')
OUT = Path(r'Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930\redraw\ppt_independent_verification.json')
SOURCE = WORK/'build/source.pptx'
FINAL = WORK/'output/MMS1_B_Vi_AE_20260720_0815_75percent.pptx'
manifest = json.loads((WORK/'build/replacement_manifest.json').read_text(encoding='utf-8'))
matches = manifest['matches']
expected_media = {row['media_part'] for row in matches}
expected_slides = {row['slide'] for row in matches}
matched_by_slide = {row['slide']:row for row in matches}
NS = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
checks = {}
report = {'verified_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
          'source':str(SOURCE),'candidate':str(FINAL),
          'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'candidate_sha256':hashlib.sha256(FINAL.read_bytes()).hexdigest(),
          'checks':checks}
checks['source_matches_manifest_hash'] = report['source_sha256']==manifest['source_sha256']
checks['manifest_has_exactly_26_unique_slides_and_media'] = len(matches)==len(expected_media)==len(expected_slides)==26

with zipfile.ZipFile(SOURCE) as before_zip, zipfile.ZipFile(FINAL) as after_zip:
    before_names,after_names = set(before_zip.namelist()),set(after_zip.namelist())
    same_names = before_names==after_names
    changed_parts = sorted(name for name in before_names&after_names if before_zip.read(name)!=after_zip.read(name))
    unchanged_parts = sorted(name for name in before_names&after_names if before_zip.read(name)==after_zip.read(name))
    checks['package_part_names_unchanged'] = same_names
    checks['only_specified_26_media_changed'] = set(changed_parts)==expected_media
    checks['other_801_parts_byte_identical'] = len(unchanged_parts)==801
    report['package'] = {'source_part_count':len(before_names),'candidate_part_count':len(after_names),
        'added_parts':sorted(after_names-before_names),'removed_parts':sorted(before_names-after_names),
        'changed_count':len(changed_parts),'unchanged_count':len(unchanged_parts),'changed_parts':changed_parts}
    slide_size = ET.fromstring(before_zip.read('ppt/presentation.xml')).find('p:sldSz',NS)
    slide_width,slide_height = int(slide_size.attrib['cx']),int(slide_size.attrib['cy'])
    geometry = {}
    image_checks = []
    for item in matches:
        slide_before = before_zip.read(item['slide_part'])
        slide_after = after_zip.read(item['slide_part'])
        old_root,new_root = ET.fromstring(slide_before),ET.fromstring(slide_after)
        old_pictures = [p for p in old_root.findall('.//p:pic',NS) if p.find('p:nvPicPr/p:cNvPr',NS).attrib['id']==item['picture_id']]
        new_pictures = [p for p in new_root.findall('.//p:pic',NS) if p.find('p:nvPicPr/p:cNvPr',NS).attrib['id']==item['picture_id']]
        assert len(old_pictures)==len(new_pictures)==1
        old_picture,new_picture = old_pictures[0],new_pictures[0]
        xfrm = old_picture.find('p:spPr/a:xfrm',NS)
        assert xfrm is not None
        off,ext = xfrm.find('a:off',NS),xfrm.find('a:ext',NS)
        frame = {key:int(value) for node in (off,ext) for key,value in node.attrib.items()}
        rotation = int(xfrm.attrib.get('rot',0))
        assert rotation==0, f'Unexpected picture rotation on slide {item["slide"]}'
        grouped = any(old_picture in list(g.iter()) for g in old_root.findall('.//p:grpSp',NS))
        assert not grouped, f'Unexpected grouped picture on slide {item["slide"]}'
        crop = old_picture.find('p:blipFill/a:srcRect',NS)
        image_checks.append({'slide':item['slide'],'media':item['media_part'],
            'slide_xml_byte_identical':slide_before==slide_after,
            'picture_xml_byte_identical':ET.tostring(old_picture)==ET.tostring(new_picture),
            'source_media_sha256_matches_manifest':hashlib.sha256(before_zip.read(item['media_part'])).hexdigest()==item['old_sha256'],
            'candidate_media_equals_exact_new_png':after_zip.read(item['media_part'])==Path(item['new_png']).read_bytes(),
            'frame_emu':frame,'crop':dict(crop.attrib) if crop is not None else None,
            'rotation':rotation,'grouped':grouped})
        geometry[item['slide']]=frame
    checks['all_26_picture_geometry_and_crop_preserved'] = all(x['picture_xml_byte_identical'] and x['slide_xml_byte_identical'] for x in image_checks)
    checks['all_26_media_match_manifest_and_new_pngs'] = all(x['source_media_sha256_matches_manifest'] and x['candidate_media_equals_exact_new_png'] for x in image_checks)
    report['replacement_images']=image_checks
    text_pages=[]
    for number in (37,50,67,87,96,110):
        part=f'ppt/slides/slide{number}.xml'
        old,new=before_zip.read(part),after_zip.read(part)
        text_pages.append({'slide':number,'slide_xml_byte_identical':old==new,
            'text_runs':[n.text or '' for n in ET.fromstring(old).findall('.//a:t',NS)]})
    checks['six_user_text_pages_xml_unchanged']=all(x['slide_xml_byte_identical'] for x in text_pages)
    report['user_text_pages']=text_pages

before_renders = list((WORK/'before').glob('*.PNG'))
after_renders = list((WORK/'after').glob('*.PNG'))
checks['both_native_render_folders_have_162_pages']=len(before_renders)==len(after_renders)==162
render_results=[]
for number in range(1,163):
    before_path=WORK/f'before/幻灯片{number}.PNG'
    after_path=WORK/f'after/幻灯片{number}.PNG'
    with Image.open(before_path) as old_image,Image.open(after_path) as new_image:
        assert old_image.size==new_image.size
        old=np.asarray(old_image.convert('RGBA'))
        new=np.asarray(new_image.convert('RGBA'))
        changed=np.any(old!=new,axis=2)
        rows,cols=np.nonzero(changed)
        count=int(changed.sum())
        result={'slide':number,'size_pixels':list(old_image.size),
                'changed_pixel_count':count,'pixel_identical':count==0}
        if number in expected_slides:
            frame=geometry[number]
            width,height=old_image.size
            exact=[frame['x']*width/slide_width,frame['y']*height/slide_height,
                (frame['x']+frame['cx'])*width/slide_width,(frame['y']+frame['cy'])*height/slide_height]
            bounds=[math.floor(exact[0]),math.floor(exact[1]),math.ceil(exact[2]),math.ceil(exact[3])]
            outside=(cols<bounds[0])|(cols>=bounds[2])|(rows<bounds[1])|(rows>=bounds[3])
            result.update({'picture_bounds_exact_pixels':exact,'picture_raster_bounds_pixels':bounds,
                'changed_pixel_bbox':([int(cols.min()),int(rows.min()),int(cols.max())+1,int(rows.max())+1] if count else None),
                'changed_pixels_outside_picture':int(outside.sum()),'all_changes_inside_picture':not bool(outside.any())})
        render_results.append(result)
unchanged=[x for x in render_results if x['slide'] not in expected_slides]
changed=[x for x in render_results if x['slide'] in expected_slides]
checks['all_136_unmodified_pages_pixel_identical']=len(unchanged)==136 and all(x['pixel_identical'] for x in unchanged)
checks['all_26_modified_pages_changed']=len(changed)==26 and all(x['changed_pixel_count']>0 for x in changed)
checks['all_26_modified_pages_differ_only_inside_original_picture']=all(x['all_changes_inside_picture'] for x in changed)
checks['six_user_text_pages_pixel_identical']=all(render_results[number-1]['pixel_identical'] for number in (37,50,67,87,96,110))
report['render_summary']={'pages':len(render_results),'unchanged_pages':len(unchanged),
    'changed_pages':len(changed),'unmodified_pages_changed_pixels':sum(x['changed_pixel_count'] for x in unchanged),
    'modified_pages_outside_picture_changed_pixels':sum(x['changed_pixels_outside_picture'] for x in changed)}
report['native_render_comparison']=render_results
report['passed']=all(checks.values())
OUT.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'passed':report['passed'],'checks':checks,'package':report['package'],
    'render_summary':report['render_summary'],'user_text_pages':report['user_text_pages'],
    'source_sha256':report['source_sha256'],'candidate_sha256':report['candidate_sha256'],
    'record':str(OUT)},ensure_ascii=False,indent=2))
assert report['passed'], 'Independent verification failed; see full record.'
