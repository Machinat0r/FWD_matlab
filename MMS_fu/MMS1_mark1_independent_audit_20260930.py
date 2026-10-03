"""只读检查用户PPT中标注为1的页面，按presentation实际顺序报告。"""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import posixpath
import re
import zipfile
import xml.etree.ElementTree as ET

SOURCE=Path(r'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1.pptx')
OUT=Path(os.environ['TEMP'])/'MMS1_mark1_independent_audit_20260930.json'
NS={'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'pr':'http://schemas.openxmlformats.org/package/2006/relationships'}

def resolve(source,target):
    return target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join(posixpath.dirname(source),target))

with zipfile.ZipFile(SOURCE) as archive:
    names=set(archive.namelist())
    rels={r.attrib['Id']:r.attrib for r in ET.fromstring(archive.read('ppt/_rels/presentation.xml.rels'))}
    presentation=ET.fromstring(archive.read('ppt/presentation.xml'))
    slides=[]
    for number,node in enumerate(presentation.find('p:sldIdLst',NS),start=1):
        rid=node.attrib['{'+NS['r']+'}id']
        part=resolve('ppt/presentation.xml',rels[rid]['Target'])
        root=ET.fromstring(archive.read(part))
        parents={child:parent for parent in root.iter() for child in parent}
        shapes=[]
        for shape in root.findall('.//p:sp',NS):
            cnv=shape.find('p:nvSpPr/p:cNvPr',NS)
            runs=[n.text or '' for n in shape.findall('.//a:t',NS)]
            paragraphs=[''.join(n.text or '' for n in p.findall('.//a:t',NS)) for p in shape.findall('p:txBody/a:p',NS)]
            text='\n'.join(paragraphs)
            if not text and not runs:
                continue
            hidden_ancestors=[]
            ancestor=shape
            while ancestor in parents:
                ancestor=parents[ancestor]
                if ancestor.tag=='{'+NS['p']+'}grpSp':
                    group_props=ancestor.find('p:nvGrpSpPr/p:cNvPr',NS)
                    if group_props is not None and group_props.attrib.get('hidden') in ('1','true'):
                        hidden_ancestors.append(dict(group_props.attrib))
            exact=text.strip()=='1'
            shapes.append({'id':cnv.attrib.get('id') if cnv is not None else None,
                'name':cnv.attrib.get('name') if cnv is not None else None,
                'hidden':cnv is not None and cnv.attrib.get('hidden') in ('1','true'),
                'hidden_parent_groups':hidden_ancestors,'text':text,'text_runs':runs,
                'exact_1':exact,'contains_1':'1' in text})
        all_text=[n.text or '' for n in root.findall('.//a:t',NS)]
        allshape_text=[t for shape in shapes for t in shape['text_runs']]
        unusual=[dict(tag=n.tag,attributes=n.attrib) for n in root.iter() if any(token in n.tag.lower() for token in ('ink','contentpart','oleobj'))]
        slide_rel_part=posixpath.join(posixpath.dirname(part),'_rels',posixpath.basename(part)+'.rels')
        extra_relationships=[]
        if slide_rel_part in names:
            for rel in ET.fromstring(archive.read(slide_rel_part)):
                rel_type=rel.attrib.get('Type','')
                if any(token in rel_type.lower() for token in ('comment','ink','notes','ole','vml')):
                    extra_relationships.append(dict(rel.attrib))
        selected=any(s['exact_1'] and not s['hidden'] and not s['hidden_parent_groups'] for s in shapes)
        slides.append({'page':number,'slide_part':part,'slide_hidden':root.attrib.get('show') in ('0','false'),
            'selected_by_visible_exact_1':selected,'all_text_runs':all_text,'text_shapes':shapes,
            'all_text_explained_by_text_shapes':all_text==allshape_text,
            'unusual_objects':unusual,'comments_ink_notes_relationships':extra_relationships})
    suspicious_parts=sorted(p for p in names if any(token in p.lower() for token in ('comments/','threadedcomments','/ink','inkml','.vml','embeddings/')))
    auxiliary_text=[]
    for part in sorted(names):
        if part.endswith('.xml') and ('notesSlides/notesSlide' in part or 'comment' in part.lower()):
            root=ET.fromstring(archive.read(part))
            text=[n.text for n in root.iter() if n.text and n.text.strip()]
            if text:
                auxiliary_text.append({'part':part,'text':text})
    report={'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'query_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'slide_count':len(slides),
        'selected_pages':[s['page'] for s in slides if s['selected_by_visible_exact_1']],
        'all_pages_with_any_text':[s for s in slides if s['all_text_runs']],
        'slides_with_unusual_objects':[s for s in slides if s['unusual_objects']],
        'slides_with_hidden_text_shapes':[s for s in slides if any(t['hidden'] or t['hidden_parent_groups'] for t in s['text_shapes'])],
        'slides_with_unexplained_text':[s for s in slides if not s['all_text_explained_by_text_shapes']],
        'suspicious_package_parts':suspicious_parts,'auxiliary_notes_or_comment_text':auxiliary_text,
        'hidden_slides':[s['page'] for s in slides if s['slide_hidden']],
        'slides':slides}
OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
brief={k:v for k,v in report.items() if k not in ('slides','all_pages_with_any_text')}
brief['all_pages_with_any_text']=[{'page':s['page'],'all_text_runs':s['all_text_runs'],
    'shapes':s['text_shapes']} for s in report['all_pages_with_any_text']]
brief['record']=str(OUT)
print(json.dumps(brief,indent=2,ensure_ascii=False))
