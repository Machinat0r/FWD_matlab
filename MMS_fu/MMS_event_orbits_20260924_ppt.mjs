import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {PresentationFile,FileBlob} from '@oai/artifact-tool';
const SKILL=String.raw`C:\Users\Administrator\.codex\plugins\cache\openai-primary-runtime\presentations\26.904.11930\skills\presentations`;
const CODE=String.raw`C:\Users\Administrator\Documents\FWD_matlab\MMS_fu`;
const OUT=String.raw`C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924`;
const AUDIT=String.raw`Z:\SPART-WORK\Data\MMS\derived\event_orbits_20260924`;
const OLD=String.raw`Z:\SPART-WORK\Data\MMS\derived\event_overview_20260924`;
const WORK=path.join(process.env.TEMP,'MMS_orbits_20260924','ppt_build');
const PYTHON=String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`;
const source=path.join(OUT,'MMS_20_events_MMS1-4_overviews_20260924.pptx');
const name='MMS_20_events_MMS1-4_overviews_with_orbits_20260924.pptx';
process.env.RUNTIME_NODE_MODULES=path.resolve(path.dirname(process.execPath),'..','node_modules');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
await fs.mkdir(path.join(WORK,'output'),{recursive:true});
const events=JSON.parse(await fs.readFile(path.join(CODE,'MMS_event_overview_20260924_events.json'),'utf8'));
const oldIndex=JSON.parse(await fs.readFile(path.join(OLD,'ppt_slide_index.json'),'utf8'));
const p=await PresentationFile.importPptx(await FileBlob.load(source));
const snapshot=await p.inspect({kind:'slide,layout',maxChars:1000000});
await fs.writeFile(path.join(WORK,'source_inspect.ndjson'),snapshot.ndjson);
const originals=[...p.slides.items];
if(originals.length!==107||oldIndex.length!==107)throw Error('Unexpected source coverage');
const font='Arial';
function text(slide,value,x,y,w,h,size=22,bold=false,color='#183047'){
 const sh=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=value;sh.text.style={typeface:font,fontSize:size,bold,color,autoFit:'none'};
}
const orbitMap=new Map();
for(const ev of [...events].reverse()){
 const record=JSON.parse(await fs.readFile(path.join(AUDIT,ev.id+'_MMS1-4_orbit_GSM.json'),'utf8'));
 if(!record.complete||record.positionKm.length!==4||record.positionKm.flat().some(x=>!Number.isFinite(x)))throw Error('Invalid orbit '+ev.id);
 const first=oldIndex.findIndex(x=>x.event===ev.id);
 if(first<1)throw Error('Missing insertion anchor '+ev.id);
 const {slide:s}=await p.slides.insert({after:originals[first-1]});
 s.background.fill='#FFFFFF';
 const title=ev.id+'  MMS1–4 spacecraft positions';
 text(s,title,38,20,1200,48,30,true);
 text(s,'GSM     '+record.snapshotUTC.replace('T',' ').replace('Z',' UTC'),38,73,1200,37,24);
 text(s,ev.label,38,115,1200,42,19,false,'#526778');
 s.images.add({blob:new Uint8Array(await fs.readFile(record.png)),contentType:'image/png',
   fit:'contain',position:{left:22,top:167,width:1236,height:715},
   alt:ev.id+' MMS1–4 GSM positions and relative configuration generated with the existing MMS_orbit program'});
 const sourceText=Number(ev.id.slice(2))<=15?'MMS MEC L2 ephemeris':'NASA SSCWeb ephemeris';
 text(s,'Position data: '+sourceText+'.  Snapshot at supplied time or interval midpoint.  R_E = 6372 km.',38,906,1200,34,18,false,'#526778');
 s.speakerNotes.textFrame.setText([
  'User request: add spacecraft positions for every event using the existing MMS_orbit program.',
  'Original script: '+path.join(String.raw`C:\Users\Administrator\Documents\FWD_matlab\新建文件夹`,'MMS_orbit.m'),
  'Entry: MMS_orbit_events_20260924.m. Compatibility copy of mms.mms4_pl_conf changes function name/callbacks and timestamp formatting only.',
  'Presentation layout: original seven panels and legend arranged into two rows. Source overview pages remain unchanged.',
  'Event UTC: '+ev.start+' to '+ev.end+'. Snapshot: '+record.snapshotUTC+'. A single supplied time is used directly.',
  'Coordinate system GSM. Absolute position km: '+JSON.stringify(record.positionKm),
  'Earth radius 6372 km follows the original plotting program.',
  'Position source: '+record.sourceType+'. '+record.sourceURL,
  'Archived original inputs: '+record.sources.join('\n'),
  'For SSC input, native GSE ephemeris is aligned only within bracketing samples using irf_resamp. Original IRFU routine performs the GSM conversion.',
  'Magnetopause and bow shock are model boundaries from the original program. '+record.modelNotes.join(' '),
  'Original SSCWeb service: https://sscweb.gsfc.nasa.gov/WebServices/REST/',
 ].join('\n'));
 orbitMap.set(ev.id,{title,kind:'orbit',event:ev.id,image:record.png,snapshotUTC:record.snapshotUTC,slideObject:s});
}
const index=[];
for(let i=0;i<oldIndex.length;i++){
 const item=oldIndex[i];
 if(item.event&&oldIndex.findIndex(x=>x.event===item.event)===i){
  const {slideObject,...orbit}=orbitMap.get(item.event);
  if(p.slides.items[index.length]!==slideObject)throw Error('Orbit slide order mismatch');
  index.push({...orbit,slide:index.length+1});
 }
 if(p.slides.items[index.length]!==originals[i])throw Error('Source slide order mismatch');
 index.push({...item,originalSlide:item.slide,slide:index.length+1});
}
if(p.slides.items.length!==127)throw Error('Expected 127 slides');
await fs.writeFile(path.join(AUDIT,'ppt_slide_index.json'),JSON.stringify(index,null,2));
console.log('EXPORT',p.slides.items.length);
const candidate=path.join(WORK,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
const finalLocal=path.join(WORK,'output',name);
await finalizePresentation({
 workspaceDir:WORK,candidatePath:candidate,finalPath:finalLocal,
 explicitTotalSlideCount:127,requiredNativeTableOwnerSlides:[],requiredNativeChartOwnerSlides:[],
 pythonExecutable:PYTHON,
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,9144000','--validate-heading-fit'],
 fontPolicy:{basis:'reference',families:[font],referencePath:source,referenceSha256:crypto.createHash('sha256').update(await fs.readFile(source)).digest('hex')},
 verifyArtifactToolImport:true,receiptPath:path.join(WORK,'validation.json')
});
await fs.copyFile(finalLocal,path.join(OUT,name));
console.log('FINALIZED',path.join(OUT,name));
