import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { Presentation, PresentationFile, FileBlob } from '@oai/artifact-tool';

const SKILL=String.raw`C:\Users\Administrator\.codex\plugins\cache\openai-primary-runtime\presentations\26.904.11930\skills\presentations`;
const CODE=String.raw`C:\Users\Administrator\Documents\FWD_matlab\MMS_fu`;
const AUDIT=String.raw`Z:\SPART-WORK\Data\MMS\derived\event_overview_20260924`;
const OLD_AUDIT=String.raw`Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923`;
const WORK=path.join(process.env.TEMP,'MMS_overview_20260924','ppt_build');
const OUT=String.raw`C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924`;
const PYTHON=String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`;
process.env.RUNTIME_NODE_MODULES=path.resolve(path.dirname(process.execPath),'..','node_modules');
const { finalizePresentation, resolvePresentationFont }=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const font=resolvePresentationFont({fontFamily:'Arial'});
const W=1280,H=960;
const PILOT=process.argv.includes('--pilot');
await fs.mkdir(WORK,{recursive:true});
await fs.mkdir(path.join(WORK,'output'),{recursive:true});
await fs.mkdir(path.join(WORK,'rendered'),{recursive:true});
await fs.mkdir(OUT,{recursive:true});
const events=JSON.parse(await fs.readFile(path.join(CODE,'MMS_event_overview_20260924_events.json'),'utf8'));
const refs=JSON.parse(await fs.readFile(path.join(OLD_AUDIT,'official_reference_index.json'),'utf8'));
const presentation=Presentation.create({slideSize:{width:W,height:H}});
const slideIndex=[];
function text(slide, value,x,y,w,h,size=22,bold=false,color='#183047'){
 const shape=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 shape.text=value;shape.text.style={typeface:font,fontSize:size,bold,color,autoFit:'none'};
 return shape;
}
function newSlide(title,kind,extra={}){
 const slide=presentation.slides.add();slide.background.fill='#FFFFFF';
 slideIndex.push({slide:slideIndex.length+1,title,kind,...extra});return slide;
}
function expanded(ev){
 const start=new Date(new Date(ev.start).getTime()-600000).toISOString().slice(0,19).replace('T',' ');
 const end=new Date(new Date(ev.end).getTime()+600000).toISOString().slice(0,19).replace('T',' ');
 return [start,end];
}
async function addImage(slide,file,x,y,w,h,alt){
 const bytes=await fs.readFile(file);
 slide.images.add({blob:new Uint8Array(bytes),contentType:'image/png',fit:'contain',position:{left:x,top:y,width:w,height:h},alt});
}
const cover=newSlide('MMS Event Overviews','cover');
text(cover,'MMS Event Overviews',65,65,1160,85,48,true);
text(cover,'20 events, MMS1–4',65,170,1160,65,32);
text(cover,'21 July – 16 September 2026',65,237,1160,55,28);
const lines=[
'Each event includes at least 10 minutes before and after the supplied interval.',
'All times are UTC. The 11 August 15:00 ±5 min event covers 14:45–15:15.',
'60 MATLAB figures read archived original CDFs through IRFU.',
'B, Vi, Ve, E, N, Ti, Te, Ei and Ee follow the labels in the original scripts.',
'B and V use GSM. MMS1–3 E uses GSM. MMS4 E retains native DSL XY.',
'Burst data take priority where available. Particle moments use available fast data.',
'11 August: only B could be redrawn from the required public L2 products.',
'28–29 August and 16 September: official quicklook references are included.',
'MMS4 electron data gaps remain blank.',
'Public-data availability was checked on 23 September 2026.',
'Event descriptions are the supplied candidate interpretations.'
];
lines.forEach((s,i)=>text(cover,s,65,337+i*43,1150,42,21,false,i>5?'#526778':'#183047'));
cover.speakerNotes.textFrame.setText('Source PDF: MMS_20_events_MMS1-4_overviews.pdf (107 pages). Updated MATLAB plots on 2026-09-24. Official data API: https://lasp.colorado.edu/mms/sdc/public/about/how-to/ . Availability evidence: '+OLD_AUDIT+'. The figures are embedded images. Regenerate and edit their plotted content using the supplied MATLAB source package.');
for(let block=0;block<2;block++){
 const s=newSlide('Event index '+(block*10+1)+'–'+(block*10+10),'index');
 text(s,'Event index '+(block*10+1)+'–'+(block*10+10),50,40,1180,60,36,true);
 text(s,'ID / date',50,123,210,35,20,true);
 text(s,'Supplied interval (UTC)',278,123,380,35,20,true);
 text(s,'Plotted context (UTC)',655,123,570,35,20,true);
 events.slice(block*10,block*10+10).forEach((e,i)=>{
  const y=180+i*72;const [a,b]=expanded(e);
  const core=e.start.slice(11,16)+(e.start===e.end?'':'–'+e.end.slice(11,16));
  const context=a.slice(0,10)===b.slice(0,10)?a.slice(11,16)+'–'+b.slice(11,16):a.slice(5,16)+' / '+b.slice(5,16);
  text(s,e.id+'  '+e.start.slice(0,10),50,y,220,29,19,true);
  text(s,core,278,y,355,29,20);
  text(s,context,655,y,570,29,20);
  text(s,e.label,278,y+29,920,39,18,false,'#526778');
 });
 s.speakerNotes.textFrame.setText('User-supplied event list. Times are UTC. At least 10 minutes of context at both ends. Original event notes retain their uncertainty. Source: MMS_event_overview_20260924_events.json.');
}
for(const ev of events){
 for(let sc=1;sc<=4;sc++){
  const id=Number(ev.id.slice(2));
  if(PILOT&&!((id===4&&sc===1)||(id===6&&sc===4)||(id===14&&sc===1)||(id===16&&sc===3)))continue;
  if(id<=15){
   const base=ev.id+'_'+ev.start.slice(0,10).replaceAll('-','')+'_MMS'+sc+'_overview';
   const record=JSON.parse(await fs.readFile(PILOT?path.join(process.env.TEMP,'MMS_overview_20260924',base+'_preview.json'):path.join(AUDIT,base+'.json'),'utf8'));
   if(!record.complete||record.readErrors.length||record.renderVersion!==4) throw new Error('Unfinished plot: '+base);
   const slide=newSlide(ev.id+' MMS'+sc,'native',{event:ev.id,spacecraft:sc,image:record.png});
   await addImage(slide,record.png,8,8,W-16,H-16,ev.id+' MMS'+sc+' nine tightly stacked scientific panels, no mode timeline');
   const coverage=record.statistics.coverage;
   slide.speakerNotes.textFrame.setText([
    ev.id+' MMS'+sc+' '+ev.label,
    'Plot interval UTC: '+record.startUTC+' to '+record.endUTC,
    'Source: original MMS L2 CDFs in '+String.raw`Z:\SPART-WORK\Data\MMS`,
    'Original plotting blocks: Overview_download.m and Overview_download_mms4.m.',
    'Run entry: run_MMS_event_overviews_20260924.m. Plot function: Overview_download_events_20260924.m.',
    'No smoothing or interpolation across gaps. Display compression retains original extrema.',
    'Coordinates: '+record.coordinateSystem,
    'Native data coverage (modes retained here, removed from the plot): '+JSON.stringify(coverage),
    'CDF sources: '+record.sources.join('\n')
   ].join('\n'));
  }
  if(id>=13){
   const ref=refs.find(r=>r.event===ev.id&&r.spacecraft===sc);
   const groups=new Map();
   for(const q of ref.quicklook){
    const key=q.startUTC+'|'+q.endUTC;
    if(!groups.has(key))groups.set(key,[]);
    groups.get(key).push(q);
   }
   for(const [key,items] of [...groups].sort((a,b)=>a[0].localeCompare(b[0]))){
    items.sort((a,b)=>a.plot.localeCompare(b.plot));
    const s=newSlide(ev.id+' MMS'+sc+' official reference','reference',{event:ev.id,spacecraft:sc,interval:key,images:items.map(x=>x.path)});
    const [a,b]=expanded(ev);
    text(s,ev.id+'  MMS'+sc+'  Official quicklook reference',38,22,1200,50,30,true);
    text(s,'Requested context: '+a+' to '+b+' UTC',38,77,1200,34,19);
    text(s,'Original image interval: '+key.replace('|',' to ').replaceAll('T',' ').replaceAll('Z','')+' UTC',38,111,1200,32,18);
    text(s,'Original coordinates and labels are retained (DMPA / DBCS / DSL).',38,145,1200,30,18,false,'#526778');
    if(items.length===1) await addImage(s,items[0].path,35,183,1210,727,items[0].plot);
    else for(let j=0;j<items.length;j++)await addImage(s,items[j].path,30+j*625,183,595,727,items[j].plot);
    text(s,'MMS Science Data Center official reference images',38,924,1200,24,16,false,'#526778');
    s.speakerNotes.textFrame.setText('Official source images, kept unchanged. These pages supplement unavailable L2 quantities. Source: https://lasp.colorado.edu/mms/sdc/public/quicklook/ .\n'+JSON.stringify(items,null,2));
   }
  }
 }
}
if(!PILOT&&presentation.slides.items.length!==107)throw new Error('Expected original PDF coverage of 107 slides, found '+presentation.slides.items.length);
await fs.writeFile(path.join(PILOT?WORK:AUDIT,PILOT?'pilot_slide_index.json':'ppt_slide_index.json'),JSON.stringify(slideIndex,null,2));
const candidate=path.join(WORK,PILOT?'pilot_candidate.pptx':'candidate.pptx');
console.log('EXPORT',presentation.slides.items.length);
await (await PresentationFile.exportPptx(presentation)).save(candidate);
console.log('EXPORTED',candidate);
const finalLocal=path.join(WORK,'output',PILOT?'pilot.pptx':'MMS_20_events_MMS1-4_overviews_20260924.pptx');
await finalizePresentation({
 workspaceDir:WORK,candidatePath:candidate,finalPath:finalLocal,
 explicitTotalSlideCount:presentation.slides.items.length,requiredNativeTableOwnerSlides:[],requiredNativeChartOwnerSlides:[],
 pythonExecutable:PYTHON,
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu',String(W*9525)+','+String(H*9525),'--validate-heading-fit'],
 fontPolicy:{basis:'design',families:[font]},
 verifyArtifactToolImport:true,receiptPath:path.join(WORK,PILOT?'pilot_validation.json':'validation.json')
});
console.log('FINALIZED');
const finalPresentation=await PresentationFile.importPptx(await FileBlob.load(finalLocal));
for(let i=0;i<finalPresentation.slides.items.length;i++){
 const slide=finalPresentation.slides.items[i];
 const png=await finalPresentation.export({slide,format:'png',scale:1});
 await fs.writeFile(path.join(WORK,'rendered',(PILOT?'pilot-':'slide-')+String(i+1).padStart(3,'0')+'.png'),new Uint8Array(await png.arrayBuffer()));
 if((i+1)%10===0) console.log('RENDERED',i+1);
}
if(!PILOT)await fs.copyFile(finalLocal,path.join(OUT,path.basename(finalLocal)));
console.log('DELIVERED',path.join(OUT,path.basename(finalLocal)));
