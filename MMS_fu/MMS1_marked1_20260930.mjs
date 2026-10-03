import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';
const work=path.join(process.env.TEMP,'MMS1_marked1_20260930');
const build=path.join(work,'build');
const manifest=JSON.parse(await fs.readFile(path.join(build,'manifest.json'),'utf8'));
const skill=String.raw`C:\Users\Administrator\.codex\plugins\cache\openai-primary-runtime\presentations\26.904.11930\skills\presentations`;
const runtime=String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies`;
process.env.RUNTIME_NODE_MODULES=path.join(runtime,'node','node_modules');
if(process.argv[2]==='finalize'){
  const {finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools','artifact_tool_utils.mjs')).href);
  await fs.mkdir(path.join(work,'output'),{recursive:true});
  const result=await finalizePresentation({workspaceDir:work,candidatePath:path.join(build,'candidate.pptx'),finalPath:path.join(work,'output','MMS1_marked1.pptx'),
    explicitTotalSlideCount:manifest.selected_pages.length,requiredNativeTableOwnerSlides:[],requiredNativeChartOwnerSlides:[],
    pythonExecutable:path.join(runtime,'python','python.exe'),integrityValidatorPath:path.join(skill,'container_tools','inspect_presentation_package_integrity.py'),
    layoutValidatorPath:path.join(skill,'container_tools','inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','16192500,8753475'],
    verifyArtifactToolImport:true,receiptPath:path.join(build,'validation.json')});
  console.log(JSON.stringify({finalPath:result.finalPath,sha256:result.finalSha256,bytes:result.byteCount}));
}else{
  const p=await PresentationFile.importPptx(await FileBlob.load(manifest.snapshot));
  const state=await p.inspect({kind:'slide,image,shape,textbox,layout',maxChars:1000000});
  await fs.writeFile(path.join(build,'artifact_before.ndjson'),state.ndjson);
  if(process.argv[2]==='inspect'){
    console.log(JSON.stringify({slidesMethods:Object.getOwnPropertyNames(Object.getPrototypeOf(p.slides)),slideMethods:Object.getOwnPropertyNames(Object.getPrototypeOf(p.slides.getItem(0)))}));
  }else if(process.argv[2]==='edit'){
    const selected=new Set(manifest.selected_pages);
    for(let index=manifest.source_slide_count-1;index>=0;index--) if(!selected.has(index+1))p.slides.getItem(index).delete();
    if(p.slides.count!==selected.size)throw new Error('Wrong selected slide count');
    await (await PresentationFile.exportPptx(p)).save(path.join(build,'artifact_candidate.pptx'));
    console.log(JSON.stringify({selected:manifest.selected_pages,slideCount:p.slides.count}));
  }
}
