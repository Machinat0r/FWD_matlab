// 仅更换26张新增离子流速图；最终包保留原PPT的其余内部文件。
import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';

const workspaceDir = path.join(process.env.TEMP, 'MMS1_BViAE_refresh_20260930');
const build = path.join(workspaceDir, 'build');
const skillDir = String.raw`C:\Users\Administrator\.codex\plugins\cache\openai-primary-runtime\presentations\26.904.11930\skills\presentations`;
const runtime = String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies`;
process.env.RUNTIME_NODE_MODULES = path.join(runtime, 'node', 'node_modules');
const manifest = JSON.parse(await fs.readFile(path.join(build, 'replacement_manifest.json'), 'utf8'));
const phase = process.argv[2] ?? 'inspect';

if (phase === 'finalize') {
  const { finalizePresentation } = await import(pathToFileURL(path.join(skillDir, 'container_tools', 'artifact_tool_utils.mjs')).href);
  const result = await finalizePresentation({
    workspaceDir,
    candidatePath: path.join(build, 'candidate.pptx'),
    finalPath: path.join(workspaceDir, 'output', 'MMS1_B_Vi_AE_20260720_0815_75percent.pptx'),
    explicitTotalSlideCount: 162,
    requiredNativeTableOwnerSlides: [],
    requiredNativeChartOwnerSlides: [],
    pythonExecutable: path.join(runtime, 'python', 'python.exe'),
    integrityValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_package_integrity.py'),
    layoutValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_layout_geometry.py'),
    layoutArgs: ['--expected-slide-size-emu', '16192500,8753475'],
    verifyArtifactToolImport: true,
    receiptPath: path.join(build, 'validation.json'),
  });
  console.log(JSON.stringify(result));
} else {
  const presentation = await PresentationFile.importPptx(await FileBlob.load(manifest.source_path));
  const inspected = await presentation.inspect({kind:'slide,image,shape,textbox,notes,layout', maxChars:2000000});
  const ndjson = typeof inspected === 'string' ? inspected : inspected.ndjson;
  if (typeof ndjson !== 'string') throw new Error('Unexpected inspect result');
  await fs.writeFile(path.join(build, 'artifact_before.ndjson'), ndjson);
  const records = ndjson.split('\n').filter(Boolean).map(line => JSON.parse(line));
  if (phase === 'inspect') {
    console.log(JSON.stringify(records.filter(r => r.slide === 128)));
  } else if (phase === 'edit') {
    const changes = [];
    for (const item of manifest.matches) {
      const anchors = records.filter(r => r.slide === item.slide && r.kind === 'image');
      if (anchors.length !== 1) throw new Error(`Slide ${item.slide}: ${anchors.length} images`);
      const image = presentation.resolve(anchors[0].id);
      const old = {};
      for (const key of ['frame','crop','fit','alt','prompt','geometry','borderRadius','rotation','flipHorizontal','flipVertical','lockAspectRatio']) old[key] = image[key];
      image.replace({blob: new Uint8Array(await fs.readFile(item.new_png)), contentType:'image/png', alt:old.alt ?? '', ...(old.fit ? {fit:old.fit} : {}), ...(old.prompt ? {prompt:old.prompt} : {})});
      for (const key of ['frame','crop','geometry','borderRadius','rotation','flipHorizontal','flipVertical','lockAspectRatio']) image[key] = old[key];
      changes.push({slide:item.slide, window_id:item.window_id, anchor:anchors[0].id, preserved:old});
    }
    await fs.writeFile(path.join(build,'artifact_replacement_geometry.json'), JSON.stringify(changes,null,2));
    await (await PresentationFile.exportPptx(presentation)).save(path.join(build,'artifact_candidate.pptx'));
    console.log(JSON.stringify({replaced:changes.length, exported:path.join(build,'artifact_candidate.pptx')}));
  } else throw new Error(`Unknown phase ${phase}`);
}
