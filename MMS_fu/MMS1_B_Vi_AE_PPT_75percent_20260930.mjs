// 将现有PPT每页图片的宽、高缩至原尺寸75%，保留比例并居中，另存新文件。
// 原PPT保持不变；构建、布局及验证记录放TEMP，交付文件放Recovery目录。
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';

const sourcePath = String.raw`C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815.pptx`;
const deliveryPath = String.raw`C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815_75percent.pptx`;
const skillDir = String.raw`C:\Users\Administrator\.codex\plugins\cache\openai-primary-runtime\presentations\26.904.11930\skills\presentations`;
const pythonExecutable = String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`;
process.env.RUNTIME_NODE_MODULES = String.raw`C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\node_modules`;
const workspaceDir = process.env.SMILE_PPT_WORK;
if (!workspaceDir || !path.isAbsolute(workspaceDir)) throw new Error('SMILE_PPT_WORK必须为绝对路径');
try { await fs.access(deliveryPath); throw new Error('目标PPT已存在，避免覆盖'); }
catch (error) { if (error.code !== 'ENOENT') throw error; }

const stagingDir = path.join(workspaceDir, 'build');
const finalDir = path.join(workspaceDir, 'output');
await fs.mkdir(stagingDir, { recursive: true });
await fs.mkdir(finalDir, { recursive: true });
const sourceCopy = path.join(stagingDir, 'source.pptx');
await fs.copyFile(sourcePath, sourceCopy);
const sourceHash = crypto.createHash('sha256').update(await fs.readFile(sourceCopy)).digest('hex');
console.log('IMPORT', sourceCopy);
const presentation = await PresentationFile.importPptx(await FileBlob.load(sourceCopy));
const before = await presentation.inspect({ kind: 'slide,image,shape,textbox,table,chart,notes', maxChars: 1000000 });
await fs.writeFile(path.join(stagingDir, 'before.ndjson'), before.ndjson);
const records = before.ndjson.split(/\r?\n/).filter(line => line.trim()).map(line => JSON.parse(line));
const slides = records.filter(record => record.kind === 'slide');
const images = records.filter(record => record.kind === 'image');
if (slides.length !== 162 || images.length !== 162) throw new Error(`预期162页、162张图，实际${slides.length}页、${images.length}图`);
const width = 1700;
const height = 919;
const ratio = 0.75;
const manifest = [];
for (const record of images) {
  if (images.filter(image => image.slide === record.slide).length !== 1) throw new Error('每页必须恰好一张图片');
  const image = presentation.resolve(record.id);
  const beforeFrame = { ...image.frame };
  if (![beforeFrame.left, beforeFrame.top, beforeFrame.width, beforeFrame.height].every(Number.isFinite)) throw new Error('图片位置含无效数值');
  const afterFrame = {
    left: (width - beforeFrame.width * ratio) / 2,
    top: (height - beforeFrame.height * ratio) / 2,
    width: beforeFrame.width * ratio,
    height: beforeFrame.height * ratio,
  };
  image.frame = afterFrame;
  manifest.push({ slide: record.slide, imageId: record.id, beforeFrame, afterFrame });
}
const after = await presentation.inspect({ kind: 'slide,image,shape,textbox,table,chart,notes', maxChars: 1000000 });
await fs.writeFile(path.join(stagingDir, 'after.ndjson'), after.ndjson);
await fs.writeFile(path.join(stagingDir, 'resize_manifest.json'), JSON.stringify({ sourcePath, sourceHash, width, height, ratio, slides: manifest }, null, 2));
const candidatePath = path.join(stagingDir, 'candidate.pptx');
console.log('EXPORT', slides.length, 'slides, scale', ratio);
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);
const { finalizePresentation } = await import(pathToFileURL(path.join(skillDir, 'container_tools', 'artifact_tool_utils.mjs')).href);
const finalPath = path.join(finalDir, path.basename(deliveryPath));
await finalizePresentation({
  workspaceDir,
  candidatePath,
  finalPath,
  explicitTotalSlideCount: slides.length,
  requiredNativeTableOwnerSlides: [],
  requiredNativeChartOwnerSlides: [],
  pythonExecutable,
  integrityValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_package_integrity.py'),
  layoutValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_layout_geometry.py'),
  layoutArgs: ['--expected-slide-size-emu', `${width * 9525},${height * 9525}`],
  verifyArtifactToolImport: true,
  receiptPath: path.join(stagingDir, 'validation.json'),
});
if (crypto.createHash('sha256').update(await fs.readFile(sourcePath)).digest('hex') !== sourceHash) throw new Error('原PPT在处理期间有变化，请重新核对');
await fs.writeFile(path.join(stagingDir, 'delivery_paths.json'), JSON.stringify({ finalPath, deliveryPath, workspaceDir }, null, 2));
console.log('FINALIZED', finalPath);
