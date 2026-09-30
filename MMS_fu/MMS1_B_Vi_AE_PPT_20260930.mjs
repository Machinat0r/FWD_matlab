// 将用户指定目录中的原始PNG依序嵌入PPT，每页仅一张图。
// 代码保存在MMS_fu，构建及验证文件放TEMP，交付PPT放Recovery-Work_SMILE-MMS。
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { Presentation, PresentationFile } from '@oai/artifact-tool';

const inputDir = String.raw`C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\B_Vi_AE`;
const deliveryPath = String.raw`C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815.pptx`;
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
const files = (await fs.readdir(inputDir, { withFileTypes: true }))
  .filter(entry => entry.isFile() && /^W\d{3}_\d{14}_MMS1_B_Vi_AE\.png$/i.test(entry.name))
  .map(entry => entry.name).sort();
if (files.length !== 162) throw new Error(`预期162张图片，实际${files.length}`);
for (let index = 0; index < files.length; index++) {
  if (!files[index].startsWith(`W${String(index + 1).padStart(3, '0')}_`)) {
    throw new Error('图片窗口序号不连续');
  }
}
const width = 1700;
const height = 919;
const presentation = Presentation.create({ slideSize: { width, height } });
const manifest = [];
for (const [index, name] of files.entries()) {
  const file = path.join(inputDir, name);
  const bytes = await fs.readFile(file);
  if (bytes.subarray(0, 8).toString('hex') !== '89504e470d0a1a0a') throw new Error(`非PNG: ${name}`);
  const imageWidth = bytes.readUInt32BE(16);
  const imageHeight = bytes.readUInt32BE(20);
  const scale = Math.min(width / imageWidth, height / imageHeight);
  const position = {
    left: (width - imageWidth * scale) / 2,
    top: (height - imageHeight * scale) / 2,
    width: imageWidth * scale,
    height: imageHeight * scale,
  };
  const slide = presentation.slides.add();
  slide.background.fill = '#FFFFFF';
  slide.images.add({ blob: new Uint8Array(bytes), contentType: 'image/png', alt: '', fit: 'contain', position });
  // 用户明确只要图片，不创建文字框、封面、页码、备注或评论。
  manifest.push({ slide: index + 1, source: file, sha256: crypto.createHash('sha256').update(bytes).digest('hex'), imageWidth, imageHeight, position });
}
await fs.writeFile(path.join(stagingDir, 'image_manifest.json'), JSON.stringify({ width, height, slides: manifest }, null, 2));
const candidatePath = path.join(stagingDir, 'candidate.pptx');
console.log('EXPORT', files.length, 'slides');
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);
const { finalizePresentation } = await import(pathToFileURL(path.join(skillDir, 'container_tools', 'artifact_tool_utils.mjs')).href);
const finalPath = path.join(finalDir, path.basename(deliveryPath));
await finalizePresentation({
  workspaceDir,
  candidatePath,
  finalPath,
  explicitTotalSlideCount: files.length,
  requiredNativeTableOwnerSlides: [],
  requiredNativeChartOwnerSlides: [],
  pythonExecutable,
  integrityValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_package_integrity.py'),
  layoutValidatorPath: path.join(skillDir, 'container_tools', 'inspect_presentation_layout_geometry.py'),
  layoutArgs: ['--expected-slide-size-emu', `${width * 9525},${height * 9525}`],
  verifyArtifactToolImport: true,
  receiptPath: path.join(stagingDir, 'validation.json'),
});
await fs.writeFile(path.join(stagingDir, 'delivery_paths.json'), JSON.stringify({ finalPath, deliveryPath, workspaceDir }, null, 2));
console.log('FINALIZED', finalPath);
