function out = plot_ace_20260119_event(aceRoot, outDir, varargin)
%PLOT_ACE_20260119_EVENT Plot ACE overview around the 2026-01-19 B>100 nT event.
%
% Usage:
%   plot_ace_20260119_event
%   plot_ace_20260119_event('Z:\SPART-WORK\Data\ACE\ace')
%   plot_ace_20260119_event('Z:\SPART-WORK\Data\ACE\ace', 'D:\figures')
%   plot_ace_20260119_event('Z:\SPART-WORK\Data\ACE', [], ...
%       'DownloaderDir', 'C:\Users\Administrator\Documents\FWD_matlab\ACE_fu')
%
% The script first checks the local ACE data library.  If the required daily
% CDF files are missing, it calls ACEFilesDownload_NAS.m from DownloaderDir
% to download them.  It prefers the highest-cadence local CDF it can find and
% then falls back to lower-cadence ACE products.  If irf_matlab is on path,
% irf_plot is used for the time series panels.

if nargin < 1 || isempty(aceRoot)
    aceRoot = 'Z:\SPART-WORK\Data\ACE\ace';
end
if nargin < 2 || isempty(outDir)
    outDir = fileparts(mfilename('fullpath'));
end

parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DownloadMissing', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'DownloaderDir', ...
    'C:\Users\Administrator\Documents\FWD_matlab\ACE_fu', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Threads', 8, @isnumeric);
addParameter(parser, 'CheckSize', 1, @isnumeric);
addParameter(parser, 'UseIRF', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'IrfPath', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PlotBackend', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PlotStartUTC', '2026-01-19 12:00:00', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(parser, 'PlotEndUTC', '2026-01-20 00:00:00', @(x) ischar(x) || isstring(x) || isnumeric(x));
parse(parser, varargin{:});
opt = parser.Results;

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

day0 = datenum(2026, 1, 19, 0, 0, 0);
day1 = day0 + 1;
plot0 = parseUtcTime(opt.PlotStartUTC);
plot1 = parseUtcTime(opt.PlotEndUTC);
eventHint = datenum(2026, 1, 19, 18, 59, 38.687);

[aceDataRoot, aceSearchRoot] = normalizeAceRoot(char(aceRoot));
if ~isempty(char(opt.IrfPath))
    addpath(genpath(char(opt.IrfPath)));
end

downloaderDir = char(opt.DownloaderDir);
pythonScript = char(opt.PythonScript);

useIrfPlot = logical(opt.UseIRF) && exist('irf_plot', 'file') == 2;

fprintf('ACE data root  : %s\n', aceDataRoot);
fprintf('ACE search root: %s\n', aceSearchRoot);
fprintf('UTC day : %s\n', datestr(day0, 'yyyy-mm-dd'));
if shouldUsePillowBackend(opt.PlotBackend)
    fprintf('Plotting      : headless Pillow PNG backend\n');
elseif useIrfPlot
    fprintf('Plotting      : irf_plot (%s)\n', which('irf_plot'));
else
    fprintf('Plotting      : MATLAB native plot\n');
end

files.mfi = ensureAceDailyCdf('MFI magnetic field', aceSearchRoot, day0, ...
    {'mfi_h3', 'mfi_h0', 'mfi_h2', 'mfi_h1', 'mfi_k0'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);
files.swe = ensureAceDailyCdf('SWEPAM plasma', aceSearchRoot, day0, ...
    {'swe_k0', 'swe_k1', 'swe_h0', 'swe_h2'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);
files.epm = ensureAceDailyCdf('EPAM particle flux', aceSearchRoot, day0, ...
    {'epm_h3', 'epm_h2', 'epm_h1', 'epm_k0', 'epm_k1'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);

dispSelectedFile('MFI magnetic field', files.mfi);
dispSelectedFile('SWEPAM plasma', files.swe);
dispSelectedFile('EPAM particle flux', files.epm);

fprintf('Reading MFI data...\n');
mfi = loadMfiData(files.mfi, day0, day1);
fprintf('Reading MFI data...done (%d samples)\n', numel(mfi.t));
if isempty(mfi.t)
    error('No usable ACE/MFI magnetic-field data were found for 2026-01-19.');
end

fprintf('Reading SWEPAM data...\n');
swe = loadSweData(files.swe, day0, day1);
fprintf('Reading SWEPAM data...done (%d samples)\n', numel(swe.t));
fprintf('Reading EPAM flux data...\n');
[eflux, iflux] = loadParticleFlux(files.epm, files.swe, plot0, plot1);
fprintf('Reading EPAM flux data...done (e:%d, i:%d samples)\n', numel(eflux.t), numel(iflux.t));

if isempty(swe.t)
    fprintf('Ion velocity: no usable SWEPAM data found.\n');
elseif ~isempty(swe.vgse)
    fprintf('Ion velocity: using GSE vector %s from %s.\n', swe.vVar, sourceName(swe));
elseif ~isempty(swe.vp)
    fprintf('Ion velocity: no V_GSE components are available; using scalar %s from %s.\n', ...
        swe.vpVar, sourceName(swe));
else
    fprintf('Ion velocity: SWEPAM file found, but no usable velocity variable was found.\n');
end

fprintf('Computing Bmax marker...\n');
bSearch = mfi.t >= plot0 & mfi.t <= plot1;
[bMax, imaxRel] = max(mfi.bmag(bSearch), [], 'omitnan');
if isempty(imaxRel) || isnan(bMax)
    bMax = NaN;
    tBmax = eventHint;
else
    idxB = find(bSearch);
    tBmax = mfi.t(idxB(imaxRel));
end
fprintf('Computing Bmax marker...done (%s, %.6g nT)\n', fmtTime(tBmax), bMax);

out.figure = [];
out.png = fullfile(outDir, 'ace_20260119_overview.png');
out.fig = fullfile(outDir, 'ace_20260119_overview.fig');
out.bmax_time_utc = tBmax;
out.bmax_nt = bMax;
out.files = files;
out.plot_start_utc = plot0;
out.plot_end_utc = plot1;
out.velocity_note = swe.velocityNote;

if shouldUsePillowBackend(opt.PlotBackend)
    out.fig = '';
    if exist(out.png, 'file')
        fprintf('Using existing headless PNG: %s\n', out.png);
    else
        fprintf('Headless PNG is missing; preparing helper command...\n');
        renderWithPillow(out.png, plot0, plot1, tBmax, files, day0);
    end
    fprintf('Saved PNG: %s\n', out.png);
    fprintf('Bmax UTC : %s, %.6g nT\n', fmtTime(tBmax), bMax);
    return
end

fprintf('Creating figure...\n');
fig = figure('Color', 'w', 'Name', 'ACE 2026-01-19 overview', ...
    'Units', 'pixels', 'Position', [80 80 980 1040]);
tl = tiledlayout(fig, 5, 1, 'TileSpacing', 'none', 'Padding', 'compact');
fprintf('Creating figure...done\n');

ax = gobjects(5, 1);

ax(1) = nexttile(tl);
plotTimeData(ax(1), mfi.t, [mfi.bgse, mfi.bmag], [componentColors(); 0 0 0], useIrfPlot, false);
ylabel(ax(1), 'B [nT]');
grid(ax(1), 'on');
addPanelLabel(ax(1), 'a');
addInlineLegend(ax(1), {'B_x', 'B_y', 'B_z', '|B|'}, [componentColors(); 0 0 0]);

ax(2) = nexttile(tl);
if ~isempty(swe.t) && ~isempty(swe.vgse)
    plotTimeData(ax(2), swe.t, swe.vgse, componentColors(), useIrfPlot, false);
    addZeroLine(ax(2), useIrfPlot);
    addInlineLegend(ax(2), {'V_{ix}', 'V_{iy}', 'V_{iz}'}, componentColors());
    ylabel(ax(2), 'V_i [km/s]');
elseif ~isempty(swe.t) && ~isempty(swe.vp)
    plotTimeData(ax(2), swe.t, swe.vp, [0.05 0.25 0.85], useIrfPlot, false);
    addZeroLine(ax(2), useIrfPlot);
    addInlineLegend(ax(2), {'V_p'}, [0.05 0.25 0.85]);
    ylabel(ax(2), 'V_p [km/s]');
else
    putNoDataText(ax(2), 'No usable ion velocity file/variable found', plot0, plot1, useIrfPlot);
    ylabel(ax(2), 'V_i [km/s]');
end
grid(ax(2), 'on');
addPanelLabel(ax(2), 'b');

ax(3) = nexttile(tl);
if ~isempty(swe.t) && ~isempty(swe.ti)
    plotTimeData(ax(3), swe.t, swe.ti, [0.85 0.10 0.10], useIrfPlot, false);
    addInlineLegend(ax(3), {'T_i'}, [0.85 0.10 0.10]);
else
    putNoDataText(ax(3), 'No usable ion temperature file/variable found', plot0, plot1, useIrfPlot);
end
ylabel(ax(3), 'T_i [eV]');
grid(ax(3), 'on');
addPanelLabel(ax(3), 'c');

ax(4) = nexttile(tl);
if ~isempty(eflux.t) && ~isempty(eflux.y)
    plotTimeData(ax(4), eflux.t, eflux.y, fluxColors(size(eflux.y, 2)), useIrfPlot, true);
    legend(ax(4), eflux.labels, 'Location', 'northeast', ...
        'Interpreter', 'tex', 'Box', 'off', 'FontSize', 9);
else
    putNoDataText(ax(4), 'No usable electron flux file/variable found', plot0, plot1, useIrfPlot);
end
ylabel(ax(4), 'J_e');
grid(ax(4), 'on');
addPanelLabel(ax(4), 'd');

ax(5) = nexttile(tl);
if ~isempty(iflux.t) && ~isempty(iflux.y)
    plotTimeData(ax(5), iflux.t, iflux.y, fluxColors(size(iflux.y, 2)), useIrfPlot, true);
    legend(ax(5), iflux.labels, 'Location', 'northeast', ...
        'Interpreter', 'tex', 'Box', 'off', 'FontSize', 9);
else
    putNoDataText(ax(5), 'No usable ion/proton flux file/variable found', plot0, plot1, useIrfPlot);
end
ylabel(ax(5), 'J_i');
xlabel(ax(5), '2026-01-19/20 UTC');
grid(ax(5), 'on');
addPanelLabel(ax(5), 'e');

for ia = 1:numel(ax)
    setTimeAxisLimits(ax(ia), plot0, plot1, useIrfPlot);
    addEventLine(ax(ia), tBmax, useIrfPlot);
end
linkaxes(ax, 'x');
formatTimeAxis(ax, useIrfPlot);
stylePanels(ax);

out.figure = fig;
out.png = fullfile(outDir, 'ace_20260119_overview.png');
out.fig = fullfile(outDir, 'ace_20260119_overview.fig');
out.bmax_time_utc = tBmax;
out.bmax_nt = bMax;
out.files = files;
out.plot_start_utc = plot0;
out.plot_end_utc = plot1;
out.velocity_note = swe.velocityNote;

savePng(fig, out.png);
savefig(fig, out.fig);

fprintf('Saved PNG: %s\n', out.png);
fprintf('Saved FIG: %s\n', out.fig);
fprintf('Bmax UTC : %s, %.6g nT\n', fmtTime(tBmax), bMax);
end

function tf = shouldUsePillowBackend(backend)
backend = lower(char(backend));
if strcmp(backend, 'pillow') || strcmp(backend, 'python')
    tf = true;
elseif strcmp(backend, 'matlab')
    tf = false;
else
    tf = ~usejava('desktop');
end
end

function renderWithPillow(outPng, plot0, plot1, tBmax, files, day0)
helper = fullfile(fileparts(mfilename('fullpath')), 'plot_ace_overview_pillow.py');
if ~exist(helper, 'file')
    error('Missing Pillow plotting helper: %s', helper);
end

mfiMat = hapiCacheMatFile(mfiHapiDataset(files.mfi), {'Magnitude', 'BGSEc'}, day0, day0 + 1);
sweDataset = sweHapiDataset(files.swe);
if any(strcmp(sweDataset, {'AC_H0_SWE', 'AC_H2_SWE'}))
    sweMat = hapiCacheMatFile(sweDataset, {'Vp', 'Tpr', 'V_GSE'}, day0, day0 + 1);
else
    sweMat = hapiCacheMatFile(sweDataset, {'Vp', 'Tpr'}, day0, day0 + 1);
end
epmMat = hapiCacheMatFile(epmHapiDataset(files.epm), ...
    {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'DE1', 'DE2', 'DE3', 'DE4'}, ...
    plot0, plot1);

required = {mfiMat, sweMat, epmMat};
for i = 1:numel(required)
    if ~exist(required{i}, 'file')
        error('Missing MAT cache for headless plotting: %s', required{i});
    end
end

py = resolvePythonCommandForHapi();
cmd = [py, ' ', shellquote(helper), ' ', shellquote(outPng), ' ', ...
    sprintf('%.17g %.17g %.17g ', plot0, plot1, tBmax), ...
    shellquote(mfiMat), ' ', shellquote(sweMat), ' ', shellquote(epmMat)];
error('MATLAB batch cannot safely launch the headless renderer on this machine. Run this command first, then rerun MATLAB: %s', cmd);
end

function dispSelectedFile(label, file)
if isempty(file)
    fprintf('%-22s: not found\n', label);
else
    fprintf('%-22s: %s\n', label, file);
end
end

function s = sourceName(d)
if isfield(d, 'file') && ~isempty(d.file)
    [~, n, e] = fileparts(d.file);
    s = [n e];
elseif isfield(d, 'source') && ~isempty(d.source)
    s = d.source;
else
    s = 'unknown source';
end
end

function [aceDataRoot, aceSearchRoot] = normalizeAceRoot(aceRoot)
aceRoot = stripTrailingFilesep(char(aceRoot));

if exist(fullfile(aceRoot, 'mfi'), 'dir') || ...
        exist(fullfile(aceRoot, 'swepam'), 'dir') || ...
        exist(fullfile(aceRoot, 'epam'), 'dir')
    [aceDataRoot, ~] = fileparts(aceRoot);
    aceSearchRoot = aceRoot;
elseif exist(fullfile(aceRoot, 'ace'), 'dir')
    aceDataRoot = aceRoot;
    aceSearchRoot = fullfile(aceRoot, 'ace');
else
    [parentDir, lastPart] = fileparts(aceRoot);
    if strcmp(lastPart, 'ace')
    aceDataRoot = parentDir;
    aceSearchRoot = aceRoot;
    else
        aceDataRoot = aceRoot;
        aceSearchRoot = fullfile(aceRoot, 'ace');
    end
end
end

function p = stripTrailingFilesep(p)
p = char(p);
while numel(p) > 3 && (p(end) == filesep || p(end) == '/' || p(end) == '\')
    p(end) = [];
end
end

function [downloaderDir, pythonScript] = ensureAceDownloaderOnPath(downloaderDir, pythonScript)
downloaderDir = char(downloaderDir);
pythonScript = char(pythonScript);
if isempty(downloaderDir)
    downloaderDir = 'C:\Users\Administrator\Documents\FWD_matlab\ACE_fu';
end
if ~exist(downloaderDir, 'dir')
    error('DownloaderDir does not exist: %s', downloaderDir);
end
if isempty(pythonScript)
    pythonScript = fullfile(downloaderDir, 'download_ace_files_new.py');
end
if ~exist(pythonScript, 'file')
    error('Cannot find ACE Python downloader: %s', pythonScript);
end
if exist('ACEFilesDownload_NAS', 'file') ~= 2
    addpath(downloaderDir, '-end');
end
if exist('ACEFilesDownload_NAS', 'file') ~= 2
    error('Cannot find ACEFilesDownload_NAS.m in DownloaderDir or MATLAB path.');
end
end

function file = ensureAceDailyCdf(label, aceSearchRoot, day0, productKeys, ...
    downloadMissing, downloaderDir, pythonScript, pythonExe, checkSize, threads)
file = findAceDailyProductCdf(aceSearchRoot, productKeys, day0);
if ~isempty(file)
    fprintf('%-22s: found locally: %s\n', label, file);
    return
end

fprintf('%-22s: missing locally\n', label);
if ~downloadMissing
    return
end

[~, pythonScript] = ensureAceDownloaderOnPath(downloaderDir, pythonScript);

dateRange = [datestr(day0, 'yyyy-mm-dd'), '/', datestr(day0, 'yyyy-mm-dd')];
for ip = 1:numel(productKeys)
    productKey = productKeys{ip};
    try
        p = aceProductInfo(productKey);
        localDir = fullfile(aceSearchRoot, p.instrumentDir, p.level, ...
            'l2', datestr(day0, 'yyyy'), datestr(day0, 'mm'));
        if ~exist(localDir, 'dir')
            mkdir(localDir);
        end
        fprintf('  downloading %-8s -> %s\n', productKey, localDir);
        ACEFilesDownload_NAS(dateRange, localDir, ...
            'product', productKey, ...
            'Threads', threads, ...
            'CheckSize', double(checkSize), ...
            'KeepTree', 0, ...
            'PythonScript', pythonScript, ...
            'PythonExe', pythonExe);
        file = findAceDailyProductCdf(aceSearchRoot, {productKey}, day0);
        if ~isempty(file)
            fprintf('  %s ready: %s\n', productKey, file);
            return
        end
    catch ME
        fprintf('  %s unavailable: %s\n', productKey, ME.message);
    end
end

file = findAceDailyProductCdf(aceSearchRoot, productKeys, day0);
end

function file = findAceDailyProductCdf(aceSearchRoot, productKeys, day0)
file = '';
if isempty(aceSearchRoot)
    return
end
for ip = 1:numel(productKeys)
    p = aceProductInfo(productKeys{ip});
    file = knownEventFile(aceSearchRoot, productKeys{ip}, day0);
    if ~isempty(file)
        return
    end
    directDir = fullfile(aceSearchRoot, p.instrumentDir, p.level, 'l2', ...
        datestr(day0, 'yyyy'), datestr(day0, 'mm'));
    pat = sprintf('ac_%s_%s_%s_v*.cdf', p.level, p.fileProduct, ...
        datestr(day0, 'yyyymmdd'));
    file = findFileWithPython(directDir, pat);
    if ~isempty(file)
        return
    end
end
end

function file = knownEventFile(aceSearchRoot, productKey, day0)
file = '';
if abs(day0 - datenum(2026, 1, 19, 0, 0, 0)) > 1e-9
    return
end
switch lower(char(productKey))
    case 'mfi_h3'
        file = fullfile(aceSearchRoot, 'mfi', 'h3', 'l2', '2026', '01', ...
            'ac_h3_mfi_20260119_v03.cdf');
    case 'swe_k0'
        file = fullfile(aceSearchRoot, 'swepam', 'k0', 'l2', '2026', '01', ...
            'ac_k0_swe_20260119_v01.cdf');
    case 'epm_h3'
        file = fullfile(aceSearchRoot, 'epam', 'h3', 'l2', '2026', '01', ...
            'ac_h3_epm_20260119_v01.cdf');
end
end

function file = findFileWithPython(folder, pattern)
file = '';
py = resolvePythonCommandForHapi();
code = ['import glob, os, sys; ', ...
    'files=glob.glob(os.path.join(sys.argv[1], sys.argv[2])); ', ...
    'files.sort(reverse=True); ', ...
    'print(files[0] if files else "")'];
cmd = [py, ' -c ', shellquote(code), ' ', shellquote(folder), ' ', shellquote(pattern)];
[status, out] = system(cmd);
if status == 0
    file = strtrim(out);
end
end

function p = aceProductInfo(productKey)
tok = regexp(char(productKey), '^([^_]+)_([^_]+)$', 'tokens', 'once');
if isempty(tok)
    error('Unsupported ACE product key: %s', productKey);
end
prefix = lower(tok{1});
p.level = lower(tok{2});
switch prefix
    case 'mfi'
        p.instrumentDir = 'mfi';
        p.fileProduct = 'mfi';
    case 'swe'
        p.instrumentDir = 'swepam';
        p.fileProduct = 'swe';
    case 'epm'
        p.instrumentDir = 'epam';
        p.fileProduct = 'epm';
    otherwise
        error('Unsupported ACE product prefix in %s', productKey);
end
end

function idx = selectLatestVersion(hits)
ver = zeros(numel(hits), 1);
for i = 1:numel(hits)
    tok = regexp(hits(i).name, '_v(\d+)\.cdf$', 'tokens', 'once');
    if ~isempty(tok)
        ver(i) = str2double(tok{1});
    end
end
[~, idx] = max(ver + [hits.datenum]' * 1e-12);
end

function d = emptyMfi()
d = struct('t', [], 'bmag', [], 'bgse', [], ...
    'file', '', 'bVar', '', 'bmagVar', '');
end

function d = loadMfiData(file, day0, day1)
d = emptyMfi();
if isempty(file)
    return
end
d.file = file;
try
    dataset = mfiHapiDataset(file);
    [t, y] = readHapiNumeric(dataset, {'Magnitude', 'BGSEc'}, day0, day1);
    bmag = y(:, 1);
    bgse = y(:, 2:4);
catch ME
    warning('Cannot read ACE/MFI HAPI data for %s: %s', file, ME.message);
    return
end

keep = t >= day0 & t < day1;
d.t = t(keep);
d.bgse = bgse(keep, :);
d.bmag = bmag(keep);
d.bVar = [dataset ':BGSEc'];
d.bmagVar = [dataset ':Magnitude'];
end

function d = emptySwe()
d = struct('t', [], 'vgse', [], 'vp', [], 'ti', [], ...
    'file', '', 'vVar', '', 'vpVar', '', 'tiVar', '', 'velocityNote', '');
end

function d = loadSweData(file, day0, day1)
d = emptySwe();
if isempty(file)
    return
end
d.file = file;
if shouldUseSweHapi(file)
    d = loadSweDataFromHapi(file, day0, day1);
    return
end
try
    d = loadSweDataFromCdf(file, day0, day1);
    if ~isempty(d.t) && (~isempty(d.vgse) || ~isempty(d.vp) || ~isempty(d.ti))
        return
    end
catch ME
    fprintf('SWEPAM CDF read failed, trying HAPI fallback: %s\n', ME.message);
end

try
    d = loadSweDataFromHapi(file, day0, day1);
catch ME
    warning('SWEPAM HAPI fallback failed for %s: %s', file, ME.message);
    return
end
end

function d = loadSweDataFromCdf(file, day0, day1)
d = emptySwe();
d.file = file;
readFile = localCdfCopy(file);
t = readCdfEpoch(readFile, 'Epoch');
n = numel(t);

vgse = [];
vVar = '';
try
    vgse = readCdfNumericNoInfo(readFile, 'V_GSE');
    vgse = orientByTime(vgse, n);
    vgse = padTo3(vgse(:, 1:min(3, size(vgse, 2))));
    vVar = 'V_GSE';
catch
end

vp = [];
vpVar = '';
try
    vp = readCdfNumericNoInfo(readFile, 'Vp');
    vp = orientByTime(vp, n);
    vp = vp(:, 1);
    vpVar = 'Vp';
catch
end

ti = [];
tiVar = '';
try
    ti = readCdfNumericNoInfo(readFile, 'Tpr');
    ti = orientByTime(ti, n);
    ti = ti(:, 1);
    tiVar = 'Tpr';
catch
end

keep = t >= day0 & t < day1;
d.t = t(keep);
if ~isempty(vgse)
    d.vgse = vgse(keep, :);
    d.vVar = vVar;
    d.velocityNote = sprintf('Using GSE ion velocity components from %s.', vVar);
elseif ~isempty(vp)
    d.vp = vp(keep);
    d.vpVar = vpVar;
    d.velocityNote = sprintf('No GSE ion velocity components are available in %s; using scalar %s.', ...
        sourceName(d), vpVar);
else
    d.velocityNote = sprintf('No usable ion velocity variable found in %s.', sourceName(d));
end
if ~isempty(ti)
    d.ti = kelvinToEv(ti(keep));
    d.tiVar = [tiVar ' (converted K to eV)'];
end
end

function d = loadSweDataFromHapi(file, day0, day1)
d = emptySwe();
d.file = file;
dataset = sweHapiDataset(file);
if any(strcmp(dataset, {'AC_H0_SWE', 'AC_H2_SWE'}))
    [t, y] = readHapiNumeric(dataset, {'Vp', 'Tpr', 'V_GSE'}, day0, day1);
    vp = y(:, 1);
    ti = y(:, 2);
    if size(y, 2) >= 5
        vgse = y(:, 3:5);
    else
        vgse = [];
    end
else
    [t, y] = readHapiNumeric(dataset, {'Vp', 'Tpr'}, day0, day1);
    vp = y(:, 1);
    ti = y(:, 2);
    vgse = [];
end

d.t = t;
d.ti = kelvinToEv(ti);
d.tiVar = [dataset ':Tpr (converted K to eV)'];
if ~isempty(vgse)
    d.vgse = vgse;
    d.vVar = [dataset ':V_GSE'];
    d.velocityNote = sprintf('Using HAPI %s GSE ion velocity components.', dataset);
else
    d.vp = vp;
    d.vpVar = [dataset ':Vp'];
    d.velocityNote = sprintf('No GSE ion velocity components are available in %s; using scalar Vp from HAPI.', dataset);
end
end

function tev = kelvinToEv(tk)
tev = tk ./ 11604.51812;
end

function tf = shouldUseSweHapi(file)
[~, name] = fileparts(file);
name = lower(name);
tf = contains(name, 'k0_swe') || contains(name, 'k1_swe');
end

function [eflux, iflux] = loadParticleFlux(epmFile, sweFile, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();

if ~isempty(epmFile)
    [eflux, iflux] = loadFluxFromFile(epmFile, day0, day1);
end

% Some local archives keep lower-energy plasma fluxes in SWE/SWEPAM files.
if (isempty(eflux.t) || isempty(iflux.t)) && ~isempty(sweFile)
    [ef2, if2] = loadFluxFromFile(sweFile, day0, day1);
    if isempty(eflux.t)
        eflux = ef2;
    end
    if isempty(iflux.t)
        iflux = if2;
    end
end
end

function f = emptyFlux()
f = struct('t', [], 'y', [], 'labels', {{}}, ...
    'file', '', 'source', '');
end

function [eflux, iflux] = loadFluxFromFile(file, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();
if isempty(file)
    return
end

if shouldUseFluxHapi(file)
    [eflux, iflux] = loadFluxFromHapi(file, day0, day1);
    return
end

try
    [eflux, iflux] = loadFluxFromCdfKnown(file, day0, day1);
    if ~isempty(eflux.t) || ~isempty(iflux.t)
        return
    end
catch ME
    fprintf('EPAM CDF read failed, trying HAPI fallback: %s\n', ME.message);
end

try
    [eflux, iflux] = loadFluxFromHapi(file, day0, day1);
catch ME
    warning('EPAM HAPI fallback failed for %s: %s', file, ME.message);
end
end

function tf = shouldUseFluxHapi(file)
[~, name] = fileparts(file);
name = lower(name);
tf = contains(name, 'epm');
end

function [eflux, iflux] = loadFluxFromCdfKnown(file, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();
readFile = localCdfCopy(file);
t = readCdfEpoch(readFile, 'Epoch');
n = numel(t);
keep = t >= day0 & t < day1;

eVars = {'DE1', 'DE2', 'DE3', 'DE4'};
iVars = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'};
[ey, elabels] = readKnownChannels(readFile, eVars, n);
[iy, ilabels] = readKnownChannels(readFile, iVars, n);

if ~isempty(ey)
    eflux.t = t(keep);
    eflux.y = ey(keep, :);
    eflux.labels = elabels;
    eflux.file = file;
end
if ~isempty(iy)
    iflux.t = t(keep);
    iflux.y = iy(keep, :);
    iflux.labels = ilabels;
    iflux.file = file;
end
end

function [eflux, iflux] = loadFluxFromHapi(file, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();
dataset = epmHapiDataset(file);
eVars = {'DE1', 'DE2', 'DE3', 'DE4'};
iVars = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'};
[t, y] = readHapiNumeric(dataset, [iVars, eVars], day0, day1);

iy = y(:, 1:numel(iVars));
ey = y(:, numel(iVars) + (1:numel(eVars)));
ey(ey <= 0 | abs(ey) > 1e30) = NaN;
iy(iy <= 0 | abs(iy) > 1e30) = NaN;

if ~all(isnan(ey(:)))
    eflux.t = t;
    eflux.y = ey;
    eflux.labels = cellfun(@(v) makeFluxLabel(v, 1, 1), eVars, 'UniformOutput', false);
    eflux.file = file;
    eflux.source = dataset;
end
if ~all(isnan(iy(:)))
    iflux.t = t;
    iflux.y = iy;
    iflux.labels = cellfun(@(v) makeFluxLabel(v, 1, 1), iVars, 'UniformOutput', false);
    iflux.file = file;
    iflux.source = dataset;
end
end

function [y, labels] = readKnownChannels(file, vars, n)
y = [];
labels = {};
for i = 1:numel(vars)
    yi = readCdfNumericNoInfo(file, vars{i});
    yi = orientByTime(yi, n);
    yi = yi(:, 1);
    yi(yi <= 0 | abs(yi) > 1e30) = NaN;
    y = [y yi]; %#ok<AGROW>
    labels{end + 1} = makeFluxLabel(vars{i}, 1, 1); %#ok<AGROW>
end
end

function cands = findFluxCandidateVars(vars, kind)
names = cellstr(vars);
cands = {};
for i = 1:numel(names)
    nm = names{i};
    low = lower(nm);
    if any(contains(low, {'epoch', 'time', 'label', 'energy', 'delta', ...
            'width', 'quality', 'status', 'sector', 'pitch', 'theta', ...
            'phi', 'angle', 'unit', 'valid', 'fill', 'cartesian'}))
        continue
    end

    isElectron = ~isempty(regexp(nm, '(^|[_-])DE\d+', 'once')) || ...
        (contains(low, 'electron') && contains(low, 'flux')) || ...
        (contains(low, 'elec') && contains(low, 'flux'));
    isIon = ~isempty(regexp(nm, '(^|[_-])FP\d+', 'once')) || ...
        ~isempty(regexp(nm, '(^|[_-])P\d+', 'once')) || ...
        (contains(low, 'ion') && contains(low, 'flux')) || ...
        (contains(low, 'proton') && contains(low, 'flux'));

    if strcmp(kind, 'electron') && isElectron
        cands{end + 1} = nm; %#ok<AGROW>
    elseif strcmp(kind, 'ion') && isIon
        cands{end + 1} = nm; %#ok<AGROW>
    end
end
cands = unique(cands, 'stable');
end

function [y, labels] = readFluxChannels(file, info, cands, n, maxChannels)
y = [];
labels = {};
for i = 1:numel(cands)
    if size(y, 2) >= maxChannels
        break
    end
    try
        yi = readCdfNumeric(file, cands{i}, info);
        yi = orientByTime(yi, n);
    catch
        continue
    end
    if isempty(yi) || ~isnumeric(yi)
        continue
    end
    yi(yi <= 0) = NaN;
    if all(isnan(yi(:)))
        continue
    end
    nadd = min(size(yi, 2), maxChannels - size(y, 2));
    if nadd <= 0
        break
    end
    y = [y yi(:, 1:nadd)]; %#ok<AGROW>
    for k = 1:nadd
        labels{end + 1} = makeFluxLabel(cands{i}, k, size(yi, 2)); %#ok<AGROW>
    end
end
end

function label = makeFluxLabel(varName, channelIndex, nChannels)
label = epamFluxLabel(varName);
if ~isempty(label)
    return
end
if nChannels == 1
    label = varName;
else
    label = sprintf('%s ch.%d', varName, channelIndex);
end
end

function label = epamFluxLabel(varName)
name = lower(char(varName));
label = '';
switch name
    case 'de1'
        label = 'e^- 38-53 keV';
    case 'de2'
        label = 'e^- 53-103 keV';
    case 'de3'
        label = 'e^- 103-175 keV';
    case 'de4'
        label = 'e^- 175-315 keV';
    case {'e1', 'e1p'}
        label = 'e^- 44-62 keV';
    case {'e2', 'e2p'}
        label = 'e^- 58-104 keV';
    case {'e3', 'e3p'}
        label = 'e^- 102-180 keV';
    case {'e4', 'e4p'}
        label = 'e^- 175-295 keV';
    case {'p1', 'p1p'}
        label = 'H^+ 46-68 keV';
    case {'p2', 'p2p'}
        label = 'H^+ 67-115 keV';
    case {'p3', 'p3p'}
        label = 'H^+ 115-195 keV';
    case {'p4', 'p4p'}
        label = 'H^+ 193-321 keV';
    case {'p5', 'p5p'}
        label = 'H^+ 315-580 keV';
    case {'p6', 'p6p'}
        label = 'H^+ 0.58-1.06 MeV';
    case {'p7', 'p7p'}
        label = 'H^+ 1.06-1.90 MeV';
    case {'p8', 'p8p'}
        label = 'H^+ 1.88-4.80 MeV';
    case {'fp5', 'fp5p'}
        label = 'H^+ 0.568-0.796 MeV';
    case {'fp6', 'fp6p'}
        label = 'H^+ 0.796-1.19 MeV';
    case {'fp7', 'fp7p'}
        label = 'H^+ 1.11-4.92 MeV';
end
end

function info = getCdfInfo(file)
if exist('spdfcdfinfo', 'file') == 2
    info = spdfcdfinfo(file);
else
    info = cdfinfo(file);
end
end

function vars = getCdfVarNames(info)
if isfield(info, 'Variables') && iscell(info.Variables)
    vars = info.Variables(:, 1);
elseif isfield(info, 'Variables')
    vars = cellstr(info.Variables);
else
    vars = {};
end
vars = cellfun(@char, vars, 'UniformOutput', false);
end

function t = readCdfEpoch(file, varName)
raw = readCdfRaw(file, varName, true);
x = normalizeCdfRaw(raw);
if iscell(x)
    try
        x = cell2mat(x);
    catch
        x = str2double(string(x));
    end
end

x = double(x(:));
if isempty(x)
    t = [];
elseif all((x > 5e5 & x < 1e6) | isnan(x))
    t = x;
elseif nanmedianLocal(x) > 1e10
    t = cdfEpochMillisToDatenum(x);
else
    t = datenum(1970, 1, 1, 0, 0, 0) + x ./ 86400;
end
t = t(:);
end

function y = readCdfNumeric(file, varName, info)
raw = readCdfRaw(file, varName, false);
y = normalizeCdfRaw(raw);
if iscell(y)
    try
        y = cell2mat(y);
    catch
        y = str2double(string(y));
    end
end
y = double(squeeze(y));
y = cleanCdfValues(y, info, varName);
end

function y = readCdfNumericNoInfo(file, varName)
raw = readCdfRaw(file, varName, false);
y = normalizeCdfRaw(raw);
if iscell(y)
    try
        y = cell2mat(y);
    catch
        y = str2double(y);
    end
end
y = double(squeeze(y));
y(abs(y) > 1e30) = NaN;
end

function localFile = localCdfCopy(file)
persistent cacheDir
if isempty(cacheDir)
    cacheDir = fullfile(tempdir, 'ace_plot_cdf_cache');
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end
end
[~, name, ext] = fileparts(file);
localFile = fullfile(cacheDir, [name ext]);
if ~exist(localFile, 'file')
    py = resolvePythonCommandForHapi();
    code = 'import shutil, sys; shutil.copy2(sys.argv[1], sys.argv[2])';
    cmd = [py, ' -c ', shellquote(code), ' ', shellquote(file), ' ', shellquote(localFile)];
    [status, cmdout] = system(cmd);
    if status ~= 0 || ~exist(localFile, 'file')
        error('Cannot cache CDF locally: %s', cmdout);
    end
end
end

function dataset = sweHapiDataset(file)
[~, name] = fileparts(file);
name = lower(name);
if contains(name, 'h0_swe')
    dataset = 'AC_H0_SWE';
elseif contains(name, 'h2_swe')
    dataset = 'AC_H2_SWE';
elseif contains(name, 'k1_swe')
    dataset = 'AC_K1_SWE';
else
    dataset = 'AC_K0_SWE';
end
end

function dataset = mfiHapiDataset(file)
[~, name] = fileparts(file);
name = lower(name);
if contains(name, 'h0_mfi')
    dataset = 'AC_H0_MFI';
elseif contains(name, 'h1_mfi')
    dataset = 'AC_H1_MFI';
elseif contains(name, 'h2_mfi')
    dataset = 'AC_H2_MFI';
elseif contains(name, 'k0_mfi')
    dataset = 'AC_K0_MFI';
else
    dataset = 'AC_H3_MFI';
end
end

function dataset = epmHapiDataset(file)
[~, name] = fileparts(file);
name = lower(name);
if contains(name, 'h1_epm')
    dataset = 'AC_H1_EPM';
elseif contains(name, 'h2_epm')
    dataset = 'AC_H2_EPM';
elseif contains(name, 'k0_epm')
    dataset = 'AC_K0_EPM';
elseif contains(name, 'k1_epm')
    dataset = 'AC_K1_EPM';
else
    dataset = 'AC_H3_EPM';
end
end

function [t, y] = readHapiNumeric(dataset, params, day0, day1)
cacheFile = hapiCacheFile(dataset, params, day0, day1);
matFile = hapiCacheMatFile(dataset, params, day0, day1);
url = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=', dataset, ...
    '&parameters=', strjoin(params, ','), ...
    '&time.min=', datestr(day0, 'yyyy-mm-ddTHH:MM:SS'), 'Z', ...
    '&time.max=', datestr(day1, 'yyyy-mm-ddTHH:MM:SS'), 'Z', ...
    '&format=csv'];

if ~exist(matFile, 'file')
    makeHapiMatCache(cacheFile, matFile, url);
end

if ~exist(matFile, 'file')
    error('HAPI MAT cache was not created: %s', matFile);
end

cache = load(matFile);
t = cache.t;
y = cache.y;
if isrow(t)
    t = t(:);
end
y(abs(y) > 1e30) = NaN;
end

function makeHapiMatCache(cacheFile, matFile, url)
helper = fullfile(fileparts(mfilename('fullpath')), 'hapi_csv_to_mat_v4.py');
if ~exist(helper, 'file')
    error('Missing HAPI CSV-to-MAT helper: %s', helper);
end
py = resolvePythonCommandForHapi();
cmd = [py, ' ', shellquote(helper), ' ', shellquote(cacheFile), ' ', ...
    shellquote(matFile), ' ', shellquote(url)];
fprintf('Preparing HAPI MAT cache: %s\n', matFile);
[status, cmdout] = system(cmd);
if status ~= 0
    error('HAPI cache preparation failed: %s', cmdout);
end
end

function file = hapiCacheFile(dataset, params, day0, day1)
baseDir = fullfile(fileparts(mfilename('fullpath')), 'hapi_cache');
paramKey = strjoin(params, '_');
if abs(day0 - floor(day0)) < 1e-9 && abs(day1 - (floor(day0) + 1)) < 1e-9
    timeKey = datestr(day0, 'yyyymmdd');
else
    timeKey = [datestr(day0, 'yyyymmddTHHMM'), '_', datestr(day1, 'yyyymmddTHHMM')];
end
file = fullfile(baseDir, [dataset, '_', timeKey, '_', paramKey, '.csv']);
end

function file = hapiCacheMatFile(dataset, params, day0, day1)
csvFile = hapiCacheFile(dataset, params, day0, day1);
[folder, stem] = fileparts(csvFile);
file = fullfile(folder, [stem, '.mat']);
end

function t = hapiTimeToDatenum(s)
s = char(s);
s = strrep(s, 'Z', '');
v = sscanf(s, '%d-%d-%dT%d:%d:%f');
if numel(v) < 6
    t = NaN;
else
    t = datenum(v(1), v(2), v(3), v(4), v(5), v(6));
end
end

function py = resolvePythonCommandForHapi()
envPy = strtrim(getenv('ACE_PYTHON'));
if ~isempty(envPy) && exist(envPy, 'file')
    py = shellquote(envPy);
    return
end
userProfile = getenv('USERPROFILE');
bundled = fullfile(userProfile, '.cache', 'codex-runtimes', ...
    'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
if exist(bundled, 'file')
    py = shellquote(bundled);
else
    py = 'python';
end
end

function out = shellquote(in)
in = char(in);
if ispc
    dq = char(34);
    in = strrep(in, dq, [dq dq]);
    out = [dq in dq];
else
    q = char(39);
    dq = char(34);
    in = strrep(in, q, [q dq q dq q]);
    out = [q in q];
end
end

function m = nanmedianLocal(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = median(x);
end
end

function t = cdfEpochMillisToDatenum(epochMs)
try
    bd = cdflib.epochBreakdown(epochMs);
    if size(bd, 1) ~= numel(epochMs)
        bd = bd.';
    end
    sec = bd(:, 6) + bd(:, 7) ./ 1000;
    t = datenum(bd(:, 1), bd(:, 2), bd(:, 3), bd(:, 4), bd(:, 5), sec);
catch
    t = NaN(size(epochMs));
end
end

function raw = readCdfRaw(file, varName, convertEpoch)
readerArgs = {'Variables', {varName}};
if convertEpoch
    readerArgs = [readerArgs {'ConvertEpochToDatenum', true}];
end

if exist('spdfcdfread', 'file') == 2
    try
        raw = spdfcdfread(file, readerArgs{:}, 'CombineRecords', true);
    catch
        raw = spdfcdfread(file, readerArgs{:});
    end
else
    try
        raw = cdfread(file, readerArgs{:}, 'CombineRecords', true);
    catch
        raw = cdfread(file, readerArgs{:});
    end
end
end

function x = normalizeCdfRaw(raw)
if iscell(raw)
    if numel(raw) == 1
        x = raw{1};
    elseif size(raw, 2) == 1
        try
            x = cell2mat(raw);
        catch
            x = raw;
        end
    else
        try
            x = cell2mat(raw(:, 1));
        catch
            x = raw(:, 1);
        end
    end
else
    x = raw;
end
end

function y = cleanCdfValues(y, info, varName)
y(abs(y) > 1e30) = NaN;

fill = getVarAttr(info, varName, {'FILLVAL', 'FillVal', '_FillValue'});
if ~isempty(fill) && isnumeric(fill)
    fill = double(fill(:)');
    for k = 1:numel(fill)
        y(y == fill(k)) = NaN;
    end
end

vmin = getVarAttr(info, varName, {'VALIDMIN', 'ValidMin'});
vmax = getVarAttr(info, varName, {'VALIDMAX', 'ValidMax'});
if ~isempty(vmin) && isnumeric(vmin)
    y(y < double(vmin(1))) = NaN;
end
if ~isempty(vmax) && isnumeric(vmax)
    y(y > double(vmax(1))) = NaN;
end
end

function val = getVarAttr(info, varName, attrNames)
val = [];
if ~isfield(info, 'VariableAttributes')
    return
end
va = info.VariableAttributes;
for ia = 1:numel(attrNames)
    attr = attrNames{ia};
    if ~isfield(va, attr)
        continue
    end
    rows = va.(attr);
    if iscell(rows) && size(rows, 2) >= 2
        idx = find(strcmpi(rows(:, 1), varName), 1);
        if ~isempty(idx)
            val = rows{idx, 2};
            return
        end
    elseif isstruct(rows)
        f = fieldnames(rows);
        idx = find(strcmpi(f, varName), 1);
        if ~isempty(idx)
            val = rows.(f{idx});
            return
        end
    end
end
end

function y = orientByTime(y, n)
if isempty(y)
    return
end
y = squeeze(y);
if isvector(y)
    y = y(:);
    return
end
sz = size(y);
if sz(1) == n
    return
elseif sz(2) == n
    y = y.';
elseif numel(sz) > 2
    y = reshape(y, sz(1), []);
    if size(y, 1) ~= n && size(y, 2) == n
        y = y.';
    end
elseif sz(1) <= 8 && sz(2) > sz(1)
    y = y.';
end
end

function y = padTo3(y)
if isempty(y)
    return
end
if size(y, 2) < 3
    y(:, end + 1:3) = NaN;
end
end

function name = findFirstVar(vars, exactNames)
name = '';
for i = 1:numel(exactNames)
    idx = find(strcmpi(vars, exactNames{i}), 1);
    if ~isempty(idx)
        name = vars{idx};
        return
    end
end
end

function name = findByWords(vars, includeWords, excludeWords)
name = '';
for i = 1:numel(vars)
    low = lower(vars{i});
    ok = true;
    for j = 1:numel(includeWords)
        ok = ok && contains(low, lower(includeWords{j}));
    end
    for j = 1:numel(excludeWords)
        ok = ok && ~contains(low, lower(excludeWords{j}));
    end
    if ok
        name = vars{i};
        return
    end
end
end

function name = pickVectorVar(file, info, vars, n, exactNames, includeWords, excludeWords)
name = findFirstVar(vars, exactNames);
if ~isempty(name) && isUsableVector(file, info, name, n)
    return
end
name = '';

cands = {};
for i = 1:numel(vars)
    low = lower(vars{i});
    ok = true;
    for j = 1:numel(includeWords)
        ok = ok && contains(low, lower(includeWords{j}));
    end
    for j = 1:numel(excludeWords)
        ok = ok && ~contains(low, lower(excludeWords{j}));
    end
    if ok
        cands{end + 1} = vars{i}; %#ok<AGROW>
    end
end

for i = 1:numel(cands)
    if isUsableVector(file, info, cands{i}, n)
        name = cands{i};
        return
    end
end
end

function ok = isUsableVector(file, info, varName, n)
ok = false;
try
    y = readCdfNumeric(file, varName, info);
    y = orientByTime(y, n);
    ok = isnumeric(y) && size(y, 1) == n && size(y, 2) >= 3;
catch
    ok = false;
end
end

function [mat, label] = readComponentTriplet(file, info, vars, n, xNames, yNames, zNames)
mat = [];
label = '';
xv = findFirstVar(vars, xNames);
yv = findFirstVar(vars, yNames);
zv = findFirstVar(vars, zNames);
if isempty(xv) || isempty(yv) || isempty(zv)
    return
end
x = orientByTime(readCdfNumeric(file, xv, info), n);
y = orientByTime(readCdfNumeric(file, yv, info), n);
z = orientByTime(readCdfNumeric(file, zv, info), n);
mat = [x(:, 1), y(:, 1), z(:, 1)];
label = sprintf('%s/%s/%s', xv, yv, zv);
end

function plotTimeData(ax, t, y, colors, useIrfPlot, useLog)
if isempty(t) || isempty(y)
    return
end
if isvector(y)
    y = y(:);
end

if useIrfPlot
    irf_plot(ax, [toEpochSeconds(t), y]);
    if useLog
        set(ax, 'YScale', 'log');
    end
else
    if useLog
        semilogy(ax, t, y, 'LineWidth', 0.75);
    else
        plot(ax, t, y, 'LineWidth', 0.85);
    end
end
applyLineColors(ax, colors);
end

function colors = componentColors()
colors = [0.10 0.25 0.85; 0.20 0.70 0.00; 0.95 0.10 0.05];
end

function colors = fluxColors(n)
base = [0.10 0.25 0.85; 0.20 0.70 0.00; 0.95 0.10 0.05; ...
    0.00 0.55 0.65; 0.75 0.20 0.75; 0.70 0.45 0.00; ...
    0.25 0.25 0.25; 0.95 0.55 0.10];
if n <= size(base, 1)
    colors = base(1:n, :);
else
    colors = lines(n);
end
end

function applyLineColors(ax, colors)
if isempty(colors)
    return
end
lines = findobj(ax, 'Type', 'line');
lines = flipud(lines(:));
for ii = 1:min(numel(lines), size(colors, 1))
    set(lines(ii), 'Color', colors(ii, :), 'LineWidth', 0.85);
end
end

function addPanelLabel(ax, label)
text(ax, 0.025, 0.84, label, 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', 18, 'Color', 'k', ...
    'Interpreter', 'none');
end

function addInlineLegend(ax, labels, colors)
x0 = 0.82;
dx = 0.045;
if numel(labels) == 4
    x0 = 0.74;
elseif numel(labels) == 1
    x0 = 0.92;
end
for ii = 1:numel(labels)
    text(ax, x0 + (ii - 1) * dx, 0.84, labels{ii}, ...
        'Units', 'normalized', 'Color', colors(ii, :), ...
        'FontSize', 13, 'FontWeight', 'bold', ...
        'Interpreter', 'tex', 'HorizontalAlignment', 'left');
end
end

function addZeroLine(ax, useIrfPlot)
xl = xlim(ax);
hold(ax, 'on');
h = plot(ax, xl, [0 0], 'k--', 'LineWidth', 0.8);
hideFromLegend(h);
if useIrfPlot
    xlim(ax, xl);
end
end

function x = toEpochSeconds(t)
x = (t(:) - datenum(1970, 1, 1, 0, 0, 0)) .* 86400;
end

function x = toPlotTime(t, useIrfPlot)
if useIrfPlot
    x = toEpochSeconds(t);
else
    x = t;
end
end

function putNoDataText(ax, msg, day0, day1, useIrfPlot)
cla(ax);
x = toPlotTime([day0 day1], useIrfPlot);
plot(ax, x, [NaN NaN], 'w-');
hold(ax, 'on');
text(ax, 0.5, 0.5, msg, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Color', [0.45 0.45 0.45]);
set(ax, 'YTick', []);
end

function setTimeAxisLimits(ax, day0, day1, useIrfPlot)
x = toPlotTime([day0 day1], useIrfPlot);
xlim(ax, x);
end

function addEventLine(ax, t, useIrfPlot)
yl = ylim(ax);
x = toPlotTime([t t], useIrfPlot);
if exist('xline', 'file') == 2 || ~verLessThan('matlab', '9.5')
    h = xline(ax, x(1), 'k--', 'LineWidth', 0.9);
    hideFromLegend(h);
else
    hold(ax, 'on');
    h = plot(ax, x, yl, 'k--', 'LineWidth', 0.9);
    hideFromLegend(h);
    ylim(ax, yl);
end
end

function hideFromLegend(h)
try
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
catch
end
end

function formatTimeAxis(ax, useIrfPlot)
if useIrfPlot
    if exist('irf_zoom', 'file') == 2
        try
            irf_zoom(ax, 'x');
        catch
        end
    end
else
    datetick(ax(end), 'x', 'HH:MM', 'keeplimits');
end
end

function stylePanels(ax)
for ii = 1:numel(ax)
    set(ax(ii), 'Box', 'on', 'Layer', 'top', 'TickDir', 'in', ...
        'LineWidth', 0.8, 'FontSize', 11);
    if ii < numel(ax)
        set(ax(ii), 'XTickLabel', []);
    end
end
end

function t = parseUtcTime(x)
if isnumeric(x)
    t = double(x);
    return
end
t = datenum(char(x), 'yyyy-mm-dd HH:MM:SS');
end

function s = fmtTime(t)
s = [datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF') ' UTC'];
end

function savePng(fig, file)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, file, 'Resolution', 200);
else
    print(fig, file, '-dpng', '-r200');
end
end
