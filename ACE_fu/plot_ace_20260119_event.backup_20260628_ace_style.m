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
parse(parser, varargin{:});
opt = parser.Results;

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

day0 = datetime(2026, 1, 19, 0, 0, 0, 'TimeZone', 'UTC');
day1 = day0 + days(1);
eventHint = datetime('2026-01-19T18:59:38.687', ...
    'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC');

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
if useIrfPlot
    fprintf('Plotting      : irf_plot (%s)\n', which('irf_plot'));
else
    fprintf('Plotting      : MATLAB native plot\n');
end

files.mfi = ensureAceDailyCdf('MFI magnetic field', aceSearchRoot, day0, ...
    {'mfi_h3', 'mfi_h0', 'mfi_h2', 'mfi_h1', 'mfi_k0'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);
files.swe = ensureAceDailyCdf('SWEPAM plasma', aceSearchRoot, day0, ...
    {'swe_h0', 'swe_h2', 'swe_k0', 'swe_k1'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);
files.epm = ensureAceDailyCdf('EPAM particle flux', aceSearchRoot, day0, ...
    {'epm_h3', 'epm_h2', 'epm_h1', 'epm_k0', 'epm_k1'}, ...
    logical(opt.DownloadMissing), downloaderDir, pythonScript, char(opt.PythonExe), ...
    opt.CheckSize, opt.Threads);

dispSelectedFile('MFI magnetic field', files.mfi);
dispSelectedFile('SWEPAM plasma', files.swe);
dispSelectedFile('EPAM particle flux', files.epm);

mfi = loadMfiData(files.mfi, day0, day1);
if isempty(mfi.t)
    error('No usable ACE/MFI magnetic-field data were found for 2026-01-19.');
end

swe = loadSweData(files.swe, day0, day1);
[eflux, iflux] = loadParticleFlux(files.epm, files.swe, day0, day1);

[bMax, imax] = max(mfi.bmag, [], 'omitnan');
if isempty(imax) || isnan(bMax)
    bMax = NaN;
    tBmax = eventHint;
else
    tBmax = mfi.t(imax);
end

fig = figure('Color', 'w', 'Name', 'ACE 2026-01-19 overview', ...
    'Units', 'pixels', 'Position', [80 80 1180 1180]);
tl = tiledlayout(fig, 6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('ACE 2026-01-19 UTC, B_{max}=%.3g nT at %s', ...
    bMax, fmtTime(tBmax)));

ax = gobjects(6, 1);

ax(1) = nexttile(tl);
plotTimeData(ax(1), mfi.t, mfi.bmag, [0 0 0], useIrfPlot, false);
ylabel(ax(1), '|B| [nT]');
grid(ax(1), 'on');
title(ax(1), sprintf('Magnetic field magnitude (%s)', sourceName(mfi)));

ax(2) = nexttile(tl);
plotTimeData(ax(2), mfi.t, mfi.bgse, componentColors(), useIrfPlot, false);
ylabel(ax(2), 'B_{GSE} [nT]');
legend(ax(2), {'B_x', 'B_y', 'B_z'}, 'Location', 'eastoutside');
grid(ax(2), 'on');
title(ax(2), sprintf('Magnetic field components (%s)', mfi.bVar));

ax(3) = nexttile(tl);
if ~isempty(swe.t) && ~isempty(swe.vgse)
    plotTimeData(ax(3), swe.t, swe.vgse, componentColors(), useIrfPlot, false);
    legend(ax(3), {'V_x', 'V_y', 'V_z'}, 'Location', 'eastoutside');
    title(ax(3), sprintf('Plasma bulk velocity (%s)', sourceName(swe)));
else
    putNoDataText(ax(3), 'No usable GSE plasma velocity file/variable found', day0, day1, useIrfPlot);
end
ylabel(ax(3), 'V_{GSE} [km/s]');
grid(ax(3), 'on');

ax(4) = nexttile(tl);
if ~isempty(swe.t) && ~isempty(swe.ti)
    plotTimeData(ax(4), swe.t, swe.ti, [0.65 0.15 0.75], useIrfPlot, false);
    title(ax(4), sprintf('Ion temperature (%s)', swe.tiVar));
else
    putNoDataText(ax(4), 'No usable ion temperature file/variable found', day0, day1, useIrfPlot);
end
ylabel(ax(4), 'T_i');
grid(ax(4), 'on');

ax(5) = nexttile(tl);
if ~isempty(eflux.t) && ~isempty(eflux.y)
    plotTimeData(ax(5), eflux.t, eflux.y, [], useIrfPlot, true);
    legend(ax(5), eflux.labels, 'Location', 'eastoutside', 'Interpreter', 'none');
    title(ax(5), sprintf('Electron flux (%s)', sourceName(eflux)));
else
    putNoDataText(ax(5), 'No usable electron flux file/variable found', day0, day1, useIrfPlot);
end
ylabel(ax(5), 'e flux');
grid(ax(5), 'on');

ax(6) = nexttile(tl);
if ~isempty(iflux.t) && ~isempty(iflux.y)
    plotTimeData(ax(6), iflux.t, iflux.y, [], useIrfPlot, true);
    legend(ax(6), iflux.labels, 'Location', 'eastoutside', 'Interpreter', 'none');
    title(ax(6), sprintf('Ion/proton flux (%s)', sourceName(iflux)));
else
    putNoDataText(ax(6), 'No usable ion/proton flux file/variable found', day0, day1, useIrfPlot);
end
ylabel(ax(6), 'ion flux');
xlabel(ax(6), 'UTC on 2026-01-19');
grid(ax(6), 'on');

for ia = 1:numel(ax)
    setTimeAxisLimits(ax(ia), day0, day1, useIrfPlot);
    addEventLine(ax(ia), tBmax, useIrfPlot);
end
linkaxes(ax, 'x');
formatTimeAxis(ax, useIrfPlot);

out.figure = fig;
out.png = fullfile(outDir, 'ace_20260119_overview.png');
out.fig = fullfile(outDir, 'ace_20260119_overview.fig');
out.bmax_time_utc = tBmax;
out.bmax_nt = bMax;
out.files = files;

savePng(fig, out.png);
savefig(fig, out.fig);

fprintf('Saved PNG: %s\n', out.png);
fprintf('Saved FIG: %s\n', out.fig);
fprintf('Bmax UTC : %s, %.6g nT\n', fmtTime(tBmax), bMax);
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
[parentDir, lastPart] = fileparts(aceRoot);
if strcmpi(lastPart, 'ace')
    aceDataRoot = parentDir;
    aceSearchRoot = aceRoot;
else
    aceDataRoot = aceRoot;
    aceSearchRoot = fullfile(aceRoot, 'ace');
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
    directDir = fullfile(aceSearchRoot, p.instrumentDir, p.level, 'l2', ...
        datestr(day0, 'yyyy'), datestr(day0, 'mm'));
    pat = sprintf('ac_%s_%s_%s_v*.cdf', p.level, p.fileProduct, ...
        datestr(day0, 'yyyymmdd'));
    hits = dir(fullfile(directDir, pat));
    if isempty(hits) && exist(aceSearchRoot, 'dir')
        hits = dir(fullfile(aceSearchRoot, '**', pat));
    end
    if ~isempty(hits)
        idx = selectLatestVersion(hits);
        file = fullfile(hits(idx).folder, hits(idx).name);
        return
    end
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
d = struct('t', datetime.empty(0, 1), 'bmag', [], 'bgse', [], ...
    'file', '', 'bVar', '', 'bmagVar', '');
end

function d = loadMfiData(file, day0, day1)
d = emptyMfi();
if isempty(file) || ~exist(file, 'file')
    return
end
d.file = file;
info = getCdfInfo(file);
vars = getCdfVarNames(info);
epochVar = findFirstVar(vars, {'Epoch', 'epoch', 'Epoch0', 'Time', 'time'});
if isempty(epochVar)
    warning('No epoch variable found in %s', file);
    return
end
t = readCdfEpoch(file, epochVar);
n = numel(t);

bVar = pickVectorVar(file, info, vars, n, ...
    {'BGSEc', 'BGSE', 'B_GSE', 'B_GSEc', 'Bvec_GSE'}, ...
    {'gse'}, {'gsm', 'sc_pos', 'pos', 'rms', 'sigma', 'label', 'epoch'});
if isempty(bVar)
    warning('No GSE magnetic-field vector variable found in %s', file);
    return
end
bgse = readCdfNumeric(file, bVar, info);
bgse = orientByTime(bgse, n);
bgse = bgse(:, 1:min(3, size(bgse, 2)));
bgse = padTo3(bgse);

bmagVar = findFirstVar(vars, ...
    {'Magnitude', 'Bmag', 'B_mag', 'B_MAG', 'Bt', 'Btot', 'B_Total'});
if ~isempty(bmagVar)
    bmag = readCdfNumeric(file, bmagVar, info);
    bmag = orientByTime(bmag, n);
    bmag = bmag(:, 1);
else
    btmp = bgse;
    btmp(isnan(btmp)) = 0;
    bmag = sqrt(sum(btmp.^2, 2));
    bmag(all(isnan(bgse), 2)) = NaN;
    bmagVar = 'norm(BGSE)';
end

keep = t >= day0 & t < day1;
d.t = t(keep);
d.bgse = bgse(keep, :);
d.bmag = bmag(keep);
d.bVar = bVar;
d.bmagVar = bmagVar;
end

function d = emptySwe()
d = struct('t', datetime.empty(0, 1), 'vgse', [], 'ti', [], ...
    'file', '', 'vVar', '', 'tiVar', '');
end

function d = loadSweData(file, day0, day1)
d = emptySwe();
if isempty(file) || ~exist(file, 'file')
    return
end
d.file = file;
info = getCdfInfo(file);
vars = getCdfVarNames(info);
epochVar = findFirstVar(vars, {'Epoch', 'epoch', 'Epoch0', 'Time', 'time'});
if isempty(epochVar)
    warning('No epoch variable found in %s', file);
    return
end
t = readCdfEpoch(file, epochVar);
n = numel(t);

vVar = pickVectorVar(file, info, vars, n, ...
    {'V_GSE', 'VGSE', 'V_GSE_p', 'Vgse', 'V_GSEc'}, ...
    {'gse'}, {'gsm', 'pos', 'label', 'epoch', 'sigma', 'quality'});
if ~isempty(vVar)
    vgse = readCdfNumeric(file, vVar, info);
    vgse = orientByTime(vgse, n);
    vgse = padTo3(vgse(:, 1:min(3, size(vgse, 2))));
else
    [vgse, vVar] = readComponentTriplet(file, info, vars, n, ...
        {'Vx_GSE', 'VX_GSE', 'V_X_GSE', 'x_dot_GSE', 'Vx'}, ...
        {'Vy_GSE', 'VY_GSE', 'V_Y_GSE', 'y_dot_GSE', 'Vy'}, ...
        {'Vz_GSE', 'VZ_GSE', 'V_Z_GSE', 'z_dot_GSE', 'Vz'});
end

tiVar = findFirstVar(vars, {'Tpr', 'Tp', 'T_p', 'Tpr_K', ...
    'Ion_Temperature', 'Ion_Temp', 'Proton_Temperature', ...
    'proton_temp', 'Tion', 'Ti'});
if isempty(tiVar)
    tiVar = findByWords(vars, {'temp'}, {'electron', 'quality', 'label', 'epoch'});
end
if ~isempty(tiVar)
    ti = readCdfNumeric(file, tiVar, info);
    ti = orientByTime(ti, n);
    ti = ti(:, 1);
else
    ti = [];
end

keep = t >= day0 & t < day1;
d.t = t(keep);
if ~isempty(vgse)
    d.vgse = vgse(keep, :);
    d.vVar = vVar;
end
if ~isempty(ti)
    d.ti = ti(keep);
    d.tiVar = tiVar;
end
end

function [eflux, iflux] = loadParticleFlux(epmFile, sweFile, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();

if ~isempty(epmFile) && exist(epmFile, 'file')
    [eflux, iflux] = loadFluxFromFile(epmFile, day0, day1);
end

% Some local archives keep lower-energy plasma fluxes in SWE/SWEPAM files.
if (isempty(eflux.t) || isempty(iflux.t)) && ~isempty(sweFile) && exist(sweFile, 'file')
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
f = struct('t', datetime.empty(0, 1), 'y', [], 'labels', {{}}, ...
    'file', '', 'source', '');
end

function [eflux, iflux] = loadFluxFromFile(file, day0, day1)
eflux = emptyFlux();
iflux = emptyFlux();
if isempty(file) || ~exist(file, 'file')
    return
end
info = getCdfInfo(file);
vars = getCdfVarNames(info);
epochVar = findFirstVar(vars, {'Epoch', 'epoch', 'Epoch0', 'Time', 'time'});
if isempty(epochVar)
    return
end
t = readCdfEpoch(file, epochVar);
n = numel(t);
keep = t >= day0 & t < day1;

eCands = findFluxCandidateVars(vars, 'electron');
iCands = findFluxCandidateVars(vars, 'ion');

[ey, elabels] = readFluxChannels(file, info, eCands, n, 6);
[iy, ilabels] = readFluxChannels(file, info, iCands, n, 6);

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
        if size(yi, 2) == 1
            labels{end + 1} = cands{i}; %#ok<AGROW>
        else
            labels{end + 1} = sprintf('%s[%d]', cands{i}, k); %#ok<AGROW>
        end
    end
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
    x = string(x);
end

if isdatetime(x)
    t = x(:);
    if isempty(t.TimeZone)
        t.TimeZone = 'UTC';
    end
elseif isnumeric(x)
    x = double(x(:));
    if all(x > 5e5 & x < 1e6 | isnan(x))
        t = datetime(x, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    else
        % Last-resort support for POSIX-like epoch seconds.
        t = datetime(x, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    end
else
    t = datetime(string(x(:)), 'TimeZone', 'UTC');
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
colors = [0.85 0.10 0.10; 0.05 0.55 0.15; 0.10 0.25 0.85];
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

function x = toEpochSeconds(t)
if isdatetime(t)
    if isempty(t.TimeZone)
        t.TimeZone = 'UTC';
    end
    x = posixtime(t(:));
else
    x = t(:);
end
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
    xline(ax, x(1), 'k--', 'B max', 'LineWidth', 0.9, ...
        'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
else
    hold(ax, 'on');
    plot(ax, x, yl, 'k--', 'LineWidth', 0.9);
    ylim(ax, yl);
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
    xtickformat(ax(end), 'HH:mm');
end
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
