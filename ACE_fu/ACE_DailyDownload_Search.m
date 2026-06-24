function ACE_DailyDownload_Search(varargin)
% ACE_DailyDownload_Search
% Download daily ACE MFI/SWEPAM high-resolution CDF files and search events.
%
% Default task:
%   1998-01-01 through today, one day at a time.
%   Data root:
%       /Volumes/SPART-WORK/Data/ACE/
%   Local tree, MMS-like:
%       ACE/ace/mfi/h0/l2/yyyy/mm/*.cdf
%       ACE/ace/swepam/h0/l2/yyyy/mm/*.cdf
%
% Event criteria:
%   1. MFI magnetic field magnitude > 100 nT
%   2. SWEPAM ion/plasma velocity GSE-x component Vix > 0 km/s
%
% Product choice:
%   MFI    : mfi_h3, about 1 s cadence, variable Magnitude.
%            Falls back to mfi_h0, about 16 s cadence, when mfi_h3 is absent.
%   SWEPAM : swe_h0, about 64 s cadence, variable V_GSE(:,1).
%   Newer swe_k0/k1 files may exist when swe_h0 is absent, but they only
%   contain scalar speed Vp, not V_GSE, so they are not used for Vix search.
%
% Output TXT files:
%   ACE_event_Bgt100nT_or_Vixgt0.txt
%   ACE_problem_dates.txt
%   ACE_processed_dates.txt
%
% Download layer:
%   Uses ACEFilenames.m + ACEFilesDownload_NAS.m + download_ace_files_new.py
%   from the same folder as this file by default.
%
% CDF reader:
%   Uses irfu_matlab dataobj/getmat first, with MATLAB native CDF reading as
%   a fallback.
%
% Example:
%   ACE_DailyDownload_Search
%
% Test a short interval:
%   ACE_DailyDownload_Search('StartDate','2024-01-01', ...
%       'EndDate','2024-01-03', 'ForceRedo', 1)

%% parameters
parser = inputParser;
parser.CaseSensitive = false;

addParameter(parser, 'StartDate', '1998-01-01', @(x) ischar(x) || isstring(x) || isdatetime(x));
addParameter(parser, 'EndDate', datetime('today'), @(x) ischar(x) || isstring(x) || isdatetime(x));
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\ACE\', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Threads', 8, @isnumeric);
addParameter(parser, 'CheckSize', 1, @isnumeric);
addParameter(parser, 'ForceRedo', 0, @isnumeric);
addParameter(parser, 'Timeout', 60, @isnumeric); % kept for old calls; Python backend manages network timeout
addParameter(parser, 'GapFactor', 3, @isnumeric);
addParameter(parser, 'DownloaderDir', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

StartDate = dayOnly(parser.Results.StartDate);
EndDate = dayOnly(parser.Results.EndDate);
DataRoot = char(parser.Results.DataRoot);
Threads = parser.Results.Threads;
CheckSize = logical(parser.Results.CheckSize);
ForceRedo = logical(parser.Results.ForceRedo);
Timeout = parser.Results.Timeout; %#ok<NASGU>
GapFactor = parser.Results.GapFactor;
DownloaderDir = char(parser.Results.DownloaderDir);
PythonScript = char(parser.Results.PythonScript);
PythonExe = char(parser.Results.PythonExe);

if EndDate < StartDate
    error('EndDate must be later than StartDate.');
end

if ~isfolder(DataRoot)
    mkdir(DataRoot);
end

[DownloaderDir, PythonScript] = ensureAceDownloaderOnPath(DownloaderDir, PythonScript);

EventFile = fullfile(DataRoot, 'ACE_event_Bgt100nT_or_Vixgt0.txt');
ProblemFile = fullfile(DataRoot, 'ACE_problem_dates.txt');
ProgressFile = fullfile(DataRoot, 'ACE_processed_dates.txt');

initTextFile(EventFile, ...
    ['date_utc', char(9), 'event_type', char(9), 'start_utc', char(9), ...
    'duration_s', char(9), 'start_value', char(9), 'max_value', char(9), ...
    'unit', char(9), 'samples', char(9), 'source_file']);
initTextFile(ProblemFile, ...
    ['date_utc', char(9), 'stage', char(9), 'product', char(9), 'message']);
initTextFile(ProgressFile, ...
    ['date_utc', char(9), 'status', char(9), 'message']);

doneDates = readProcessedDates(ProgressFile);

fprintf('========== ACE daily download/search started ==========\n');
fprintf('Date range : %s / %s\n', ymd(StartDate), ymd(EndDate));
fprintf('Data root  : %s\n', DataRoot);
fprintf('Downloader : %s\n', DownloaderDir);
fprintf('Event file : %s\n', EventFile);
fprintf('Problem log: %s\n', ProblemFile);

products = aceProducts();
allDays = StartDate:caldays(1):EndDate;

for ii = 1:numel(allDays)
    thisDate = allDays(ii);
    dateStr = ymd(thisDate);

    if ~ForceRedo && isKey(doneDates, dateStr)
        fprintf('[%d/%d] %s already processed, skip.\n', ii, numel(allDays), dateStr);
        continue
    end

    fprintf('\n[%d/%d] %s\n', ii, numel(allDays), dateStr);

    mfiFile = '';
    mfiProduct = [];
    sweFile = '';
    dayHadProblem = false;

    try
        [mfiFile, mfiProduct] = downloadAceDailyFirstAvailable(thisDate, ...
            [products.mfi_h3, products.mfi_h0], DataRoot, CheckSize, Threads, ...
            PythonScript, PythonExe);
    catch ME
        dayHadProblem = true;
        logProblem(ProblemFile, thisDate, 'download', 'mfi_h3_or_mfi_h0', ME.message);
    end

    try
        sweFile = downloadAceDailyProduct(thisDate, products.swe_h0, DataRoot, CheckSize, Threads, PythonScript, PythonExe);
    catch ME
        dayHadProblem = true;
        logProblem(ProblemFile, thisDate, 'download', products.swe_h0.key, ME.message);
    end

    try
        if ~isempty(mfiFile)
            [tB, bMag] = readAceBmag(mfiFile);
            bEvents = findContinuousEvents(tB, bMag, bMag > 100, GapFactor);
            appendEvents(EventFile, thisDate, 'Bmag_gt_100nT', bEvents, 'nT', mfiFile);
            fprintf('  B events: %d (%s)\n', numel(bEvents), mfiProduct.key);
        else
            fprintf('  B data: missing\n');
        end
    catch ME
        dayHadProblem = true;
        if isempty(mfiProduct)
            mfiKey = 'mfi_h3_or_mfi_h0';
        else
            mfiKey = mfiProduct.key;
        end
        logProblem(ProblemFile, thisDate, 'read/search', mfiKey, ME.message);
    end

    try
        if ~isempty(sweFile)
            [tV, vix] = readAceVix(sweFile);
            vEvents = findContinuousEvents(tV, vix, vix > 0, GapFactor);
            appendEvents(EventFile, thisDate, 'Vix_gt_0_km_s', vEvents, 'km/s', sweFile);
            fprintf('  Vix events: %d\n', numel(vEvents));
        else
            fprintf('  Vix data: missing\n');
        end
    catch ME
        dayHadProblem = true;
        logProblem(ProblemFile, thisDate, 'read/search', products.swe_h0.key, ME.message);
    end

    if dayHadProblem
        appendLine(ProgressFile, sprintf('%s\tDONE_WITH_PROBLEM\tsee problem log', dateStr));
    else
        appendLine(ProgressFile, sprintf('%s\tDONE\tok', dateStr));
    end
end

fprintf('\n========== ACE daily download/search finished ==========\n');
end

function [downloaderDir, pythonScript] = ensureAceDownloaderOnPath(downloaderDir, pythonScript)
if isempty(downloaderDir)
    downloaderDir = fileparts(mfilename('fullpath'));
end
if isempty(downloaderDir)
    downloaderDir = pwd;
end
if ~isfolder(downloaderDir)
    error('DownloaderDir does not exist: %s', downloaderDir);
end

if exist('ACEFilenames', 'file') ~= 2 || exist('ACEFilesDownload_NAS', 'file') ~= 2
    addpath(downloaderDir);
end

if exist('ACEFilenames', 'file') ~= 2
    error('Cannot find ACEFilenames.m. Put it in DownloaderDir or on MATLAB path.');
end
if exist('ACEFilesDownload_NAS', 'file') ~= 2
    error('Cannot find ACEFilesDownload_NAS.m. Put it in DownloaderDir or on MATLAB path.');
end

if isempty(pythonScript)
    pythonScript = fullfile(downloaderDir, 'download_ace_files_new.py');
end
if ~isfile(pythonScript)
    error('Cannot find Python backend: %s', pythonScript);
end
end

%% product definitions
function products = aceProducts()
products.mfi_h3 = struct( ...
    'key', 'mfi_h3', ...
    'localParts', {{'ace', 'mfi', 'h3', 'l2'}});

products.mfi_h0 = struct( ...
    'key', 'mfi_h0', ...
    'localParts', {{'ace', 'mfi', 'h0', 'l2'}});

products.swe_h0 = struct( ...
    'key', 'swe_h0', ...
    'localParts', {{'ace', 'swepam', 'h0', 'l2'}});
end

%% download
function [localFile, productUsed] = downloadAceDailyFirstAvailable(thisDate, productList, dataRoot, checkSize, threads, pythonScript, pythonExe)
if isempty(productList)
    error('Product preference list is empty.');
end
messages = {};

for ii = 1:numel(productList)
    product = productList(ii);
    try
        localFile = downloadAceDailyProduct(thisDate, product, dataRoot, ...
            checkSize, threads, pythonScript, pythonExe);
        productUsed = product;
        if ii > 1
            fprintf('  %s fallback succeeded: %s\n', product.key, getFileName(localFile));
        end
        return
    catch ME
        messages{end + 1} = sprintf('%s: %s', product.key, ME.message); %#ok<AGROW>
        if ii < numel(productList)
            fprintf('  %s unavailable, trying %s.\n', product.key, productList(ii + 1).key);
        end
    end
end

error('No ACE CDF found for preferred product list on %s. %s', ...
    ymd(thisDate), strjoin(messages, ' | '));
end

function localFile = downloadAceDailyProduct(thisDate, product, dataRoot, checkSize, threads, pythonScript, pythonExe)
dateRange = [ymd(thisDate), '/', ymd(thisDate)];
localDir = localProductDir(thisDate, product, dataRoot);
if ~isfolder(localDir)
    mkdir(localDir);
end

[filenames, ~] = ACEFilenames(dateRange, 'product', product.key, ...
    'PythonScript', pythonScript, ...
    'PythonExe', pythonExe);
if isempty(filenames)
    error('No remote CDF found for %s on %s.', product.key, ymd(thisDate));
end

remoteFilename = filenames{1};
localFile = fullfile(localDir, remoteFilename);

if isfile(localFile) && ~checkSize
    fprintf('  %s exists, skip size check: %s\n', product.key, remoteFilename);
    return
end

fprintf('  %s using ACEFilesDownload_NAS: %s\n', product.key, remoteFilename);
ACEFilesDownload_NAS(dateRange, localDir, ...
    'product', product.key, ...
    'Threads', threads, ...
    'CheckSize', double(checkSize), ...
    'KeepTree', 0, ...
    'PythonScript', pythonScript, ...
    'PythonExe', pythonExe);

if ~isfile(localFile)
    localFile = findLocalDownloadedFile(localDir, thisDate);
end
if isempty(localFile) || ~isfile(localFile)
    error('ACEFilesDownload_NAS finished but local CDF was not found for %s on %s.', ...
        product.key, ymd(thisDate));
end
end

function name = getFileName(filePath)
[~, stem, ext] = fileparts(filePath);
name = [stem, ext];
end

function localFile = findLocalDownloadedFile(localDir, thisDate)
localFile = '';
dateToken = datestr(thisDate, 'yyyymmdd');
files = dir(fullfile(localDir, ['*', dateToken, '_v*.cdf']));
if isempty(files)
    return
end

versions = -inf(numel(files), 1);
for ii = 1:numel(files)
    tok = regexp(files(ii).name, '_v(\d+)\.cdf$', 'tokens', 'once');
    if ~isempty(tok)
        versions(ii) = str2double(tok{1});
    end
end
[~, idx] = max(versions);
localFile = fullfile(localDir, files(idx).name);
end

function d = localProductDir(thisDate, product, dataRoot)
parts = product.localParts;
d = dataRoot;
for ii = 1:numel(parts)
    d = fullfile(d, parts{ii});
end
d = fullfile(d, datestr(thisDate, 'yyyy'), datestr(thisDate, 'mm'));
end

%% CDF readers
function [t, bMag] = readAceBmag(cdfFile)
irfuMsg = '';
try
    [t, data] = readIrfuTimeSeries(cdfFile, 'Magnitude');
    bMag = data(:, 1);
catch ME
    irfuMsg = ME.message;
    try
        epoch = readCdfVar(cdfFile, 'Epoch');
        bMag = double(readCdfVar(cdfFile, 'Magnitude'));
        bMag = bMag(:);
        t = cdfEpochToDatetime(epoch);
    catch ME2
        error('Cannot read Magnitude from %s. IRFU: %s | Native CDF: %s', ...
            cdfFile, irfuMsg, ME2.message);
    end
end

valid = isfinite(bMag) & abs(bMag) < 1e30 & ~isnat(t);
t = t(valid);
bMag = bMag(valid);
end

function [t, vix] = readAceVix(cdfFile)
irfuMsg = '';
try
    [t, vGse] = readIrfuTimeSeries(cdfFile, 'V_GSE');
catch ME
    irfuMsg = ME.message;
    try
        epoch = readCdfVar(cdfFile, 'Epoch');
        vGse = double(readCdfVar(cdfFile, 'V_GSE'));
        if size(vGse, 2) < 3 && size(vGse, 1) == 3
            vGse = vGse.';
        end
        t = cdfEpochToDatetime(epoch);
    catch ME2
        error('Cannot read V_GSE from %s. IRFU: %s | Native CDF: %s', ...
            cdfFile, irfuMsg, ME2.message);
    end
end

if size(vGse, 2) < 1
    error('V_GSE has no data columns in %s.', cdfFile);
end
vGse = squeeze(vGse);
if size(vGse, 1) ~= numel(t) && size(vGse, 2) == numel(t)
    vGse = vGse.';
end
if size(vGse, 1) ~= numel(t)
    error('V_GSE time/data length mismatch in %s.', cdfFile);
end
vix = vGse(:, 1);

valid = isfinite(vix) & abs(vix) < 1e30 & ~isnat(t);
t = t(valid);
vix = vix(valid);
end

function [t, data] = readIrfuTimeSeries(cdfFile, varName)
dataobjFile = which('dataobj');
if isempty(dataobjFile)
    error('IRFU dataobj is not on MATLAB path.');
end
getmatFile = fullfile(fileparts(dataobjFile), 'getmat.m');
if ~isfile(getmatFile)
    error('IRFU getmat.m was not found next to dataobj.m.');
end

dobj = dataobj(cdfFile);
mat = getmat(dobj, varName);
if isempty(mat)
    error('Variable %s was not found by IRFU getmat.', varName);
end

if isstruct(mat)
    if ~isfield(mat, 't') || ~isfield(mat, 'data')
        error('IRFU getmat returned an unsupported struct for %s.', varName);
    end
    epochSec = double(mat.t(:));
    data = double(mat.data);
    if size(data, 1) ~= numel(epochSec) && size(data, 2) == numel(epochSec)
        data = data.';
    end
    if size(data, 1) ~= numel(epochSec)
        data = reshape(data, numel(epochSec), []);
    end
else
    mat = double(mat);
    if size(mat, 2) < 2
        error('IRFU getmat returned no data columns for %s.', varName);
    end
    epochSec = mat(:, 1);
    data = mat(:, 2:end);
end

t = datetime(epochSec, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
t = t(:);
if size(data, 1) ~= numel(t)
    error('IRFU time/data length mismatch for %s.', varName);
end
end

function data = readCdfVar(cdfFile, varName)
try
    raw = spdfcdfread(cdfFile, 'Variables', {varName}, 'ConvertEpochToDatenum', true);
catch
    raw = cdfread(cdfFile, 'Variables', {varName}, ...
        'CombineRecords', true, 'ConvertEpochToDatenum', true);
end
data = unwrapCdf(raw);
end

function data = unwrapCdf(raw)
if iscell(raw)
    if isempty(raw)
        data = [];
        return
    end
    if numel(raw) == 1
        data = raw{1};
        return
    end
    try
        data = cell2mat(raw);
    catch
        data = raw;
    end
else
    data = raw;
end
end

function t = cdfEpochToDatetime(epoch)
epoch = unwrapCdf(epoch);

if isdatetime(epoch)
    t = epoch(:);
    return
end

if iscell(epoch)
    epoch = cell2mat(epoch);
end

epoch = double(epoch(:));
if isempty(epoch)
    t = datetime.empty(0, 1);
    return
end

% With ConvertEpochToDatenum=true, MATLAB usually returns datenum values
% near 730000. Raw CDF_EPOCH values are milliseconds near 6e13.
if median(epoch, 'omitnan') > 1e10
    try
        bd = cdflib.epochBreakdown(epoch);
        if size(bd, 1) ~= numel(epoch)
            bd = bd.';
        end
        t = datetime(bd(:, 1), bd(:, 2), bd(:, 3), bd(:, 4), bd(:, 5), ...
            bd(:, 6) + bd(:, 7) ./ 1000, 'TimeZone', 'UTC');
    catch
        t = NaT(size(epoch));
    end
else
    t = datetime(epoch, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
end
t = t(:);
end

%% event finding
function events = findContinuousEvents(t, x, mask, gapFactor)
events = struct('startTime', {}, 'durationSec', {}, 'startValue', {}, ...
    'maxValue', {}, 'nSamples', {});
if isempty(t) || isempty(x) || isempty(mask)
    return
end

t = t(:);
x = x(:);
mask = logical(mask(:)) & isfinite(x) & ~isnat(t);

idx = find(mask);
if isempty(idx)
    return
end

cadenceSec = estimateCadenceSeconds(t);
if isnan(cadenceSec) || cadenceSec <= 0
    cadenceSec = 0;
end

breakHere = true(size(idx));
breakHere(1) = true;
for ii = 2:numel(idx)
    adjacentRecord = (idx(ii) == idx(ii - 1) + 1);
    if cadenceSec > 0
        gapSec = seconds(t(idx(ii)) - t(idx(ii - 1)));
        adjacentTime = gapSec <= gapFactor * cadenceSec;
    else
        adjacentTime = adjacentRecord;
    end
    breakHere(ii) = ~(adjacentRecord && adjacentTime);
end

runStarts = find(breakHere);
runEnds = [runStarts(2:end) - 1; numel(idx)];

for iRun = 1:numel(runStarts)
    runIdx = idx(runStarts(iRun):runEnds(iRun));
    t0 = t(runIdx(1));
    t1 = t(runIdx(end));
    if numel(runIdx) == 1
        durationSec = cadenceSec;
    else
        durationSec = seconds(t1 - t0) + cadenceSec;
    end
    if isnan(durationSec) || durationSec < 0
        durationSec = 0;
    end

    events(end + 1).startTime = t0; %#ok<AGROW>
    events(end).durationSec = durationSec;
    events(end).startValue = x(runIdx(1));
    events(end).maxValue = max(x(runIdx), [], 'omitnan');
    events(end).nSamples = numel(runIdx);
end
end

function cadenceSec = estimateCadenceSeconds(t)
cadenceSec = NaN;
if numel(t) < 2
    return
end
dtSec = seconds(diff(t));
dtSec = dtSec(isfinite(dtSec) & dtSec > 0);
if isempty(dtSec)
    return
end
cadenceSec = median(dtSec, 'omitnan');
end

function appendEvents(eventFile, thisDate, eventType, events, unit, sourceFile)
for ii = 1:numel(events)
    line = sprintf('%s\t%s\t%s\t%.3f\t%.6g\t%.6g\t%s\t%d\t%s', ...
        ymd(thisDate), eventType, fmtTime(events(ii).startTime), ...
        events(ii).durationSec, events(ii).startValue, events(ii).maxValue, ...
        unit, events(ii).nSamples, sourceFile);
    appendLine(eventFile, line);
end
end

%% logs and helpers
function initTextFile(filePath, header)
if ~isfile(filePath)
    fid = fopen(filePath, 'w');
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', header);
    clear cleaner
end
end

function appendLine(filePath, line)
fid = fopen(filePath, 'a');
if fid < 0
    error('Cannot open text file for append: %s', filePath);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', line);
clear cleaner
end

function logProblem(problemFile, thisDate, stage, productKey, msg)
msg = regexprep(msg, '[\r\n\t]+', ' ');
appendLine(problemFile, sprintf('%s\t%s\t%s\t%s', ymd(thisDate), stage, productKey, msg));
fprintf('  problem [%s/%s]: %s\n', stage, productKey, msg);
end

function doneDates = readProcessedDates(progressFile)
doneDates = containers.Map('KeyType', 'char', 'ValueType', 'logical');
if ~isfile(progressFile)
    return
end
txt = fileread(progressFile);
lines = regexp(txt, '\r?\n', 'split');
for ii = 1:numel(lines)
    tok = regexp(lines{ii}, '^(\d{4}-\d{2}-\d{2})\tDONE\t', 'tokens', 'once');
    if ~isempty(tok)
        doneDates(tok{1}) = true;
    end
end
end

function d = dayOnly(x)
if isdatetime(x)
    d = dateshift(x, 'start', 'day');
else
    d = datetime(char(x), 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
end
d.TimeZone = 'UTC';
d = dateshift(d, 'start', 'day');
end

function s = ymd(t)
s = datestr(t, 'yyyy-mm-dd');
end

function s = fmtTime(t)
s = datestr(t, 'yyyy-mm-ddTHH:MM:SS.FFF');
end
