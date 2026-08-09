function report = Voyager_Plot_Interstellar_Monthly(varargin)
% Voyager_Plot_Interstellar_Monthly Plot all locally archived VLISM data.
%
% The function creates one comprehensive figure per spacecraft per UTC
% calendar month after the accepted heliopause crossings:
%   Voyager 1: 2012-08-25
%   Voyager 2: 2018-11-05
%
% Science-product priority:
%   1) reviewed VIM MAG, 48 s (highest quantitative field quality)
%   2) PLS Level-3 fitted moments, 12-192 s, where available
%   3) COHO one-hour merged plasma, LECP, CRS and position products
%
% NASA marks post-1991 1.92 s MAG as generally not science quality, so it
% is intentionally excluded from the quantitative monthly overview.
%
% Example
%   report = Voyager_Plot_Interstellar_Monthly( ...
%       'DataRoot', 'Z:\SPART-WORK\Data\Voyager', ...
%       'OutputRoot', 'D:\Voyager_monthly', ...
%       'DownloadMissing', true);

parser = inputParser;
parser.CaseSensitive = false;
programRoot = fileparts(mfilename('fullpath'));

addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'HelperRoot', ...
    'C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu', @isTextScalar);
addParameter(parser, 'OutputRoot', fullfile(pwd, 'Voyager_monthly_plots'), @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, 'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'StartTime', [], @isTimeInput);
addParameter(parser, 'EndTime', 'auto', @isEndTimeInput);
addParameter(parser, 'DownloadMissing', false, @isLogicalScalar);
addParameter(parser, 'Overwrite', false, @isLogicalScalar);
addParameter(parser, 'Visible', false, @isLogicalScalar);
addParameter(parser, 'MakeEmptyFigures', true, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 180, @(x) isnumeric(x) && isscalar(x) && x >= 72);
addParameter(parser, 'SaveFormats', {'png'}, @isFormatInput);
addParameter(parser, 'MaxPlotPoints', 25000, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1000);
addParameter(parser, 'FillGaps', true, @isLogicalScalar);
addParameter(parser, 'MagCadenceSeconds', 48, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(parser, 'PythonExe', '', @isTextScalar);
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.HelperRoot = char(opts.HelperRoot);
opts.OutputRoot = char(opts.OutputRoot);
opts.CacheRoot = char(opts.CacheRoot);
opts.SaveFormats = normalizeFormats(opts.SaveFormats);
opts.Spacecraft = unique(double(opts.Spacecraft(:).'), 'stable');
if isempty(opts.Spacecraft) || any(~ismember(opts.Spacecraft, [1 2]))
    error('VoyagerMonthly:InvalidSpacecraft', ...
        'Spacecraft must contain 1 and/or 2.');
end
if ~isfolder(opts.DataRoot)
    error('VoyagerMonthly:DataRootMissing', ...
        'Voyager data root is unavailable: %s', opts.DataRoot);
end
if ~isfolder(opts.OutputRoot)
    [ok, message] = mkdir(opts.OutputRoot);
    if ~ok
        error('VoyagerMonthly:OutputCreateFailed', '%s', message);
    end
end
if ~isfolder(opts.CacheRoot)
    [ok, message] = mkdir(opts.CacheRoot);
    if ~ok
        error('VoyagerMonthly:CacheCreateFailed', '%s', message);
    end
end

bridgeFile = fullfile(programRoot, 'voyager_cdf_bridge.py');
if ~isfile(bridgeFile)
    error('VoyagerMonthly:BridgeMissing', 'Missing CDF bridge: %s', bridgeFile);
end
pythonCommand = resolvePython(opts.PythonExe, bridgeFile);

if logical(opts.DownloadMissing)
    downloadMissingProducts(opts);
end

crossing = [datetime(2012, 8, 25, 0, 0, 0, 'TimeZone', 'UTC'), ...
    datetime(2018, 11, 5, 0, 0, 0, 'TimeZone', 'UTC')];
requestedStart = parseOptionalUtc(opts.StartTime);
requestedEnd = parseEndUtc(opts.EndTime);

rows = cell(0, 20);
memoryCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

for spacecraft = opts.Spacecraft
    spacecraftStart = crossing(spacecraft);
    if ~isempty(requestedStart)
        spacecraftStart = max(spacecraftStart, requestedStart);
    end
    if isempty(requestedEnd)
        spacecraftEnd = latestAvailableEnd(opts.DataRoot, spacecraft);
    else
        spacecraftEnd = requestedEnd;
    end
    spacecraftEnd = min(spacecraftEnd, datetime('now', 'TimeZone', 'UTC'));
    if spacecraftEnd < spacecraftStart
        warning('VoyagerMonthly:NoInterval', ...
            'No plot interval remains for Voyager %d.', spacecraft);
        continue
    end

    monthStart = dateshift(spacecraftStart, 'start', 'month');
    finalMonth = dateshift(spacecraftEnd, 'start', 'month');
    while monthStart <= finalMonth
        monthEnd = monthStart + calmonths(1);
        dataStart = max(monthStart, spacecraftStart);
        dataEnd = min(monthEnd, spacecraftEnd + seconds(1));
        monthText = char(string(monthStart, 'yyyy-MM'));
        fprintf('[Voyager %d] %s\n', spacecraft, monthText);

        outputDir = fullfile(opts.OutputRoot, sprintf('Voyager%d', spacecraft), ...
            char(string(monthStart, 'yyyy')));
        if ~isfolder(outputDir)
            mkdir(outputDir);
        end
        outputBase = fullfile(outputDir, sprintf( ...
            'Voyager%d_VLISM_%s', spacecraft, monthText));

        magFile = findMag48sFile(opts.DataRoot, spacecraft, year(monthStart));
        cohoFile = findCohoFile(opts.DataRoot, spacecraft, monthStart);
        plsFile = findPlsFile(opts.DataRoot, spacecraft, year(monthStart));

        [magFull, magError] = tryReadBridgeProduct(magFile, 'mag48s', ...
            pythonCommand, bridgeFile, opts.CacheRoot, memoryCache, true);
        [coho, cohoError] = tryReadBridgeProduct(cohoFile, 'coho', ...
            pythonCommand, bridgeFile, opts.CacheRoot, memoryCache, false);
        [plsFull, plsError] = tryReadBridgeProduct(plsFile, 'pls', ...
            pythonCommand, bridgeFile, opts.CacheRoot, memoryCache, true);
        if coho.available
            coho = subsetByTime(coho, 'Epoch', dataStart, dataEnd);
        end
        errors = {magError, cohoError, plsError};
        monthError = strjoin(errors(~cellfun('isempty', errors)), ' | ');
        if ~isempty(monthError)
            warning('VoyagerMonthly:ReadFailed', ...
                'Voyager %d %s: %s', spacecraft, monthText, monthError);
        end

        mag = subsetByTime(magFull, 'Epoch', dataStart, dataEnd);
        pls = subsetByTime(plsFull, 'Epoch', dataStart, dataEnd);
        position = extractPosition(magFull, coho, dataStart, dataEnd);
        rawSummary = summarizeMonth(mag, pls, coho, position);
        fillStats = emptyFillStats();
        if logical(opts.FillGaps)
            [mag, coho, position, fillStats] = makeMonthContinuous( ...
                mag, coho, position, dataStart, dataEnd, ...
                double(opts.MagCadenceSeconds));
        end
        summary = summarizeMonth(mag, pls, coho, position);
        summary.raw_mag_records = rawSummary.mag_records;
        summary.raw_pls_records = rawSummary.pls_records;
        summary.raw_plasma_records = rawSummary.plasma_records;
        summary.raw_lecp_valid = rawSummary.lecp_valid;
        summary.raw_crs_valid = rawSummary.crs_valid;
        summary.mag_low_precision_filled = fillStats.mag_low_precision;
        summary.mag_interpolated = fillStats.mag_interpolated;
        summary.coho_interpolated = fillStats.coho_interpolated;
        figureFile = '';
        status = 'ok';
        if ~summary.has_any_data
            status = 'empty';
        end
        if ~isempty(monthError) && summary.has_any_data
            status = 'partial_error';
        elseif ~isempty(monthError)
            status = 'error';
        end

        shouldPlot = summary.has_any_data || logical(opts.MakeEmptyFigures);
        if shouldPlot
            try
                figureFile = exportMonthlyFigure(spacecraft, monthStart, monthEnd, ...
                    crossing(spacecraft), mag, pls, coho, position, ...
                    outputBase, opts, summary);
            catch ME
                status = 'plot_error';
                monthError = sprintf('%s: %s', ME.identifier, ME.message);
                warning('VoyagerMonthly:PlotFailed', ...
                    'Voyager %d %s: %s', spacecraft, monthText, monthError);
            end
        end

        rows(end + 1, :) = {spacecraft, monthStart, dataStart, dataEnd, ...
            status, magFile, rawSummary.mag_records, plsFile, ...
            rawSummary.pls_records, cohoFile, rawSummary.plasma_records, ...
            summary.lecp_channels, rawSummary.lecp_valid, ...
            summary.crs_channels, rawSummary.crs_valid, ...
            figureFile, monthError, fillStats.mag_low_precision, ...
            fillStats.mag_interpolated, fillStats.coho_interpolated}; %#ok<AGROW>
        monthStart = monthEnd;
    end
end

variableNames = {'Spacecraft', 'MonthUTC', 'DataStartUTC', 'DataEndUTC', ...
    'Status', 'MAGFile', 'MAGRecords', 'PLSFile', 'PLSRecords', ...
    'COHOFile', 'PlasmaRecords', 'LECPChannels', 'LECPValidValues', ...
    'CRSChannels', 'CRSValidValues', 'FigureFile', 'Notes', ...
    'MAGLowPrecisionFilled', 'MAGInterpolated', 'COHOInterpolated'};
report = cell2table(rows, 'VariableNames', variableNames);
coverageFile = fullfile(opts.OutputRoot, 'Voyager_VLISM_monthly_coverage.csv');
writetable(report, coverageFile);
save(fullfile(opts.OutputRoot, 'Voyager_VLISM_monthly_report.mat'), ...
    'report', 'opts', '-v7');
fprintf('Coverage report: %s\n', coverageFile);
end

function downloadMissingProducts(opts)
if ~isfolder(opts.HelperRoot)
    error('VoyagerMonthly:HelperRootMissing', ...
        'Downloader helper directory is unavailable: %s', opts.HelperRoot);
end
addpath(opts.HelperRoot);
if exist('Voyager_Download', 'file') ~= 2
    error('VoyagerMonthly:DownloaderMissing', ...
        'Voyager_Download.m was not found in %s.', opts.HelperRoot);
end
crossing = [datetime(2012, 8, 25, 'TimeZone', 'UTC'), ...
    datetime(2018, 11, 5, 'TimeZone', 'UTC')];
requestedStart = parseOptionalUtc(opts.StartTime);
requestedEnd = parseEndUtc(opts.EndTime);
if isempty(requestedEnd)
    requestedEnd = datetime('today', 'TimeZone', 'UTC');
end
for spacecraft = opts.Spacecraft
    startTime = crossing(spacecraft);
    if ~isempty(requestedStart), startTime = max(startTime, requestedStart); end
    if requestedEnd < startTime, continue, end
    dateRange = sprintf('%s/%s', ...
        char(string(startTime, 'yyyy-MM-dd')), ...
        char(string(requestedEnd, 'yyyy-MM-dd')));
    Voyager_Download( ...
        'Date', dateRange, ...
        'Spacecraft', spacecraft, ...
        'Products', 'mag48s_vim,coho1hr,plasma_fine,position1hr', ...
        'DataRoot', opts.DataRoot, ...
        'StageDir', fullfile(tempdir, 'Voyager_monthly_gapfill'), ...
        'ManifestName', sprintf('voyager%d_monthly_gapfill_manifest.json', spacecraft), ...
        'Threads', 5, 'CheckSize', true, 'Force', false);
end
end

function product = readBridgeProduct(sourceFile, profile, pythonCommand, ...
        bridgeFile, cacheRoot, memoryCache, retainInMemory)
cacheKey = [profile, '|', sourceFile];
if retainInMemory && isKey(memoryCache, cacheKey)
    product = memoryCache(cacheKey);
    return
end
command = sprintf('%s %s --source %s --cache-root %s --profile %s', ...
    pythonCommand, quoteArgument(bridgeFile), quoteArgument(sourceFile), ...
    quoteArgument(cacheRoot), quoteArgument(profile));
[status, output] = system(command);
if status ~= 0
    error('VoyagerMonthly:CDFBridgeFailed', ...
        'CDF bridge failed for %s.\n%s', sourceFile, output);
end
match = regexp(output, 'CACHE_JSON\|([^\r\n]+)', 'tokens', 'once');
if isempty(match)
    error('VoyagerMonthly:CDFBridgeOutput', ...
        'CDF bridge did not return a cache path. Output:\n%s', output);
end
product = readBinaryCache(strtrim(match{1}));
if retainInMemory
    memoryCache(cacheKey) = product; %#ok<NASGU>
end
end

function [product, errorText] = tryReadBridgeProduct(sourceFile, profile, ...
        pythonCommand, bridgeFile, cacheRoot, memoryCache, retainInMemory)
product = emptyProduct();
errorText = '';
if isempty(sourceFile), return, end
try
    product = readBridgeProduct(sourceFile, profile, pythonCommand, ...
        bridgeFile, cacheRoot, memoryCache, retainInMemory);
catch ME
    errorText = sprintf('%s [%s]: %s', profile, ME.identifier, ME.message);
end
end

function product = readBinaryCache(metadataFile)
metadata = jsondecode(fileread(metadataFile));
binaryFile = fullfile(fileparts(metadataFile), metadata.binary_file);
fid = fopen(binaryFile, 'rb', 'ieee-le');
if fid < 0
    error('VoyagerMonthly:BinaryCacheOpen', ...
        'Cannot open converted binary cache: %s', binaryFile);
end
cleanup = onCleanup(@() fclose(fid));
product = struct('available', true, 'source_file', metadata.source_file, ...
    'profile', metadata.profile, 'variable_meta', struct);
variables = metadata.variables;
for ii = 1:numel(variables)
    item = variables(ii);
    field = matlab.lang.makeValidName(item.name);
    fseek(fid, double(item.offset_bytes), 'bof');
    value = fread(fid, double(item.count), 'double=>double');
    shape = double(item.shape(:).');
    if isempty(shape)
        shape = [1 1];
    elseif isscalar(shape)
        shape = [shape 1]; %#ok<AGROW>
    end
    value = reshape(value, shape);
    if logical(item.is_time)
        value = datetime(value, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    end
    product.(field) = value;
    product.variable_meta.(field) = item;
end
clear cleanup
end

function output = subsetByTime(input, timeName, startTime, endTime)
output = input;
if ~isfield(input, 'available') || ~input.available || ~isfield(input, timeName)
    return
end
time = input.(timeName);
mask = time >= startTime & time < endTime;
recordCount = numel(time);
fields = fieldnames(input);
for ii = 1:numel(fields)
    field = fields{ii};
    value = input.(field);
    if (isnumeric(value) || isdatetime(value)) && ~isempty(value) && ...
            size(value, 1) == recordCount
        output.(field) = value(mask, :);
    end
end
end

function position = extractPosition(mag, coho, startTime, endTime)
position = struct('available', false, 'time', NaT(0, 1, 'TimeZone', 'UTC'), ...
    'radius_au', zeros(0, 1), 'latitude_deg', zeros(0, 1), ...
    'longitude_deg', zeros(0, 1), 'source', '');
if isfield(mag, 'Epoch_ephem') && isfield(mag, 'Radius')
    mask = mag.Epoch_ephem >= startTime & mag.Epoch_ephem < endTime;
    position.available = any(mask);
    position.time = mag.Epoch_ephem(mask);
    position.radius_au = mag.Radius(mask, :);
    if isfield(mag, 'hg_lat'), position.latitude_deg = mag.hg_lat(mask, :); end
    if isfield(mag, 'hgi_lon'), position.longitude_deg = mag.hgi_lon(mask, :); end
    position.source = mag.source_file;
elseif isfield(coho, 'Epoch') && isfield(coho, 'heliocentricDistance')
    position.available = ~isempty(coho.Epoch);
    position.time = coho.Epoch;
    position.radius_au = coho.heliocentricDistance;
    if isfield(coho, 'heliographicLatitude')
        position.latitude_deg = coho.heliographicLatitude;
    end
    if isfield(coho, 'heliographicLongitude')
        position.longitude_deg = coho.heliographicLongitude;
    end
    position.source = coho.source_file;
end
end

function stats = emptyFillStats()
stats = struct('mag_low_precision', 0, 'mag_interpolated', 0, ...
    'coho_interpolated', 0);
end

function [mag, coho, position, stats] = makeMonthContinuous( ...
        mag, coho, position, startTime, endTime, magCadenceSeconds)
stats = emptyFillStats();
[mag, stats.mag_low_precision, stats.mag_interpolated] = ...
    makeContinuousMag(mag, coho, startTime, endTime, magCadenceSeconds);
[coho, stats.coho_interpolated] = regularizeCoho( ...
    coho, startTime, endTime);
position = regularizePosition(position, startTime, endTime);
end

function [output, lowCount, interpolationCount] = makeContinuousMag( ...
        mag, coho, startTime, endTime, cadenceSeconds)
output = mag;
lowCount = 0;
interpolationCount = 0;
hasHighTime = isfield(mag, 'Epoch') && ~isempty(mag.Epoch);
hasLowTime = isfield(coho, 'Epoch') && ~isempty(coho.Epoch);
if ~hasHighTime && ~hasLowTime, return, end

grid = (startTime:seconds(cadenceSeconds): ...
    endTime - seconds(cadenceSeconds)).';
if isempty(grid), return, end
if ~isfield(output, 'available') || ~output.available
    output = emptyProduct();
    output.available = true;
    if isfield(coho, 'source_file'), output.source_file = coho.source_file; end
end
output.profile = 'continuous_mag';
output.Epoch = grid;

mapping = { ...
    'F1', 'ABS_B'; ...
    'BR', 'BR'; ...
    'BT', 'BT'; ...
    'BN', 'BN'};
for ii = 1:size(mapping, 1)
    highName = mapping{ii, 1};
    lowName = mapping{ii, 2};
    high = nan(numel(grid), 1);
    low = nan(numel(grid), 1);
    if hasHighTime && isfield(mag, highName)
        high = sampleNearest(mag.Epoch, mag.(highName), grid, ...
            max(90, 2 * cadenceSeconds));
    end
    if hasLowTime && isfield(coho, lowName)
        low = interpolateSeries(coho.Epoch, coho.(lowName), grid, false);
    end
    combined = high;
    useLow = ~isfinite(combined) & isfinite(low);
    combined(useLow) = low(useLow);
    lowCount = lowCount + sum(useLow);
    beforeInterpolation = combined;
    combined = fillContinuousVector(combined, false);
    interpolationCount = interpolationCount + sum( ...
        ~isfinite(beforeInterpolation) & isfinite(combined));
    output.(highName) = combined;
end
end

function [output, interpolationCount] = regularizeCoho( ...
        coho, startTime, endTime)
output = coho;
interpolationCount = 0;
if ~isfield(coho, 'available') || ~coho.available || ...
        ~isfield(coho, 'Epoch') || isempty(coho.Epoch)
    return
end
grid = (dateshift(startTime, 'start', 'hour'):hours(1): ...
    endTime - hours(1)).';
if isempty(grid), return, end
recordCount = numel(coho.Epoch);
fields = fieldnames(coho);
for ii = 1:numel(fields)
    field = fields{ii};
    value = coho.(field);
    if ~isnumeric(value) || isempty(value) || size(value, 1) ~= recordCount
        continue
    end
    positiveLog = ~isempty(regexp(field, '^protonFlux\d+_', 'once'));
    filledValue = nan(numel(grid), size(value, 2));
    for column = 1:size(value, 2)
        sampled = sampleNearest(coho.Epoch, value(:, column), grid, 31 * 60);
        filled = fillContinuousVector(sampled, positiveLog);
        interpolationCount = interpolationCount + sum( ...
            ~isfinite(sampled) & isfinite(filled));
        filledValue(:, column) = filled;
    end
    output.(field) = filledValue;
end
output.Epoch = grid;
end

function output = regularizePosition(position, startTime, endTime)
output = position;
if ~position.available || isempty(position.time), return, end
grid = (dateshift(startTime, 'start', 'day'):days(1): ...
    endTime - days(1)).';
if isempty(grid), return, end
names = {'radius_au', 'latitude_deg', 'longitude_deg'};
hasPosition = false;
for ii = 1:numel(names)
    name = names{ii};
    if ~isfield(position, name), continue, end
    value = interpolateSeries(position.time, position.(name), grid, false);
    output.(name) = value;
    hasPosition = hasPosition || any(isfinite(value));
end
output.time = grid;
output.available = hasPosition;
end

function sampled = sampleNearest(sourceTime, sourceValue, targetTime, toleranceSeconds)
sampled = nan(numel(targetTime), 1);
sourceTime = sourceTime(:);
sourceValue = sourceValue(:);
valid = ~isnat(sourceTime) & isfinite(sourceValue);
if ~any(valid), return, end
x = posixtime(sourceTime(valid));
y = sourceValue(valid);
[x, order] = sort(x);
y = y(order);
[x, uniqueIndex] = unique(x, 'stable');
y = y(uniqueIndex);
target = posixtime(targetTime(:));
if isscalar(x)
    distance = abs(target - x);
    sampled(distance <= toleranceSeconds) = y;
    return
end
nearestIndex = interp1(x, (1:numel(x)).', target, 'nearest', NaN);
keep = isfinite(nearestIndex);
nearestIndex(keep) = round(nearestIndex(keep));
distance = inf(size(target));
distance(keep) = abs(target(keep) - x(nearestIndex(keep)));
keep = keep & distance <= toleranceSeconds;
sampled(keep) = y(nearestIndex(keep));
end

function output = interpolateSeries(sourceTime, sourceValue, targetTime, positiveLog)
output = nan(numel(targetTime), 1);
sourceTime = sourceTime(:);
sourceValue = sourceValue(:);
valid = ~isnat(sourceTime) & isfinite(sourceValue);
if positiveLog, valid = valid & sourceValue > 0; end
if ~any(valid), return, end
x = posixtime(sourceTime(valid));
y = sourceValue(valid);
[x, order] = sort(x);
y = y(order);
[x, uniqueIndex] = unique(x, 'stable');
y = y(uniqueIndex);
target = posixtime(targetTime(:));
if isscalar(x)
    [distance, index] = min(abs(target - x));
    if distance <= 3600, output(index) = y; end
    return
end
if positiveLog, y = log10(y); end
output = interp1(x, y, target, 'linear', NaN);
output(target < x(1)) = y(1);
output(target > x(end)) = y(end);
if positiveLog, output = 10 .^ output; end
end

function output = fillContinuousVector(input, positiveLog)
output = input;
valid = isfinite(input);
if positiveLog, valid = valid & input > 0; end
if sum(valid) < 2, return, end
if positiveLog
    work = nan(size(input));
    work(valid) = log10(input(valid));
    work = fillmissing(work, 'linear', 'EndValues', 'nearest');
    output = 10 .^ work;
else
    output = fillmissing(input, 'linear', 'EndValues', 'nearest');
end
end

function summary = summarizeMonth(mag, pls, coho, position)
summary = struct;
summary.mag_records = validRowCount(mag, {'F1', 'BR', 'BT', 'BN'});
summary.pls_records = validRowCount(pls, {'V', 'dens', 'T', 'V_rtn'});
summary.plasma_records = validRowCount(coho, {'V', 'protonDensity', 'protonTemp'});
lecpNames = productVariableNames(coho, 'LECP');
crsNames = productVariableNames(coho, 'CRS');
summary.lecp_channels = numel(lecpNames);
summary.crs_channels = numel(crsNames);
summary.lecp_valid = validValueCount(coho, lecpNames);
summary.crs_valid = validValueCount(coho, crsNames);
summary.position_records = sum(isfinite(position.radius_au));
summary.has_any_data = summary.mag_records > 0 || summary.pls_records > 0 || ...
    summary.plasma_records > 0 || summary.lecp_valid > 0 || ...
    summary.crs_valid > 0 || summary.position_records > 0;
end

function outputFile = exportMonthlyFigure(spacecraft, monthStart, monthEnd, ...
        crossing, mag, pls, coho, position, outputBase, opts, summary)
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [50 50 1600 2300]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 8, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axesList = gobjects(8, 1);
for ii = 1:8
    axesList(ii) = nexttile(layout, ii);
    hold(axesList(ii), 'on');
end

plotMagMagnitude(axesList(1), mag, coho, opts.MaxPlotPoints);
plotMagComponents(axesList(2), mag, coho, opts.MaxPlotPoints);
plotPlasmaSpeed(axesList(3), pls, coho, opts.MaxPlotPoints);
plotPlasmaDensity(axesList(4), pls, coho, opts.MaxPlotPoints);
plotPlasmaTemperature(axesList(5), pls, coho, opts.MaxPlotPoints);
plotLecp(axesList(6), coho, opts.MaxPlotPoints);
plotCrs(axesList(7), coho);
plotPosition(axesList(8), position, opts.MaxPlotPoints);

xLimits = datenum([monthStart monthEnd]); %#ok<DATNM>
for ii = 1:numel(axesList)
    xlim(axesList(ii), xLimits);
    axesList(ii).FontSize = 9;
    axesList(ii).LineWidth = 0.75;
    if ii < numel(axesList)
        axesList(ii).XTickLabel = [];
    end
    if crossing >= monthStart && crossing < monthEnd
        xline(axesList(ii), datenum(crossing), '--', 'Heliopause', ...
            'Color', [0.65 0.1 0.65], 'LineWidth', 1.1, ...
            'LabelVerticalAlignment', 'bottom', 'FontSize', 8, ...
            'HandleVisibility', 'off'); %#ok<DATNM>
    end
end
linkaxes(axesList, 'x');
tickTimes = monthStart:caldays(5):monthEnd;
if tickTimes(end) < monthEnd, tickTimes(end + 1) = monthEnd; end
axesList(end).XTick = datenum(tickTimes); %#ok<DATNM>
datetick(axesList(end), 'x', 'dd', 'keepticks', 'keeplimits'); %#ok<DATIC>
xlabel(axesList(end), sprintf('Day of %s (UTC)', char(string(monthStart, 'yyyy-MM'))));

titleText = sprintf(['Voyager %d VLISM monthly overview — %s UTC\n' ...
    'reviewed MAG 48 s; PLS L3 where present; COHO LECP/CRS/plasma'], ...
    spacecraft, char(string(monthStart, 'yyyy-MM')));
sgtitle(layout, titleText, 'FontWeight', 'bold', 'FontSize', 14);
annotation(fig, 'textbox', [0.01 0.001 0.98 0.018], ...
    'String', sprintf(['Raw: MAG %d, PLS %d, plasma %d | filled: ' ...
    'MAG low-res %d, MAG interpolated %d, COHO interpolated %d'], ...
    summary.raw_mag_records, summary.raw_pls_records, ...
    summary.raw_plasma_records, summary.mag_low_precision_filled, ...
    summary.mag_interpolated, summary.coho_interpolated), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'Interpreter', 'none');

outputFile = '';
for ii = 1:numel(opts.SaveFormats)
    extension = opts.SaveFormats{ii};
    file = [outputBase, '.', extension];
    if isfile(file) && ~logical(opts.Overwrite)
        if isempty(outputFile), outputFile = file; end
        continue
    end
    switch extension
        case 'png'
            exportgraphics(fig, file, 'Resolution', opts.ExportDPI);
        case 'pdf'
            exportgraphics(fig, file, 'ContentType', 'image', ...
                'Resolution', opts.ExportDPI);
        otherwise
            exportgraphics(fig, file, 'Resolution', opts.ExportDPI);
    end
    if isempty(outputFile), outputFile = file; end
end
clear cleanup
end

function plotMagMagnitude(ax, mag, coho, maxPoints)
if isfield(mag, 'Epoch') && isfield(mag, 'F1') && any(isfinite(mag.F1))
    plotReduced(ax, mag.Epoch, mag.F1, maxPoints, [0.05 0.05 0.05], 0.7);
    uncertainty = scalarFinite(mag, 'dF');
    if isfield(mag, 'profile') && strcmp(mag.profile, 'continuous_mag')
        title(ax, sprintf(['|B| — reviewed 48 s + COHO 1 h/interpolation; ' ...
            'dF = %.3g nT'], uncertainty));
    else
        title(ax, sprintf('|B| — reviewed VIM MAG, 48 s; dF = %.3g nT', uncertainty));
    end
elseif isfield(coho, 'Epoch') && isfield(coho, 'ABS_B') && any(isfinite(coho.ABS_B))
    plotReduced(ax, coho.Epoch, coho.ABS_B, maxPoints, [0.05 0.05 0.05], 0.9);
    title(ax, '|B| — COHO hourly fallback');
else
    emptyPanel(ax, '|B| — no valid measurement');
end
ylabel(ax, 'nT'); grid(ax, 'on');
end

function plotMagComponents(ax, mag, coho, maxPoints)
source = mag;
sourceTitle = 'reviewed VIM MAG, 48 s';
if isfield(mag, 'profile') && strcmp(mag.profile, 'continuous_mag')
    sourceTitle = 'reviewed 48 s + COHO 1 h/interpolation';
end
if ~isfield(source, 'Epoch') || ~all(isfield(source, {'BR', 'BT', 'BN'}))
    source = coho;
    sourceTitle = 'COHO hourly fallback';
end
if isfield(source, 'Epoch') && all(isfield(source, {'BR', 'BT', 'BN'})) && ...
        any(isfinite([source.BR; source.BT; source.BN]))
    colors = lines(3);
    plotReduced(ax, source.Epoch, source.BR, maxPoints, colors(1, :), 0.65);
    plotReduced(ax, source.Epoch, source.BT, maxPoints, colors(2, :), 0.65);
    plotReduced(ax, source.Epoch, source.BN, maxPoints, colors(3, :), 0.65);
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
    legend(ax, {'B_R', 'B_T', 'B_N'}, 'Location', 'eastoutside');
    title(ax, ['RTN magnetic components — ', sourceTitle]);
else
    emptyPanel(ax, 'RTN magnetic components — no valid measurement');
end
ylabel(ax, 'nT'); grid(ax, 'on');
end

function plotPlasmaSpeed(ax, pls, coho, maxPoints)
if isfield(pls, 'Epoch') && isfield(pls, 'V') && any(isfinite(pls.V))
    plotReduced(ax, pls.Epoch, pls.V, maxPoints, [0.1 0.35 0.8], 0.9);
    title(ax, 'Plasma bulk speed — PLS Level 3');
elseif isfield(coho, 'Epoch') && isfield(coho, 'V') && any(isfinite(coho.V))
    plotReduced(ax, coho.Epoch, coho.V, maxPoints, [0.1 0.35 0.8], 0.9);
    title(ax, 'Plasma bulk speed — COHO hourly, interpolated gaps');
else
    emptyPanel(ax, 'Plasma bulk speed — no valid measurement');
end
ylabel(ax, 'km s^{-1}'); grid(ax, 'on');
end

function plotPlasmaDensity(ax, pls, coho, maxPoints)
if isfield(pls, 'Epoch') && isfield(pls, 'dens') && any(isfinite(pls.dens))
    plotReduced(ax, pls.Epoch, pls.dens, maxPoints, [0.1 0.55 0.25], 0.9);
    title(ax, 'Proton density — PLS Level 3');
elseif isfield(coho, 'Epoch') && isfield(coho, 'protonDensity') && ...
        any(isfinite(coho.protonDensity))
    plotReduced(ax, coho.Epoch, coho.protonDensity, maxPoints, [0.1 0.55 0.25], 0.9);
    title(ax, 'Proton density — COHO hourly, interpolated gaps');
else
    emptyPanel(ax, 'Proton density — no valid measurement');
end
ylabel(ax, 'cm^{-3}'); grid(ax, 'on');
end

function plotPlasmaTemperature(ax, pls, coho, maxPoints)
if isfield(pls, 'Epoch') && isfield(pls, 'T') && any(isfinite(pls.T))
    temperatureK = pls.T .* 11604.51812;
    plotReduced(ax, pls.Epoch, temperatureK, maxPoints, [0.85 0.25 0.1], 0.9);
    title(ax, 'Proton temperature — PLS Level 3');
elseif isfield(pls, 'Epoch') && isfield(pls, 'w') && any(isfinite(pls.w))
    temperatureK = 60.5 .* pls.w .^ 2;
    plotReduced(ax, pls.Epoch, temperatureK, maxPoints, [0.85 0.25 0.1], 0.9);
    title(ax, 'Proton temperature derived from PLS thermal speed');
elseif isfield(coho, 'Epoch') && isfield(coho, 'protonTemp') && ...
        any(isfinite(coho.protonTemp))
    plotReduced(ax, coho.Epoch, coho.protonTemp, maxPoints, [0.85 0.25 0.1], 0.9);
    title(ax, 'Proton temperature — COHO hourly, interpolated gaps');
else
    emptyPanel(ax, 'Proton temperature — no valid measurement');
end
ylabel(ax, 'K'); grid(ax, 'on'); set(ax, 'YScale', 'log');
end

function plotLecp(ax, coho, maxPoints)
names = productVariableNames(coho, 'LECP');
if isempty(names) || ~isfield(coho, 'Epoch')
    emptyPanel(ax, 'LECP proton flux — no valid measurement');
    ylabel(ax, 'flux'); return
end
colors = lines(numel(names));
labels = cell(size(names));
hasData = false;
for ii = 1:numel(names)
    value = coho.(names{ii});
    value(value <= 0) = NaN;
    hasData = hasData || any(isfinite(value));
    plotReduced(ax, coho.Epoch, value, maxPoints, colors(ii, :), 0.85);
    labels{ii} = variableLabel(coho, names{ii});
end
if hasData
    set(ax, 'YScale', 'log');
    legend(ax, labels, 'Location', 'eastoutside', 'Interpreter', 'none');
    title(ax, 'LECP proton differential flux — log-interpolated gaps');
else
    emptyPanel(ax, 'LECP proton flux — no valid measurement');
end
ylabel(ax, 'cm^{-2} s^{-1} sr^{-1} MeV^{-1}'); grid(ax, 'on');
end

function plotCrs(ax, coho)
names = productVariableNames(coho, 'CRS');
if isempty(names) || ~isfield(coho, 'Epoch')
    emptyPanel(ax, 'CRS proton flux — no valid measurement');
    ylabel(ax, 'MeV'); return
end
flux = nan(numel(coho.Epoch), numel(names));
energy = nan(numel(names), 1);
for ii = 1:numel(names)
    flux(:, ii) = coho.(names{ii});
    energy(ii) = energyCenter(variableLabel(coho, names{ii}), ii);
end
flux(flux <= 0) = NaN;
if ~any(isfinite(flux(:)))
    emptyPanel(ax, 'CRS proton flux — no valid measurement');
    ylabel(ax, 'MeV'); return
end
x = datenum(coho.Epoch); %#ok<DATNM>
[xGrid, yGrid] = meshgrid(x, energy);
surface(ax, xGrid, yGrid, zeros(size(xGrid)), log10(flux.'), ...
    'EdgeColor', 'none', 'FaceColor', 'flat');
view(ax, 2); set(ax, 'YScale', 'log');
colormap(ax, turbo(256));
cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'log_{10} flux';
title(ax, sprintf(['CRS proton differential flux — %d energy channels, ' ...
    'log-interpolated gaps'], numel(names)));
ylabel(ax, 'MeV');
end

function plotPosition(ax, position, maxPoints)
if ~position.available || isempty(position.time) || ...
        ~any(isfinite(position.radius_au))
    emptyPanel(ax, 'Heliocentric position — no valid measurement');
    ylabel(ax, 'AU'); return
end
plotReduced(ax, position.time, position.radius_au, maxPoints, ...
    [0.35 0.15 0.65], 1.0);
latitude = finiteMedian(position.latitude_deg);
longitude = finiteMedian(position.longitude_deg);
title(ax, sprintf(['Heliocentric position — daily interpolation; ' ...
    'median HGI lat %.2f deg, lon %.2f deg'], ...
    latitude, longitude));
ylabel(ax, 'radius (AU)'); grid(ax, 'on');
end

function plotReduced(ax, time, value, maxPoints, color, lineWidth)
time = time(:);
value = value(:);
validTime = ~isnat(time);
time = time(validTime);
value = value(validTime);
if isempty(time), return, end
x = datenum(time); %#ok<DATNM>
finiteX = x(isfinite(x));
cadence = NaN;
if numel(finiteX) >= 3
    cadence = median(diff(finiteX), 'omitnan');
end
if numel(x) > maxPoints
    [x, value] = minMaxEnvelope(x, value, maxPoints);
end
[x, value] = breakTimeGaps(x, value, cadence);
if sum(isfinite(value)) <= 3
    plot(ax, x, value, 'o-', 'Color', color, 'LineWidth', lineWidth, ...
        'MarkerSize', 4, 'MarkerFaceColor', color);
else
    plot(ax, x, value, 'Color', color, 'LineWidth', lineWidth);
end
end

function [x, y] = minMaxEnvelope(x, y, maxPoints)
halfBudget = max(1, floor(maxPoints / 2));
binSize = ceil(numel(x) / halfBudget);
binCount = ceil(numel(x) / binSize);
padCount = binCount * binSize - numel(x);
yBlock = reshape([y; nan(padCount, 1)], binSize, binCount);
[~, minAt] = min(yBlock, [], 1, 'omitnan');
[~, maxAt] = max(yBlock, [], 1, 'omitnan');
base = (0:(binCount - 1)) .* binSize;
index = [base + minAt; base + maxAt];
index = sort(index, 1);
index = index(:);
index = index(index <= numel(x));
x = x(index);
y = y(index);
end

function [x, y] = breakTimeGaps(x, y, cadence)
if ~isfinite(cadence) || cadence <= 0, return, end
gap = find(diff(x) > 5 * cadence);
if isempty(gap), return, end
insertAt = gap + (1:numel(gap)).';
for ii = 1:numel(insertAt)
    index = insertAt(ii);
    x = [x(1:index-1); NaN; x(index:end)];
    y = [y(1:index-1); NaN; y(index:end)];
end
end

function emptyPanel(ax, message)
text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45], ...
    'FontAngle', 'italic', 'Interpreter', 'none');
title(ax, message, 'Interpreter', 'none');
end

function names = productVariableNames(product, suffix)
names = {};
if ~isfield(product, 'available') || ~product.available, return, end
fields = fieldnames(product);
pattern = ['^protonFlux(\d+)_', suffix, '$'];
number = [];
for ii = 1:numel(fields)
    token = regexp(fields{ii}, pattern, 'tokens', 'once');
    if ~isempty(token)
        names{end + 1, 1} = fields{ii}; %#ok<AGROW>
        number(end + 1, 1) = str2double(token{1}); %#ok<AGROW>
    end
end
[~, order] = sort(number);
names = names(order);
end

function label = variableLabel(product, name)
label = name;
if ~isfield(product, 'variable_meta') || ~isfield(product.variable_meta, name)
    return
end
item = product.variable_meta.(name);
if isfield(item, 'attributes')
    attributes = item.attributes;
    candidates = {'FIELDNAM', 'LABLAXIS', 'CATDESC'};
    for ii = 1:numel(candidates)
        if isfield(attributes, candidates{ii})
            raw = attributes.(candidates{ii});
            if ischar(raw) || (isstring(raw) && isscalar(raw))
                label = char(raw);
                return
            end
        end
    end
end
end

function output = energyCenter(label, fallback)
token = regexp(label, '([0-9.]+)\s*-\s*([0-9.]+)\s*MeV', ...
    'tokens', 'once', 'ignorecase');
if isempty(token)
    output = fallback;
else
    low = str2double(token{1});
    high = str2double(token{2});
    output = sqrt(low * high);
end
end

function count = validRowCount(product, names)
count = 0;
if ~isfield(product, 'available') || ~product.available, return, end
matrix = [];
for ii = 1:numel(names)
    if isfield(product, names{ii})
        value = product.(names{ii});
        if isvector(value), value = value(:); end
        if isempty(matrix), matrix = false(size(value, 1), 1); end
        if size(value, 1) == numel(matrix)
            matrix = matrix | any(isfinite(value), 2);
        end
    end
end
count = sum(matrix);
end

function count = validValueCount(product, names)
count = 0;
for ii = 1:numel(names)
    count = count + sum(isfinite(product.(names{ii})), 'all');
end
end

function value = scalarFinite(product, name)
value = NaN;
if isfield(product, name)
    raw = product.(name);
    finite = raw(isfinite(raw));
    if ~isempty(finite), value = finite(1); end
end
end

function value = finiteMedian(input)
input = input(isfinite(input));
if isempty(input), value = NaN; else, value = median(input); end
end

function file = findMag48sFile(dataRoot, spacecraft, dataYear)
folder = fullfile(dataRoot, sprintf('voyager%d', spacecraft), ...
    'mag', '48s', 'reviewed_vim', sprintf('%04d', dataYear));
file = latestVersionFile(folder);
end

function file = findCohoFile(dataRoot, spacecraft, monthStart)
folder = fullfile(dataRoot, sprintf('voyager%d', spacecraft), ...
    'coho', '1hr', 'l2', 'merged_mag_plasma', ...
    char(string(monthStart, 'yyyy')), char(string(monthStart, 'MM')));
file = latestVersionFile(folder);
end

function file = findPlsFile(dataRoot, spacecraft, dataYear)
file = '';
if spacecraft == 2
    folder = fullfile(dataRoot, 'voyager2', 'pls', 'hires', 'l3', ...
        'heliosheath', sprintf('%04d', dataYear));
    file = latestVersionFile(folder);
end
end

function file = latestVersionFile(folder)
file = '';
if ~isfolder(folder), return, end
files = dir(fullfile(folder, '*.cdf'));
if isempty(files), return, end
version = zeros(numel(files), 1);
for ii = 1:numel(files)
    token = regexp(files(ii).name, '_v(\d+)\.cdf$', 'tokens', 'once');
    if ~isempty(token), version(ii) = str2double(token{1}); end
end
[~, index] = max(version + (1:numel(files)).' .* eps);
file = fullfile(files(index).folder, files(index).name);
end

function endTime = latestAvailableEnd(dataRoot, spacecraft)
endTime = datetime(1970, 1, 1, 'TimeZone', 'UTC');
magRoot = fullfile(dataRoot, sprintf('voyager%d', spacecraft), ...
    'mag', '48s', 'reviewed_vim');
if isfolder(magRoot)
    years = dir(magRoot);
    years = years([years.isdir] & ~ismember({years.name}, {'.', '..'}));
    for ii = 1:numel(years)
        dataYear = str2double(years(ii).name);
        if isfinite(dataYear)
            candidate = datetime(dataYear, 12, 31, 23, 59, 59, 'TimeZone', 'UTC');
            endTime = max(endTime, candidate);
        end
    end
end
cohoRoot = fullfile(dataRoot, sprintf('voyager%d', spacecraft), ...
    'coho', '1hr', 'l2', 'merged_mag_plasma');
if isfolder(cohoRoot)
    files = dir(fullfile(cohoRoot, '*', '*', '*.cdf'));
    for ii = 1:numel(files)
        token = regexp(files(ii).name, '_(\d{6})01_v\d+\.cdf$', 'tokens', 'once');
        if ~isempty(token)
            candidate = datetime(token{1}, 'InputFormat', 'yyyyMM', 'TimeZone', 'UTC');
            candidate = dateshift(candidate, 'end', 'month');
            endTime = max(endTime, candidate);
        end
    end
end
if year(endTime) == 1970
    error('VoyagerMonthly:NoLocalData', ...
        'No local post-heliopause products found for Voyager %d.', spacecraft);
end
end

function product = emptyProduct()
product = struct('available', false, 'source_file', '', 'profile', '', ...
    'variable_meta', struct);
end

function command = resolvePython(explicit, bridgeFile)
candidates = {};
explicit = strtrim(char(explicit));
if ~isempty(explicit), candidates{end + 1} = quoteArgument(explicit); end
environment = strtrim(getenv('VOYAGER_PYTHON'));
if ~isempty(environment), candidates{end + 1} = quoteArgument(environment); end
if ispc
    bundled = fullfile(getenv('USERPROFILE'), '.cache', 'codex-runtimes', ...
        'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
    if isfile(bundled), candidates{end + 1} = quoteArgument(bundled); end
    candidates = [candidates, {'py -3', 'python', 'python3'}];
else
    candidates = [candidates, {'python3', 'python'}];
end
for ii = 1:numel(candidates)
    test = sprintf('%s %s --help', candidates{ii}, quoteArgument(bridgeFile));
    [status, ~] = system(test);
    if status == 0
        command = candidates{ii};
        return
    end
end
error('VoyagerMonthly:PythonUnavailable', ...
    ['A Python interpreter with numpy and cdflib is required. Retain the ' ...
    'bundled python_packages directory or pass PythonExe explicitly.']);
end

function output = quoteArgument(input)
input = char(input);
dq = char(34);
output = [dq, strrep(input, dq, [dq, dq]), dq];
end

function value = parseOptionalUtc(input)
if isempty(input), value = []; return, end
value = parseUtc(input);
end

function value = parseEndUtc(input)
if (ischar(input) || isstring(input)) && strcmpi(strtrim(char(input)), 'auto')
    value = [];
else
    value = parseUtc(input);
end
end

function value = parseUtc(input)
if isdatetime(input)
    value = input;
else
    value = datetime(input, 'TimeZone', 'UTC');
end
if isempty(value) || ~isscalar(value) || isnat(value)
    error('VoyagerMonthly:InvalidTime', 'Time inputs must be valid scalars.');
end
value.TimeZone = 'UTC';
end

function formats = normalizeFormats(input)
if ischar(input) || isstring(input)
    formats = cellstr(lower(string(input)));
else
    formats = cellfun(@(x) lower(char(x)), input, 'UniformOutput', false);
end
formats = cellfun(@(x) regexprep(x, '^\.', ''), formats, 'UniformOutput', false);
allowed = {'png', 'pdf', 'jpg', 'jpeg', 'tif', 'tiff'};
if isempty(formats) || any(~ismember(formats, allowed))
    error('VoyagerMonthly:InvalidFormat', ...
        'SaveFormats must contain PNG, PDF, JPG, or TIFF formats.');
end
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = isLogicalScalar(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isfinite(double(value)) && ismember(double(value), [0 1]);
end

function tf = isTimeInput(value)
tf = isempty(value) || isdatetime(value) || isTextScalar(value);
end

function tf = isEndTimeInput(value)
tf = isdatetime(value) || isTextScalar(value);
end

function tf = isFormatInput(value)
tf = ischar(value) || isstring(value) || iscellstr(value);
end
