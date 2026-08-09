function report = Voyager_Plot_HighestPrecision_Segments(varargin)
%Voyager_Plot_HighestPrecision_Segments Plot reviewed 48 s MAG segments.
%   This workflow uses only reviewed VIM magnetic-field measurements. A new
%   availability segment starts whenever adjacent records are more than
%   two hours apart. This threshold reproduces the 59 visible high-cadence
%   clumps in Voyager 1 during 2013-01. Every plotted value is an original
%   instrument record; no interpolation is performed.

parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'OutputFolder', fullfile(pwd, ...
    'HighestPrecision_Segments'), @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'StartTime', [], @isTimeInput);
addParameter(parser, 'EndTime', [], @isTimeInput);
addParameter(parser, 'SegmentGapSeconds', 7200, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 48);
addParameter(parser, 'LineBreakSeconds', 240, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 48);
addParameter(parser, 'MinimumRecords', 2, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);
addParameter(parser, 'Overwrite', false, @isLogicalScalar);
addParameter(parser, 'Visible', false, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 120, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 72);
addParameter(parser, 'PythonExe', '', @isTextScalar);
addParameter(parser, 'ReportTag', '', @isTextScalar);
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.OutputFolder = char(opts.OutputFolder);
opts.CacheRoot = char(opts.CacheRoot);
opts.Spacecraft = unique(double(opts.Spacecraft(:).'), 'stable');
if ~isfolder(opts.DataRoot)
    error('VoyagerHigh:DataRootMissing', 'Data root is unavailable: %s', opts.DataRoot);
end
if isempty(opts.Spacecraft) || any(~ismember(opts.Spacecraft, [1 2]))
    error('VoyagerHigh:Spacecraft', 'Spacecraft must contain 1 and/or 2.');
end
if ~isfolder(opts.OutputFolder)
    [ok, message] = mkdir(opts.OutputFolder);
    if ~ok, error('VoyagerHigh:OutputCreateFailed', '%s', message); end
end

crossing = [datetime(2012, 8, 25, 'TimeZone', 'UTC'), ...
    datetime(2018, 11, 5, 'TimeZone', 'UTC')];
requestedStart = parseOptionalUtc(opts.StartTime);
requestedEnd = parseOptionalUtc(opts.EndTime);
rows = cell(0, 14);

for spacecraft = opts.Spacecraft
    startTime = crossing(spacecraft);
    if ~isempty(requestedStart), startTime = max(startTime, requestedStart); end
    endTime = requestedEnd;
    fprintf('\n[Highest precision] Voyager %d: loading reviewed VIM MAG 48 s\n', spacecraft);
    [data, sourceFiles, sourceIndex] = loadReviewedMag( ...
        opts, spacecraft, startTime, endTime);
    if isempty(data.Epoch), continue, end

    mask = data.Epoch >= startTime;
    if ~isempty(endTime), mask = mask & data.Epoch < endTime; end
    data = subsetRecords(data, mask);
    sourceIndex = sourceIndex(mask);
    [data.Epoch, order] = sort(data.Epoch);
    data = reorderRecords(data, order, 'Epoch');
    sourceIndex = sourceIndex(order);
    [~, uniqueAt] = unique(data.Epoch, 'stable');
    data = subsetRecords(data, ismember((1:numel(data.Epoch)).', uniqueAt));
    sourceIndex = sourceIndex(uniqueAt);
    if numel(data.Epoch) < opts.MinimumRecords, continue, end

    deltaSeconds = seconds(diff(data.Epoch));
    boundaries = [1; find(deltaSeconds > opts.SegmentGapSeconds) + 1; ...
        numel(data.Epoch) + 1];
    keptNumber = 0;
    for segmentNumber = 1:(numel(boundaries) - 1)
        first = boundaries(segmentNumber);
        last = boundaries(segmentNumber + 1) - 1;
        recordCount = last - first + 1;
        if recordCount < opts.MinimumRecords, continue, end
        keptNumber = keptNumber + 1;
        segment = subsetRecords(data, false(numel(data.Epoch), 1));
        fields = dataFields();
        for ff = 1:numel(fields)
            segment.(fields{ff}) = data.(fields{ff})(first:last, :);
        end
        segmentSources = unique(sourceIndex(first:last), 'stable');
        sources = strjoin(sourceFiles(segmentSources), '; ');
        segmentStart = segment.Epoch(1);
        segmentEnd = segment.Epoch(end);
        segmentDelta = seconds(diff(segment.Epoch));
        medianCadence = median(segmentDelta, 'omitnan');
        if isempty(segmentDelta), maximumGap = 0; else, maximumGap = max(segmentDelta); end
        durationHours = hours(segmentEnd - segmentStart);

        startText = compactUtc(segmentStart);
        endText = compactUtc(segmentEnd);
        outputFile = fullfile(opts.OutputFolder, sprintf( ...
            'V%d_MAG48s_segment_%04d_%s_%s.png', spacecraft, keptNumber, ...
            startText, endText));
        status = 'ok';
        note = '';
        try
            if logical(opts.Overwrite) || ~isfile(outputFile)
                exportSegmentFigure(spacecraft, keptNumber, segment, ...
                    outputFile, opts, durationHours, medianCadence, maximumGap);
            else
                status = 'existing';
            end
        catch ME
            status = 'plot_error';
            note = sprintf('%s: %s', ME.identifier, ME.message);
            warning('VoyagerHigh:PlotFailed', 'V%d segment %d: %s', ...
                spacecraft, keptNumber, note);
        end
        rows(end + 1, :) = {spacecraft, keptNumber, segmentStart, ...
            segmentEnd, durationHours, recordCount, medianCadence, ...
            maximumGap, opts.SegmentGapSeconds, opts.LineBreakSeconds, ...
            sources, outputFile, status, note}; %#ok<AGROW>
        fprintf('[Highest precision] V%d segment %d: %s -- %s (%d records)\n', ...
            spacecraft, keptNumber, char(string(segmentStart)), ...
            char(string(segmentEnd)), recordCount);
    end
end

names = {'Spacecraft', 'Segment', 'StartUTC', 'EndUTC', 'DurationHours', ...
    'Records', 'MedianCadenceSeconds', 'MaximumRecordedGapSeconds', ...
    'SegmentSplitThresholdSeconds', 'LineBreakThresholdSeconds', ...
    'SourceCDFs', 'FigureFile', 'Status', 'Notes'};
report = cell2table(rows, 'VariableNames', names);
reportSuffix = char(opts.ReportTag);
if ~isempty(reportSuffix) && reportSuffix(1) ~= '_'
    reportSuffix = ['_', reportSuffix];
end
writetable(report, fullfile(opts.OutputFolder, ...
    ['HighestPrecision_Segments_manifest', reportSuffix, '.csv']));
save(fullfile(opts.OutputFolder, ...
    ['HighestPrecision_Segments_report', reportSuffix, '.mat']), ...
    'report', 'opts', '-v7.3');
end

function [data, sourceFiles, sourceIndex] = loadReviewedMag( ...
        opts, spacecraft, startTime, endTime)
data = emptyMagData();
sourceFiles = {};
sourceIndex = zeros(0, 1);
root = fullfile(opts.DataRoot, sprintf('voyager%d', spacecraft), ...
    'mag', '48s', 'reviewed_vim');
if ~isfolder(root), return, end
yearFolders = dir(root);
yearFolders = yearFolders([yearFolders.isdir] & ...
    ~ismember({yearFolders.name}, {'.', '..'}));
yearValues = cellfun(@str2double, {yearFolders.name});
[~, order] = sort(yearValues);
yearFolders = yearFolders(order);
for ii = 1:numel(yearFolders)
    dataYear = str2double(yearFolders(ii).name);
    if ~isfinite(dataYear), continue, end
    if dataYear < year(startTime), continue, end
    if ~isempty(endTime) && dataYear > year(endTime - seconds(1)), continue, end
    file = latestVersionFile(fullfile(yearFolders(ii).folder, yearFolders(ii).name));
    if isempty(file), continue, end
    fprintf('  reading %s\n', file);
    product = Voyager_Read_CDF_Product(file, 'mag48s', ...
        'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
    if ~isfield(product, 'Epoch') || isempty(product.Epoch), continue, end
    sourceFiles{end + 1, 1} = file; %#ok<AGROW>
    n = numel(product.Epoch);
    data.Epoch = [data.Epoch; product.Epoch(:)]; %#ok<AGROW>
    numericFields = setdiff(dataFields(), {'Epoch'}, 'stable');
    for ff = 1:numel(numericFields)
        field = numericFields{ff};
        data.(field) = [data.(field); alignedColumn(product, field, n)]; %#ok<AGROW>
    end
    sourceIndex = [sourceIndex; repmat(numel(sourceFiles), n, 1)]; %#ok<AGROW>
end
end

function output = alignedColumn(product, field, recordCount)
output = nan(recordCount, 1);
if ~isfield(product, field), return, end
value = product.(field);
if isscalar(value)
    output = repmat(double(value), recordCount, 1);
elseif size(value, 1) == recordCount
    output = value(:, 1);
elseif numel(value) == recordCount
    output = value(:);
end
end

function data = emptyMagData()
data = struct('Epoch', NaT(0, 1, 'TimeZone', 'UTC'), 'F1', zeros(0, 1), ...
    'BR', zeros(0, 1), 'BT', zeros(0, 1), 'BN', zeros(0, 1), ...
    'dF', zeros(0, 1), 'dBR', zeros(0, 1), 'dBT', zeros(0, 1), ...
    'dBN', zeros(0, 1));
end

function fields = dataFields()
fields = {'Epoch', 'F1', 'BR', 'BT', 'BN', 'dF', 'dBR', 'dBT', 'dBN'};
end

function output = subsetRecords(input, mask)
output = input;
fields = dataFields();
for ii = 1:numel(fields)
    output.(fields{ii}) = input.(fields{ii})(mask, :);
end
end

function output = reorderRecords(input, order, skipField)
output = input;
fields = dataFields();
for ii = 1:numel(fields)
    if strcmp(fields{ii}, skipField), continue, end
    output.(fields{ii}) = input.(fields{ii})(order, :);
end
end

function exportSegmentFigure(spacecraft, number, data, outputFile, opts, ...
        durationHours, medianCadence, maximumGap)
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [40 40 1700 1050]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout); hold(ax1, 'on');
plotMeasuredWithBreaks(ax1, data.Epoch, data.F1, ...
                    opts.LineBreakSeconds, [0.05 0.05 0.05], 0.7);
ylabel(ax1, '|B| (nT)'); grid(ax1, 'on');
title(ax1, 'Reviewed VIM MAG magnitude, nominal 48 s records');

ax2 = nexttile(layout); hold(ax2, 'on');
colors = lines(3);
plotMeasuredWithBreaks(ax2, data.Epoch, data.BR, ...
    opts.LineBreakSeconds, colors(1, :), 0.65);
plotMeasuredWithBreaks(ax2, data.Epoch, data.BT, ...
    opts.LineBreakSeconds, colors(2, :), 0.65);
plotMeasuredWithBreaks(ax2, data.Epoch, data.BN, ...
    opts.LineBreakSeconds, colors(3, :), 0.65);
yline(ax2, 0, ':', 'Color', [0.5 0.5 0.5]);
legend(ax2, {'B_R', 'B_T', 'B_N'}, 'Location', 'eastoutside');
ylabel(ax2, 'B_{RTN} (nT)'); grid(ax2, 'on');
title(ax2, 'Reviewed VIM MAG RTN components');

ax3 = nexttile(layout); hold(ax3, 'on');
uncertaintyFields = {'dF', 'dBR', 'dBT', 'dBN'};
uncertaintyLabels = {'d|B|', 'dB_R', 'dB_T', 'dB_N'};
uncertaintyColors = lines(4);
for ii = 1:numel(uncertaintyFields)
    plotMeasuredWithBreaks(ax3, data.Epoch, data.(uncertaintyFields{ii}), ...
        opts.LineBreakSeconds, uncertaintyColors(ii, :), 0.65);
end
legend(ax3, uncertaintyLabels, 'Location', 'eastoutside');
ylabel(ax3, 'uncertainty (nT)'); grid(ax3, 'on');
title(ax3, 'Uncertainty fields supplied with the reviewed product');
xlabel(ax3, 'UTC');

axesList = [ax1 ax2 ax3];
linkaxes(axesList, 'x');
left = datenum(data.Epoch(1)); %#ok<DATNM>
right = datenum(data.Epoch(end)); %#ok<DATNM>
if right <= left, right = left + 1 / 1440; end
for ii = 1:numel(axesList)
    xlim(axesList(ii), [left right]);
    datetick(axesList(ii), 'x', 'keeplimits'); %#ok<DATIC>
end

sgtitle(layout, sprintf(['Voyager %d highest-precision segment %04d | ' ...
    '%s to %s UTC'], spacecraft, number, ...
    char(string(data.Epoch(1), 'yyyy-MM-dd HH:mm:ss')), ...
    char(string(data.Epoch(end), 'yyyy-MM-dd HH:mm:ss'))), ...
    'FontWeight', 'bold', 'FontSize', 13);
annotation(fig, 'textbox', [0.01 0.001 0.98 0.025], ...
    'String', sprintf(['Original reviewed 48 s instrument records only; ' ...
    'no interpolation | N=%d | duration=%.2f h | median cadence=%.2f s | ' ...
    'maximum recorded gap=%.1f s | line breaks > %.0f s | ' ...
    'new segment after gaps > %.0f s'], ...
    numel(data.Epoch), durationHours, medianCadence, maximumGap, ...
    opts.LineBreakSeconds, opts.SegmentGapSeconds), 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'Interpreter', 'none');
exportgraphics(fig, outputFile, 'Resolution', opts.ExportDPI);
clear cleanup
end

function plotMeasuredWithBreaks(ax, time, value, breakSeconds, color, width)
time = time(:); value = value(:);
mask = ~isnat(time) & isfinite(value);
time = time(mask); value = value(mask);
if isempty(time), return, end
[time, order] = sort(time); value = value(order);
x = datenum(time); %#ok<DATNM>
gaps = find(seconds(diff(time)) > breakSeconds);
for ii = numel(gaps):-1:1
    at = gaps(ii) + 1;
    x = [x(1:at-1); NaN; x(at:end)]; %#ok<AGROW>
    value = [value(1:at-1); NaN; value(at:end)]; %#ok<AGROW>
end
plot(ax, x, value, 'Color', color, 'LineWidth', width);
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

function value = parseOptionalUtc(input)
if isempty(input), value = []; return, end
if isdatetime(input), value = input; else, value = datetime(input, 'TimeZone', 'UTC'); end
if isempty(value) || ~isscalar(value) || isnat(value)
    error('VoyagerHigh:Time', 'Time inputs must be valid scalar values.');
end
value.TimeZone = 'UTC';
end

function output = compactUtc(value)
output = char(string(value, 'yyyyMMdd''T''HHmmss'));
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = isTimeInput(value)
tf = isempty(value) || isdatetime(value) || isTextScalar(value);
end

function tf = isLogicalScalar(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isfinite(double(value)) && ismember(double(value), [0 1]);
end
