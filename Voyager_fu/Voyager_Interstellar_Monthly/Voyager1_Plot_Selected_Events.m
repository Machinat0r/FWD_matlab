function report = Voyager1_Plot_Selected_Events(varargin)
%Voyager1_Plot_Selected_Events Plot the V1 intervals circled in 0809.pptx.
%   The source is the existing Voyager 1 COHO merged one-hour product.
%   Recorded values are retained at their native cadence. No smoothing,
%   interpolation, resampling, or gap filling is applied. The CRS variables
%   in this COHO product are documented as 6-hour values; they remain on the
%   one-hour master time axis with missing records left empty.

programRoot = fileparts(mfilename('fullpath'));
parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'OutputFolder', fullfile(programRoot, ...
    'Voyager1_Selected_Events_1h'), @isTextScalar);
addParameter(parser, 'DatabaseFile', '', @isTextScalar);
addParameter(parser, 'SourcePptx', ...
    'C:\Users\Administrator\.codex\attachments\9277234b-55db-45aa-bd78-7a2eeb200d14\0809.pptx', ...
    @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'PythonExe', '', @isTextScalar);
addParameter(parser, 'OverwriteFigures', false, @isLogicalScalar);
addParameter(parser, 'OverwriteDatabase', false, @isLogicalScalar);
addParameter(parser, 'Visible', false, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 200, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 72);
addParameter(parser, 'GapBreakHours', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
addParameter(parser, 'CRSColorLimits', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && x(1) < x(2)));
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.OutputFolder = char(opts.OutputFolder);
opts.DatabaseFile = char(opts.DatabaseFile);
opts.SourcePptx = char(opts.SourcePptx);
opts.CacheRoot = char(opts.CacheRoot);
opts.PythonExe = char(opts.PythonExe);
if isempty(opts.DatabaseFile)
    opts.DatabaseFile = fullfile(opts.DataRoot, 'voyager1', 'events', ...
        'Voyager1_Selected_Events_1h.sqlite');
end
if ~isfolder(opts.DataRoot)
    error('VoyagerEvent:DataRootMissing', ...
        'Voyager data root is unavailable: %s', opts.DataRoot);
end
ensureFolder(opts.OutputFolder);
ensureFolder(fileparts(opts.DatabaseFile));

catalog = selectedEventCatalog(opts.SourcePptx);
catalogFile = fullfile(opts.OutputFolder, ...
    'Voyager1_Selected_Events_1h_catalog.csv');
writetable(catalog, catalogFile);

monthlyFiles = listMonthlyCoho(opts.DataRoot);
eventData = cell(height(catalog), 1);
sourceFiles = cell(height(catalog), 1);
loadNotes = strings(height(catalog), 1);
allCRSLogFlux = zeros(0, 1);

for ii = 1:height(catalog)
    fprintf('[V1 selected event] %s  %s to %s UTC\n', ...
        catalog.EventID{ii}, char(catalog.StartUTC(ii)), ...
        char(catalog.EndUTCInclusive(ii)));
    try
        [eventData{ii}, sourceFiles{ii}] = loadEventProduct( ...
            monthlyFiles, catalog.StartUTC(ii), ...
            catalog.EndUTCExclusive(ii), opts);
        crsNames = productVariableNames(eventData{ii}, 'CRS');
        for jj = 1:numel(crsNames)
            value = eventData{ii}.(crsNames{jj});
            value = value(isfinite(value) & value > 0);
            allCRSLogFlux = [allCRSLogFlux; log10(value(:))]; %#ok<AGROW>
        end
    catch ME
        loadNotes(ii) = string(sprintf('%s: %s', ME.identifier, ME.message));
        warning('VoyagerEvent:LoadFailed', '%s: %s', ...
            catalog.EventID{ii}, loadNotes(ii));
        eventData{ii} = emptyProduct();
        sourceFiles{ii} = cell(0, 1);
    end
end

if isempty(opts.CRSColorLimits)
    colorLimits = robustColorLimits(allCRSLogFlux);
else
    colorLimits = double(opts.CRSColorLimits(:).');
end

figureFiles = cell(height(catalog), 1);
status = repmat({'ok'}, height(catalog), 1);
notes = cellstr(loadNotes);
epochRecords = zeros(height(catalog), 1);
magneticRecords = zeros(height(catalog), 1);
lecpValidValues = zeros(height(catalog), 1);
crsValidValues = zeros(height(catalog), 1);

for ii = 1:height(catalog)
    eventID = catalog.EventID{ii};
    fileName = sprintf('V1_%s_%s_%s_COHO1h_raw.png', ...
        eventID, datestr(catalog.StartUTC(ii), 'yyyymmdd'), ...
        datestr(catalog.EndUTCInclusive(ii), 'yyyymmdd')); %#ok<DATST>
    figureFiles{ii} = fullfile(opts.OutputFolder, fileName);
    data = eventData{ii};
    epochRecords(ii) = variableRecordCount(data, 'Epoch');
    magneticRecords(ii) = anyFiniteRecordCount(data, ...
        {'ABS_B', 'F', 'BR', 'BT', 'BN'});
    lecpValidValues(ii) = validValueCount(data, ...
        productVariableNames(data, 'LECP'));
    crsValidValues(ii) = validValueCount(data, ...
        productVariableNames(data, 'CRS'));
    if strlength(loadNotes(ii)) > 0
        status{ii} = 'load_error';
    end
    try
        if logical(opts.OverwriteFigures) || ~isfile(figureFiles{ii})
            exportEventFigure(catalog(ii, :), data, figureFiles{ii}, ...
                opts, colorLimits);
        else
            status{ii} = 'existing';
        end
    catch ME
        status{ii} = 'plot_error';
        if isempty(notes{ii})
            notes{ii} = sprintf('%s: %s', ME.identifier, ME.message);
        else
            notes{ii} = sprintf('%s | %s: %s', notes{ii}, ...
                ME.identifier, ME.message);
        end
        warning('VoyagerEvent:PlotFailed', '%s: %s', eventID, notes{ii});
    end
end

report = catalog(:, {'EventID', 'SourceSlide', 'PPTCircleColor', ...
    'StartUTC', 'EndUTCInclusive', 'EndUTCExclusive'});
report.SourceFileCount = cellfun(@numel, sourceFiles);
report.EpochRecords = epochRecords;
report.MagneticRecords = magneticRecords;
report.LECPValidValues = lecpValidValues;
report.CRSValidValues = crsValidValues;
report.OriginalInstrumentValuesOnly = true(height(report), 1);
report.GeneratedOrInterpolatedValues = zeros(height(report), 1);
report.FigureFile = figureFiles;
report.Status = status;
report.Notes = notes;

manifestFile = fullfile(opts.OutputFolder, ...
    'Voyager1_Selected_Events_1h_manifest.csv');
writetable(report, manifestFile);
save(fullfile(opts.OutputFolder, ...
    'Voyager1_Selected_Events_1h_report.mat'), ...
    'report', 'catalog', 'colorLimits', 'opts', '-v7.3');

writeEventDatabase(opts.DatabaseFile, catalog, report, eventData, ...
    sourceFiles, colorLimits, opts);

fprintf('Completed %d events. Figures: %s\n', height(report), opts.OutputFolder);
fprintf('SQLite database: %s\n', opts.DatabaseFile);
end

function catalog = selectedEventCatalog(sourcePptx)
% Ellipses are ordered chronologically within each source slide.
years = [2013; 2013; 2014; 2015; 2016; 2016; 2017; ...
    2018; 2018; 2018; 2019; 2020; 2020; 2021; 2021; 2022];
slides = [4; 4; 5; 6; 7; 7; 8; 9; 9; 9; 10; 11; 11; 12; 12; 13];
orders = [1; 2; 1; 1; 1; 2; 1; 1; 2; 3; 1; 1; 2; 1; 2; 1];
xEMU = [4083801; 6052086; 3835828; 1983782; 6827002; 9043259; ...
    4014061; 2247254; 2905931; 8888277; 3246893; 6168325; ...
    8880528; 3316637; 5114441; 2468104];
widthEMU = [325465; 325465; 325465; 325465; 325465; 325465; ...
    658678; 565688; 154983; 154983; 325465; 325465; 325465; ...
    286719; 286719; 286719];
colors = {'accent1'; 'accent1'; 'accent1'; 'accent1'; 'accent1'; ...
    'accent1'; 'accent1'; 'accent1'; 'accent1'; 'accent1'; ...
    'accent1'; 'accent2'; 'accent1'; 'accent1'; 'accent2'; 'accent2'};

% The annual PNG is 2783 px wide and is placed identically on slides 4-13.
% Pixel columns 118 and 2443 are the left and right plot-axis boundaries.
pictureXEMU = 1459019;
pictureWidthEMU = 9273962;
sourceWidthPixels = 2783;
plotLeftEMU = pictureXEMU + ...
    118 / sourceWidthPixels * pictureWidthEMU;
plotRightEMU = pictureXEMU + ...
    2443 / sourceWidthPixels * pictureWidthEMU;

rawStart = NaT(size(years), 'TimeZone', 'UTC');
rawEnd = NaT(size(years), 'TimeZone', 'UTC');
for ii = 1:numel(years)
    yearStart = datetime(years(ii), 1, 1, 'TimeZone', 'UTC');
    yearEnd = datetime(years(ii) + 1, 1, 1, 'TimeZone', 'UTC');
    yearDays = days(yearEnd - yearStart);
    fractionStart = (xEMU(ii) - plotLeftEMU) / ...
        (plotRightEMU - plotLeftEMU);
    fractionEnd = (xEMU(ii) + widthEMU(ii) - plotLeftEMU) / ...
        (plotRightEMU - plotLeftEMU);
    rawStart(ii) = yearStart + days(fractionStart * yearDays);
    rawEnd(ii) = yearStart + days(fractionEnd * yearDays);
end
startUTC = dateshift(rawStart, 'start', 'day');
endUTCExclusive = dateshift(rawEnd, 'start', 'day');
needsNextDay = rawEnd > endUTCExclusive;
endUTCExclusive(needsNextDay) = ...
    endUTCExclusive(needsNextDay) + caldays(1);
endUTCInclusive = endUTCExclusive - caldays(1);

eventID = cell(numel(years), 1);
for ii = 1:numel(years)
    eventID{ii} = sprintf('%04d-E%02d', years(ii), orders(ii));
end
sourcePptxColumn = repmat({sourcePptx}, numel(years), 1);
extractionMethod = repmat({['PPTX vector ellipse x bounds mapped ', ...
    'linearly to the annual UTC axis; rounded outward to whole days']}, ...
    numel(years), 1);
catalog = table(eventID, years, slides, orders, colors, xEMU, widthEMU, ...
    rawStart, rawEnd, startUTC, endUTCInclusive, endUTCExclusive, ...
    sourcePptxColumn, extractionMethod, ...
    'VariableNames', {'EventID', 'YearUTC', 'SourceSlide', ...
    'CircleOrderOnSlide', 'PPTCircleColor', 'PPTCircleX_EMU', ...
    'PPTCircleWidth_EMU', 'PPTMappedStartUTC', 'PPTMappedEndUTC', ...
    'StartUTC', 'EndUTCInclusive', 'EndUTCExclusive', ...
    'SourcePPTX', 'ExtractionMethod'});
end

function files = listMonthlyCoho(dataRoot)
root = fullfile(dataRoot, 'voyager1', 'coho', '1hr', 'l2', ...
    'merged_mag_plasma');
if ~isfolder(root)
    error('VoyagerEvent:COHORootMissing', ...
        'Voyager 1 COHO one-hour folder is unavailable: %s', root);
end
candidates = dir(fullfile(root, '*', '*', '*.cdf'));
months = NaT(0, 1, 'TimeZone', 'UTC');
paths = cell(0, 1);
versions = zeros(0, 1);
for ii = 1:numel(candidates)
    token = regexp(candidates(ii).name, ...
        '_(\d{8})_v(\d+)\.cdf$', 'tokens', 'once');
    if isempty(token), continue, end
    month = datetime(token{1}(1:6), 'InputFormat', 'yyyyMM', ...
        'TimeZone', 'UTC');
    fullPath = fullfile(candidates(ii).folder, candidates(ii).name);
    version = str2double(token{2});
    existing = find(months == month, 1);
    if isempty(existing)
        months(end + 1, 1) = month; %#ok<AGROW>
        paths{end + 1, 1} = fullPath; %#ok<AGROW>
        versions(end + 1, 1) = version; %#ok<AGROW>
    elseif version >= versions(existing)
        paths{existing} = fullPath;
        versions(existing) = version;
    end
end
[months, order] = sort(months);
paths = paths(order);
files = table(months, paths, ...
    'VariableNames', {'MonthUTC', 'SourceCDF'});
end

function [product, sources] = loadEventProduct(files, startTime, endTime, opts)
firstMonth = dateshift(startTime, 'start', 'month');
lastMonth = dateshift(endTime - seconds(1), 'start', 'month');
use = files.MonthUTC >= firstMonth & files.MonthUTC <= lastMonth;
sources = files.SourceCDF(use);
if isempty(sources)
    error('VoyagerEvent:NoSourceCDF', ...
        'No COHO CDF overlaps %s to %s.', char(startTime), char(endTime));
end
product = struct;
for ii = 1:numel(sources)
    current = Voyager_Read_CDF_Product(sources{ii}, 'coho', ...
        'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
    current = subsetByTime(current, startTime, endTime);
    product = appendProduct(product, current);
end
product = sortAndDeduplicate(product);
if ~isfield(product, 'Epoch') || isempty(product.Epoch)
    error('VoyagerEvent:NoRecords', ...
        'The selected CDF files contain no records in this event interval.');
end
end

function output = subsetByTime(input, startTime, endTime)
output = input;
if ~isfield(input, 'Epoch'), return, end
mask = input.Epoch >= startTime & input.Epoch < endTime;
recordCount = numel(input.Epoch);
fields = fieldnames(input);
for ii = 1:numel(fields)
    value = input.(fields{ii});
    if (isnumeric(value) || isdatetime(value)) && ~isempty(value) && ...
            size(value, 1) == recordCount
        output.(fields{ii}) = value(mask, :);
    end
end
end

function output = appendProduct(output, input)
if ~isfield(input, 'Epoch') || isempty(input.Epoch), return, end
if ~isfield(output, 'Epoch') || isempty(output.Epoch)
    output = input;
    return
end
oldCount = numel(output.Epoch);
newCount = numel(input.Epoch);
oldFields = numericRecordFields(output, oldCount);
newFields = numericRecordFields(input, newCount);
allFields = union(oldFields, newFields, 'stable');
for ii = 1:numel(allFields)
    field = allFields{ii};
    hasOld = isfield(output, field) && size(output.(field), 1) == oldCount;
    hasNew = isfield(input, field) && size(input.(field), 1) == newCount;
    if hasOld && hasNew
        if size(output.(field), 2) ~= size(input.(field), 2)
            error('VoyagerEvent:ColumnMismatch', ...
                'Column count changed for variable %s.', field);
        end
        output.(field) = [output.(field); input.(field)]; %#ok<AGROW>
    elseif hasOld
        output.(field) = [output.(field); ...
            nan(newCount, size(output.(field), 2))]; %#ok<AGROW>
    elseif hasNew
        output.(field) = [nan(oldCount, size(input.(field), 2)); ...
            input.(field)];
    end
end
output.Epoch = [output.Epoch; input.Epoch];
if isfield(input, 'variable_meta')
    metaFields = fieldnames(input.variable_meta);
    for ii = 1:numel(metaFields)
        output.variable_meta.(metaFields{ii}) = ...
            input.variable_meta.(metaFields{ii});
    end
end
end

function fields = numericRecordFields(product, recordCount)
fields = {};
names = fieldnames(product);
for ii = 1:numel(names)
    value = product.(names{ii});
    if isnumeric(value) && ~isempty(value) && size(value, 1) == recordCount
        fields{end + 1, 1} = names{ii}; %#ok<AGROW>
    end
end
end

function output = sortAndDeduplicate(input)
output = input;
if ~isfield(input, 'Epoch') || isempty(input.Epoch), return, end
[sortedTime, order] = sort(input.Epoch);
[sortedTime, uniqueAt] = unique(sortedTime, 'stable');
order = order(uniqueAt);
recordCount = numel(input.Epoch);
fields = numericRecordFields(input, recordCount);
for ii = 1:numel(fields)
    output.(fields{ii}) = input.(fields{ii})(order, :);
end
output.Epoch = sortedTime;
end

function exportEventFigure(eventRow, data, outputFile, opts, colorLimits)
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [30 30 1800 1200]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 6, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout, 1); hold(ax1, 'on');
field = firstExistingField(data, {'ABS_B', 'F'});
if isempty(field)
    emptyPanel(ax1, 'No recorded |B| values');
else
    plotRecorded(ax1, data.Epoch, data.(field), [0.05 0.05 0.05], ...
        0.9, false, opts.GapBreakHours);
end
ylabel(ax1, {'|B|', '(nT)'}, 'FontSize', 12);
panelLabel(ax1, '(a)');

ax2 = nexttile(layout, 2); hold(ax2, 'on');
componentFields = {'BR', 'BT', 'BN'};
componentLabels = {'B_R', 'B_T', 'B_N'};
componentColors = [0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
labels = cell(0, 1);
for ii = 1:3
    if isfield(data, componentFields{ii}) && ...
            any(isfinite(data.(componentFields{ii})), 'all')
        plotRecorded(ax2, data.Epoch, data.(componentFields{ii}), ...
            componentColors(ii, :), 0.8, false, opts.GapBreakHours);
        labels{end + 1, 1} = componentLabels{ii}; %#ok<AGROW>
    end
end
if isempty(labels)
    emptyPanel(ax2, 'No recorded RTN magnetic components');
else
    yline(ax2, 0, ':', 'Color', [0.45 0.45 0.45], ...
        'HandleVisibility', 'off');
    legend(ax2, labels, 'Location', 'northeast', 'Box', 'off', ...
        'FontSize', 10, 'Interpreter', 'tex', 'Orientation', 'horizontal');
end
ylabel(ax2, {'B_R, B_T, B_N', '(nT)'}, 'FontSize', 12);
panelLabel(ax2, '(b)');

ax3 = nexttile(layout, 3, [2 1]); hold(ax3, 'on');
lecpNames = productVariableNames(data, 'LECP');
plotLECPChannels(ax3, data, lecpNames, opts.GapBreakHours);
ylabel(ax3, {'LECP proton J', ...
    '(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, 'FontSize', 12);
panelLabel(ax3, '(c)');

ax4 = nexttile(layout, 5, [2 1]); hold(ax4, 'on');
plotCRSSpectrogram(ax4, data, eventRow.StartUTC, ...
    eventRow.EndUTCExclusive, colorLimits);
ylabel(ax4, {'CRS proton', 'energy (MeV)'}, 'FontSize', 12);
xlabel(ax4, 'UTC', 'FontSize', 12);
panelLabel(ax4, '(d)');

axesList = [ax1 ax2 ax3 ax4];
startTime = eventRow.StartUTC;
endTime = eventRow.EndUTCExclusive;
left = datenum(startTime); right = datenum(endTime); %#ok<DATNM>
linkaxes(axesList, 'x');
tickTimes = eventTimeTicks(startTime, endTime);
tickValues = datenum(tickTimes); %#ok<DATNM>
tickLabels = cellstr(datestr(tickTimes, 'dd-mmm')); %#ok<DATST>
for ii = 1:numel(axesList)
    xlim(axesList(ii), [left right]);
    set(axesList(ii), 'XTick', tickValues, 'XTickLabel', tickLabels, ...
        'FontSize', 10, 'LineWidth', 0.8, 'TickDir', 'out', ...
        'Box', 'on', 'Layer', 'top');
    grid(axesList(ii), 'on');
    axesList(ii).GridAlpha = 0.18;
    axesList(ii).MinorGridAlpha = 0.08;
end
set([ax1 ax2 ax3], 'XTickLabel', []);
sgtitle(layout, sprintf('Voyager 1  %s   %s--%s UTC', ...
    eventRow.EventID{1}, datestr(startTime, 'yyyy-mm-dd'), ...
    datestr(eventRow.EndUTCInclusive, 'yyyy-mm-dd')), ...
    'FontWeight', 'bold', 'FontSize', 15, 'Interpreter', 'none'); %#ok<DATST>
exportgraphics(fig, outputFile, 'Resolution', opts.ExportDPI);
clear cleanup
end

function success = plotRecorded(ax, time, value, color, width, ...
        positiveOnly, gapBreakHours)
success = false;
if isempty(time) || isempty(value), return, end
time = time(:); value = value(:);
mask = ~isnat(time) & isfinite(value);
if positiveOnly, mask = mask & value > 0; end
time = time(mask); value = value(mask);
if isempty(time), return, end
[time, order] = sort(time); value = value(order);
[time, uniqueAt] = unique(time, 'stable'); value = value(uniqueAt);
x = datenum(time); %#ok<DATNM>
gapBefore = [false; seconds(diff(time)) > gapBreakHours * 3600];
positions = (1:numel(x)).' + cumsum(gapBefore);
xPlot = nan(numel(x) + sum(gapBefore), 1);
yPlot = nan(size(xPlot));
xPlot(positions) = x;
yPlot(positions) = value;
plot(ax, xPlot, yPlot, '-', 'Color', color, 'LineWidth', width);
success = true;
end

function plotLECPChannels(ax, data, names, gapBreakHours)
if isempty(names) || ~isfield(data, 'Epoch')
    emptyPanel(ax, 'No recorded LECP proton flux');
    return
end
colors = [0.2706 0.0000 0.3294; 0.1276 0.5669 0.5506; ...
    0.9932 0.9062 0.1439];
labels = cell(0, 1);
hasData = false;
for ii = 1:numel(names)
    plotted = plotRecorded(ax, data.Epoch, data.(names{ii}), ...
        colors(min(ii, size(colors, 1)), :), 1.0, true, ...
        gapBreakHours);
    if plotted
        labels{end + 1, 1} = conciseEnergyLabel(data, names{ii}); %#ok<AGROW>
        hasData = true;
    end
end
if hasData
    set(ax, 'YScale', 'log');
    legend(ax, labels, 'Location', 'northeast', 'Box', 'off', ...
        'FontSize', 9, 'Interpreter', 'none');
else
    emptyPanel(ax, 'No recorded LECP proton flux');
end
end

function plotCRSSpectrogram(ax, data, startTime, endTime, colorLimits)
names = productVariableNames(data, 'CRS');
if isempty(names) || ~isfield(data, 'Epoch') || isempty(data.Epoch)
    emptyPanel(ax, 'No recorded CRS proton flux');
    return
end
[energyLow, energyHigh] = energyBounds(data, names);
validChannels = isfinite(energyLow) & isfinite(energyHigh) & ...
    energyLow > 0 & energyHigh > energyLow;
names = names(validChannels);
energyLow = energyLow(validChannels);
energyHigh = energyHigh(validChannels);
if isempty(names)
    emptyPanel(ax, 'CRS energy metadata are unavailable');
    return
end

timeGrid = (startTime:hours(1):(endTime - hours(1))).';
flux = nan(numel(timeGrid), numel(names));
[isRecorded, location] = ismember(data.Epoch, timeGrid);
for ii = 1:numel(names)
    values = data.(names{ii});
    flux(location(isRecorded), ii) = values(isRecorded);
end
flux(~isfinite(flux) | flux <= 0) = nan;
logFlux = log10(flux);
timeEdges = datenum([timeGrid; endTime]); %#ok<DATNM>
for ii = 1:numel(names)
    xData = [timeEdges.'; timeEdges.'];
    yData = [repmat(energyLow(ii), 1, numel(timeEdges)); ...
        repmat(energyHigh(ii), 1, numel(timeEdges))];
    cRow = [logFlux(:, ii).', nan];
    cData = [cRow; cRow];
    surface(ax, xData, yData, zeros(size(xData)), cData, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'CDataMapping', 'scaled');
end
set(ax, 'YScale', 'log', 'YDir', 'normal');
ylim(ax, [min(energyLow) max(energyHigh)]);
yticks(ax, [3 10 30 100 300]);
colormap(ax, turbo(256));
caxis(ax, colorLimits);
cb = colorbar(ax, 'eastoutside');
cb.FontSize = 11;
cb.Label.String = {'log_{10} J', ...
    'cm^{-2} s^{-1} sr^{-1} MeV^{-1}'};
cb.Label.Interpreter = 'tex';
cb.Label.FontSize = 12;
end

function ticks = eventTimeTicks(startTime, endTime)
spanDays = days(endTime - startTime);
if spanDays <= 10
    step = 1;
elseif spanDays <= 22
    step = 2;
else
    step = 3;
end
ticks = startTime:caldays(step):endTime;
if ticks(end) < endTime
    ticks(end + 1) = endTime; %#ok<AGROW>
end
end

function panelLabel(ax, label)
text(ax, 0.009, 0.93, label, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.1 0.1 0.1]);
end

function writeEventDatabase(databaseFile, catalog, report, eventData, ...
        sourceFiles, colorLimits, opts)
if isfile(databaseFile)
    if ~logical(opts.OverwriteDatabase)
        error('VoyagerEvent:DatabaseExists', ...
            ['Database already exists: %s. Set OverwriteDatabase=true ', ...
            'to replace it.'], databaseFile);
    end
    conn = sqlite(databaseFile, 'connect');
else
    conn = sqlite(databaseFile, 'create');
end
cleanup = onCleanup(@() close(conn));
exec(conn, 'PRAGMA foreign_keys = OFF');
exec(conn, 'DROP TABLE IF EXISTS particle_flux');
exec(conn, 'DROP TABLE IF EXISTS magnetic_field');
exec(conn, 'DROP TABLE IF EXISTS event_sources');
exec(conn, 'DROP TABLE IF EXISTS events');
exec(conn, 'DROP TABLE IF EXISTS provenance');
exec(conn, 'DROP TABLE IF EXISTS literature_style_references');
exec(conn, 'PRAGMA foreign_keys = ON');
exec(conn, ['CREATE TABLE events (' ...
    'event_id TEXT PRIMARY KEY, spacecraft INTEGER NOT NULL, ' ...
    'year_utc INTEGER NOT NULL, source_slide INTEGER NOT NULL, ' ...
    'circle_order INTEGER NOT NULL, ppt_circle_color TEXT, ' ...
    'ppt_circle_x_emu INTEGER, ppt_circle_width_emu INTEGER, ' ...
    'ppt_mapped_start_utc TEXT, ppt_mapped_end_utc TEXT, ' ...
    'start_utc TEXT NOT NULL, end_utc_inclusive TEXT NOT NULL, ' ...
    'end_utc_exclusive TEXT NOT NULL, source_pptx TEXT, ' ...
    'extraction_method TEXT, figure_file TEXT, status TEXT, notes TEXT)']);
exec(conn, ['CREATE TABLE event_sources (' ...
    'event_id TEXT NOT NULL, source_file TEXT NOT NULL, ' ...
    'source_product TEXT NOT NULL, PRIMARY KEY(event_id, source_file), ' ...
    'FOREIGN KEY(event_id) REFERENCES events(event_id))']);
exec(conn, ['CREATE TABLE magnetic_field (' ...
    'event_id TEXT NOT NULL, epoch_utc TEXT NOT NULL, ' ...
    'abs_b_nt REAL, br_nt REAL, bt_nt REAL, bn_nt REAL, ' ...
    'heliocentric_distance_au REAL, ' ...
    'PRIMARY KEY(event_id, epoch_utc), ' ...
    'FOREIGN KEY(event_id) REFERENCES events(event_id))']);
exec(conn, ['CREATE TABLE particle_flux (' ...
    'event_id TEXT NOT NULL, epoch_utc TEXT NOT NULL, ' ...
    'instrument TEXT NOT NULL, channel INTEGER NOT NULL, species TEXT, ' ...
    'energy_low_mev REAL, energy_high_mev REAL, channel_label TEXT, ' ...
    'differential_flux REAL, unit TEXT, native_cadence_hours REAL, ' ...
    'is_missing INTEGER NOT NULL, ' ...
    'PRIMARY KEY(event_id, epoch_utc, instrument, channel), ' ...
    'FOREIGN KEY(event_id) REFERENCES events(event_id))']);
exec(conn, ['CREATE TABLE provenance (' ...
    'key TEXT PRIMARY KEY, value TEXT NOT NULL)']);
exec(conn, ['CREATE TABLE literature_style_references (' ...
    'reference_id TEXT PRIMARY KEY, citation TEXT NOT NULL, ' ...
    'doi TEXT, url TEXT, plotting_feature_used TEXT)']);

eventsTable = table(catalog.EventID, ones(height(catalog), 1), ...
    catalog.YearUTC, catalog.SourceSlide, catalog.CircleOrderOnSlide, ...
    catalog.PPTCircleColor, catalog.PPTCircleX_EMU, ...
    catalog.PPTCircleWidth_EMU, isoUtc(catalog.PPTMappedStartUTC), ...
    isoUtc(catalog.PPTMappedEndUTC), isoUtc(catalog.StartUTC), ...
    isoUtc(catalog.EndUTCInclusive), isoUtc(catalog.EndUTCExclusive), ...
    catalog.SourcePPTX, catalog.ExtractionMethod, report.FigureFile, ...
    report.Status, report.Notes, 'VariableNames', ...
    {'event_id', 'spacecraft', 'year_utc', 'source_slide', ...
    'circle_order', 'ppt_circle_color', 'ppt_circle_x_emu', ...
    'ppt_circle_width_emu', 'ppt_mapped_start_utc', ...
    'ppt_mapped_end_utc', 'start_utc', 'end_utc_inclusive', ...
    'end_utc_exclusive', 'source_pptx', 'extraction_method', ...
    'figure_file', 'status', 'notes'});
sqlwrite(conn, 'events', eventsTable);

try
    for ii = 1:height(catalog)
        eventID = catalog.EventID{ii};
        if ~isempty(sourceFiles{ii})
            count = numel(sourceFiles{ii});
            sourceTable = table(repmat({eventID}, count, 1), ...
                sourceFiles{ii}(:), repmat( ...
                {'Voyager 1 COHO one-hour merged MAG/particle CDF'}, ...
                count, 1), 'VariableNames', ...
                {'event_id', 'source_file', 'source_product'});
            sqlwrite(conn, 'event_sources', sourceTable);
        end
        data = eventData{ii};
        if ~isfield(data, 'Epoch') || isempty(data.Epoch), continue, end
        n = numel(data.Epoch);
        magneticTable = table(repmat({eventID}, n, 1), ...
            isoUtc(data.Epoch), getRecordField(data, {'ABS_B', 'F'}, n), ...
            getRecordField(data, {'BR'}, n), ...
            getRecordField(data, {'BT'}, n), ...
            getRecordField(data, {'BN'}, n), ...
            getRecordField(data, {'heliocentricDistance'}, n), ...
            'VariableNames', {'event_id', 'epoch_utc', 'abs_b_nt', ...
            'br_nt', 'bt_nt', 'bn_nt', 'heliocentric_distance_au'});
        sqlwrite(conn, 'magnetic_field', magneticTable);
        particleTable = particleLongTable(eventID, data);
        if ~isempty(particleTable)
            sqlwrite(conn, 'particle_flux', particleTable);
        end
    end
catch ME
    rethrow(ME)
end

created = datetime('now', 'TimeZone', 'UTC');
provenanceKeys = {'created_utc'; 'source_pptx'; 'source_data_product'; ...
    'mag_master_cadence'; 'crs_documented_cadence'; ...
    'value_processing'; 'gap_handling'; 'event_mapping'; ...
    'crs_color_limits_log10'; 'program_file'};
provenanceValues = {isoUtcScalar(created); opts.SourcePptx; ...
    'Voyager 1 COHO one-hour merged MAG/particle CDF'; ...
    '1 hour'; '6 hours in the CDF CATDESC metadata'; ...
    'Original recorded values; no smoothing, interpolation, averaging, resampling, or gap filling'; ...
    'Missing values are SQL NULL and plot as line breaks or blank spectrogram cells'; ...
    ['Ellipse x bounds on slides 4-13 mapped to PNG axis pixels ', ...
    '118-2443 of a 2783-pixel annual figure; bounds rounded outward to full UTC days']; ...
    sprintf('[%.3f, %.3f]', colorLimits(1), colorLimits(2)); ...
    mfilename('fullpath')};
provenanceTable = table(provenanceKeys, provenanceValues, ...
    'VariableNames', {'key', 'value'});
sqlwrite(conn, 'provenance', provenanceTable);

referenceID = {'Gurnett2015'; 'Gurnett2021'; 'Hill2020'};
citation = {['Gurnett et al. (2015), Precursors to Interstellar ', ...
    'Shocks of Solar Origin, ApJ 809:121']; ...
    ['Gurnett et al. (2021), A Foreshock Model for Interstellar ', ...
    'Shocks of Solar Origin, AJ 161:11']; ...
    ['Hill et al. (2020), Influence of Solar Disturbances on ', ...
    'Galactic Cosmic Rays..., ApJ 905:69']};
doi = {'10.1088/0004-637X/809/2/121'; ...
    '10.3847/1538-3881/abc337'; '10.3847/1538-4357/abb408'};
url = {'https://doi.org/10.1088/0004-637X/809/2/121'; ...
    'https://doi.org/10.3847/1538-3881/abc337'; ...
    'https://doi.org/10.3847/1538-4357/abb408'};
feature = {['Compact event-aligned magnetic-field and selected ', ...
    'particle time profiles']; ...
    'Common UTC axis, panel labels, uncluttered line profiles'; ...
    'Logarithmic energetic-particle intensity presentation across locations'};
referenceTable = table(referenceID, citation, doi, url, feature, ...
    'VariableNames', {'reference_id', 'citation', 'doi', 'url', ...
    'plotting_feature_used'});
sqlwrite(conn, 'literature_style_references', referenceTable);
exec(conn, 'CREATE INDEX idx_mag_epoch ON magnetic_field(epoch_utc)');
exec(conn, ['CREATE INDEX idx_flux_epoch_instrument ON ', ...
    'particle_flux(epoch_utc, instrument)']);
clear cleanup
end

function output = particleLongTable(eventID, data)
output = table;
if ~isfield(data, 'Epoch') || isempty(data.Epoch), return, end
epochText = isoUtc(data.Epoch);
n = numel(data.Epoch);
instruments = {'LECP', 'CRS'};
for kk = 1:numel(instruments)
    instrument = instruments{kk};
    names = productVariableNames(data, instrument);
    for ii = 1:numel(names)
        values = data.(names{ii});
        [low, high] = energyBounds(data, names(ii));
        label = variableLabel(data, names{ii});
        unit = variableUnit(data, names{ii});
        cadence = variableCadenceHours(data, names{ii});
        block = table(repmat({eventID}, n, 1), epochText, ...
            repmat({instrument}, n, 1), repmat(ii, n, 1), ...
            repmat({'H'}, n, 1), repmat(low, n, 1), ...
            repmat(high, n, 1), repmat({label}, n, 1), values(:), ...
            repmat({unit}, n, 1), repmat(cadence, n, 1), ...
            double(~isfinite(values(:)) | values(:) <= 0), ...
            'VariableNames', {'event_id', 'epoch_utc', 'instrument', ...
            'channel', 'species', 'energy_low_mev', 'energy_high_mev', ...
            'channel_label', 'differential_flux', 'unit', ...
            'native_cadence_hours', 'is_missing'});
        block.differential_flux(block.is_missing == 1) = nan;
        output = [output; block]; %#ok<AGROW>
    end
end
end

function value = getRecordField(data, candidates, count)
value = nan(count, 1);
for ii = 1:numel(candidates)
    if isfield(data, candidates{ii}) && ...
            size(data.(candidates{ii}), 1) == count
        value = data.(candidates{ii})(:, 1);
        return
    end
end
end

function textValues = isoUtc(time)
if isempty(time)
    textValues = cell(0, 1);
    return
end
time = time(:);
time.TimeZone = 'UTC';
textValues = cellstr(datestr(time, 'yyyy-mm-ddTHH:MM:SSZ')); %#ok<DATST>
end

function value = isoUtcScalar(time)
values = isoUtc(time);
value = values{1};
end

function names = productVariableNames(product, suffix)
names = {};
if isempty(product), return, end
fields = fieldnames(product);
number = [];
pattern = ['^protonFlux(\d+)_', suffix, '$'];
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

function [low, high] = energyBounds(product, names)
low = nan(numel(names), 1);
high = nan(numel(names), 1);
for ii = 1:numel(names)
    label = variableLabel(product, names{ii});
    token = regexp(label, ...
        'H\s+([0-9.]+)\s*-\s*([0-9.]+)\s*MeV', 'tokens', 'once');
    if ~isempty(token)
        low(ii) = str2double(token{1});
        high(ii) = str2double(token{2});
    end
end
end

function label = conciseEnergyLabel(product, name)
[low, high] = energyBounds(product, {name});
if isfinite(low) && isfinite(high)
    label = sprintf('%g-%g MeV', low, high);
else
    label = variableLabel(product, name);
end
end

function label = variableLabel(product, name)
label = name;
attrs = variableAttributes(product, name);
candidates = {'FIELDNAM', 'LABLAXIS', 'CATDESC'};
for ii = 1:numel(candidates)
    if isfield(attrs, candidates{ii})
        raw = attrs.(candidates{ii});
        if ischar(raw) || (isstring(raw) && isscalar(raw))
            label = char(raw);
            return
        end
    end
end
end

function unit = variableUnit(product, name)
unit = '1/(cm^2 sec ster MeV)';
attrs = variableAttributes(product, name);
if isfield(attrs, 'UNITS')
    raw = attrs.UNITS;
    if ischar(raw) || (isstring(raw) && isscalar(raw))
        unit = char(raw);
    end
end
end

function cadence = variableCadenceHours(product, name)
cadence = 1;
attrs = variableAttributes(product, name);
if isfield(attrs, 'CATDESC')
    raw = char(string(attrs.CATDESC));
    token = regexp(raw, '(\d+(?:\.\d+)?)\s*-\s*hr', ...
        'tokens', 'once', 'ignorecase');
    if ~isempty(token)
        cadence = str2double(token{1});
    end
end
end

function attrs = variableAttributes(product, name)
attrs = struct;
if ~isfield(product, 'variable_meta') || ...
        ~isfield(product.variable_meta, name)
    return
end
item = product.variable_meta.(name);
if isfield(item, 'attributes')
    attrs = item.attributes;
end
end

function field = firstExistingField(product, candidates)
field = '';
for ii = 1:numel(candidates)
    if isfield(product, candidates{ii}) && ...
            any(isfinite(product.(candidates{ii})), 'all')
        field = candidates{ii};
        return
    end
end
end

function count = variableRecordCount(product, field)
count = 0;
if isfield(product, field), count = size(product.(field), 1); end
end

function count = anyFiniteRecordCount(product, names)
mask = false(variableRecordCount(product, 'Epoch'), 1);
for ii = 1:numel(names)
    if isfield(product, names{ii}) && ...
            size(product.(names{ii}), 1) == numel(mask)
        mask = mask | any(isfinite(product.(names{ii})), 2);
    end
end
count = sum(mask);
end

function count = validValueCount(product, names)
count = 0;
for ii = 1:numel(names)
    value = product.(names{ii});
    count = count + nnz(isfinite(value) & value > 0);
end
end

function limits = robustColorLimits(values)
values = sort(values(isfinite(values)));
if isempty(values)
    limits = [-5 -1];
    return
end
n = numel(values);
lo = values(max(1, round(0.01 * n)));
hi = values(min(n, max(1, round(0.99 * n))));
lo = floor(lo * 4) / 4;
hi = ceil(hi * 4) / 4;
if hi - lo < 1
    middle = (lo + hi) / 2;
    lo = middle - 0.5;
    hi = middle + 0.5;
end
limits = [lo hi];
end

function product = emptyProduct()
product = struct;
product.Epoch = NaT(0, 1, 'TimeZone', 'UTC');
product.variable_meta = struct;
end

function emptyPanel(ax, message)
text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Color', [0.45 0.45 0.45], 'FontAngle', 'italic', ...
    'Interpreter', 'none');
end

function ensureFolder(folder)
if isfolder(folder), return, end
[ok, message] = mkdir(folder);
if ~ok
    error('VoyagerEvent:FolderCreateFailed', ...
        'Unable to create %s: %s', folder, message);
end
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = isLogicalScalar(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isfinite(double(value)) && ismember(double(value), [0 1]);
end
