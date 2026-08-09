function report = Voyager_Plot_COHO1h_Raw_Monthly(varargin)
%Voyager_Plot_COHO1h_Raw_Monthly Plot uniform raw COHO one-hour products.
%   One flat PNG is written for each available post-heliopause month. All
%   panels use values contained in the same COHO one-hour instrument file.
%   Missing samples are omitted before plotting, so adjacent measured points
%   form a continuous trace. No interpolation or resampling is performed.

parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'OutputFolder', fullfile(pwd, ...
    'UniformLowPrecision_Monthly_Raw'), @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'StartTime', [], @isTimeInput);
addParameter(parser, 'EndTime', [], @isTimeInput);
addParameter(parser, 'Overwrite', false, @isLogicalScalar);
addParameter(parser, 'Visible', false, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 150, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 72);
addParameter(parser, 'PythonExe', '', @isTextScalar);
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.OutputFolder = char(opts.OutputFolder);
opts.CacheRoot = char(opts.CacheRoot);
opts.Spacecraft = unique(double(opts.Spacecraft(:).'), 'stable');
if ~isfolder(opts.DataRoot)
    error('VoyagerLow:DataRootMissing', 'Data root is unavailable: %s', opts.DataRoot);
end
if isempty(opts.Spacecraft) || any(~ismember(opts.Spacecraft, [1 2]))
    error('VoyagerLow:Spacecraft', 'Spacecraft must contain 1 and/or 2.');
end
if ~isfolder(opts.OutputFolder)
    [ok, message] = mkdir(opts.OutputFolder);
    if ~ok, error('VoyagerLow:OutputCreateFailed', '%s', message); end
end

crossing = [datetime(2012, 8, 25, 'TimeZone', 'UTC'), ...
    datetime(2018, 11, 5, 'TimeZone', 'UTC')];
requestedStart = parseOptionalUtc(opts.StartTime);
requestedEnd = parseOptionalUtc(opts.EndTime);
rows = cell(0, 17);

for spacecraft = opts.Spacecraft
    firstTime = crossing(spacecraft);
    if ~isempty(requestedStart), firstTime = max(firstTime, requestedStart); end
    files = listMonthlyCoho(opts.DataRoot, spacecraft);
    for ii = 1:height(files)
        monthStart = files.MonthUTC(ii);
        monthEnd = monthStart + calmonths(1);
        if monthEnd <= firstTime, continue, end
        if ~isempty(requestedEnd) && monthStart >= requestedEnd, continue, end
        dataStart = max(monthStart, firstTime);
        dataEnd = monthEnd;
        if ~isempty(requestedEnd), dataEnd = min(dataEnd, requestedEnd); end
        if dataEnd <= dataStart, continue, end

        sourceFile = files.SourceCDF{ii};
        monthText = char(string(monthStart, 'yyyy-MM'));
        outputFile = fullfile(opts.OutputFolder, sprintf( ...
            'V%d_COHO1h_raw_%s.png', spacecraft, monthText));
        fprintf('[Uniform low precision] Voyager %d %s\n', spacecraft, monthText);
        status = 'ok';
        note = '';
        product = struct;
        try
            product = Voyager_Read_CDF_Product(sourceFile, 'coho', ...
                'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
            product = subsetByTime(product, dataStart, dataEnd);
            if logical(opts.Overwrite) || ~isfile(outputFile)
                exportRawMonthFigure(spacecraft, monthStart, monthEnd, ...
                    product, outputFile, opts);
            else
                status = 'existing';
            end
        catch ME
            status = 'error';
            note = sprintf('%s: %s', ME.identifier, ME.message);
            warning('VoyagerLow:MonthFailed', 'V%d %s: %s', ...
                spacecraft, monthText, note);
        end

        epochRecords = variableRecordCount(product, 'Epoch');
        magRecords = anyFiniteRecordCount(product, {'ABS_B', 'F', 'BR', 'BT', 'BN'});
        plasmaRecords = anyFiniteRecordCount(product, ...
            {'V', 'protonDensity', 'protonTemp'});
        lecpNames = productVariableNames(product, 'LECP');
        crsNames = productVariableNames(product, 'CRS');
        lecpValues = validValueCount(product, lecpNames);
        crsValues = validValueCount(product, crsNames);
        rows(end + 1, :) = {spacecraft, monthStart, dataStart, dataEnd, ...
            sourceFile, epochRecords, magRecords, plasmaRecords, ...
            numel(lecpNames), lecpValues, numel(crsNames), crsValues, ...
            true, 0, outputFile, status, note}; %#ok<AGROW>
    end
end

names = {'Spacecraft', 'MonthUTC', 'DataStartUTC', 'DataEndUTC', ...
    'SourceCDF', 'EpochRecords', 'MagneticRecords', 'PlasmaRecords', ...
    'LECPChannels', 'LECPValidValues', 'CRSChannels', 'CRSValidValues', ...
    'OriginalInstrumentValuesOnly', 'GeneratedOrInterpolatedValues', ...
    'FigureFile', 'Status', 'Notes'};
report = cell2table(rows, 'VariableNames', names);
writetable(report, fullfile(opts.OutputFolder, ...
    'UniformLowPrecision_Monthly_Raw_manifest.csv'));
save(fullfile(opts.OutputFolder, 'UniformLowPrecision_Monthly_Raw_report.mat'), ...
    'report', 'opts', '-v7.3');
end

function files = listMonthlyCoho(dataRoot, spacecraft)
root = fullfile(dataRoot, sprintf('voyager%d', spacecraft), ...
    'coho', '1hr', 'l2', 'merged_mag_plasma');
if ~isfolder(root)
    files = table(NaT(0, 1, 'TimeZone', 'UTC'), cell(0, 1), ...
        'VariableNames', {'MonthUTC', 'SourceCDF'});
    return
end
candidates = dir(fullfile(root, '*', '*', '*.cdf'));
months = NaT(0, 1, 'TimeZone', 'UTC');
paths = cell(0, 1);
versions = zeros(0, 1);
for ii = 1:numel(candidates)
    token = regexp(candidates(ii).name, '_(\d{8})_v(\d+)\.cdf$', ...
        'tokens', 'once');
    if isempty(token), continue, end
    month = datetime(token{1}(1:6), 'InputFormat', 'yyyyMM', ...
        'TimeZone', 'UTC');
    fullPath = fullfile(candidates(ii).folder, candidates(ii).name);
    existing = find(months == month, 1);
    version = str2double(token{2});
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
files = table(months, paths, 'VariableNames', {'MonthUTC', 'SourceCDF'});
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

function exportRawMonthFigure(spacecraft, monthStart, monthEnd, data, ...
        outputFile, opts)
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [30 30 1800 1650]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 8, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout); hold(ax1, 'on');
field = firstExistingField(data, {'ABS_B', 'F'});
if isempty(field)
    emptyPanel(ax1, '|B| -- no recorded value');
else
    plotRecorded(ax1, data.Epoch, data.(field), [0.05 0.05 0.05], 0.8, false);
    title(ax1, '|B| -- raw COHO one-hour records');
end
ylabel(ax1, 'nT'); grid(ax1, 'on');

ax2 = nexttile(layout); hold(ax2, 'on');
colors = lines(3);
hasComponents = false;
componentFields = {'BR', 'BT', 'BN'};
for ii = 1:3
    if isfield(data, componentFields{ii})
        hasComponents = plotRecorded(ax2, data.Epoch, ...
            data.(componentFields{ii}), colors(ii, :), 0.7, false) || hasComponents;
    end
end
if hasComponents
    legend(ax2, {'B_R', 'B_T', 'B_N'}, 'Location', 'eastoutside');
    title(ax2, 'RTN magnetic components -- raw COHO one-hour records');
else
    emptyPanel(ax2, 'RTN components -- no recorded value');
end
ylabel(ax2, 'nT'); grid(ax2, 'on');

ax3 = nexttile(layout); hold(ax3, 'on');
if isfield(data, 'V') && plotRecorded(ax3, data.Epoch, data.V, ...
        [0.1 0.35 0.8], 0.8, false)
    title(ax3, 'Plasma bulk speed -- raw COHO one-hour records');
else
    emptyPanel(ax3, 'Plasma bulk speed -- no recorded value');
end
ylabel(ax3, 'km s^{-1}'); grid(ax3, 'on');

ax4 = nexttile(layout); hold(ax4, 'on');
if isfield(data, 'protonDensity') && plotRecorded(ax4, data.Epoch, ...
        data.protonDensity, [0.1 0.55 0.25], 0.8, true)
    title(ax4, 'Proton density -- raw COHO one-hour records');
else
    emptyPanel(ax4, 'Proton density -- no recorded value');
end
ylabel(ax4, 'cm^{-3}'); grid(ax4, 'on');

ax5 = nexttile(layout); hold(ax5, 'on');
if isfield(data, 'protonTemp') && plotRecorded(ax5, data.Epoch, ...
        data.protonTemp, [0.85 0.25 0.1], 0.8, true)
    title(ax5, 'Proton temperature -- raw COHO one-hour records');
else
    emptyPanel(ax5, 'Proton temperature -- no recorded value');
end
ylabel(ax5, 'K'); grid(ax5, 'on');

ax6 = nexttile(layout); hold(ax6, 'on');
lecpNames = productVariableNames(data, 'LECP');
plotParticleChannels(ax6, data, lecpNames, 'LECP');
ylabel(ax6, 'flux'); grid(ax6, 'on');

ax7 = nexttile(layout); hold(ax7, 'on');
crsNames = productVariableNames(data, 'CRS');
plotParticleChannels(ax7, data, crsNames, 'CRS');
ylabel(ax7, 'flux'); grid(ax7, 'on');

ax8 = nexttile(layout); hold(ax8, 'on');
if isfield(data, 'heliocentricDistance') && plotRecorded(ax8, data.Epoch, ...
        data.heliocentricDistance, [0.35 0.15 0.65], 0.9, false)
    title(ax8, 'Heliocentric distance -- raw COHO one-hour records');
else
    emptyPanel(ax8, 'Heliocentric distance -- no recorded value');
end
ylabel(ax8, 'AU'); grid(ax8, 'on'); xlabel(ax8, 'Day of month (UTC)');

axesList = [ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8];
left = datenum(monthStart); right = datenum(monthEnd); %#ok<DATNM>
for ii = 1:numel(axesList)
    xlim(axesList(ii), [left right]);
    datetick(axesList(ii), 'x', 'dd', 'keeplimits'); %#ok<DATIC>
end
linkaxes(axesList, 'x');
sgtitle(layout, sprintf(['Voyager %d uniform low-precision monthly data | ' ...
    '%s UTC | COHO one-hour instrument product'], spacecraft, ...
    char(string(monthStart, 'yyyy-MM'))), 'FontWeight', 'bold', 'FontSize', 14);
annotation(fig, 'textbox', [0.01 0.001 0.98 0.022], ...
    'String', ['Original instrument samples only; no fillmissing, interp1, ' ...
    'retime, resampling, or generated values. Dots mark recorded samples; ' ...
    'lines connect consecutive available records.'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'Interpreter', 'none');
exportgraphics(fig, outputFile, 'Resolution', opts.ExportDPI);
clear cleanup
end

function success = plotRecorded(ax, time, value, color, width, positiveOnly)
success = false;
if isempty(time) || isempty(value), return, end
time = time(:); value = value(:);
mask = ~isnat(time) & isfinite(value);
if positiveOnly, mask = mask & value > 0; end
time = time(mask); value = value(mask);
if isempty(time), return, end
[time, order] = sort(time); value = value(order);
[time, uniqueAt] = unique(time, 'stable'); value = value(uniqueAt);
plot(ax, datenum(time), value, '-', 'Color', color, ... %#ok<DATNM>
    'LineWidth', width, 'Marker', '.', 'MarkerSize', 4);
success = true;
end

function plotParticleChannels(ax, data, names, instrument)
if isempty(names) || ~isfield(data, 'Epoch')
    emptyPanel(ax, sprintf('%s proton flux -- no recorded value', instrument));
    return
end
colors = turbo(max(numel(names), 2));
labels = cell(numel(names), 1);
hasData = false;
for ii = 1:numel(names)
    labels{ii} = variableLabel(data, names{ii});
    hasData = plotRecorded(ax, data.Epoch, data.(names{ii}), ...
        colors(ii, :), 0.65, true) || hasData;
end
if hasData
    set(ax, 'YScale', 'log');
    legend(ax, labels, 'Location', 'eastoutside', 'Interpreter', 'none', ...
        'FontSize', 6);
    title(ax, sprintf(['%s proton flux -- %d raw one-hour channels; ' ...
        'recorded samples only'], instrument, numel(names)));
else
    emptyPanel(ax, sprintf('%s proton flux -- no recorded value', instrument));
end
end

function names = productVariableNames(product, suffix)
names = {};
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
[~, order] = sort(number); names = names(order);
end

function label = variableLabel(product, name)
label = name;
if ~isfield(product, 'variable_meta') || ...
        ~isfield(product.variable_meta, name), return, end
item = product.variable_meta.(name);
if ~isfield(item, 'attributes'), return, end
candidates = {'FIELDNAM', 'LABLAXIS', 'CATDESC'};
for ii = 1:numel(candidates)
    if isfield(item.attributes, candidates{ii})
        raw = item.attributes.(candidates{ii});
        if ischar(raw) || (isstring(raw) && isscalar(raw))
            label = char(raw); return
        end
    end
end
end

function field = firstExistingField(product, candidates)
field = '';
for ii = 1:numel(candidates)
    if isfield(product, candidates{ii}) && ...
            any(isfinite(product.(candidates{ii})), 'all')
        field = candidates{ii}; return
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
    if isfield(product, names{ii}) && size(product.(names{ii}), 1) == numel(mask)
        mask = mask | any(isfinite(product.(names{ii})), 2);
    end
end
count = sum(mask);
end

function count = validValueCount(product, names)
count = 0;
for ii = 1:numel(names)
    count = count + sum(isfinite(product.(names{ii})), 'all');
end
end

function emptyPanel(ax, message)
text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45], ...
    'FontAngle', 'italic', 'Interpreter', 'none');
title(ax, message, 'Interpreter', 'none');
end

function value = parseOptionalUtc(input)
if isempty(input), value = []; return, end
if isdatetime(input), value = input; else, value = datetime(input, 'TimeZone', 'UTC'); end
if isempty(value) || ~isscalar(value) || isnat(value)
    error('VoyagerLow:Time', 'Time inputs must be valid scalar values.');
end
value.TimeZone = 'UTC';
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
