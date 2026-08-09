function report = Voyager_Plot_COHO1h_Raw_Yearly(varargin)
%Voyager_Plot_COHO1h_Raw_Yearly Plot annual raw COHO one-hour products.
%   One flat PNG is written for each available post-heliopause year. All
%   Values come from monthly COHO files, with Voyager 2 MAG after 2020
%   supplied by the NASA COHOWeb one-hour listing. No averaging,
%   interpolation, resampling, or gap filling is performed. Lines break
%   across missing one-hour records.

programRoot = fileparts(mfilename('fullpath'));
parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'OutputFolder', fullfile(pwd, ...
    'UniformLowPrecision_Yearly_Raw'), @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'StartTime', [], @isTimeInput);
addParameter(parser, 'EndTime', [], @isTimeInput);
addParameter(parser, 'Overwrite', false, @isLogicalScalar);
addParameter(parser, 'Visible', false, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 150, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 72);
addParameter(parser, 'GapBreakHours', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
addParameter(parser, 'PythonExe', '', @isTextScalar);
addParameter(parser, 'NasaV2MagFile', fullfile(programRoot, ...
    'NASA_COHO_Voyager2_MAG_1h_2021_2024.html'), @isTextScalar);
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.OutputFolder = char(opts.OutputFolder);
opts.CacheRoot = char(opts.CacheRoot);
opts.NasaV2MagFile = char(opts.NasaV2MagFile);
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
rows = cell(0, 19);

for spacecraft = opts.Spacecraft
    firstTime = crossing(spacecraft);
    if ~isempty(requestedStart), firstTime = max(firstTime, requestedStart); end
    files = listMonthlyCoho(opts.DataRoot, spacecraft);
    nasaProduct = struct;
    if spacecraft == 2 && isfile(opts.NasaV2MagFile)
        nasaProduct = readNasaV2CohoMag(opts.NasaV2MagFile);
    end
    availableYears = zeros(0, 1);
    if ~isempty(files), availableYears = unique(year(files.MonthUTC)); end
    if isfield(nasaProduct, 'Epoch') && ~isempty(nasaProduct.Epoch)
        availableYears = unique([availableYears; year(nasaProduct.Epoch)]);
    end
    if isempty(availableYears), continue, end
    for ii = 1:numel(availableYears)
        yearNumber = availableYears(ii);
        yearStart = datetime(yearNumber, 1, 1, 'TimeZone', 'UTC');
        yearEnd = datetime(yearNumber + 1, 1, 1, 'TimeZone', 'UTC');
        if yearEnd <= firstTime, continue, end
        if ~isempty(requestedEnd) && yearStart >= requestedEnd, continue, end
        dataStart = max(yearStart, firstTime);
        dataEnd = yearEnd;
        if ~isempty(requestedEnd), dataEnd = min(dataEnd, requestedEnd); end
        if dataEnd <= dataStart, continue, end

        useNasaMag = spacecraft == 2 && yearNumber > 2020 && ...
            isfield(nasaProduct, 'Epoch') && ...
            any(year(nasaProduct.Epoch) == yearNumber);
        if useNasaMag
            sourceFiles = {opts.NasaV2MagFile};
            sourceProduct = 'NASA COHOWeb one-hour MAG';
            outputFile = fullfile(opts.OutputFolder, sprintf( ...
                'V%d_COHOweb1h_MAG_raw_%04d.png', spacecraft, yearNumber));
        else
            inYear = year(files.MonthUTC) == yearNumber;
            sourceFiles = files.SourceCDF(inYear);
            sourceProduct = 'COHO one-hour merged instrument product';
            outputFile = fullfile(opts.OutputFolder, sprintf( ...
                'V%d_COHO1h_raw_%04d.png', spacecraft, yearNumber));
        end
        fprintf('[Annual one-hour raw] Voyager %d %04d\n', spacecraft, yearNumber);
        status = 'ok';
        note = '';
        product = struct;
        try
            if useNasaMag
                product = subsetByTime(nasaProduct, dataStart, dataEnd);
            else
                for jj = 1:numel(sourceFiles)
                    current = Voyager_Read_CDF_Product(sourceFiles{jj}, 'coho', ...
                        'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
                    current = subsetByTime(current, dataStart, dataEnd);
                    product = appendProduct(product, current);
                end
            end
            product = sortAndDeduplicate(product);
            if logical(opts.Overwrite) || ~isfile(outputFile)
                exportRawYearFigure(spacecraft, yearStart, yearEnd, ...
                    product, outputFile, opts, sourceProduct);
            else
                status = 'existing';
            end
        catch ME
            status = 'error';
            note = sprintf('%s: %s', ME.identifier, ME.message);
            warning('VoyagerLow:YearFailed', 'V%d %04d: %s', ...
                spacecraft, yearNumber, note);
        end

        epochRecords = variableRecordCount(product, 'Epoch');
        magRecords = anyFiniteRecordCount(product, {'ABS_B', 'F', 'BR', 'BT', 'BN'});
        plasmaRecords = anyFiniteRecordCount(product, ...
            {'V', 'protonDensity', 'protonTemp'});
        lecpNames = productVariableNames(product, 'LECP');
        crsNames = productVariableNames(product, 'CRS');
        lecpValues = validValueCount(product, lecpNames);
        crsValues = validValueCount(product, crsNames);
        rows(end + 1, :) = {spacecraft, yearStart, dataStart, dataEnd, ...
            sourceProduct, numel(sourceFiles), strjoin(sourceFiles, ';'), ...
            epochRecords, magRecords, plasmaRecords, ...
            numel(lecpNames), lecpValues, numel(crsNames), crsValues, ...
            true, 0, outputFile, status, note}; %#ok<AGROW>
    end
end

names = {'Spacecraft', 'YearUTC', 'DataStartUTC', 'DataEndUTC', ...
    'SourceProduct', 'SourceFileCount', 'SourceFiles', 'EpochRecords', ...
    'MagneticRecords', 'PlasmaRecords', ...
    'LECPChannels', 'LECPValidValues', 'CRSChannels', 'CRSValidValues', ...
    'OriginalInstrumentValuesOnly', 'GeneratedOrInterpolatedValues', ...
    'FigureFile', 'Status', 'Notes'};
report = cell2table(rows, 'VariableNames', names);
writetable(report, fullfile(opts.OutputFolder, ...
    'UniformLowPrecision_Yearly_Raw_manifest.csv'));
save(fullfile(opts.OutputFolder, 'UniformLowPrecision_Yearly_Raw_report.mat'), ...
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

function product = readNasaV2CohoMag(sourceFile)
% Parse the numeric listing returned by NASA COHOWeb nx1.cgi.
text = fileread(sourceFile);
number = '([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)';
pattern = ['(?m)^\s*(\d{4})\s+(\d{1,3})\s+(\d{1,2})\s+', ...
    number, '\s+', number, '\s+', number, '\s+', number, '\s+', number];
tokens = regexp(text, pattern, 'tokens');
if isempty(tokens)
    error('VoyagerLow:NasaParse', ...
        'No NASA COHOWeb numeric records found in %s.', sourceFile);
end
matrix = nan(numel(tokens), 8);
for ii = 1:numel(tokens)
    matrix(ii, :) = cellfun(@str2double, tokens{ii});
end
epoch = datetime(matrix(:, 1), ones(size(matrix, 1), 1), ...
    ones(size(matrix, 1), 1), 'TimeZone', 'UTC') + ...
    days(matrix(:, 2) - 1) + hours(matrix(:, 3));
values = matrix(:, 4:8);
values(abs(values) >= 900) = nan;
product = struct;
product.Epoch = epoch;
product.heliocentricDistance = values(:, 1);
product.ABS_B = values(:, 2);
product.BR = values(:, 3);
product.BT = values(:, 4);
product.BN = values(:, 5);
product.variable_meta = struct;
product.source_file = sourceFile;
product.available = true;
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
    output.source_files = {input.source_file};
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
            error('VoyagerLow:ColumnMismatch', ...
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
output.source_files{end + 1, 1} = input.source_file;
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

function exportRawYearFigure(spacecraft, yearStart, yearEnd, data, ...
        outputFile, opts, sourceProduct)
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [30 30 1800 1350]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile(layout); hold(ax1, 'on');
field = firstExistingField(data, {'ABS_B', 'F'});
if isempty(field)
    emptyPanel(ax1, 'No recorded |B| value');
else
    plotRecorded(ax1, data.Epoch, data.(field), [0.05 0.05 0.05], ...
        0.7, false, opts.GapBreakHours);
end
ylabel(ax1, 'Magnetic field |B| (nT)', 'FontSize', 11);
grid(ax1, 'on');

ax2 = nexttile(layout); hold(ax2, 'on');
colors = lines(3);
hasComponents = false;
componentFields = {'BR', 'BT', 'BN'};
componentLabels = {'B_R', 'B_T', 'B_N'};
labels = cell(0, 1);
for ii = 1:3
    if isfield(data, componentFields{ii}) && ...
            any(isfinite(data.(componentFields{ii})), 'all')
        plotRecorded(ax2, data.Epoch, ...
            data.(componentFields{ii}), colors(ii, :), 0.6, false, ...
            opts.GapBreakHours);
        labels{end + 1, 1} = componentLabels{ii}; %#ok<AGROW>
        hasComponents = true;
    end
end
if hasComponents
    legend(ax2, labels, 'Location', 'eastoutside', 'FontSize', 10);
else
    emptyPanel(ax2, 'No recorded RTN magnetic components');
end
ylabel(ax2, 'RTN field B_R, B_T, B_N (nT)', 'FontSize', 11);
grid(ax2, 'on');

ax3 = nexttile(layout); hold(ax3, 'on');
lecpNames = productVariableNames(data, 'LECP');
plotParticleChannels(ax3, data, lecpNames, 'LECP', opts.GapBreakHours);
ylabel(ax3, {'LECP proton flux', 'cm^{-2} s^{-1} sr^{-1} MeV^{-1}'}, ...
    'FontSize', 11);
grid(ax3, 'on');

ax4 = nexttile(layout); hold(ax4, 'on');
crsNames = productVariableNames(data, 'CRS');
plotParticleChannels(ax4, data, crsNames, 'CRS', opts.GapBreakHours);
ylabel(ax4, {'CRS proton flux', 'cm^{-2} s^{-1} sr^{-1} MeV^{-1}'}, ...
    'FontSize', 11);
grid(ax4, 'on');

ax5 = nexttile(layout); hold(ax5, 'on');
if isfield(data, 'heliocentricDistance') && plotRecorded(ax5, data.Epoch, ...
        data.heliocentricDistance, [0.35 0.15 0.65], 0.8, false, ...
        opts.GapBreakHours)
else
    emptyPanel(ax5, 'No recorded heliocentric distance');
end
ylabel(ax5, 'Heliocentric distance (AU)', 'FontSize', 11);
grid(ax5, 'on'); xlabel(ax5, 'Month (UTC)', 'FontSize', 11);

axesList = [ax1 ax2 ax3 ax4 ax5];
left = datenum(yearStart); right = datenum(yearEnd); %#ok<DATNM>
monthTicks = yearStart:calmonths(1):yearEnd;
monthTickValues = datenum(monthTicks); %#ok<DATNM>
monthTickLabels = cellstr(datestr(monthTicks, 'mmm')); %#ok<DATST>
linkaxes(axesList, 'x');
for ii = 1:numel(axesList)
    xlim(axesList(ii), [left right]);
    set(axesList(ii), 'XTick', monthTickValues, ...
        'XTickLabel', monthTickLabels);
end
sgtitle(layout, sprintf('%04d', year(yearStart)), ...
    'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
exportgraphics(fig, outputFile, 'Resolution', opts.ExportDPI);
clear cleanup
end

function success = plotRecorded(ax, time, value, color, width, positiveOnly, gapBreakHours)
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
plot(ax, xPlot, yPlot, '-', 'Color', color, ...
    'LineWidth', width, 'Marker', '.', 'MarkerSize', 4);
success = true;
end

function plotParticleChannels(ax, data, names, instrument, gapBreakHours)
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
        colors(ii, :), 0.55, true, gapBreakHours) || hasData;
end
if hasData
    set(ax, 'YScale', 'log');
    legend(ax, labels, 'Location', 'eastoutside', 'Interpreter', 'none', ...
        'FontSize', 8);
else
    emptyPanel(ax, sprintf('No recorded %s proton flux', instrument));
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
