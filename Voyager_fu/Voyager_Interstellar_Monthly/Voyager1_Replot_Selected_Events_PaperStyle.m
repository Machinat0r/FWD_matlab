function report = Voyager1_Replot_Selected_Events_PaperStyle(varargin)
%Voyager1_Replot_Selected_Events_PaperStyle Replot V1 selected events.
%   Each figure covers one calendar month before the event start through
%   one calendar month after the event end. COHO one-hour MAG observations
%   are black and daily averages of those one-hour values are gray.
%   Days with no observation remain gaps. Direct COHO speed, density, and
%   temperature variables are plotted only when valid values exist. Pitch
%   angles are derived only when requested; no missing value is interpolated.

programRoot = fileparts(mfilename('fullpath'));
parser = inputParser;
parser.CaseSensitive = false;
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', @isTextScalar);
addParameter(parser, 'OutputFolder', fullfile(programRoot, ...
    'Voyager1_Selected_Events_1h'), @isTextScalar);
addParameter(parser, 'CatalogFile', '', @isTextScalar);
addParameter(parser, 'EventIDs', {}, ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(parser, 'ReportFolder', '', @isTextScalar);
addParameter(parser, 'ReportTag', '', @isTextScalar);
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'PythonExe', '', @isTextScalar);
addParameter(parser, 'ContextMonths', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && fix(x) == x);
addParameter(parser, 'ContextDays', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0 && fix(x) == x));
addParameter(parser, 'Overwrite', true, @isLogicalScalar);
addParameter(parser, 'Visible', true, @isLogicalScalar);
addParameter(parser, 'ShowEventBoundaries', true, @isLogicalScalar);
addParameter(parser, 'ExportDPI', 200, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 72);
addParameter(parser, 'MAGGapHours', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
addParameter(parser, 'FluxGapHours', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 1);
addParameter(parser, 'LECPDailyAverage', false, @isLogicalScalar);
addParameter(parser, 'CRSDisplay', 'spectrogram', @isTextScalar);
addParameter(parser, 'LECPSectoredDailyCDF', fullfile(programRoot, ...
    'Voyager1_LECP_Sectored_Daily', ...
    'voyager1_lecp_hydrogen_sectored_daily_20130301_20220401.cdf'), ...
    @isTextScalar);
% Zero keeps the cadence of the supplied sectored CDF.  A positive odd
% value requests a centered running mean of daily records.
addParameter(parser, 'LECPSectorAverageDays', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 0 && fix(x) == x && (x == 0 || mod(x, 2) == 1));
addParameter(parser, 'PitchAngleMethod', 'florinski2008', @isTextScalar);
% Public CDAWeb LECP Level-2 sector fluxes carry no event-specific
% penetrating-particle background correction.  "none" preserves those
% published values.  "s8_daily_approx" applies the blocked-sector daily
% value as an explicitly labelled approximation; Florinski et al. (2008)
% state that the background was removed, but do not publish that correction
% series or an automatic algorithm for deriving it.
addParameter(parser, 'LECPBackgroundMode', 'none', @isTextScalar);
addParameter(parser, 'MinimumDailyCoverageFraction', 0.50, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
addParameter(parser, 'MinimumDailySectorSamples', 6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && fix(x) == x);
addParameter(parser, 'MinimumPitchAngleSpanDeg', 45, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 180);
addParameter(parser, 'PitchAngleDataFolder', '', @isTextScalar);
addParameter(parser, 'LineYMarginFraction', 0.06, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(parser, 'CRSColorLimits', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && x(1) < x(2)));
parse(parser, varargin{:});
opts = parser.Results;

opts.DataRoot = char(opts.DataRoot);
opts.OutputFolder = char(opts.OutputFolder);
opts.CatalogFile = char(opts.CatalogFile);
opts.ReportFolder = char(opts.ReportFolder);
opts.ReportTag = char(opts.ReportTag);
opts.CacheRoot = char(opts.CacheRoot);
opts.PythonExe = char(opts.PythonExe);
opts.CRSDisplay = lower(char(opts.CRSDisplay));
opts.LECPSectoredDailyCDF = char(opts.LECPSectoredDailyCDF);
opts.PitchAngleMethod = lower(char(opts.PitchAngleMethod));
opts.LECPBackgroundMode = lower(char(opts.LECPBackgroundMode));
opts.PitchAngleDataFolder = char(opts.PitchAngleDataFolder);
if ~ismember(opts.CRSDisplay, ...
        {'spectrogram', 'channel_panels', 'lecp_p1_sectors', ...
        'lecp_p1_pitch_angle'})
    error('VoyagerPaper:CRSDisplay', ...
        ['CRSDisplay must be spectrogram, channel_panels, ', ...
        'lecp_p1_sectors, or lecp_p1_pitch_angle.']);
end
if ~ismember(opts.PitchAngleMethod, ...
        {'florinski2008', 'fullvector_native'})
    error('VoyagerPaper:PitchAngleMethod', ...
        ['PitchAngleMethod must be florinski2008 or ', ...
        'fullvector_native.']);
end
if ~ismember(opts.LECPBackgroundMode, {'none', 's8_daily_approx'})
    error('VoyagerPaper:LECPBackgroundMode', ...
        ['LECPBackgroundMode must be none or s8_daily_approx. ', ...
        'The latter is an approximation, not the unpublished LECP-team ', ...
        'background correction used by Florinski et al. (2008).']);
end
if strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle')
    if strcmp(opts.PitchAngleMethod, 'florinski2008') && ...
            opts.LECPSectorAverageDays ~= 1
        error('VoyagerPaper:FlorinskiCadence', ...
            ['Florinski et al. (2008) processing requires one UTC-day ', ...
            'averages: set LECPSectorAverageDays to 1.']);
    elseif strcmp(opts.PitchAngleMethod, 'fullvector_native') && ...
            opts.LECPSectorAverageDays ~= 0
        error('VoyagerPaper:NativeCadence', ...
            ['fullvector_native keeps the supplied cadence: set ', ...
            'LECPSectorAverageDays to 0.']);
    end
end
if isempty(opts.CatalogFile)
    opts.CatalogFile = fullfile(opts.OutputFolder, ...
        'Voyager1_Selected_Events_1h_catalog.csv');
end
if isempty(opts.ReportFolder)
    opts.ReportFolder = opts.OutputFolder;
end
if isempty(opts.PitchAngleDataFolder)
    opts.PitchAngleDataFolder = fullfile(programRoot, ...
        'Voyager1_Selected_Event_Data', 'derived', ...
        'lecp_pitch_angle_daily');
end
if ~isfolder(opts.DataRoot)
    error('VoyagerPaper:DataRootMissing', ...
        'Voyager data root is unavailable: %s', opts.DataRoot);
end
if ~isfile(opts.CatalogFile)
    error('VoyagerPaper:CatalogMissing', ...
        'Event catalog is unavailable: %s', opts.CatalogFile);
end
ensureFolder(opts.OutputFolder);
ensureFolder(opts.ReportFolder);
if strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle')
    ensureFolder(opts.PitchAngleDataFolder);
end

catalog = readEventCatalog(opts.CatalogFile);
if ~isempty(opts.EventIDs)
    requestedIDs = string(opts.EventIDs(:));
    catalog = catalog(ismember(string(catalog.EventID), requestedIDs), :);
    if isempty(catalog)
        error('VoyagerPaper:EventSelection', ...
            'None of the requested EventIDs occur in the event catalog.');
    end
end
if isempty(opts.ContextDays)
    catalog.PlotStartUTC = catalog.StartUTC - calmonths(opts.ContextMonths);
    catalog.PlotEndUTCExclusive = ...
        catalog.EndUTCExclusive + calmonths(opts.ContextMonths);
else
    catalog.PlotStartUTC = catalog.StartUTC - days(opts.ContextDays);
    catalog.PlotEndUTCExclusive = ...
        catalog.EndUTCExclusive + days(opts.ContextDays);
end

lecpSectored = emptyProduct();
lecpSectorLoadNote = "";
if ismember(opts.CRSDisplay, ...
        {'lecp_p1_sectors', 'lecp_p1_pitch_angle'})
    try
        if ~isfile(opts.LECPSectoredDailyCDF)
            error('VoyagerPaper:LECPSectorSourceMissing', ...
                'LECP sectored CDF is unavailable: %s', ...
                opts.LECPSectoredDailyCDF);
        end
        lecpSectored = Voyager_Read_CDF_Product( ...
            opts.LECPSectoredDailyCDF, 'lecp_sector_daily', ...
            'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
    catch ME
        lecpSectorLoadNote = appendNote(lecpSectorLoadNote, ME);
    end
end

cohoFiles = listMonthlyCoho(opts.DataRoot);
memoryCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
cohoData = cell(height(catalog), 1);
cohoSources = cell(height(catalog), 1);
loadNotes = strings(height(catalog), 1);
allCRSLogFlux = zeros(0, 1);

for ii = 1:height(catalog)
    fprintf('[Paper-style event] %s  plot %s to %s UTC\n', ...
        catalog.EventID{ii}, char(catalog.PlotStartUTC(ii)), ...
        char(catalog.PlotEndUTCExclusive(ii)));
    try
        [cohoData{ii}, cohoSources{ii}] = loadPeriodProduct( ...
            cohoFiles, 'MonthUTC', catalog.PlotStartUTC(ii), ...
            catalog.PlotEndUTCExclusive(ii), 'coho', opts, memoryCache);
        cohoData{ii} = normalizeCohoMagProduct(cohoData{ii});
    catch ME
        cohoData{ii} = emptyProduct();
        cohoSources{ii} = cell(0, 1);
        loadNotes(ii) = appendNote(loadNotes(ii), ME);
    end
    crsNames = productVariableNames(cohoData{ii}, 'CRS');
    for jj = 1:numel(crsNames)
        value = cohoData{ii}.(crsNames{jj});
        value = value(isfinite(value) & value > 0);
        allCRSLogFlux = [allCRSLogFlux; log10(value(:))]; %#ok<AGROW>
    end
end

if isempty(opts.CRSColorLimits)
    colorLimits = robustColorLimits(allCRSLogFlux);
else
    colorLimits = double(opts.CRSColorLimits(:).');
end

figureFiles = cell(height(catalog), 1);
pitchAngleFiles = repmat({''}, height(catalog), 1);
status = repmat({'ok'}, height(catalog), 1);
notes = cellstr(loadNotes);
mag1hRecords = zeros(height(catalog), 1);
magDailyDays = zeros(height(catalog), 1);
cohoRecords = zeros(height(catalog), 1);
directPlasmaVariables = cell(height(catalog), 1);

for ii = 1:height(catalog)
    productTag = 'COHO1h_raw';
    if logical(opts.LECPDailyAverage)
        productTag = 'MAG1h_LECPdaily';
    end
    if strcmp(opts.CRSDisplay, 'channel_panels')
        productTag = [productTag, '_CRSchannels'];
    elseif strcmp(opts.CRSDisplay, 'lecp_p1_sectors')
        if opts.LECPSectorAverageDays == 0
            productTag = [productTag, '_LECP_P1_sectors_1h'];
        else
            productTag = [productTag, sprintf('_LECP_P1_sectors_%dd', ...
                opts.LECPSectorAverageDays)];
        end
    elseif strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle')
        if strcmp(opts.PitchAngleMethod, 'florinski2008')
            productTag = [productTag, ...
                '_LECP_P1_pitch_angle_Florinski2008_1d'];
        else
            productTag = [productTag, ...
                '_LECP_P1_pitch_angle_fullVector_native'];
        end
    end
    fileName = sprintf('V1_%s_%s_%s_%s.png', ...
        catalog.EventID{ii}, datestr(catalog.StartUTC(ii), 'yyyymmdd'), ...
        datestr(catalog.EndUTCInclusive(ii), 'yyyymmdd'), productTag); %#ok<DATST>
    figureFiles{ii} = fullfile(opts.OutputFolder, fileName);
    if strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle')
        if strcmp(opts.PitchAngleMethod, 'florinski2008')
            pitchAngleFiles{ii} = fullfile(opts.PitchAngleDataFolder, ...
                sprintf(['V1_%s_%s_%s_LECP_P1_Florinski2008_', ...
                'pitch_angle_1d.csv'], catalog.EventID{ii}, ...
                datestr(catalog.StartUTC(ii), 'yyyymmdd'), ...
                datestr(catalog.EndUTCInclusive(ii), 'yyyymmdd'))); %#ok<DATST>
        else
            pitchAngleFiles{ii} = fullfile(opts.PitchAngleDataFolder, ...
                sprintf(['V1_%s_%s_%s_LECP_P1_fullVector_', ...
                'pitch_angle_native.csv'], catalog.EventID{ii}, ...
                datestr(catalog.StartUTC(ii), 'yyyymmdd'), ...
                datestr(catalog.EndUTCInclusive(ii), 'yyyymmdd'))); %#ok<DATST>
        end
    end
    mag1hRecords(ii) = finiteMagRecordCount(cohoData{ii});
    cohoRecords(ii) = variableRecordCount(cohoData{ii}, 'Epoch');
    daily = dailyMagneticMeans(cohoData{ii}, catalog.PlotStartUTC(ii), ...
        catalog.PlotEndUTCExclusive(ii));
    magDailyDays(ii) = nnz(isfinite(daily.F1));
    directPlasmaVariables{ii} = strjoin( ...
        availableDirectPlasmaFields(cohoData{ii}), ';');
    if strlength(loadNotes(ii)) > 0
        status{ii} = 'partial_error';
    end
    if ismember(opts.CRSDisplay, ...
            {'lecp_p1_sectors', 'lecp_p1_pitch_angle'}) && ...
            strlength(lecpSectorLoadNote) > 0
        status{ii} = 'partial_error';
        notes{ii} = char(appendNote(string(notes{ii}), ...
            MException('VoyagerPaper:LECPSectorLoad', ...
            '%s', lecpSectorLoadNote)));
    end
    try
        if logical(opts.Overwrite) || ~isfile(figureFiles{ii})
            pitchAngleTable = exportPaperStyleFigure( ...
                catalog(ii, :), cohoData{ii}, ...
                cohoData{ii}, lecpSectored, figureFiles{ii}, ...
                colorLimits, opts);
            if strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle') && ...
                    ~isempty(pitchAngleTable)
                writetable(pitchAngleTable, pitchAngleFiles{ii});
            end
        else
            status{ii} = 'existing';
        end
    catch ME
        status{ii} = 'plot_error';
        notes{ii} = char(appendNote(string(notes{ii}), ME));
        warning('VoyagerPaper:PlotFailed', '%s: %s', ...
            catalog.EventID{ii}, notes{ii});
    end
end

report = table(catalog.EventID, catalog.StartUTC, ...
    catalog.EndUTCInclusive, catalog.PlotStartUTC, ...
    catalog.PlotEndUTCExclusive, cellfun(@numel, cohoSources), ...
    mag1hRecords, magDailyDays, cellfun(@numel, cohoSources), ...
    cohoRecords, directPlasmaVariables, figureFiles, pitchAngleFiles, ...
    status, notes, ...
    'VariableNames', {'EventID', 'EventStartUTC', ...
    'EventEndUTCInclusive', 'PlotStartUTC', 'PlotEndUTCExclusive', ...
    'MAG1hSourceFileCount', 'MAG1hRecords', 'MAGDailyMeanDays', ...
    'COHOSourceFileCount', 'COHORecords', 'DirectPlasmaVariables', ...
    'FigureFile', 'PitchAngleDataFile', 'Status', 'Notes'});
reportSuffix = '';
if ~isempty(opts.ReportTag)
    reportSuffix = ['_', regexprep(opts.ReportTag, '[^A-Za-z0-9_-]', '_')];
end
writetable(report, fullfile(opts.ReportFolder, ...
    ['Voyager1_Selected_Events_PaperStyle_manifest', reportSuffix, '.csv']));
save(fullfile(opts.ReportFolder, ...
    ['Voyager1_Selected_Events_PaperStyle_report', reportSuffix, '.mat']), ...
    'report', 'catalog', 'colorLimits', 'opts', '-v7.3');
fprintf('Replotted %d events in %s\n', height(report), opts.OutputFolder);
end

function catalog = readEventCatalog(catalogFile)
% R2026a can crash inside readtable on some Windows installations with
% third-party virtual display drivers.  This catalog is a simple quoted CSV,
% so read it with textscan and construct the table explicitly.
fid = fopen(catalogFile, 'rt');
if fid < 0
    error('VoyagerPaper:CatalogOpen', ...
        'Unable to open event catalog: %s', catalogFile);
end
cleanup = onCleanup(@() fclose(fid));
headerLine = fgetl(fid);
if ~ischar(headerLine)
    error('VoyagerPaper:CatalogEmpty', ...
        'Event catalog is empty: %s', catalogFile);
end
variableNames = cellstr(split(string(headerLine), ','));
formatSpec = repmat('%q', 1, numel(variableNames));
columns = textscan(fid, formatSpec, 'Delimiter', ',', ...
    'EndOfLine', '\n', 'ReturnOnError', false, ...
    'Whitespace', '');
clear cleanup
catalog = table;
for ii = 1:numel(variableNames)
    field = matlab.lang.makeValidName(strtrim(variableNames{ii}));
    catalog.(field) = string(columns{ii});
end
required = {'EventID', 'StartUTC', 'EndUTCInclusive', 'EndUTCExclusive'};
if ~all(ismember(required, catalog.Properties.VariableNames))
    error('VoyagerPaper:CatalogColumns', ...
        'Catalog is missing one or more required columns.');
end
catalog.EventID = cellstr(catalog.EventID);
catalog.StartUTC = parseUtcColumn(catalog.StartUTC);
catalog.EndUTCInclusive = parseUtcColumn(catalog.EndUTCInclusive);
catalog.EndUTCExclusive = parseUtcColumn(catalog.EndUTCExclusive);
end

function value = parseUtcColumn(input)
if isdatetime(input)
    value = input;
    value.TimeZone = 'UTC';
    return
end
value = datetime(string(input), 'TimeZone', 'UTC');
end

function files = listMonthlyCoho(dataRoot)
root = fullfile(dataRoot, 'voyager1', 'coho', '1hr', 'l2', ...
    'merged_mag_plasma');
files = latestPeriodFiles(root, '*', '*', ...
    '_(\d{8})_v(\d+)\.cdf$', 'month');
files.Properties.VariableNames{1} = 'MonthUTC';
end

function product = normalizeCohoMagProduct(product)
% Use one consistent magnitude field name for the magnetic plotting code.
if isfield(product, 'F1') && any(isfinite(product.F1), 'all')
    return
end
if isfield(product, 'ABS_B') && any(isfinite(product.ABS_B), 'all')
    product.F1 = product.ABS_B;
elseif isfield(product, 'F') && any(isfinite(product.F), 'all')
    product.F1 = product.F;
end
end

function files = latestPeriodFiles(root, level1, level2, pattern, periodType)
if ~isfolder(root)
    error('VoyagerPaper:ProductRootMissing', ...
        'Product folder is unavailable: %s', root);
end
if isempty(level2)
    candidates = dir(fullfile(root, level1, '*.cdf'));
else
    candidates = dir(fullfile(root, level1, level2, '*.cdf'));
end
periods = NaT(0, 1, 'TimeZone', 'UTC');
paths = cell(0, 1);
versions = zeros(0, 1);
for ii = 1:numel(candidates)
    token = regexp(candidates(ii).name, pattern, 'tokens', 'once');
    if isempty(token), continue, end
    if strcmp(periodType, 'month')
        period = datetime(token{1}(1:6), 'InputFormat', 'yyyyMM', ...
            'TimeZone', 'UTC');
    else
        period = datetime(token{1}(1:4), 'InputFormat', 'yyyy', ...
            'TimeZone', 'UTC');
    end
    version = str2double(token{2});
    fullPath = fullfile(candidates(ii).folder, candidates(ii).name);
    existing = find(periods == period, 1);
    if isempty(existing)
        periods(end + 1, 1) = period; %#ok<AGROW>
        paths{end + 1, 1} = fullPath; %#ok<AGROW>
        versions(end + 1, 1) = version; %#ok<AGROW>
    elseif version >= versions(existing)
        paths{existing} = fullPath;
        versions(existing) = version;
    end
end
[periods, order] = sort(periods);
files = table(periods, paths(order), ...
    'VariableNames', {'PeriodUTC', 'SourceCDF'});
end

function [product, sources] = loadPeriodProduct(files, periodName, ...
        startTime, endTime, profile, opts, memoryCache)
period = files.(periodName);
if strcmp(periodName, 'MonthUTC')
    firstPeriod = dateshift(startTime, 'start', 'month');
    lastPeriod = dateshift(endTime - seconds(1), 'start', 'month');
else
    firstPeriod = dateshift(startTime, 'start', 'year');
    lastPeriod = dateshift(endTime - seconds(1), 'start', 'year');
end
use = period >= firstPeriod & period <= lastPeriod;
sources = files.SourceCDF(use);
if isempty(sources)
    error('VoyagerPaper:NoSource', ...
        'No %s source file overlaps %s to %s.', ...
        profile, char(startTime), char(endTime));
end
product = struct;
for ii = 1:numel(sources)
    key = [profile, '|', sources{ii}];
    if isKey(memoryCache, key)
        current = memoryCache(key);
    else
        current = Voyager_Read_CDF_Product(sources{ii}, profile, ...
            'CacheRoot', opts.CacheRoot, 'PythonExe', opts.PythonExe);
        memoryCache(key) = current;
    end
    current = subsetByTime(current, startTime, endTime);
    product = appendProduct(product, current);
end
product = sortAndDeduplicate(product);
if ~isfield(product, 'Epoch') || isempty(product.Epoch)
    error('VoyagerPaper:NoRecords', ...
        'No %s record falls in the requested plot interval.', profile);
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

function pitchAngleTable = exportPaperStyleFigure( ...
        eventRow, mag, coho, lecpSectored, ...
        outputFile, colorLimits, opts)
pitchAngleTable = table;
plasmaFields = availableDirectPlasmaFields(coho);
crsNames = productVariableNames(coho, 'CRS');
smallPanels = 4 + numel(plasmaFields);
if strcmp(opts.CRSDisplay, 'channel_panels')
    crsRows = max(1, numel(crsNames));
    totalRows = smallPanels + 2 + crsRows;
    figureHeight = max(1600, 170 * totalRows);
else
    totalRows = smallPanels + 4;
    figureHeight = 1600 + 170 * numel(plasmaFields);
end
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [30 30 1800 figureHeight]);
cleanup = [];
if ~logical(opts.Visible)
    cleanup = onCleanup(@() close(fig));
end
layout = tiledlayout(fig, totalRows, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
axesList = gobjects(0, 1);
axSector = gobjects(0, 1);
tile = 1;
panelIndex = 0;
dailyMag = dailyMagneticMeans(mag, eventRow.PlotStartUTC, ...
    eventRow.PlotEndUTCExclusive);

magFields = {'F1', 'BR', 'BT', 'BN'};
magLabels = {{'|B|', '(nT)'}, {'B_R', '(nT)'}, ...
    {'B_T', '(nT)'}, {'B_N', '(nT)'}};
for ii = 1:numel(magFields)
    ax = nexttile(layout, tile); tile = tile + 1;
    axesList(end + 1, 1) = ax; %#ok<AGROW>
    hold(ax, 'on');
    if isfield(mag, 'Epoch') && isfield(mag, magFields{ii}) && ...
            any(isfinite(mag.(magFields{ii})))
        rawHandle = plotGapLine(ax, mag.Epoch, mag.(magFields{ii}), ...
            [0.03 0.03 0.03], 0.62, opts.MAGGapHours, false);
        dailyHandle = plot(ax, datenum(dailyMag.Epoch), ... %#ok<DATNM>
            dailyMag.(magFields{ii}), '-', 'Color', [0.78 0.78 0.78], ...
            'LineWidth', 1.45);
        if ii == 1
            legend(ax, [rawHandle dailyHandle], ...
                {'COHO MAG 1 h', 'daily mean from 1 h'}, ...
                'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
        end
        rawValues = mag.(magFields{ii});
        dailyValues = dailyMag.(magFields{ii});
        applyPaddedLinearYLim(ax, [rawValues(:); dailyValues(:)], ...
            opts.LineYMarginFraction);
    else
        emptyPanel(ax, sprintf('No recorded %s', magFields{ii}));
    end
    ylabel(ax, magLabels{ii}, 'FontSize', 11);
    panelIndex = panelIndex + 1; panelLabel(ax, panelIndex);
end

for ii = 1:numel(plasmaFields)
    ax = nexttile(layout, tile); tile = tile + 1;
    axesList(end + 1, 1) = ax; %#ok<AGROW>
    hold(ax, 'on');
    field = plasmaFields{ii};
    rawHandle = plotGapLine(ax, coho.Epoch, coho.(field), ...
        [0.72 0.72 0.72], 0.55, opts.FluxGapHours, false);
    daily = dailyMeanSeries(coho.Epoch, coho.(field), ...
        eventRow.PlotStartUTC, eventRow.PlotEndUTCExclusive);
    dailyHandle = plot(ax, datenum(daily.Epoch), daily.Value, '-', ... %#ok<DATNM>
        'Color', [0.03 0.03 0.03], 'LineWidth', 1.35);
    if ii == 1
        legend(ax, [rawHandle dailyHandle], {'COHO 1 h', 'daily mean'}, ...
            'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
    end
    ylabel(ax, plasmaLabel(field), 'FontSize', 11);
    panelIndex = panelIndex + 1; panelLabel(ax, panelIndex);
end

axLECP = nexttile(layout, tile, [2 1]); tile = tile + 2;
axesList(end + 1, 1) = axLECP; hold(axLECP, 'on'); %#ok<AGROW>
plotLECPChannels(axLECP, coho, ...
    productVariableNames(coho, 'LECP'), opts.FluxGapHours, ...
    eventRow.PlotStartUTC, eventRow.PlotEndUTCExclusive, ...
    logical(opts.LECPDailyAverage));
ylabel(axLECP, {'LECP proton J', ...
    '(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, 'FontSize', 11);
panelIndex = panelIndex + 1; panelLabel(axLECP, panelIndex);

if strcmp(opts.CRSDisplay, 'channel_panels')
    if isempty(crsNames)
        axCRS = nexttile(layout, tile); tile = tile + 1;
        axesList(end + 1, 1) = axCRS; hold(axCRS, 'on'); %#ok<AGROW>
        emptyPanel(axCRS, 'No recorded CRS proton flux');
        ylabel(axCRS, {'CRS proton J', '(no channel)'}, 'FontSize', 9);
        panelIndex = panelIndex + 1; panelLabel(axCRS, panelIndex);
    else
        channelColors = turbo(numel(crsNames));
        for jj = 1:numel(crsNames)
            axCRS = nexttile(layout, tile); tile = tile + 1;
            axesList(end + 1, 1) = axCRS; hold(axCRS, 'on'); %#ok<AGROW>
            plotCRSChannelPointLine(axCRS, coho, crsNames{jj}, ...
                opts.FluxGapHours, channelColors(jj, :), ...
                opts.LineYMarginFraction);
            ylabel(axCRS, {conciseEnergyLabel(coho, crsNames{jj}), ...
                'J (cm^{-2} s^{-1}', 'sr^{-1} MeV^{-1})'}, ...
                'FontSize', 8);
            panelIndex = panelIndex + 1; panelLabel(axCRS, panelIndex);
        end
    end
    xlabel(axesList(end), 'UTC', 'FontSize', 11);
elseif strcmp(opts.CRSDisplay, 'lecp_p1_sectors')
    axSector = nexttile(layout, tile, [2 1]);
    axesList(end + 1, 1) = axSector; hold(axSector, 'on'); %#ok<AGROW>
    plotLECPP1SectorPairs(axSector, lecpSectored, ...
        eventRow.PlotStartUTC, eventRow.PlotEndUTCExclusive, ...
        opts.LECPSectorAverageDays);
    ylabel(axSector, {'LECP H 0.57-1.78 MeV', ...
        'J (cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, 'FontSize', 10);
    xlabel(axSector, 'UTC', 'FontSize', 11);
    panelIndex = panelIndex + 1; panelLabel(axSector, panelIndex);
elseif strcmp(opts.CRSDisplay, 'lecp_p1_pitch_angle')
    axPitch = nexttile(layout, tile, [2 1]);
    axesList(end + 1, 1) = axPitch; hold(axPitch, 'on'); %#ok<AGROW>
    pitchAngleTable = plotLECPP1PitchAngleFlux(axPitch, ...
        lecpSectored, mag, eventRow.PlotStartUTC, ...
        eventRow.PlotEndUTCExclusive, opts);
    set(axPitch, 'YAxisLocation', 'left');
    ylabel(axPitch, 'PA (deg)', 'FontSize', 11);
    xlabel(axPitch, 'UTC', 'FontSize', 11);
    panelIndex = panelIndex + 1; panelLabel(axPitch, panelIndex);
else
    axCRS = nexttile(layout, tile, [2 1]);
    axesList(end + 1, 1) = axCRS; hold(axCRS, 'on'); %#ok<AGROW>
    plotCRSSpectrogram(axCRS, coho, eventRow.PlotStartUTC, ...
        eventRow.PlotEndUTCExclusive, colorLimits);
    ylabel(axCRS, {'CRS proton', 'energy (MeV)'}, 'FontSize', 11);
    xlabel(axCRS, 'UTC', 'FontSize', 11);
    panelIndex = panelIndex + 1; panelLabel(axCRS, panelIndex);
end

left = datenum(eventRow.PlotStartUTC); %#ok<DATNM>
right = datenum(eventRow.PlotEndUTCExclusive); %#ok<DATNM>
tickTimes = eventTimeTicks(eventRow.PlotStartUTC, ...
    eventRow.PlotEndUTCExclusive);
tickValues = datenum(tickTimes); %#ok<DATNM>
tickLabels = cellstr(datestr(tickTimes, 'dd-mmm-yyyy')); %#ok<DATST>
linkaxes(axesList, 'x');
for ii = 1:numel(axesList)
    ax = axesList(ii);
    xlim(ax, [left right]);
    set(ax, 'XTick', tickValues, 'XTickLabel', tickLabels, ...
        'FontSize', 9, 'LineWidth', 0.8, 'TickDir', 'out', ...
        'Box', 'on', 'Layer', 'top');
    grid(ax, 'off');
    if logical(opts.ShowEventBoundaries)
        xline(ax, datenum(eventRow.StartUTC), '--', ... %#ok<DATNM>
            'Color', [0.15 0.35 0.75], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        xline(ax, datenum(eventRow.EndUTCExclusive), '--', ... %#ok<DATNM>
            'Color', [0.15 0.35 0.75], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end
end
if numel(axesList) > 1
    set(axesList(1:end-1), 'XTickLabel', []);
end
if isempty(opts.ContextDays)
    contextText = sprintf('%d-month context on each side', ...
        opts.ContextMonths);
else
    contextText = sprintf('%d-day context on each side', opts.ContextDays);
end
sgtitle(layout, sprintf(['Voyager 1  %s   event %s--%s UTC  ', ...
    '(%s)'], eventRow.EventID{1}, ...
    datestr(eventRow.StartUTC, 'yyyy-mm-dd'), ...
    datestr(eventRow.EndUTCInclusive, 'yyyy-mm-dd'), contextText), ...
    'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none'); %#ok<DATST>
if isgraphics(axSector)
    addLECPSectorDisk(axSector, mag, eventRow.StartUTC, ...
        eventRow.EndUTCExclusive);
end
if logical(opts.Visible)
    figure(fig);
    drawnow;
    enableVoyagerDataCursor(fig);
end
exportgraphics(fig, outputFile, 'Resolution', opts.ExportDPI);
clear cleanup
end

function enableVoyagerDataCursor(fig)
% Keep the exported figure open and enable point-by-point inspection.
try
    cursorMode = datacursormode(fig);
    cursorMode.Enable = 'on';
    cursorMode.SnapToDataVertex = 'on';
    cursorMode.UpdateFcn = @voyagerDataTipText;
catch ME
    warning('VoyagerPaper:DataCursor', ...
        'Could not enable data-cursor mode: %s', ME.message);
end
end

function textOutput = voyagerDataTipText(~, eventObject)
% Format datenum x coordinates as UTC and expose heat-map color values.
position = eventObject.Position;
try
    timeText = datestr(position(1), 'yyyy-mm-dd HH:MM'); %#ok<DATST>
catch
    timeText = sprintf('%.8g', position(1));
end
textOutput = {sprintf('UTC: %s', timeText), ...
    sprintf('Y: %.8g', position(2))};

target = eventObject.Target;
if isgraphics(target)
    targetType = lower(get(target, 'Type'));
    if ismember(targetType, {'surface', 'image'})
        colorData = get(target, 'CData');
        colorData = colorData(isfinite(colorData));
        if ~isempty(colorData)
            logFlux = colorData(1);
            textOutput{end + 1} = sprintf('log10 J: %.6g', logFlux); %#ok<AGROW>
            textOutput{end + 1} = sprintf('J: %.6g', 10.^logFlux); %#ok<AGROW>
        end
    end
end
end

function daily = dailyMagneticMeans(mag, startTime, endTime)
daily = struct;
daysGrid = (dateshift(startTime, 'start', 'day'):caldays(1): ...
    dateshift(endTime - seconds(1), 'start', 'day')).';
daily.Epoch = daysGrid + hours(12);
fields = {'F1', 'BR', 'BT', 'BN'};
for ii = 1:numel(fields)
    if isfield(mag, 'Epoch') && isfield(mag, fields{ii})
        series = dailyMeanSeries(mag.Epoch, mag.(fields{ii}), ...
            startTime, endTime);
        daily.(fields{ii}) = series.Value;
    else
        daily.(fields{ii}) = nan(numel(daysGrid), 1);
    end
end
end

function hourly = hourlyMagneticMeans(mag, startTime, endTime)
hourGrid = (dateshift(startTime, 'start', 'hour'):hours(1): ...
    dateshift(endTime - seconds(1), 'start', 'hour')).';
hourly = struct;
hourly.Epoch = hourGrid + minutes(30);
fields = {'F1', 'BR', 'BT', 'BN'};
for ii = 1:numel(fields)
    if isfield(mag, 'Epoch') && isfield(mag, fields{ii})
        values = nan(numel(hourGrid), 1);
        sourceHour = dateshift(mag.Epoch(:), 'start', 'hour');
        sourceValue = mag.(fields{ii});
        sourceValue = sourceValue(:);
        for jj = 1:numel(hourGrid)
            use = sourceHour == hourGrid(jj) & isfinite(sourceValue);
            if any(use)
                values(jj) = mean(sourceValue(use), 'omitnan');
            end
        end
        hourly.(fields{ii}) = values;
    else
        hourly.(fields{ii}) = nan(numel(hourGrid), 1);
    end
end
end

function series = dailyMeanSeries(time, value, startTime, endTime)
dayGrid = (dateshift(startTime, 'start', 'day'):caldays(1): ...
    dateshift(endTime - seconds(1), 'start', 'day')).';
means = nan(numel(dayGrid), 1);
time = time(:); value = value(:);
for ii = 1:numel(dayGrid)
    use = time >= dayGrid(ii) & time < dayGrid(ii) + caldays(1) & ...
        isfinite(value);
    if any(use)
        means(ii) = mean(value(use), 'omitnan');
    end
end
series = struct('Epoch', dayGrid + hours(12), 'Value', means);
end

function handle = plotGapLine(ax, time, value, color, width, ...
        gapBreakHours, positiveOnly)
handle = gobjects(1);
if isempty(time) || isempty(value)
    return
end
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
xPlot(positions) = x; yPlot(positions) = value;
handle = plot(ax, xPlot, yPlot, '-', 'Color', color, 'LineWidth', width);
end

function fields = availableDirectPlasmaFields(coho)
fields = {};
candidates = {'V', 'protonDensity', 'protonTemp'};
for ii = 1:numel(candidates)
    field = candidates{ii};
    if isfield(coho, field) && any(isfinite(coho.(field)))
        fields{end + 1, 1} = field; %#ok<AGROW>
    end
end
end

function label = plasmaLabel(field)
switch field
    case 'V'
        label = {'Proton speed', '(km s^{-1})'};
    case 'protonDensity'
        label = {'Proton density', '(cm^{-3})'};
    case 'protonTemp'
        label = {'Proton temperature', '(K)'};
    otherwise
        label = field;
end
end

function plotLECPChannels(ax, data, names, gapBreakHours, startTime, ...
        endTime, dailyAverage)
if isempty(names) || ~isfield(data, 'Epoch')
    emptyPanel(ax, 'No recorded LECP proton flux');
    return
end
colors = [0.2706 0.0000 0.3294; 0.1276 0.5669 0.5506; ...
    0.9932 0.9062 0.1439];
handles = gobjects(0, 1); labels = cell(0, 1);
for ii = 1:numel(names)
    value = data.(names{ii});
    if ~any(isfinite(value) & value > 0), continue, end
    time = data.Epoch;
    lineWidth = 0.9;
    label = conciseEnergyLabel(data, names{ii});
    if dailyAverage
        daily = dailyMeanSeries(time, value, startTime, endTime);
        time = daily.Epoch;
        value = daily.Value;
        gapBreakHours = 36;
        lineWidth = 1.25;
        label = sprintf('%s daily mean', label);
    end
    handles(end + 1, 1) = plotGapLine(ax, time, value, ... %#ok<AGROW>
        colors(min(ii, size(colors, 1)), :), lineWidth, ...
        gapBreakHours, true);
    labels{end + 1, 1} = label; %#ok<AGROW>
end

if isempty(handles)
    emptyPanel(ax, 'No recorded LECP proton flux');
else
    set(ax, 'YScale', 'log');
    legend(ax, handles, labels, 'Location', 'northeast', 'Box', 'off', ...
        'FontSize', 8, 'Interpreter', 'none');
end
end

function plotLECPP1SectorPairs(ax, data, startTime, endTime, averageDays)
required = {'Epoch', 'FHDU_SectoredFluxes', 'FHDU_Energy'};
if ~all(isfield(data, required)) || isempty(data.Epoch)
    emptyPanel(ax, 'No recorded LECP P1 sectored flux');
    return
end

energy = data.FHDU_Energy;
if size(energy, 1) == numel(data.Epoch)
    energy = energy(find(any(isfinite(energy), 2), 1), :);
end
energy = energy(:);
[~, p1Index] = min(abs(energy - 0.73));
flux = data.FHDU_SectoredFluxes;
if ndims(flux) ~= 3 || size(flux, 1) ~= numel(data.Epoch) || ...
        size(flux, 2) < p1Index || size(flux, 3) < 7
    emptyPanel(ax, 'LECP P1 sectored array has an unexpected shape');
    return
end

use = data.Epoch >= startTime & data.Epoch < endTime;
time = data.Epoch(use);
values = squeeze(flux(use, p1Index, :));
if isempty(time) || isempty(values)
    emptyPanel(ax, 'No LECP P1 sectored flux in this interval');
    return
end
values(~isfinite(values) | values <= 0) = nan;

sectorPairs = {[1 5], [3 7], [2 6]};
pairFlux = nan(numel(time), numel(sectorPairs));
for ii = 1:numel(sectorPairs)
    current = values(:, sectorPairs{ii});
    valid = all(isfinite(current), 2);
    pairFlux(valid, ii) = mean(current(valid, :), 2);
end
    if averageDays > 0
        [time, pairFlux] = centeredDailyAverageNoInterpolation( ...
            time, pairFlux, averageDays);
        gapBreakHours = 36;
    else
        [time, order] = sort(time(:));
        pairFlux = pairFlux(order, :);
        [time, uniqueAt] = unique(time, 'stable');
        pairFlux = pairFlux(uniqueAt, :);
        gapBreakHours = 1.5;
    end

colors = [0.00 0.28 0.92; 0.90 0.08 0.08; 0.08 0.52 0.10];
labels = {'Sectors 1 & 5', 'Sectors 3 & 7', 'Sectors 2 & 6'};
handles = gobjects(0, 1);
for ii = 1:size(pairFlux, 2)
    handle = plotGapLine(ax, time, pairFlux(:, ii), colors(ii, :), ...
        1.15, gapBreakHours, true);
    if isgraphics(handle)
        set(handle, 'Marker', '.', 'MarkerSize', 5.5, ...
            'MarkerEdgeColor', colors(ii, :));
        handles(end + 1, 1) = handle; %#ok<AGROW>
    end
end
if isempty(handles)
    emptyPanel(ax, 'No valid LECP P1 sector-pair averages');
    return
end
set(ax, 'YScale', 'log');
validValues = pairFlux(isfinite(pairFlux) & pairFlux > 0);
if ~isempty(validValues)
    applyPaddedLogYLim(ax, validValues, 0.08);
end
legend(ax, handles, labels(1:numel(handles)), ...
    'Location', 'northeast', 'Box', 'on', 'FontSize', 9, ...
    'Interpreter', 'none');
if averageDays > 0
    processingText = sprintf('%d-day average; no interpolation', averageDays);
else
    processingText = 'native source cadence; no added averaging/interpolation';
end
text(ax, 0.01, 0.08, processingText, 'Units', 'normalized', ...
    'FontSize', 8, 'Color', [0.25 0.25 0.25], 'Interpreter', 'none');
end

function output = plotLECPP1PitchAngleFlux( ...
        ax, data, mag, startTime, endTime, opts)
% Construct LECP P1 pitch angles following Florinski et al. (2008).
% The seven active sectors are averaged separately over a UTC day and their
% entry vectors are dotted with the complete daily mean MAG vector.  The
% nominal LECP scan plane is taken to be parallel to RT, as in the paper.
% Background handling is explicit because the paper's LECP-team correction
% series is not included in the public CDAWeb products.  No interpolation,
% fit, or time-dependent SEDR attitude correction is used.
output = table;
required = {'Epoch', 'FHDU_Energy'};
if ~all(isfield(data, required)) || isempty(data.Epoch)
    emptyPanel(ax, 'No recorded LECP P1 sectored measurement');
    return
end

if isfield(data, 'FHDU_SectoredFluxes')
    sectoredValues = data.FHDU_SectoredFluxes;
    sourceValueName = 'Flux';
    valueName = 'Flux';
    sourceToIntensityFactor = 1;
    uncertaintyField = 'FHDU_SectoredFluxUncertainties';
    colorBarLabel = 'log_{10} J';
elseif isfield(data, 'FHDU_SectoredRates')
    sectoredValues = data.FHDU_SectoredRates;
    sourceValueName = 'Rate';
    valueName = 'Flux';
    % Florinski et al. give G=0.44 cm^2 sr for Ch1.  The channel width is
    % 1.78-0.57=1.21 MeV, so J=(rate-background)/(G*dE).
    sourceToIntensityFactor = 1 / (0.44 * (1.78 - 0.57));
    uncertaintyField = 'FHDU_SectoredRateUncertainties';
    colorBarLabel = 'log_{10} J';
else
    emptyPanel(ax, 'No recorded LECP P1 sectored measurement');
    return
end

energy = data.FHDU_Energy;
if size(energy, 1) == numel(data.Epoch)
    firstValid = find(any(isfinite(energy), 2), 1);
    if isempty(firstValid)
        emptyPanel(ax, 'No valid LECP hydrogen energy table');
        return
    end
    energy = energy(firstValid, :);
end
energy = energy(:);
[~, p1Index] = min(abs(energy - 0.73));
if ndims(sectoredValues) ~= 3 || ...
        size(sectoredValues, 1) ~= numel(data.Epoch) || ...
        size(sectoredValues, 2) < p1Index || ...
        size(sectoredValues, 3) < 8
    emptyPanel(ax, 'LECP P1 sectored array has an unexpected shape');
    return
end

use = data.Epoch >= startTime & data.Epoch < endTime;
sourceTime = data.Epoch(use);
rawSector = squeeze(sectoredValues(use, p1Index, 1:8));
if isempty(sourceTime) || isempty(rawSector)
    emptyPanel(ax, 'No LECP P1 sectored measurement in this interval');
    return
end
rawSector = reshape(double(rawSector), numel(sourceTime), 8);
rawSector(~isfinite(rawSector) | rawSector < 0) = nan;

sourceSigma = nan(size(rawSector));
if isfield(data, uncertaintyField)
    allSigma = data.(uncertaintyField);
    if ndims(allSigma) == 3 && size(allSigma, 1) == numel(data.Epoch) && ...
            size(allSigma, 2) >= p1Index && size(allSigma, 3) >= 8
        sourceSigma = squeeze(allSigma(use, p1Index, 1:8));
        sourceSigma = reshape(double(sourceSigma), numel(sourceTime), 8);
        sourceSigma(~isfinite(sourceSigma) | sourceSigma < 0) = nan;
    end
end

[sourceTime, order] = sort(sourceTime(:));
rawSector = rawSector(order, :);
sourceSigma = sourceSigma(order, :);
[sourceTime, uniqueAt] = unique(sourceTime, 'stable');
rawSector = rawSector(uniqueAt, :);
sourceSigma = sourceSigma(uniqueAt, :);

if strcmp(opts.PitchAngleMethod, 'florinski2008')
    [time, rawSector, rawSigma, minimumSectorSamples, ...
        coverageFraction, requiredSamples, sectorSampleCount] = ...
        florinskiDailySectorAverage( ...
        sourceTime, rawSector, sourceSigma, startTime, endTime, ...
        opts.MinimumDailyCoverageFraction, opts.MinimumDailySectorSamples);
    matchedMag = magneticVectorMeans(mag, startTime, endTime, 'day');
    lecpBin = dateshift(time, 'start', 'day');
    magBin = dateshift(matchedMag.Epoch, 'start', 'day');
    magVariableNames = {'BR_daily_nT', 'BT_daily_nT', 'BN_daily_nT'};
    valueCadenceSuffix = '1d';
    cellHalfWidthDays = 0.5;
    methodText = ["Florinski et al. (2008) angle construction: " + ...
        "per-sector 24 h mean; mu=vhat dot <B>/|<B>|; " + ...
        "nominal scan plane parallel to RT; no SEDR correction"];
else
    time = sourceTime;
    rawSigma = sourceSigma;
    complete = all(isfinite(rawSector(:, 1:7)), 2);
    sectorSampleCount = double(isfinite(rawSector));
    minimumSectorSamples = min(sectorSampleCount(:, 1:7), [], 2);
    coverageFraction = double(complete);
    requiredSamples = ones(numel(time), 1);
    matchedMag = magneticVectorMeans(mag, startTime, endTime, 'hour');
    lecpBin = dateshift(time, 'start', 'hour');
    magBin = dateshift(matchedMag.Epoch, 'start', 'hour');
    magVariableNames = {'BR_hourly_nT', 'BT_hourly_nT', 'BN_hourly_nT'};
    valueCadenceSuffix = 'native';
    if numel(time) > 1
        cellHalfWidthDays = max(1 / (24 * 60), ...
            median(days(diff(time)), 'omitnan') / 2);
    else
        cellHalfWidthDays = 1 / 48;
    end
    methodText = ["native-cadence diagnostic using Florinski geometry; " + ...
        "mu=vhat dot B/|B|; " + ...
        "nominal scan plane parallel to RT; no SEDR correction"];
end

% S8 is blocked by the sun shield and is retained as a background diagnostic.
% Instrument documentation says that the background must be determined from
% context with help from S8.  It does not support blindly equating every S8
% daily value with the correction.  Keep the public Level-2 values by default;
% the optional S8 mode is explicitly labelled as an approximation.
sectorValue = rawSector;
sectorSigma = rawSigma;
background = rawSector(:, 8);
backgroundSigma = rawSigma(:, 8);
backgroundApplied = strcmp(opts.LECPBackgroundMode, 's8_daily_approx');
if backgroundApplied
    for sector = 1:7
        sectorValue(:, sector) = rawSector(:, sector) - background;
        hasSigma = isfinite(rawSigma(:, sector)) & ...
            isfinite(backgroundSigma);
        sectorSigma(hasSigma, sector) = sqrt( ...
            rawSigma(hasSigma, sector).^2 + backgroundSigma(hasSigma).^2);
    end
    methodText = methodText + ...
        "; approximate background: same-day blocked S8 value subtracted";
else
    methodText = methodText + ...
        "; public source values retained without background correction";
end
sectorValue(~isfinite(sectorValue) | sectorValue <= 0) = nan;
sectorValue(:, 1:7) = sectorValue(:, 1:7) * sourceToIntensityFactor;
sectorSigma(:, 1:7) = sectorSigma(:, 1:7) * sourceToIntensityFactor;
sectorValue(:, 8) = background * sourceToIntensityFactor;
sectorSigma(:, 8) = backgroundSigma * sourceToIntensityFactor;
if sourceToIntensityFactor ~= 1
    methodText = methodText + ...
        "; Ch1 rate converted with G=0.44 cm^2 sr and dE=1.21 MeV";
end

[present, magIndex] = ismember(lecpBin, magBin);
br = nan(numel(time), 1); bt = br; bn = br;
magSamples = zeros(numel(time), 1);
magDirectionRMSDeg = nan(numel(time), 1);
br(present) = matchedMag.BR(magIndex(present));
bt(present) = matchedMag.BT(magIndex(present));
bn(present) = matchedMag.BN(magIndex(present));
magSamples(present) = matchedMag.SampleCount(magIndex(present));
magDirectionRMSDeg(present) = ...
    matchedMag.DirectionRMSDeg(magIndex(present));
bMagnitude = sqrt(br.^2 + bt.^2 + bn.^2);
bAzimuthDeg = mod(atan2d(bt, br), 360);
bElevationDeg = atan2d(bn, hypot(br, bt));

% Published disk orientation: +T is right and +R is down.  The disk vectors
% below are telescope look/boresight directions.  A detected particle enters
% along the opposite vector.  Its nominal N component is zero because the
% paper treats the LECP scan plane as nearly parallel to RT.
sectorPlotAngle = 112.5 + (0:7) * 45;
boresightUT = cosd(sectorPlotAngle);
boresightUR = -sind(sectorPlotAngle);
particleUT = -boresightUT;
particleUR = -boresightUR;
particleUN = zeros(size(particleUR));

pitchCosine = nan(numel(time), 8);
pitchAngle = nan(numel(time), 8);
for sector = 1:8
    pitchCosine(:, sector) = ...
        (br * particleUR(sector) + bt * particleUT(sector) + ...
        bn * particleUN(sector)) ./ bMagnitude;
    pitchCosine(:, sector) = max(-1, min(1, ...
        pitchCosine(:, sector)));
    pitchAngle(:, sector) = acosd(pitchCosine(:, sector));
end
invalidB = ~isfinite(bMagnitude) | bMagnitude <= 0;
pitchCosine(invalidB, :) = nan;
pitchAngle(invalidB, :) = nan;

activeSectors = 1:7;
pitchSpanDeg = nan(numel(time), 1);
activeSectorCount = zeros(numel(time), 1);
relativeUncertaintyMedian = nan(numel(time), 1);
padUsable = false(numel(time), 1);
quality = strings(numel(time), 1);
for ii = 1:numel(time)
    validActive = isfinite(sectorValue(ii, activeSectors)) & ...
        sectorValue(ii, activeSectors) > 0 & ...
        isfinite(pitchAngle(ii, activeSectors));
    activeSectorCount(ii) = nnz(validActive);
    currentPitch = pitchAngle(ii, activeSectors(validActive));
    if ~isempty(currentPitch)
        pitchSpanDeg(ii) = max(currentPitch) - min(currentPitch);
    end
    rel = sectorSigma(ii, activeSectors) ./ ...
        sectorValue(ii, activeSectors);
    rel = rel(isfinite(rel) & rel >= 0);
    if ~isempty(rel)
        relativeUncertaintyMedian(ii) = median(rel);
    end

    reasons = strings(0, 1);
    if minimumSectorSamples(ii) < requiredSamples(ii)
        reasons(end + 1) = sprintf( ...
            'insufficient per-sector samples (%d<%d)', ...
            minimumSectorSamples(ii), requiredSamples(ii)); %#ok<AGROW>
    end
    if magSamples(ii) < 1 || invalidB(ii)
        reasons(end + 1) = 'no valid mean vector magnetic field'; %#ok<AGROW>
    end
    if activeSectorCount(ii) < 7
        reasons(end + 1) = sprintf( ...
            'only %d/7 active sectors valid after selected background handling', ...
            activeSectorCount(ii)); %#ok<AGROW>
    end
    if isfinite(pitchSpanDeg(ii)) && ...
            pitchSpanDeg(ii) <= opts.MinimumPitchAngleSpanDeg
        reasons(end + 1) = sprintf( ...
            'pitch-angle span %.1f deg <= %.1f deg', ...
            pitchSpanDeg(ii), opts.MinimumPitchAngleSpanDeg); %#ok<AGROW>
    end
    if isfinite(magDirectionRMSDeg(ii)) && magDirectionRMSDeg(ii) > 20
        reasons(end + 1) = sprintf( ...
            'within-bin MAG direction RMS %.1f deg > 20 deg', ...
            magDirectionRMSDeg(ii)); %#ok<AGROW>
    end
    if isfinite(relativeUncertaintyMedian(ii)) && ...
            relativeUncertaintyMedian(ii) > 1
        reasons(end + 1) = sprintf( ...
            'median relative flux uncertainty %.2f > 1', ...
            relativeUncertaintyMedian(ii)); %#ok<AGROW>
    end
    if isempty(reasons)
        padUsable(ii) = true;
        quality(ii) = 'usable';
    else
        quality(ii) = "excluded: " + strjoin(reasons, '; ');
    end
end

output = table(time, br, bt, bn, bMagnitude, bAzimuthDeg, ...
    bElevationDeg, magSamples, magDirectionRMSDeg, minimumSectorSamples, ...
    requiredSamples, coverageFraction, activeSectorCount, pitchSpanDeg, ...
    relativeUncertaintyMedian, padUsable, quality, ...
    'VariableNames', [{'EpochUTC'}, magVariableNames, ...
    {'BMeanMagnitude_nT', 'BAzimuth_R_to_T_deg', 'BElevation_deg', ...
    'MAGVectorSampleCount', 'MAGDirectionRMS_deg', ...
    'LECPMinimumSectorSampleCount', 'LECPRequiredSampleCount', ...
    'LECPCoverageFraction', 'ValidActiveSectorCount', ...
    'PitchAngleSpan_deg', 'MedianRelativeUncertainty', ...
    'PADUsable', 'PADQuality'}]);
for sector = 1:8
    output.(sprintf('Mu_S%d', sector)) = pitchCosine(:, sector);
    output.(sprintf('PA_S%d_deg', sector)) = pitchAngle(:, sector);
    output.(sprintf('Raw%s_S%d_%s', sourceValueName, sector, ...
        valueCadenceSuffix)) = rawSector(:, sector);
    output.(sprintf('%s_S%d_%s', valueName, sector, ...
        valueCadenceSuffix)) = sectorValue(:, sector);
    output.(sprintf('%sUncertainty_S%d_%s', valueName, sector, ...
        valueCadenceSuffix)) = sectorSigma(:, sector);
    output.(sprintf('Samples_S%d_%s', sector, valueCadenceSuffix)) = ...
        sectorSampleCount(:, sector);
end
output.(sprintf('S8_Diagnostic_%s_%s', sourceValueName, ...
    valueCadenceSuffix)) = background;
output.SourceToDifferentialFluxFactor = repmat( ...
    sourceToIntensityFactor, height(output), 1);
output.BackgroundCorrectionMode = repmat( ...
    string(opts.LECPBackgroundMode), height(output), 1);
output.BackgroundCorrectionApplied = repmat( ...
    backgroundApplied, height(output), 1);
output.PitchAngleMethod = repmat(methodText, height(output), 1);
output.SectorAxisAssumption = repmat( ...
    "eight 45-deg nominal look axes; scan plane parallel to RT", ...
    height(output), 1);

displayFlux = sectorValue(:, activeSectors);
displayPitch = pitchAngle(:, activeSectors);
usableRows = padUsable;
validFlux = displayFlux(repmat(usableRows, 1, 7) & ...
    isfinite(displayFlux) & displayFlux > 0 & isfinite(displayPitch));
if isempty(validFlux)
    emptyPanel(ax, ['No daily PAD passed the Florinski-style ', ...
        'coverage/geometry checks']);
    return
end
colorLimits = robustColorLimits(log10(validFlux));

for ii = 1:numel(time)
    if ~padUsable(ii), continue, end
    valid = isfinite(displayPitch(ii, :)) & ...
        isfinite(displayFlux(ii, :)) & displayFlux(ii, :) > 0;
    currentPitch = displayPitch(ii, valid);
    currentFlux = displayFlux(ii, valid);
    [currentPitch, order] = sort(currentPitch);
    currentFlux = currentFlux(order);

    % A spectrogram cannot show two sector values at exactly the same PA.
    % Only for display, sector centers within 2 deg are represented by their
    % arithmetic mean.  Every original sector value remains in the CSV.
    group = cumsum([1, diff(currentPitch) > 2]);
    groupIDs = unique(group, 'stable');
    mergedPitch = nan(1, numel(groupIDs));
    mergedFlux = nan(1, numel(groupIDs));
    for jj = 1:numel(groupIDs)
        inGroup = group == groupIDs(jj);
        mergedPitch(jj) = mean(currentPitch(inGroup), 'omitnan');
        mergedFlux(jj) = mean(currentFlux(inGroup), 'omitnan');
    end

    if numel(mergedPitch) == 1
        pitchEdges = [max(0, mergedPitch - 22.5), ...
            min(180, mergedPitch + 22.5)];
    else
        midpoints = (mergedPitch(1:end-1) + mergedPitch(2:end)) / 2;
        pitchEdges = [max(0, mergedPitch(1) - ...
            (midpoints(1) - mergedPitch(1))), midpoints, ...
            min(180, mergedPitch(end) + ...
            (mergedPitch(end) - midpoints(end)))];
    end
    xCenter = datenum(time(ii)); %#ok<DATNM>
    for jj = 1:numel(mergedPitch)
        xData = [xCenter - cellHalfWidthDays, ...
            xCenter + cellHalfWidthDays; ...
            xCenter - cellHalfWidthDays, ...
            xCenter + cellHalfWidthDays];
        yData = [pitchEdges(jj), pitchEdges(jj); ...
            pitchEdges(jj + 1), pitchEdges(jj + 1)];
        colorData = log10(mergedFlux(jj)) * ones(2);
        surface(ax, xData, yData, zeros(2), colorData, ...
            'FaceColor', 'flat', 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

view(ax, 2);
set(ax, 'YDir', 'normal', 'YLim', [0 180], ...
    'YTick', [0 45 90 135 180]);
colormap(ax, turbo(256));
caxis(ax, colorLimits);
colorBar = colorbar(ax, 'Location', 'eastoutside');
colorBar.Label.String = colorBarLabel;
colorBar.Label.Interpreter = 'tex';
colorBar.FontSize = 8;
text(ax, 0.99, 0.94, '0.57-1.78 MeV', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'Color', [0 0 0], ...
    'FontSize', 9, 'Interpreter', 'none');
if strcmp(opts.PitchAngleMethod, 'florinski2008')
    processingLabel = '24 h sector means; full vector B';
else
    processingLabel = 'native cadence; full vector B';
end
if backgroundApplied
    processingLabel = [processingLabel, '; S8 daily background approximation'];
else
    processingLabel = [processingLabel, '; public flux, background uncorrected'];
end
text(ax, 0.01, 0.05, processingLabel, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'Color', [0.15 0.15 0.15], 'FontSize', 8, 'Interpreter', 'none');
end

function [outputTime, outputValues, outputSigma, minimumSectorSamples, ...
        coverageFraction, requiredSamples, sectorSampleCount] = ...
        florinskiDailySectorAverage(time, values, sigma, startTime, ...
        endTime, minimumCoverageFraction, minimumDailySectorSamples)
% Average each active motor sector over the same UTC day.  A Level-1 CDF
% record can contain only a subset of sectors, so coverage is assessed from
% S1-S7 only.  S8 is a blocked diagnostic direction and is not required for
% a PAD unless the caller explicitly requests the S8 approximation.
dayGrid = (dateshift(startTime, 'start', 'day'):caldays(1): ...
    dateshift(endTime - seconds(1), 'start', 'day')).';
outputTime = dayGrid + hours(12);
outputValues = nan(numel(dayGrid), 8);
outputSigma = nan(numel(dayGrid), 8);
minimumSectorSamples = zeros(numel(dayGrid), 1);
coverageFraction = zeros(numel(dayGrid), 1);
sectorSampleCount = zeros(numel(dayGrid), 8);

uniqueTime = unique(sort(time(:)));
if numel(uniqueTime) > 1
    cadenceHours = median(hours(diff(uniqueTime)), 'omitnan');
else
    cadenceHours = 24;
end
if ~isfinite(cadenceHours) || cadenceHours <= 0
    cadenceHours = 24;
end
sourceIsDailyAverage = cadenceHours >= 12;
requiredSamples = ones(numel(dayGrid), 1);

for ii = 1:numel(dayGrid)
    inDay = time >= dayGrid(ii) & time < dayGrid(ii) + caldays(1);
    for sector = 1:8
        sectorSampleCount(ii, sector) = ...
            nnz(inDay & isfinite(values(:, sector)));
    end
    minimumSectorSamples(ii) = min(sectorSampleCount(ii, 1:7));
    maximumSectorSamples = max(sectorSampleCount(ii, 1:7));
    if maximumSectorSamples > 0
        coverageFraction(ii) = ...
            minimumSectorSamples(ii) / maximumSectorSamples;
    end
    if sourceIsDailyAverage
        requiredSamples(ii) = 1;
    else
        requiredSamples(ii) = max(minimumDailySectorSamples, ...
            ceil(minimumCoverageFraction * maximumSectorSamples));
    end
    if minimumSectorSamples(ii) < requiredSamples(ii)
        continue
    end
    for sector = 1:8
        useSector = inDay & isfinite(values(:, sector));
        outputValues(ii, sector) = ...
            mean(values(useSector, sector), 'omitnan');
        currentSigma = sigma(useSector, sector);
        if all(isfinite(currentSigma))
            outputSigma(ii, sector) = ...
                sqrt(sum(currentSigma.^2)) / sectorSampleCount(ii, sector);
        end
    end
end
end

function output = magneticVectorMeans(mag, startTime, endTime, binType)
% Average BR, BT, and BN over identical records and retain a direction RMS.
if strcmp(binType, 'day')
    binStart = (dateshift(startTime, 'start', 'day'):caldays(1): ...
        dateshift(endTime - seconds(1), 'start', 'day')).';
    binEnd = binStart + caldays(1);
    outputTime = binStart + hours(12);
else
    binStart = (dateshift(startTime, 'start', 'hour'):hours(1): ...
        dateshift(endTime - seconds(1), 'start', 'hour')).';
    binEnd = binStart + hours(1);
    outputTime = binStart + minutes(30);
end

number = numel(binStart);
output = struct('Epoch', outputTime, 'BR', nan(number, 1), ...
    'BT', nan(number, 1), 'BN', nan(number, 1), ...
    'SampleCount', zeros(number, 1), ...
    'DirectionRMSDeg', nan(number, 1));
if ~isfield(mag, 'Epoch') || ~all(isfield(mag, {'BR', 'BT', 'BN'}))
    return
end

br = mag.BR(:); bt = mag.BT(:); bn = mag.BN(:);
for ii = 1:number
    use = mag.Epoch >= binStart(ii) & mag.Epoch < binEnd(ii) & ...
        isfinite(br) & isfinite(bt) & isfinite(bn);
    output.SampleCount(ii) = nnz(use);
    if ~any(use), continue, end
    output.BR(ii) = mean(br(use), 'omitnan');
    output.BT(ii) = mean(bt(use), 'omitnan');
    output.BN(ii) = mean(bn(use), 'omitnan');
    meanVector = [output.BR(ii), output.BT(ii), output.BN(ii)];
    meanMagnitude = norm(meanVector);
    sampleVectors = [br(use), bt(use), bn(use)];
    sampleMagnitude = sqrt(sum(sampleVectors.^2, 2));
    validDirection = sampleMagnitude > 0 & isfinite(sampleMagnitude) & ...
        isfinite(meanMagnitude) & meanMagnitude > 0;
    if any(validDirection)
        directionCosine = sampleVectors(validDirection, :) * ...
            meanVector.' ./ (sampleMagnitude(validDirection) * ...
            meanMagnitude);
        directionCosine = max(-1, min(1, directionCosine));
        directionDifference = acosd(directionCosine);
        output.DirectionRMSDeg(ii) = ...
            sqrt(mean(directionDifference.^2, 'omitnan'));
    end
end
end

function addLECPSectorDisk(parentAx, mag, eventStart, eventEnd)
% Draw the LECP eight-sector viewing disk in the same T-R orientation used
% by the reference paper.  The orange arrow is the event-mean COHO 1 h
% magnetic-field direction projected into the RT plane.
drawnow;
fig = ancestor(parentAx, 'figure');
parentPosition = getpixelposition(parentAx, true);
diameter = min(180, max(140, 0.70 * parentPosition(4)));
insetPosition = [parentPosition(1) + 34, ...
    parentPosition(2) + parentPosition(4) - diameter - 12, ...
    diameter, diameter];
inset = axes(fig, 'Units', 'pixels', 'Position', insetPosition, ...
    'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none', ...
    'HitTest', 'off', 'PickableParts', 'none');
hold(inset, 'on');
axis(inset, 'equal');
xlim(inset, [-1.48 1.48]);
ylim(inset, [-1.48 1.48]);
set(inset, 'Clipping', 'off', 'Visible', 'off');

structureColor = [0.12 0.12 0.12];
angle = linspace(0, 2 * pi, 361);
plot(inset, cos(angle), sin(angle), '-', 'Color', structureColor, ...
    'LineWidth', 1.15, 'HandleVisibility', 'off');
for boundary = 0:45:315
    x = cosd(boundary); y = sind(boundary);
    plot(inset, [-x x], [-y y], '-', 'Color', structureColor, ...
        'LineWidth', 0.85, 'HandleVisibility', 'off');
end

blue = [0.00 0.28 0.92];
red = [0.90 0.08 0.08];
green = [0.08 0.52 0.10];
neutral = [0.32 0.32 0.32];
sectorColors = [blue; green; red; neutral; blue; green; red; neutral];
sectorAngles = 112.5 + (0:7) * 45;
for sector = 1:8
    text(inset, 0.64 * cosd(sectorAngles(sector)), ...
        0.64 * sind(sectorAngles(sector)), sprintf('%d', sector), ...
        'Color', sectorColors(sector, :), 'FontSize', 10, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none', ...
        'Clipping', 'off');
end

quiver(inset, 0, 0, 1.22, 0, 0, 'Color', structureColor, ...
    'LineWidth', 1.0, 'MaxHeadSize', 0.22, 'HandleVisibility', 'off');
text(inset, 1.35, 0.02, 'T', 'Color', structureColor, ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
quiver(inset, 0, 0, 0, -1.22, 0, 'Color', structureColor, ...
    'LineWidth', 1.0, 'MaxHeadSize', 0.22, 'HandleVisibility', 'off');
text(inset, 0, -1.38, 'R', 'Color', structureColor, ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

if isfield(mag, 'Epoch') && isfield(mag, 'BR') && isfield(mag, 'BT')
    valid = mag.Epoch >= eventStart & mag.Epoch < eventEnd & ...
        isfinite(mag.BR) & isfinite(mag.BT);
    if any(valid)
        meanBR = mean(mag.BR(valid));
        meanBT = mean(mag.BT(valid));
        b = [meanBT, -meanBR];
        magnitude = hypot(b(1), b(2));
        if isfinite(magnitude) && magnitude > 0
            b = b / magnitude;
            orange = [0.95 0.43 0.05];
            quiver(inset, -1.12 * b(1), -1.12 * b(2), ...
                2.24 * b(1), 2.24 * b(2), 0, 'Color', orange, ...
                'LineWidth', 1.45, 'MaxHeadSize', 0.19, ...
                'HandleVisibility', 'off');
            text(inset, 1.34 * b(1), 1.34 * b(2), 'B', ...
                'Color', orange, 'FontSize', 11, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Interpreter', 'none', ...
                'Clipping', 'off');
        end
    end
end
end

function [outputTime, outputValues] = centeredDailyAverageNoInterpolation( ...
        time, values, windowDays)
[time, order] = sort(time(:));
values = values(order, :);
sampleDay = dateshift(time, 'start', 'day');
[day, ~, dayGroup] = unique(sampleDay, 'stable');
dailyValues = nan(numel(day), size(values, 2));
for ii = 1:numel(day)
    inDay = dayGroup == ii;
    for jj = 1:size(values, 2)
        current = values(inDay, jj);
        current = current(isfinite(current));
        if ~isempty(current)
            dailyValues(ii, jj) = mean(current);
        end
    end
end
values = dailyValues;
outputTime = day + hours(12);
outputValues = nan(size(values));
halfWindow = (windowDays - 1) / 2;
for ii = (1 + halfWindow):(numel(day) - halfWindow)
    index = (ii - halfWindow):(ii + halfWindow);
    if any(diff(day(index)) ~= days(1))
        continue
    end
    current = values(index, :);
    validColumns = all(isfinite(current), 1);
    outputValues(ii, validColumns) = ...
        mean(current(:, validColumns), 1);
end
end

function applyPaddedLinearYLim(ax, values, marginFraction)
values = values(isfinite(values));
if isempty(values), return, end
low = min(values); high = max(values);
span = high - low;
if span <= 0
    span = max(abs([low high]));
    if span <= 0, span = 1; end
end
padding = marginFraction * span;
ylim(ax, [low - padding, high + padding]);
end

function plotCRSChannelPointLine(ax, data, name, gapBreakHours, color, ...
        marginFraction)
if ~isfield(data, 'Epoch') || ~isfield(data, name)
    emptyPanel(ax, sprintf('No recorded %s', name));
    return
end
values = data.(name);
valid = isfinite(values) & values > 0;
if ~any(valid)
    emptyPanel(ax, sprintf('No positive records in %s', name));
    return
end
handle = plotGapLine(ax, data.Epoch, values, color, 0.8, ...
    gapBreakHours, true);
if isgraphics(handle)
    set(handle, 'Marker', 'o', 'MarkerSize', 2.6, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', color);
end
set(ax, 'YScale', 'log');
applyPaddedLogYLim(ax, values(valid), marginFraction);
end

function applyPaddedLogYLim(ax, values, marginFraction)
values = values(isfinite(values) & values > 0);
if isempty(values), return, end
logValues = log10(values);
low = min(logValues); high = max(logValues);
span = high - low;
if span <= 0, span = 1; end
padding = marginFraction * span;
ylim(ax, 10 .^ [low - padding, high + padding]);
end

function plotCRSSpectrogram(ax, data, startTime, endTime, colorLimits)
names = productVariableNames(data, 'CRS');
if isempty(names) || ~isfield(data, 'Epoch') || isempty(data.Epoch)
    emptyPanel(ax, 'No recorded CRS proton flux');
    return
end
[energyLow, energyHigh] = energyBounds(data, names);
valid = isfinite(energyLow) & isfinite(energyHigh) & ...
    energyLow > 0 & energyHigh > energyLow;
names = names(valid); energyLow = energyLow(valid); energyHigh = energyHigh(valid);
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
    surface(ax, xData, yData, zeros(size(xData)), [cRow; cRow], ...
        'FaceColor', 'flat', 'EdgeColor', 'none', 'CDataMapping', 'scaled');
end
set(ax, 'YScale', 'log', 'YDir', 'normal');
ylim(ax, [min(energyLow) max(energyHigh)]);
yticks(ax, [3 10 30 100 300]);
colormap(ax, turbo(256));
caxis(ax, colorLimits);
cb = colorbar(ax, 'eastoutside');
cb.FontSize = 10;
cb.Label.String = {'log_{10} J', ...
    'cm^{-2} s^{-1} sr^{-1} MeV^{-1}'};
cb.Label.Interpreter = 'tex'; cb.Label.FontSize = 10;
end

function ticks = eventTimeTicks(startTime, endTime)
spanDays = days(endTime - startTime);
if spanDays <= 50
    step = 7;
elseif spanDays <= 100
    step = 14;
else
    step = 21;
end
ticks = startTime:caldays(step):endTime;
if ticks(end) < endTime
    ticks(end + 1) = endTime; %#ok<AGROW>
end
end

function panelLabel(ax, index, xPosition)
if nargin < 3, xPosition = 0.009; end
label = sprintf('(%c)', char(double('a') + index - 1));
text(ax, xPosition, 0.91, label, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'FontSize', 10, 'Color', [0.08 0.08 0.08]);
end

function names = productVariableNames(product, suffix)
names = {};
if isempty(product), return, end
fields = fieldnames(product); number = [];
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

function [low, high] = energyBounds(product, names)
low = nan(numel(names), 1); high = nan(numel(names), 1);
for ii = 1:numel(names)
    label = variableLabel(product, names{ii});
    token = regexp(label, ...
        'H\s+([0-9.]+)\s*-\s*([0-9.]+)\s*MeV', 'tokens', 'once');
    if ~isempty(token)
        low(ii) = str2double(token{1}); high(ii) = str2double(token{2});
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
if ~isfield(product, 'variable_meta') || ...
        ~isfield(product.variable_meta, name), return, end
item = product.variable_meta.(name);
if ~isfield(item, 'attributes'), return, end
attrs = item.attributes;
for candidate = {'FIELDNAM', 'LABLAXIS', 'CATDESC'}
    if isfield(attrs, candidate{1})
        raw = attrs.(candidate{1});
        if ischar(raw) || (isstring(raw) && isscalar(raw))
            label = char(raw); return
        end
    end
end
end

function count = variableRecordCount(product, field)
count = 0;
if isfield(product, field), count = size(product.(field), 1); end
end

function count = finiteMagRecordCount(product)
recordCount = variableRecordCount(product, 'Epoch');
valid = false(recordCount, 1);
fields = {'F1', 'BR', 'BT', 'BN'};
for ii = 1:numel(fields)
    if isfield(product, fields{ii}) && ...
            size(product.(fields{ii}), 1) == recordCount
        valid = valid | any(isfinite(product.(fields{ii})), 2);
    end
end
count = nnz(valid);
end

function limits = robustColorLimits(values)
values = sort(values(isfinite(values)));
if isempty(values), limits = [-5 -1]; return, end
n = numel(values);
lo = values(max(1, round(0.01 * n)));
hi = values(min(n, max(1, round(0.99 * n))));
lo = floor(lo * 4) / 4; hi = ceil(hi * 4) / 4;
if hi - lo < 1
    middle = (lo + hi) / 2; lo = middle - 0.5; hi = middle + 0.5;
end
limits = [lo hi];
end

function note = appendNote(note, exception)
newText = string(sprintf('%s: %s', exception.identifier, exception.message));
if strlength(note) == 0
    note = newText;
else
    note = note + " | " + newText;
end
end

function product = emptyProduct()
product = struct('Epoch', NaT(0, 1, 'TimeZone', 'UTC'), ...
    'variable_meta', struct);
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
    error('VoyagerPaper:FolderCreateFailed', ...
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
