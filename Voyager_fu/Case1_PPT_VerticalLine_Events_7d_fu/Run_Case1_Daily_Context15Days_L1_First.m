function result = Run_Case1_Daily_Context15Days_L1_First(baselineFile, eventIDs)
%Run_Case1_Daily_Context15Days_L1_First Update all 47 daily event overviews.
%   Window [D-15,D+16) contains 31 complete UTC days for a one-day event.
%   Existing L1-first processing is unchanged. Hourly and standalone PAD
%   figures are protected; only existing daily overviews are overwritten.
%   Optional EventIDs runs a documented partial batch while a source is
%   unavailable. A later full run reuses the original baselineFile.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-04

%% original catalog plus the two established supplementary events
if nargin < 1, baselineFile = ''; end
if nargin < 2, eventIDs = strings(0, 1); end
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
catalog = Case1_Event_Catalog;
catalog = catalog(catalog.Spacecraft == 1, :);
assert(height(catalog) == 45, 'Unexpected original V1 event count.');
supplementaryFiles = strings(2, 1);
dates = {'20200730', '20200822'};
for ii = 1:numel(dates)
    supplementaryFiles(ii) = string(fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
        'validation', ['supplemental_', dates{ii}], ...
        ['Case1_Added_', dates{ii}, '_selection.mat']));
    added = load(supplementaryFiles(ii), 'catalog');
    assert(height(added.catalog) == 1 && added.catalog.Spacecraft == 1 && ...
        string(added.catalog.EventID) == "Case1-Added-" + string(dates{ii}), ...
        'Unexpected supplementary catalog.');
    catalog = [catalog; added.catalog]; %#ok<AGROW>
end
assert(height(catalog) == 47 && numel(unique(string(catalog.EventID))) == 47, ...
    'Expected 47 unique established V1 events.');
plotCatalog = catalog;
if ~isempty(eventIDs)
    assert(all(ismember(string(eventIDs), string(catalog.EventID))), 'Unknown requested EventID.');
    plotCatalog = catalog(ismember(string(catalog.EventID), string(eventIDs)), :);
    assert(~isempty(plotCatalog), 'No events selected.');
end
cfg.ContextDays = 15;
cfg.PADCadence = 'day';
cfg.LECPSourcePriority = 'l1_first';
cfg.ExportPeakPAD = false;

%% explicitly require each annual source; a wholly absent year must not slip by
requiredYears = [];
for ii = 1:height(plotCatalog)
    begin = plotCatalog.StartUTC(ii)-days(15);
    stop = plotCatalog.EndUTCExclusive(ii)+days(15)-seconds(1);
    requiredYears = [requiredYears, year(begin):year(stop)]; %#ok<AGROW>
end
requiredYears = unique(requiredYears);
requireYears(cfg.LECPLevel1CDFs, requiredYears, 'L1');
requireYears(cfg.LECPNativeDailyCDFs, requiredYears, 'L2 daily');
cfg.LECPLevel1CDFs = Case1_Restrict_Source_Years(cfg.LECPLevel1CDFs, plotCatalog, 15);
cfg.LECPNativeDailyCDFs = Case1_Restrict_Source_Years(cfg.LECPNativeDailyCDFs, plotCatalog, 15);
sourceManifest = Case1_Check_Data(cfg, plotCatalog);

%% save the old daily numerical payloads and hashes before any overwrite
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
auditFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', ...
    'context15d_daily_l1_first', stamp);
if ~isfolder(auditFolder), mkdir(auditFolder); end
if strlength(string(baselineFile)) == 0
    baseline = captureBaseline(cfg, catalog);
    baseline.AuditFile = string(fullfile(auditFolder, 'baseline.mat'));
    save(baseline.AuditFile, 'baseline', '-v7.3');
else
    prior = load(baselineFile, 'baseline');
    baseline = prior.baseline;
    assert(isequal(string(baseline.Catalog.EventID), string(catalog.EventID)), ...
        'The preserved baseline has a different event catalog.');
end
fprintf('Saved baseline for %d daily images and %d protected artifacts.\n', ...
    height(baseline.DailyFiles), height(baseline.ProtectedFiles));

%% existing renderer; all current non-window scientific settings are retained
result = struct('CreatedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'Cadence', "day", 'ContextDays', 15, 'SourcePriority', "l1_first", ...
    'Catalog', catalog, 'SupplementaryCatalogFiles', supplementaryFiles, ...
    'PlotCatalog', plotCatalog, 'CompleteBatch', height(plotCatalog) == 47, ...
    'RequiredSourceYears', requiredYears, 'PreflightSourceManifest', sourceManifest, ...
    'BaselineFile', baseline.AuditFile, 'AuditFolder', string(auditFolder));
result.AuditFile = string(fullfile(auditFolder, 'replot_daily_15days.mat'));
result.Run = Run_Case1_PPT_VerticalLine_Events_7d( ...
    'RunPlots', true, 'Overwrite', true, 'Visible', false, ...
    'PADCadence', 'day', 'ExportPeakPAD', false, ...
    'ContextDays', 15, 'LECPSourcePriority', 'l1_first', ...
    'EventCatalog', plotCatalog);
result.CodeFiles = string(fullfile(cfg.CodeRoot, ...
    {'Run_Case1_Daily_Context15Days_L1_First.m', ...
    'Run_Case1_PPT_VerticalLine_Events_7d.m', 'Case1_Config.m', ...
    'Case1_Check_Data.m', 'Voyager_Case1_Plot_Events.m'})).';
result.CodeSHA256 = strings(size(result.CodeFiles));
for ii = 1:numel(result.CodeFiles)
    result.CodeSHA256(ii) = Case1_File_SHA256(result.CodeFiles(ii));
end
save(result.AuditFile, 'result', '-v7.3');
assert(height(result.Run.ReportV1) == height(plotCatalog) && ...
    ~any(contains(string(result.Run.ReportV1.Status), 'error')), ...
    'A daily production plot failed; inspect the saved report.');
if ~result.CompleteBatch
    fprintf('Partial batch completed: %d/47 daily images. Original baseline: %s\n', ...
        height(plotCatalog), baseline.AuditFile);
    fprintf('Resume the full batch with this baseline after missing sources are available.\n');
    return
end
result.Validation = Case1_Validate_Daily15_Plots(result, baseline);
save(result.AuditFile, 'result', '-v7.3');
copyfile(fullfile(cfg.CodeRoot, 'README_日均前后15天_L1优先.md'), ...
    fullfile(cfg.OutputRoot, 'Epoch_daily', 'README_日均前后15天_L1优先.md'), 'f');
assert(result.Validation.Passed, 'Daily artifact verification failed; inspect the audit.');
fprintf('Completed 47 daily +/-15-day figures. Run audit: %s\n', result.AuditFile);
end

function requireYears(files, years, label)
files = string(files(:));
for yy = years
    assert(any(contains(files, sprintf('_%d0101_v', yy))), ...
        'Missing entire %s source year %d. Download the official raw CDF first.', label, yy);
end
end

function baseline = captureBaseline(cfg, catalog)
baseline = struct('CreatedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'Catalog', catalog, 'Items', struct([]));
dailyFolder = fullfile(cfg.OutputRoot, 'Epoch_daily');
dataFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1d', ...
    'derived', 'pitch_angle', '2013-2021', 'predicted_ck');
dailyFiles = strings(height(catalog), 1);
for ii = 1:height(catalog)
    entries = dir(fullfile(dataFolder, sprintf('V1_%s_*_1d_nativeCDF_Epoch.mat', ...
        string(catalog.EventID(ii)))));
    assert(numel(entries) == 1, 'Expected one existing daily numerical audit.');
    auditFile = fullfile(entries.folder, entries.name);
    [~, name] = fileparts(auditFile);
    dailyFiles(ii) = string(fullfile(dailyFolder, [name, '.png']));
    assert(isfile(dailyFiles(ii)), 'An expected existing daily figure is missing.');
    baseline.Items(ii).EventID = string(catalog.EventID(ii));
    baseline.Items(ii).AuditFile = string(auditFile);
    baseline.Items(ii).AuditSHA256 = Case1_File_SHA256(auditFile);
    baseline.Items(ii).Saved = load(auditFile);
end
baseline.DailyFiles = hashFiles(dailyFiles);
entries = dir(fullfile(cfg.OutputRoot, '**', '*.png'));
protected = string(fullfile({entries.folder}, {entries.name})).';
protected = protected(~ismember(protected, dailyFiles));
hourlyRoot = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1h', 'derived');
entries = dir(fullfile(hourlyRoot, '**', '*.mat'));
protected = unique([protected; string(fullfile({entries.folder}, {entries.name})).']);
baseline.ProtectedFiles = hashFiles(protected);
end

function manifest = hashFiles(files)
files = string(files(:));
hashes = strings(size(files));
for ii = 1:numel(files), hashes(ii) = Case1_File_SHA256(files(ii)); end
manifest = table(files, hashes, 'VariableNames', {'SourceFile', 'SHA256'});
end
