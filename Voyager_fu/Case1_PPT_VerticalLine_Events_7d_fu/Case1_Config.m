function cfg = Case1_Config()
%Case1_Config Central configuration for the Case 1 event figures.
%   All formal Voyager data paths point to the classified archive on Z:.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-02

%% project paths
cfg.CodeRoot = fileparts(mfilename('fullpath'));
cfg.DataRoot = 'Z:\SPART-WORK\Data\Voyager';
cfg.IRFURoot = 'C:\Users\Administrator\Documents\irfu-matlab-master';
cfg.OutputRoot = ['C:\Users\Administrator\Documents\', ...
    'Recovery-Work-Voyager_betatron\Case1_PPT_VerticalLine_Events_7d'];

cfg.LECPSectoredDailyCDF = fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
    '1d', 'l2', 'sectored_flux', '2013-2021', ...
    'voyager1_lecp_p1_sectored_daily_case1_20130102_20210525.cdf');
% The historical server-generated subset is retained for comparison only.
% New runs use original annual CDFs, with their website filenames intact.
cfg.LECPNativeDailyCDFs = nativeFiles(cfg.DataRoot, '1d', 'daily');
cfg.LECPNativeHourlyCDFs = nativeFiles(cfg.DataRoot, '1h', 'hourly');
% User-approved restoration of the historical rate path, 2026-09-03.
cfg.LECPLevel1Fallback = true;
% Latest event request: complete L1 UTC means take priority over L2.
cfg.LECPSourcePriority = 'l1_first';
entries = dir(fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'native', 'l1', ...
    'sectored_rates', '*', 'voyager-1_lecp_lev-1-rates_*_v*.cdf'));
cfg.LECPLevel1CDFs = sort(string(fullfile({entries.folder}, {entries.name})).');
% Optional legacy external interface; the approved default reads raw CK.
cfg.LECPSectorPointingFile = '';
cfg.PitchAngleDataFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
    '1d', 'derived', 'pitch_angle', '2013-2021');

%% plotting window and display
cfg.ContextDays = 3; % Complete UTC window [event day-3, event day+4).
cfg.ShowEventBoundaries = false;
cfg.ExportDPI = 220;
cfg.LineYMarginFraction = 0.08;
cfg.MAGGapHours = 1.5;
cfg.FluxGapHours = 1.5;
cfg.LECPDailyAverage = false;
cfg.CRSDisplay = 'spectrogram';
cfg.CRSColorLimits = [];
cfg.ColorPercentiles = [1 99];

%% Voyager 1 LECP sector treatment
% Explicit user approval on 2026-09-03; geometry and transforms tested.
cfg.PredictedAttitudeApproved = true;
cfg.NominalLECPGeometryApproved = true;
cfg.AttitudeDailyHourUTC = 12; % Legacy compatibility only; Epoch mode ignores it.
cfg.PADCadence = 'day';
% Latest user instruction: anchor retained daily/hourly records at Epoch.
cfg.HourlyAttitudeApproved = true;
cfg.AccumulationPolicy = 'epoch_drop_negative_deltat';
cfg.LECPSectorAverageDays = 0; % Keep official product values without re-averaging.
cfg.PitchAngleMethod = 'predicted_ck';
cfg.LECPBackgroundMode = 'none';
cfg.PitchMergeToleranceDeg = 2;
cfg.ExportPitchAngleTable = false;
% User display convention (2026-09-03); original CDF bounds stay unchanged.
cfg.P1DisplayEnergyMeV = [0.57 1.78];
cfg.ExportPeakPAD = true; % Hourly: L2 Epoch or identified L1 bin center.
cfg.PeakPADFolder = fullfile(cfg.OutputRoot, 'Peak_PAD_5times_hourly');
cfg.PeakPADDataFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
    '1h', 'derived', 'peak_pitch_angle', '2013-2021');

%% output names
cfg.CatalogFile = fullfile(cfg.OutputRoot, ...
    'Case1_event_catalog_refactored.csv');
cfg.SourceManifestFile = fullfile(cfg.OutputRoot, ...
    'Case1_source_manifest_refactored.csv');
cfg.RunAuditFile = fullfile(cfg.OutputRoot, ...
    'Case1_data_processing_assumptions_refactored.txt');
cfg.ReadmeSource = fullfile(cfg.CodeRoot, 'README_数据处理与假设.md');
cfg.ReadmeCopy = fullfile(cfg.OutputRoot, ...
    'README_数据处理与假设_refactored.md');
end

function files = nativeFiles(dataRoot, cadence, product)
folder = fullfile(dataRoot, 'voyager1', 'lecp', cadence, 'l2', 'sectored_flux');
entries = dir(fullfile(folder, '*', ['voyager-1_lecp_lev-2-', ...
    product, '-avg_*_v*.cdf']));
files = string(fullfile({entries.folder}, {entries.name})).';
files = sort(files);
end
