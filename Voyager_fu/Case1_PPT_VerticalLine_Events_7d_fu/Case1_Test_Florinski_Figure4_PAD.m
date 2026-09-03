function report = Case1_Test_Florinski_Figure4_PAD()
%Case1_Test_Florinski_Figure4_PAD Synthetic four-day selection regression.
%   No satellite files are read or changed. Test output remains in the
%   classified Voyager validation folder; no figure is created.

%% synthetic daily-product rows and their unchanged source metadata
cfg = Case1_Config;
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.Checks = table('Size', [0 3], ...
    'VariableTypes', {'string', 'logical', 'string'}, ...
    'VariableNames', {'Check', 'Passed', 'Detail'});
opts = struct('RenderFigure', false, 'Visible', false);
startTime = datetime(2004, 8, 8, 'TimeZone', 'UTC');
time = startTime+days((0:3).')+hours([1; 2; 3; 4])+seconds([0.125; 0.250; 0.375; 0.500]);
data = table(time, true(4, 1), 'VariableNames', {'EpochUTC', 'PADUsable'});
maximum = [0.4; 0.2; 0.3; 0.6];
for iSector = 1:7
    data.(sprintf('RawFlux_S%d_1d', iSector)) = maximum*iSector/7;
    data.(sprintf('FluxUncertainty_S%d_1d', iSector)) = maximum*iSector/70;
    data.(sprintf('PA_S%d_deg', iSector)) = repmat((iSector-1)*30, 4, 1);
end
data.RawFlux_S8_1d = repmat(1e8, 4, 1);
data.SourceCDF = repmat("synthetic-v1-lecp-daily.cdf", 4, 1);
data.SourceCDFRecord = (221:224).';
data.P1EnergyLower_MeV = repmat(0.57, 4, 1);
data.P1EnergyUpper_MeV = repmat(0.89, 4, 1);
data.DeltaT_s = [0; NaN; 864000; 86400];
data.Properties.UserData = struct('Synthetic', true, 'Method', "preserved upstream pointing audit");
original = data;

%% complete record selection, fixed normalization and provenance
audit = Case1_Plot_Florinski_Figure4_PAD(data, '', opts);
addCheck('four_days_in_order', isequal(audit.SelectedTableRows, (1:4).') && ...
    isequal(audit.DayOfYear, (221:224).') && all(audit.RecordsPerDay == 1), ...
    'One unchanged original record maps to each Figure 4 calendar day.');
addCheck('complete_status', audit.PlottedCount == 4 && audit.SelectedCount == 4 && ...
    audit.Status == "complete" && all(audit.StatusPerDay == "complete"), ...
    'All four complete PAD records are available.');
addCheck('epoch_not_relabelled', isequal(audit.SelectedEpochUTC, time), ...
    'Original subsecond Epochs are retained without midday labels or averaging.');
addCheck('own_maximum_normalization', max(abs(audit.NormalizationFlux-maximum)) < 1e-15 && ...
    all(max(audit.NormalizedFlux, [], 2) == 1), ...
    'Each day has its own S1-S7 maximum and S8 cannot change the denominator.');
addCheck('sigma_fixed_denominator', isequaln(audit.NormalizedSigma, ...
    audit.RawSigma./audit.NormalizationFlux), ...
    'Uncertainty uses the same fixed daily denominator without additional propagation.');
addCheck('seven_separate_angles', isequal(audit.PA_deg, repmat(0:30:180, 4, 1)) && ...
    isequal(size(audit.RawFlux), [4 7]) && isequal(audit.Sectors, 1:7), ...
    'All seven original sector values are retained without angle merging.');
addCheck('source_energy_and_display', isequal(audit.SourceEnergyMeV, repmat([0.57 0.89], 4, 1)) && ...
    isequal(audit.DisplayEnergyMeV, [0.57 1.78]), ...
    'CDF source bounds and the approved display convention are separately retained.');
addCheck('source_provenance', isequal(audit.SourceCDF, data.SourceCDF) && ...
    isequal(audit.SourceCDFRecord, data.SourceCDFRecord) && ...
    audit.SelectedRows{3}.SourceCDFRecord == 223 && ...
    isequaln(audit.SourceTableUserData, data.Properties.UserData), ...
    'The full selected rows, file records and upstream scientific audit are preserved.');
addCheck('input_unchanged', isequaln(data, original), 'No input row or metadata is mutated.');
addCheck('no_figures_during_test', ~audit.FigureCreated && audit.FigureFile == "" && ...
    isequal(audit.DisplayXLimits_deg, [0 180]), ...
    'RenderFigure=false skips graphics and the angle display is explicitly 0-180 degrees.');

%% missing dates, missing B/PA, incomplete sectors and existing gates
missingDay = Case1_Plot_Florinski_Figure4_PAD(data([1 3 4], :), '', opts);
addCheck('missing_date_stays_blank', isequaln(missingDay.SelectedTableRows, [1; NaN; 2; 3]) && ...
    missingDay.PlottedCount == 3 && isnat(missingDay.SelectedEpochUTC(2)) && ...
    all(isnan(missingDay.NormalizedFlux(2, :))), ...
    'A missing reference day remains in its correct slot without a substitute record.');
missingB = data; missingB.PADUsable(2) = false;
for iSector = 1:7, missingB.(sprintf('PA_S%d_deg', iSector))(2) = NaN; end
missingBAudit = Case1_Plot_Florinski_Figure4_PAD(missingB, '', opts);
addCheck('missing_B_PA_stays_blank', ~missingBAudit.SelectedPADUsable(2) && ...
    missingBAudit.SelectedTableRows(2) == 2 && missingBAudit.PlottedCount == 3 && ...
    all(isfinite(missingBAudit.RawFlux(2, :))), ...
    'The original flux record remains audited while the missing-B PAD remains blank.');
missingFlux = data; missingFlux.RawFlux_S3_1d(3) = NaN;
missingFluxAudit = Case1_Plot_Florinski_Figure4_PAD(missingFlux, '', opts);
addCheck('incomplete_sector_stays_blank', ~missingFluxAudit.SelectedPADUsable(3) && ...
    missingFluxAudit.PlottedCount == 3, ...
    'Even a stale true PADUsable flag cannot bypass incomplete S1-S7 flux.');
zeroFlux = data; zeroFlux.RawFlux_S4_1d(4) = 0;
zeroAudit = Case1_Plot_Florinski_Figure4_PAD(zeroFlux, '', opts);
addCheck('positive_flux_gate_preserved', ~zeroAudit.SelectedPADUsable(4) && ...
    zeroAudit.RawFlux(4, 4) == 0, 'The existing positive-flux gate is preserved and the raw zero is audited.');
existingGate = data; existingGate.PADUsable(1) = false;
existingAudit = Case1_Plot_Florinski_Figure4_PAD(existingGate, '', opts);
addCheck('upstream_PADUsable_preserved', ~existingAudit.SelectedPADUsable(1) && ...
    existingAudit.PlottedCount == 3, 'The renderer cannot re-enable a rejected upstream PAD.');

%% one-day uniqueness, exact UTC boundaries and unordered source rows
reordered = data([4 2 1 3], :);
reorderedAudit = Case1_Plot_Florinski_Figure4_PAD(reordered, '', opts);
addCheck('unordered_input', isequal(reorderedAudit.SourceCDFRecord, (221:224).') && ...
    isequal(reorderedAudit.SelectedTableRows, [3; 2; 4; 1]), ...
    'Calendar-day selection remains independent of source table order.');
duplicate = [data; data(2, :)];
duplicate.EpochUTC(5) = duplicate.EpochUTC(5)+hours(1);
duplicateError = false;
try
    Case1_Plot_Florinski_Figure4_PAD(duplicate, '', opts);
catch exception
    duplicateError = strcmp(exception.identifier, 'Case1:FlorinskiMultipleDailyRecords');
end
addCheck('duplicate_daily_record_explicit_error', duplicateError, ...
    'Multiple retained daily rows raise an explicit error instead of arbitrary selection or averaging.');
boundary = data;
boundary.EpochUTC = [startTime; startTime+days(2)-milliseconds(1); ...
    startTime+days(2); startTime+days(4)-milliseconds(1)];
outside = data([1 4], :);
outside.EpochUTC = [startTime-milliseconds(1); startTime+days(4)];
boundary = [outside(1, :); boundary; outside(2, :)];
boundaryAudit = Case1_Plot_Florinski_Figure4_PAD(boundary, '', opts);
addCheck('UTC_half_open_boundaries', isequal(boundaryAudit.SelectedTableRows, (2:5).') && ...
    boundaryAudit.PlottedCount == 4, ...
    'Start midnight is included, end midnight is excluded, and original within-day Epochs are retained.');
zoned = data; zoned.EpochUTC.TimeZone = 'Asia/Shanghai';
zonedAudit = Case1_Plot_Florinski_Figure4_PAD(zoned, '', opts);
addCheck('timezone_converted_to_UTC', isequal(zonedAudit.SelectedEpochUTC, time) && ...
    strcmp(zonedAudit.SelectedEpochUTC.TimeZone, 'UTC'), ...
    'Calendar-day membership uses UTC even when input datetime has another display timezone.');

%% uncertain points are displayed without adding a scientific threshold
sigma = data; sigma.FluxUncertainty_S2_1d(1) = NaN;
sigma.FluxUncertainty_S3_1d(2) = -1;
sigma.FluxUncertainty_S4_1d(3) = 3;
sigmaAudit = Case1_Plot_Florinski_Figure4_PAD(sigma, '', opts);
lower = sigmaAudit.NormalizedFlux(3, 4)-sigmaAudit.NormalizedSigma(3, 4);
upper = sigmaAudit.NormalizedFlux(3, 4)+sigmaAudit.NormalizedSigma(3, 4);
addCheck('sigma_no_new_quality_gate', sigmaAudit.PlottedCount == 4 && ...
    isnan(sigmaAudit.NormalizedSigma(1, 2)) && sigmaAudit.NormalizedSigma(2, 3) < 0, ...
    'Missing or negative sigma suppresses that error bar only; source uncertainty is retained.');
addCheck('large_sigma_in_common_limits', sigmaAudit.DisplayYLimits(1) <= lower && ...
    sigmaAudit.DisplayYLimits(2) >= upper, 'Shared display limits contain all finite nonnegative error bars.');
emptyAudit = Case1_Plot_Florinski_Figure4_PAD(table, '', opts);
addCheck('empty_input_four_blank_days', emptyAudit.SelectedCount == 0 && ...
    emptyAudit.PlottedCount == 0 && numel(emptyAudit.DayStartUTC) == 4 && ...
    all(isnat(emptyAudit.SelectedEpochUTC)), 'Empty input preserves the four requested dates without invented data.');

%% save synthetic validation under the satellite data root
report.Passed = all(report.Checks.Passed);
report.PassedCount = sum(report.Checks.Passed);
report.CheckCount = height(report.Checks);
report.SyntheticInput = data;
report.SelectionAudit = audit;
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'florinski_figure4');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = fullfile(folder, ['florinski_figure4_synthetic_', stamp, '.mat']);
save(report.AuditFile, 'report', '-v7.3');
fprintf('Florinski Figure 4 PAD tests: %d/%d passed.\n', report.PassedCount, report.CheckCount);
if ~report.Passed
    disp(report.Checks(~report.Checks.Passed, :));
    error('Case1:FlorinskiPADTestFailure', 'Florinski Figure 4 PAD regression failed.');
end

    function addCheck(name, passed, detail)
        report.Checks(end+1, :) = {string(name), logical(passed), string(detail)};
    end
end
