function report = Case1_Test_Peak_PAD()
%Case1_Test_Peak_PAD Synthetic tests of peak selection and normalization.
%   Tests do not read, alter, average or interpolate satellite data.
%   Generated test audit is saved in the classified Voyager validation area.

%% synthetic original-record table
cfg = Case1_Config;
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.Checks = table('Size', [0 3], ...
    'VariableTypes', {'string', 'logical', 'string'}, ...
    'VariableNames', {'Check', 'Passed', 'Detail'});
startTime = datetime(2020, 5, 22, 0, 0, 0, 'TimeZone', 'UTC');
event = struct('EventID', 'synthetic-peak', 'PlotStartUTC', startTime, ...
    'PlotEndUTCExclusive', startTime+days(1));
opts = struct('PADCadence', 'hour', 'RenderFigure', false, ...
    'Visible', false, 'P1DisplayEnergyMeV', [0.57 1.78]);
time = startTime+hours((1:10).')+seconds((1:10).'/8);
data = table(time, true(10, 1), 'VariableNames', {'EpochUTC', 'PADUsable'});
maximum = [1; 2; 3; 4; 50; 6; 7; 8; 9; 10];
for iSector = 1:7
    data.(sprintf('RawFlux_S%d_1h', iSector)) = maximum*iSector/7;
    data.(sprintf('FluxUncertainty_S%d_1h', iSector)) = maximum*iSector/70;
    data.(sprintf('PA_S%d_deg', iSector)) = repmat((iSector-1)*30, 10, 1);
end
data.RawFlux_S8_1h = repmat(1e8, 10, 1);
data.PADUsable([4 6]) = false;
data.PA_S2_deg(4) = NaN;
data.RawFlux_S1_1h(6) = NaN;
data.FluxUncertainty_S2_1h(7) = NaN;
data.SourceRow = (1:10).';
data.SourceCDF = repmat("synthetic.cdf", 10, 1);
data.SourceCDFRecord = (101:110).';
data.P1EnergyLower_MeV = repmat(0.57, 10, 1);
data.P1EnergyUpper_MeV = repmat(0.89, 10, 1);
data.DeltaT_s = [0; NaN; 864000; repmat(3600, 7, 1)];
before = data;

%% normal peak with missing adjacent rows
audit = Case1_Plot_Peak_PAD(data, event, '', opts);
expected = [2; 3; 5; 7; 8];
addCheck('nearest_complete_records', isequal(audit.SelectedTableRows, expected), ...
    'Both missing immediate neighbors are skipped, retaining chronological slots.');
addCheck('peak_metric_and_time', audit.PeakTableRow == 5 && ...
    audit.PeakFlux == 50 && audit.PeakEpochUTC == time(5), ...
    'Peak is raw maximum S1-S7 flux and preserves the fractional original Epoch.');
addCheck('skipped_missing_audit', audit.SkippedMissingBeforeCount == 1 && ...
    audit.SkippedMissingAfterCount == 1 && ...
    isequal(audit.SkippedMissingBeforeTableRows, 4) && ...
    isequal(audit.SkippedMissingAfterTableRows, 6), ...
    'Both skipped original table rows are recorded.');
addCheck('normalization_own_maximum', ...
    isequal(audit.NormalizationFlux, maximum(expected)) && ...
    all(max(audit.NormalizedFlux, [], 2) == 1), ...
    'Each panel independently reaches one; the peak denominator is not shared.');
addCheck('sigma_fixed_scale', isequaln(audit.NormalizedSigma, ...
    audit.RawSigma./audit.NormalizationFlux), ...
    'Sigma has exactly the same fixed denominator, without denominator propagation.');
addCheck('seven_original_points', isequal(size(audit.RawFlux), [5 7]) && ...
    isequal(audit.PA_deg, repmat(0:30:180, 5, 1)), ...
    'All seven sector centers remain separate and S8 cannot select the peak.');
addCheck('missing_sigma_not_a_PAD_gate', audit.SelectedPADUsable(4) && ...
    isnan(audit.NormalizedSigma(4, 2)) && isfinite(audit.NormalizedFlux(4, 2)), ...
    'The point with missing uncertainty remains selected.');
addCheck('source_metadata_and_display_separated', ...
    isequal(audit.SourceEnergyMeV, repmat([0.57 0.89], 5, 1)) && ...
    isequal(audit.DisplayEnergyMeV, [0.57 1.78]), ...
    'Raw CDF energy metadata is preserved separately from the approved label.');
addCheck('source_row_provenance', audit.SelectedRows{2}.SourceCDFRecord == 103 && ...
    audit.SelectedRows{2}.SourceCDF == "synthetic.cdf", ...
    'Selected full rows preserve original file and record mappings.');
addCheck('all_five_complete', audit.SelectedCount == 5 && audit.PlottedCount == 5 && ...
    audit.Status == "complete", 'Expected selected/plotted counts and status.');
addCheck('input_unchanged', isequaln(data, before), 'The input table is not mutated.');

%% incomplete actual peak remains the center
missingPeak = data; missingPeak.PADUsable(5) = false;
missingPeak.PA_S3_deg(5) = NaN;
incomplete = Case1_Plot_Peak_PAD(missingPeak, event, '', opts);
addCheck('true_peak_not_replaced', incomplete.PeakTableRow == 5 && ...
    incomplete.SelectedTableRows(3) == 5 && ~incomplete.SelectedPADUsable(3) && ...
    incomplete.PlottedCount == 4 && incomplete.Status == "peak_pad_unavailable", ...
    'The maximum-flux record remains central when its PA is unavailable.');

%% peak ties, unordered rows, window boundaries and limited neighbors
ties = data;
for iSector = 1:7
    ties.(sprintf('RawFlux_S%d_1h', iSector))(9) = 50*iSector/7;
end
tieAudit = Case1_Plot_Peak_PAD(ties, event, '', opts);
addCheck('ties_earliest_epoch', tieAudit.PeakTableRow == 5, ...
    'Equal-flux later records do not displace the earliest peak.');
unordered = data([10 5 3 9 6 1 7 4 8 2], :);
unorderedAudit = Case1_Plot_Peak_PAD(unordered, event, '', opts);
selectedSources = cellfun(@(t) t.SourceRow, unorderedAudit.SelectedRows);
addCheck('unordered_input_time_selection', isequal(selectedSources, expected), ...
    'Neighbors are selected by original Epoch, independently of input row order.');
limitedEvent = event;
limitedEvent.PlotStartUTC = time(5);
limitedEvent.PlotEndUTCExclusive = time(9);
limited = Case1_Plot_Peak_PAD(data, limitedEvent, '', opts);
addCheck('window_and_missing_slots', isequaln(limited.SelectedTableRows, ...
    [NaN; NaN; 5; 7; 8]) && limited.PlottedCount == 3 && ...
    limited.Status == "partial_neighbors", ...
    'No outside-window records are used and unavailable slots are not fabricated.');
single = Case1_Plot_Peak_PAD(data(5, :), event, '', opts);
addCheck('single_record', isequaln(single.SelectedTableRows, [NaN; NaN; 1; NaN; NaN]) && ...
    single.SelectedCount == 1 && single.PlottedCount == 1, ...
    'A one-record window retains a centered peak with four empty slots.');
oneBefore = data; oneBefore.PADUsable([1 2]) = false;
oneBeforeAudit = Case1_Plot_Peak_PAD(oneBefore, event, '', opts);
addCheck('missing_neighbor_search_reaches_boundary', ...
    isequaln(oneBeforeAudit.SelectedTableRows, [NaN; 3; 5; 7; 8]) && ...
    oneBeforeAudit.SkippedMissingBeforeCount == 3, ...
    'When fewer than two good neighbors exist, all searched missing rows are audited.');

%% no flux, no records, nonfinite flux and large uncertainty
noFlux = data;
for iSector = 1:7, noFlux.(sprintf('RawFlux_S%d_1h', iSector))(:) = NaN; end
none = Case1_Plot_Peak_PAD(noFlux, event, '', opts);
addCheck('no_positive_flux', none.Status == "no_positive_sector_flux" && ...
    none.SelectedCount == 0 && isnat(none.PeakEpochUTC), ...
    'S8 cannot provide a peak when all S1-S7 fluxes are absent.');
empty = Case1_Plot_Peak_PAD(table, event, '', opts);
addCheck('empty_input', empty.Status == "no_retained_records" && ...
    all(isnan(empty.SelectedTableRows)) && isnat(empty.PeakEpochUTC), ...
    'An empty source table has no invented times or fluxes.');
nonfinite = data;
nonfinite.RawFlux_S1_1h(1) = Inf;
nonfinite.RawFlux_S7_1h(10) = -1e10;
nonfiniteAudit = Case1_Plot_Peak_PAD(nonfinite, event, '', opts);
addCheck('peak_positive_finite_only', nonfiniteAudit.PeakTableRow == 5, ...
    'Nonfinite and negative flux cannot win the peak metric.');
largeSigma = data;
largeSigma.FluxUncertainty_S2_1h(3) = 300;
large = Case1_Plot_Peak_PAD(largeSigma, event, '', opts);
upper = large.NormalizedFlux(2, 2)+large.NormalizedSigma(2, 2);
lower = large.NormalizedFlux(2, 2)-large.NormalizedSigma(2, 2);
addCheck('large_sigma_kept_and_visible', large.PlottedCount == 5 && ...
    large.DisplayYLimits(2) >= upper && large.DisplayYLimits(1) <= lower, ...
    'Large uncertainty is retained and accommodated by shared display limits.');

%% save reproducible test results
report.Passed = all(report.Checks.Passed);
report.PassedCount = sum(report.Checks.Passed);
report.CheckCount = height(report.Checks);
report.SyntheticInput = data;
report.SelectionAudit = audit;
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'peak_pad');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = fullfile(folder, ['peak_pad_synthetic_', stamp, '.mat']);
save(report.AuditFile, 'report', '-v7.3');
fprintf('Peak PAD synthetic tests: %d/%d passed.\n', report.PassedCount, report.CheckCount);
if ~report.Passed
    disp(report.Checks(~report.Checks.Passed, :));
    error('Case1:PeakPADTestFailure', 'Peak PAD synthetic regression failed.');
end

    function addCheck(name, passed, detail)
        report.Checks(end+1, :) = {string(name), logical(passed), string(detail)};
    end
end
