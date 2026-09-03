function report = Case1_Test_L1_Priority()
%Case1_Test_L1_Priority L1-first regression, including the unchanged 51 tests.
%   Synthetic tests only. L1 rate means, error propagation, source time
%   selection and whole-record provenance are held fixed across priorities.

%% existing 51 tests establish backward-compatible default behavior
legacy = Case1_Test_L1_Fallback;
l2 = legacy.L2Input; rates = legacy.L1Input;
t0 = datetime(2004, 8, 8, 'TimeZone', 'UTC');
% Hour 3 has an incomplete L1 mean and a complete L2 observation.
l2.FHDU_SectoredFluxes(4, 10, 7) = 77;
beforeL2 = l2; beforeRates = rates;
default = Case1_Apply_L1_Fallback(l2, rates, t0, t0+hours(5), 'hour');
l2First = Case1_Apply_L1_Fallback(l2, rates, t0, t0+hours(5), 'hour', 'l2_first');
l1First = Case1_Apply_L1_Fallback(l2, rates, t0, t0+hours(5), 'hour', 'l1_first');
nCheck = 0;
check(isequaln(default, l2First));
check(isequaln(l2, beforeL2) && isequaln(rates, beforeRates));
check(l1First.L1FallbackAudit.SourcePriority == "l1_first");
check(l2First.L1FallbackAudit.SourcePriority == "l2_first");
check(contains(l1First.L1FallbackAudit.Method, 'L1-first'));
check(contains(l1First.L1FallbackAudit.MethodDescription, 'including complete L2'));
check(numel(l1First.Epoch) == 4);
check(isequal(l1First.Epoch, t0+minutes([30; 90; 150; 190])));
check(isequal(l1First.SourceProduct, ["L1_UTC_mean"; "L1_UTC_mean"; "L1_UTC_mean"; "L2"]));
check(isequaln(l1First.OriginalL2SourceRow, [NaN; NaN; NaN; 4]));
check(l1First.L1FallbackAudit.L2ReplacedRecords == 3);
check(isequal(l1First.L1FallbackAudit.ReplacedL2.OriginalRows, [1; 2; 3]));
check(isequaln(l1First.L1FallbackAudit.ReplacedL2.FHDU_SectoredFluxes, l2.FHDU_SectoredFluxes(1:3, :, :)));
check(isequaln(l1First.L1FallbackAudit.ReplacedL2.FHDU_SectoredFluxUncertainties, l2.FHDU_SectoredFluxUncertainties(1:3, :, :)));
check(isequaln(l1First.L1FallbackAudit.ReplacedL2.Epoch, l2.Epoch(1:3)));
check(isequaln(l1First.L1FallbackAudit.ReplacedL2.DeltaT, l2.DeltaT(1:3)));
check(isequal(l1First.L1FallbackAudit.ReplacedL2.SourceRecordNumber, [101; 102; 103]));
check(isequaln(l1First.FHDU_SectoredFluxes(4, :, :), l2.FHDU_SectoredFluxes(4, :, :)));
check(l1First.L1FallbackAudit.Candidates.Decision(1) == "l1_first_replace_all_L2_with_complete_L1_bin_mean");
check(l1First.L1FallbackAudit.Candidates.Decision(4) == "l1_first_incomplete_L1_mean_keep_original_L2");
check(l1First.L1FallbackAudit.Candidates.Decision(3) == "l1_first_add_complete_L1_bin_mean");

%% both priorities use exactly the same source means, error and timestamps
sameTimes = ismember(l2First.Epoch, l1First.Epoch(2:3));
check(isequaln(l2First.FHDU_SectoredFluxes(sameTimes, :, :), l1First.FHDU_SectoredFluxes(2:3, :, :)));
check(isequaln(l2First.FHDU_SectoredFluxUncertainties(sameTimes, :, :), l1First.FHDU_SectoredFluxUncertainties(2:3, :, :)));
check(isequal(l2First.Epoch(sameTimes), l1First.Epoch(2:3)));
check(isequaln(l2First.L1FallbackAudit.Candidates.MeanRate, l1First.L1FallbackAudit.Candidates.MeanRate));
check(isequaln(l2First.L1FallbackAudit.Candidates.MeanRateSigma, l1First.L1FallbackAudit.Candidates.MeanRateSigma));
check(isequal(l2First.L1FallbackAudit.Candidates.SectorSampleCount, l1First.L1FallbackAudit.Candidates.SectorSampleCount));
check(isequal(l1First.L1SourceRecords{2}.SourceCDFRecord, [102; 103]));
check(isequal(l1First.L1SourceRecords{2}.DeltaT_s, [5000; 90000]));
check(all(l1First.SectorSampleCount(2, :) == 2));
check(isnan(l1First.FHDU_SectoredFluxUncertainties(2, 10, 2)));
check(l1First.L1FallbackAudit.L1RecordSelection.NegativeDeltaTRejected == 1);
check(l1First.L1FallbackAudit.L2RecordSelection.NegativeDeltaTRejected == 1);
check(~isfield(l1First, 'ParticleUN') && ~isfield(l1First, 'BR'));
noRates = rates; noRates.Epoch = rates.Epoch+days(100);
noL1 = Case1_Apply_L1_Fallback(l2, noRates, t0, t0+hours(5), 'hour', 'l1_first');
check(isequaln(noL1.FHDU_SectoredFluxes, l2.FHDU_SectoredFluxes(1:4, :, :)));
check(all(noL1.L1FallbackAudit.Candidates.Decision == "l1_first_no_L1_Epoch_keep_original_L2"));
check(noL1.L1FallbackAudit.L1AddedRecords == 0 && noL1.L1FallbackAudit.L2ReplacedRecords == 0);

%% seven-day half-open event window: event day minus 3 through plus 3
beginTime = t0-days(3); stopTime = t0+days(4);
windowL2 = l2; windowRates = rates;
windowL2.Epoch = [beginTime-seconds(1); beginTime; beginTime+hours(1); ...
    stopTime-seconds(1); stopTime; stopTime+seconds(1)];
windowL2.DeltaT(:) = 3600;
windowRates.Epoch = [beginTime-seconds(1); beginTime; beginTime+minutes(20); ...
    beginTime+hours(1); t0; stopTime-hours(1); stopTime-seconds(1); ...
    stopTime; stopTime+seconds(1); stopTime+hours(1)];
windowRates.DeltaT(:) = 3600;
windowRates.FHDU_SectoredRates(:, 10, :) = 2;
windowRates.FHDU_SectoredRateUncertainties(:, 10, :) = 0.1;
sevenDay = Case1_Apply_L1_Fallback(windowL2, windowRates, beginTime, stopTime, 'hour', 'l1_first');
check(stopTime-beginTime == days(7));
check(all(sevenDay.Epoch >= beginTime & sevenDay.Epoch < stopTime));
check(sevenDay.L1FallbackAudit.StartUTC == beginTime && sevenDay.L1FallbackAudit.EndUTC == stopTime);
check(isequal(sevenDay.L1FallbackAudit.L1RecordSelection.KeptSourceRows, (2:7).'));
check(isequal(sevenDay.L1FallbackAudit.L2RecordSelection.KeptSourceRows, (2:4).'));
check(isequal(sevenDay.L1SourceRecords{1}.SourceCDFRecord, [102; 103]));
check(isequal(sevenDay.L1SourceRecords{end}.SourceCDFRecord, [106; 107]));
check(sevenDay.Epoch(1) == beginTime+minutes(30));
check(sevenDay.Epoch(end) == stopTime-minutes(30));
check(all(sevenDay.SourceProduct == "L1_UTC_mean"));
check(isequal(sevenDay.L1FallbackAudit.ReplacedL2.OriginalRows, (2:4).'));
dayL1 = Case1_Apply_L1_Fallback(windowL2, windowRates, beginTime, stopTime, 'day', 'l1_first');
check(isequal(dayL1.Epoch, [beginTime; t0; stopTime-days(1)]+hours(12)));
check(all(dayL1.SourceProduct == "L1_UTC_mean"));
check(dayL1.L1FallbackAudit.SourcePriority == "l1_first");

%% invalid priority cannot silently switch scientific products
invalidCaught = false;
try
    Case1_Apply_L1_Fallback(l2, rates, t0, t0+hours(5), 'hour', 'automatic');
catch ME
    invalidCaught = strcmp(ME.identifier, 'VoyagerCase1:L1SourcePriority');
end
check(invalidCaught);

%% standalone validation record, all synthetic data under Voyager archive
report = struct('Passed', true, 'AssertionCount', nCheck, ...
    'LegacyAssertionCount', legacy.AssertionCount, ...
    'LegacyAuditFile', legacy.AuditFile, 'L1First', l1First, ...
    'L2First', l2First, 'SevenDay', sevenDay, 'SevenDayDaily', dayL1);
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.CodeFile = string(fullfile(fileparts(mfilename('fullpath')), 'Case1_Apply_L1_Fallback.m'));
report.CodeSHA256 = string(Case1_File_SHA256(char(report.CodeFile)));
report.TestFile = string([mfilename('fullpath'), '.m']);
report.TestSHA256 = string(Case1_File_SHA256(char(report.TestFile)));
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(fileparts(legacy.AuditFile), ['l1_priority_synthetic_', stamp, '.mat']));
save(report.AuditFile, 'report', '-v7.3');
fprintf('L1 priority regression: %d/%d new tests; %d/%d legacy tests passed.\n', ...
    nCheck, nCheck, legacy.AssertionCount, legacy.AssertionCount);
disp(report.AuditFile);

    function check(value)
        nCheck = nCheck+1;
        assert(isscalar(value) && value, 'L1 priority assertion %d failed.', nCheck);
    end
end
