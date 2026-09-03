function report = Case1_Test_Nearest_MAG()
%Case1_Test_Nearest_MAG Synthetic regression of nearest original MAG rows.
%   No source CDF or scientific figure is changed. The audit is stored under
%   the project's Voyager data root. Test gap limits exercise the interface;
%   they do not select the limit used in a scientific processing run.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% configuration
cfg = Case1_Config;
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
Case1_Add_IRFU_Path(cfg.IRFURoot);
report = struct;
report.StartedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.TestSHA256 = Case1_File_SHA256(report.TestFile);
report.HelperFile = fullfile(fileparts(report.TestFile), 'Case1_Match_Nearest_MAG.m');
report.HelperSHA256 = Case1_File_SHA256(report.HelperFile);
report.MATLABVersion = version;
report.Checks = table('Size', [0 3], ...
    'VariableTypes', {'string', 'logical', 'string'}, ...
    'VariableNames', {'Check', 'Passed', 'Detail'});
report.FatalError = "";
report.Cases = struct;

try
    %% raw row matching and time direction
    startTime = datetime(2020, 1, 1, 'TimeZone', 'UTC');
    candidates = makeCandidates([4; 0; 2], [3 4 0; 1 2 3; -4 0 0]);
    query = startTime + hours([-1; 0; 0.25; 1; 1.75; 2; 3; 4; 5]);
    result = Case1_Match_Nearest_MAG(query, candidates, Inf);
    expectedRows = [2; 2; 2; 2; 3; 3; 3; 1; 1];
    rawB = [candidates.BR, candidates.BT, candidates.BN];
    addCheck('nearest_unsorted_original_rows', isequal(result.SourceRow, expectedRows), ...
        'Nearest points on both sides and endpoint extrapolation retain original row numbers.');
    addCheck('copy_original_vector_exactly', isequal(result.B_RTN, rawB(expectedRows, :)), ...
        'Raw vectors copied unchanged; no direction normalization or component interpolation.');
    addCheck('tie_earlier_epoch', result.SourceRow(4) == 2 && result.SourceRow(7) == 3, ...
        'Equal distances select the earlier Epoch.');
    addCheck('signed_delta_seconds', isequal(result.DeltaSeconds, ...
        [3600; 0; -900; -3600; 900; 0; -3600; 0; -3600]), ...
        'Positive delta means magnetic source follows particle Epoch.');
    addCheck('absolute_delta_hours', isequal(result.AbsDeltaHours, abs(result.DeltaSeconds)/3600), ...
        'Absolute gap agrees with the signed raw Epoch difference.');
    addCheck('source_provenance', isequal(result.SourceCDF, string(candidates.SourceCDF(expectedRows))) && ...
        isequal(result.SourceCDFRecord, candidates.SourceCDFRecord(expectedRows)) && ...
        isequal(result.SourceEpochUTC, candidates.EpochUTC(expectedRows)), ...
        'Selected source CDF, record and Epoch preserved.');
    addCheck('unlimited_acceptance', all(result.Accepted) && all(result.Status == "nearest_valid_source"), ...
        'Inf is accepted only as an explicit caller choice.');
    report.Cases.RawMatching = result;

    %% inclusive cap and rejected-candidate audit
    query = startTime + hours([1; 1.0001; 0]);
    one = candidates(2, :);
    limited = Case1_Match_Nearest_MAG(query, one, 1);
    addCheck('gap_inclusive_boundary', isequal(limited.Accepted, [true; false; true]), ...
        'Exactly one hour is accepted; a larger gap is rejected.');
    addCheck('rejected_B_remains_nan', all(isnan(limited.B_RTN(2, :))), ...
        'A rejected candidate cannot be consumed as a valid B vector.');
    addCheck('rejected_candidate_audit', isequal(limited.CandidateB_RTN(2, :), [1 2 3]) && ...
        limited.SourceRow(2) == 1 && ~isnat(limited.SourceEpochUTC(2)) && ...
        limited.SourceCDF(2) == one.SourceCDF && limited.Status(2) == "gap_exceeds_limit", ...
        'Rejected candidate vector, provenance and gap remain visible to the caller.');
    exact = Case1_Match_Nearest_MAG(startTime + seconds([0; 0.125]), one, 0);
    addCheck('zero_gap_exact_only', isequal(exact.Accepted, [true; false]), ...
        'A zero-hour cap accepts exact source timestamps only.');
    report.Cases.GapLimit = limited;

    %% missing, partial, zero and infinite source values
    bad = makeCandidates([0; 1; 2; 3; 4], ...
        [NaN 1 2; 0 0 0; 1 Inf 2; 2 3 4; -2 0 1]);
    bad.EpochUTC(4) = NaT(1, 1, 'TimeZone', 'UTC');
    badQuery = [startTime; NaT(1, 1, 'TimeZone', 'UTC')];
    selected = Case1_Match_Nearest_MAG(badQuery, bad, Inf);
    addCheck('reject_partial_zero_infinite_and_NaT_source', ...
        isequal(selected.ValidSourceRows, 5) && selected.SourceRow(1) == 5, ...
        'Only complete finite nonzero vectors with real source Epochs enter the match.');
    addCheck('NaT_query_stays_missing', isnan(selected.SourceRow(2)) && ...
        all(isnan(selected.B_RTN(2, :))) && ~selected.Accepted(2) && ...
        selected.Status(2) == "invalid_query_epoch", ...
        'No artificial timestamp is assigned to an invalid query.');
    noValid = Case1_Match_Nearest_MAG(startTime, bad(1:4, :), Inf);
    addCheck('no_valid_sources', ~noValid.Accepted && isnan(noValid.SourceRow) && ...
        noValid.Status == "no_valid_mag_candidate", 'All-invalid candidates leave the query unavailable.');
    empty = Case1_Match_Nearest_MAG(startTime, candidates([], :), Inf);
    addCheck('empty_sources', ~empty.Accepted && isnat(empty.SourceEpochUTC) && ...
        isempty(empty.ValidSourceRows), 'An empty source table returns an explicit unavailable record.');
    noQuery = Case1_Match_Nearest_MAG(NaT(0, 1, 'TimeZone', 'UTC'), candidates, Inf);
    addCheck('empty_query_shape', isequal(size(noQuery.B_RTN), [0 3]) && isempty(noQuery.SourceRow), ...
        'Empty query output retains its column dimensions.');
    report.Cases.InvalidSources = selected;

    %% duplicate Epochs and source ambiguity
    duplicates = makeCandidates([2; 0; 2; 4], [1 2 3; 4 5 6; 1 2 3; 7 8 9]);
    duplicateResult = Case1_Match_Nearest_MAG(startTime + hours(2), duplicates, Inf);
    addCheck('identical_duplicate_stable_original_row', duplicateResult.SourceRow == 1 && ...
        isequal(duplicateResult.DuplicateIdenticalSourceRows, 3), ...
        'Identical duplicate vectors retain the earliest input row.');
    duplicates.BN(3) = 30;
    addCheck('conflicting_duplicate_error', throwsError( ...
        @() Case1_Match_Nearest_MAG(startTime + hours(2), duplicates, Inf), ...
        'VoyagerCase1:MAGConflictingDuplicate'), ...
        'Conflicting duplicate Epoch vectors raise an explicit error; no farther source is substituted.');
    addCheck('conflicting_duplicates_only_error', throwsError( ...
        @() Case1_Match_Nearest_MAG(startTime, duplicates([1 3], :), Inf), ...
        'VoyagerCase1:MAGConflictingDuplicate'), ...
        'An entirely ambiguous source set raises the same explicit error.');
    report.Cases.IdenticalDuplicates = duplicateResult;

    %% fractional timestamps, query ordering and zones
    fraction = makeCandidates([0; 1], [1 2 3; 4 5 6]);
    fraction.EpochUTC = startTime + seconds([0.125; 0.875]);
    fractionalQuery = startTime + seconds([0.75; 0.5; 0.125]);
    frac = Case1_Match_Nearest_MAG(fractionalQuery, fraction, Inf);
    addCheck('fractional_epoch_and_unsorted_queries', isequal(frac.SourceRow, [2; 1; 1]) && ...
        isequal(frac.DeltaSeconds, [0.125; -0.375; 0]), ...
        'Fractional seconds and the original query order are retained.');
    fraction.EpochUTC.TimeZone = 'Asia/Shanghai';
    fractionalQuery.TimeZone = 'America/New_York';
    zoned = Case1_Match_Nearest_MAG(fractionalQuery, fraction, Inf);
    addCheck('timezone_same_instant', isequal(zoned.SourceRow, frac.SourceRow) && ...
        isequal(zoned.DeltaSeconds, frac.DeltaSeconds) && ...
        strcmp(zoned.SourceEpochUTC.TimeZone, 'UTC'), ...
        'Different zoned representations are matched as the same physical instant.');
    addCheck('policy_audit', contains(string(result.Policy), 'earlier Epoch') && ...
        result.MaxGapHours == Inf, 'The tie rule, exact-copy policy and explicit cap are auditable.');
    addCheck('negative_cap_rejected', throwsError(@() Case1_Match_Nearest_MAG(startTime, candidates, -1)), ...
        'A negative gap cap fails input validation.');
    addCheck('NaN_cap_rejected', throwsError(@() Case1_Match_Nearest_MAG(startTime, candidates, NaN)), ...
        'An undefined gap cap fails input validation.');
    noZone = startTime;
    noZone.TimeZone = '';
    addCheck('unqualified_timezone_rejected', throwsError(@() Case1_Match_Nearest_MAG(noZone, candidates, Inf)), ...
        'Unzoned datetime values are never silently interpreted as UTC.');
catch ME
    report.FatalError = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
end

%% save the reproducible validation audit
report.CompletedUTC = datetime('now', 'TimeZone', 'UTC');
report.AllPassed = strlength(report.FatalError) == 0 && all(report.Checks.Passed);
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'nearest_mag');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = fullfile(folder, ['nearest_mag_synthetic_', stamp, '.mat']);
save(report.AuditFile, 'report', '-v7.3');
disp(report.Checks);
fprintf('Nearest MAG synthetic checks: %d/%d passed.\n', ...
    sum(report.Checks.Passed), height(report.Checks));
fprintf('Audit: %s\n', report.AuditFile);
if ~report.AllPassed
    error('VoyagerCase1:NearestMAGTestFailed', '%s', char(report.FatalError));
end

    function data = makeCandidates(timeHours, B)
        EpochUTC = startTime + hours(timeHours(:));
        BR = B(:, 1); BT = B(:, 2); BN = B(:, 3);
        SourceCDF = "synthetic_source_" + string((1:numel(timeHours)).') + ".cdf";
        SourceCDFRecord = (101:100 + numel(timeHours)).';
        data = table(EpochUTC, BR, BT, BN, SourceCDF, SourceCDFRecord);
    end

    function addCheck(name, passed, detail)
        report.Checks(end + 1, :) = {string(name), logical(passed), string(detail)};
    end

    function failed = throwsError(action, identifier)
        failed = false;
        try
            action();
        catch ME
            if nargin == 1
                failed = true;
            else
                failed = strcmp(ME.identifier, identifier);
            end
        end
    end
end
