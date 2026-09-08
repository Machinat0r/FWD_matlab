function report = Case1_Test_LECP_Single_Year()
%Case1_Test_LECP_Single_Year Check single-file 2020 source-index bookkeeping.
%   Read only the original V1 2020 L1 and daily/hourly L2 CDFs. Compare every
%   retained payload with its source record, then exercise the mixed L1/L2
%   supplemental event window. No figure, source CDF or old audit is changed.

%% fixed three-source scope and current original-CDF readers
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
files = strings(3, 1);
files(1) = one2020File(cfg.LECPLevel1CDFs);
files(2) = one2020File(cfg.LECPNativeDailyCDFs);
files(3) = one2020File(cfg.LECPNativeHourlyCDFs);
roles = ["L1_native"; "L2_daily"; "L2_hourly"];
report = struct('CreatedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'Scope', 'Only three original V1 2020 LECP CDFs and [2020-07-27,2020-08-03) UTC; source-index orientation regression.', ...
    'Checks', table('Size', [0 4], ...
    'VariableTypes', {'string', 'string', 'logical', 'string'}, ...
    'VariableNames', {'Product', 'Check', 'Passed', 'Detail'}));
products = cell(3, 1);
beforeHash = strings(3, 1); afterHash = beforeHash;
sourceCount = zeros(3, 1); retainedCount = sourceCount;
currentRole = "";

%% each index must identify an unchanged original CDF record
for ii = 1:3
    currentRole = roles(ii);
    beforeHash(ii) = string(Case1_File_SHA256(char(files(ii))));
    raw = Voyager_Read_CDF_Product(files(ii), 'lecp_sector_daily');
    if ii == 1
        p = Case1_Read_LECP_Rates(files(ii));
        fields = {'Epoch', 'DeltaT', 'FHDU_SectoredRates', ...
            'FHDU_SectoredRateUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
    else
        p = Case1_Read_LECP_CDFs(files(ii));
        fields = {'Epoch', 'DeltaT', 'FHDU_SectoredFluxes', ...
            'FHDU_SectoredFluxUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
    end
    n = numel(p.Epoch); records = p.SourceRecordNumber;
    sourceCount(ii) = numel(raw.Epoch); retainedCount(ii) = n;
    addCheck('file_index_column', isequal(size(p.SourceFileIndex), [n, 1]), ...
        sprintf('SourceFileIndex size %s, expected [%d 1].', mat2str(size(p.SourceFileIndex)), n));
    addCheck('record_index_column', isequal(size(records), [n, 1]), ...
        sprintf('SourceRecordNumber size %s.', mat2str(size(records))));
    addCheck('valid_original_record_indices', all(isfinite(records) & ...
        records >= 1 & records <= sourceCount(ii) & records == fix(records)), ...
        'One-based integer indices into the original CDF.');
    addCheck('single_source_identity', height(p.SourceManifest) == 1 && ...
        string(p.SourceManifest.SourceFile) == files(ii) && ...
        all(p.SourceFileIndex(:) == 1), files(ii));
    addCheck('source_hash_recorded', string(p.SourceManifest.SHA256) == beforeHash(ii), ...
        'Reader manifest hash equals the directly hashed original CDF.');
    addCheck('sorted_unique_epochs', issorted(p.Epoch) && ...
        numel(unique(p.Epoch)) == n, 'Existing identical-Epoch deduplication retained.');
    for jj = 1:numel(fields)
        field = fields{jj}; reference = raw.(field);
        reference = reference(records, :, :);
        addCheck(['unchanged_', field], isequaln(p.(field), reference), ...
            'Exact equality with the indicated original record after the common CDF quality screening.');
    end
    addCheck('energy_metadata_unchanged', ...
        isequaln(p.variable_meta.FHDU_EnergyRange, raw.variable_meta.FHDU_EnergyRange), ...
        'The shape fix does not alter energy metadata or labels.');
    addCheck('duplicate_count_consistent', ...
        p.DuplicateIdenticalRecordsRemoved == sourceCount(ii)-n, ...
        'Existing duplicate handling and source record accounting retained.');
    afterHash(ii) = string(Case1_File_SHA256(char(files(ii))));
    addCheck('source_CDF_unchanged', beforeHash(ii) == afterHash(ii), ...
        'Original file hash is identical before and after this test.');
    products{ii} = p;
end

%% the exact supplementary window previously exposed a row/column mismatch
startUTC = datetime(2020, 7, 27, 'TimeZone', 'UTC');
endUTC = datetime(2020, 8, 3, 'TimeZone', 'UTC');
report.WindowStartUTC = startUTC; report.WindowEndUTCExclusive = endUTC;
for ii = 2:3
    cadence = 'day';
    if ii == 3, cadence = 'hour'; end
    currentRole = "supplemental_"+string(cadence);
    out = Case1_Apply_L1_Fallback(products{ii}, products{1}, ...
        startUTC, endUTC, cadence, 'l1_first');
    n = numel(out.Epoch);
    addCheck('combined_source_index_column', isequal(size(out.SourceFileIndex), [n, 1]), ...
        'Whole-record L1/L2 assembly completes with a column source index.');
    addCheck('all_epochs_in_event_window', all(out.Epoch >= startUTC & out.Epoch < endUTC), ...
        'Only the new seven-day event window is exercised.');
    addCheck('combined_source_indices_valid', all(out.SourceFileIndex >= 1 & ...
        out.SourceFileIndex <= height(out.SourceManifest)), ...
        'Every output row maps to the original CDF source manifest.');
    if ii == 3
        nL1 = nnz(out.SourceProduct == "L1_UTC_mean");
        nL2 = nnz(out.SourceProduct == "L2");
        addCheck('mixed_hourly_regression', n == 60 && nL1 == 54 && nL2 == 6, ...
            sprintf('%d output rows: %d L1 UTC means and %d retained original L2.', n, nL1, nL2));
        report.HourlyRows = n; report.HourlyL1Rows = nL1; report.HourlyL2Rows = nL2;
    else
        report.DailyRows = n;
    end
end

%% classified regression record
report.Sources = table(roles, files, beforeHash, afterHash, sourceCount, retainedCount, ...
    'VariableNames', {'Product', 'SourceFile', 'SHA256Before', 'SHA256After', ...
    'OriginalRecords', 'RetainedRecords'});
report.CodeFiles = string(fullfile(cfg.CodeRoot, ...
    {'Case1_Test_LECP_Single_Year.m', 'Case1_Read_LECP_CDFs.m', ...
    'Case1_Read_LECP_Rates.m', 'Case1_Apply_L1_Fallback.m', ...
    'Voyager_Read_CDF_Product.m'})).';
report.CodeSHA256 = strings(size(report.CodeFiles));
for ii = 1:numel(report.CodeFiles)
    report.CodeSHA256(ii) = string(Case1_File_SHA256(char(report.CodeFiles(ii))));
end
report.Passed = all(report.Checks.Passed);
report.PassedCount = nnz(report.Checks.Passed);
report.CheckCount = height(report.Checks);
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'supplemental_20200730');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(folder, ['single_year_source_index_', stamp, '.mat']));
save(report.AuditFile, 'report', '-v7.3');
fprintf('Single-year source-index regression: %d/%d checks passed.\n%s\n', ...
    report.PassedCount, report.CheckCount, report.AuditFile);
assert(report.Passed, 'Single-year source-index regression failed; inspect saved checks.');

    function addCheck(name, passed, detail)
        report.Checks(end+1, :) = {currentRole, string(name), logical(passed), string(detail)};
    end
end

function file = one2020File(files)
files = string(files);
selected = files(contains(files, '_20200101_'));
assert(isscalar(selected) && isfile(selected), 'Expected one original annual 2020 CDF per product.');
file = selected;
end
