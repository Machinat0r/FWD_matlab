function report = Case1_Test_Epoch_Anchoring()
%Case1_Test_Epoch_Anchoring Regression of approved original-Epoch selection.
%   Official daily/hourly CDF records are selected without averaging.
%   Within [startTime,endTime), only DeltaT < 0 is rejected by this policy.
%   NaN, zero and positive long/cross-boundary DeltaT remain eligible.
%   Numerical tolerances below test geometry; they do not screen PAD.
%   Run in an isolated MATLAB process with an empty SPICE kernel pool.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% configuration and original source records
cfg = Case1_Config;
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
Case1_Add_IRFU_Path(cfg.IRFURoot);
report = struct;
report.StartedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.TestSHA256 = Case1_File_SHA256(report.TestFile);
report.MATLABVersion = version;
report.Scope = ['Original-Epoch selection and geometry regression. ', ...
    'No magnetic-field matching, scientific averaging or figure generation.'];
report.Checks = table('Size', [0 5], ...
    'VariableTypes', {'string', 'logical', 'double', 'double', 'string'}, ...
    'VariableNames', {'Check', 'Passed', 'Measured', 'Limit', 'Detail'});
report.FatalError = "";
report.SourceFiles = table('Size', [0 7], ...
    'VariableTypes', {'string', 'string', 'string', 'double', ...
    'double', 'double', 'double'}, ...
    'VariableNames', {'Cadence', 'SourceFile', 'SHA256', 'InputRecords', ...
    'NegativeDeltaT', 'RetainedRecords', 'NegativeDeltaTInEventWindows'});
report.ExcludedRecords = table('Size', [0 5], ...
    'VariableTypes', {'string', 'double', 'datetime', 'double', 'double'}, ...
    'VariableNames', {'SourceFile', 'SourceRow', 'EpochUTC', 'DeltaT_s', ...
    'RawDeltaT_s'});
report.ExcludedRecords.EpochUTC.TimeZone = 'UTC';
report.SelectionAudits = cell(0, 1);

try
    %% synthetic boundary cases
    startTime = datetime(2020, 5, 22, 0, 0, 0, 'TimeZone', 'UTC');
    endTime = startTime + days(1);
    data = struct;
    data.Epoch = startTime + seconds([-1; 0; 17.125; 59.75; 86399.5; 123.25; 86400]);
    data.DeltaT = [3600; 0; NaN; 153*86400; 3600; -1; 3600];
    data.Flux = reshape(1:56, 7, 8);
    data.FluxUncertainty = data.Flux/100;
    data.Flux(3, 4) = NaN;
    data.FluxUncertainty(3, 4) = NaN;
    data.SourceFileIndex = ones(7, 1);
    data.SourceRecordNumber = (11:17).';
    data.SourceManifest = table("synthetic_source.cdf", 'VariableNames', {'SourceFile'});
    before = data;
    [keep, audit] = Case1_Select_LECP_Epoch(data, startTime, endTime);
    expected = [false; true; true; true; true; false; false];
    addCheck('synthetic_exact_keep_mask', isequal(keep(:), expected), ...
        sum(keep), 4, 'Start inclusive, end exclusive; negative DeltaT only.');
    addCheck('synthetic_input_not_mutated', isequaln(data, before), ...
        double(isequaln(data, before)), 1, 'All original fields are unchanged.');
    addCheck('synthetic_epoch_flux_sigma_exact', ...
        isequal(data.Epoch(keep), before.Epoch(expected)) && ...
        isequaln(data.Flux(keep, :), before.Flux(expected, :)) && ...
        isequaln(data.FluxUncertainty(keep, :), before.FluxUncertainty(expected, :)), ...
        sum(keep), 4, 'Direct source rows retained, including NaN.');
    addCheck('synthetic_special_cases_retained', all(keep(2:5)), sum(keep(2:5)), 4, ...
        'Zero, NaN, 153-day positive and midnight-crossing DeltaT retained.');
    addCheck('synthetic_negative_audit', audit.NegativeDeltaTRejected == 1 && ...
        height(audit.ExcludedRecords) == 1 && audit.ExcludedRecords.SourceRow(1) == 6 && ...
        audit.ExcludedRecords.EpochUTC(1) == data.Epoch(6) && ...
        audit.ExcludedRecords.DeltaT_s(1) == -1, audit.NegativeDeltaTRejected, 1, ...
        'Rejected source row, original Epoch and signed DeltaT preserved.');
    addCheck('synthetic_audit_counts', audit.InputRecords == 5 && ...
        audit.RetainedRecords == 4 && isfield(audit, 'Policy') && ...
        strlength(string(audit.Policy)) > 0, audit.RetainedRecords, 4, ...
        'InputRecords counts five in-window rows; four retained.');
    addCheck('synthetic_original_source_mapping', ...
        audit.ExcludedRecords.SourceCDF(1) == "synthetic_source.cdf" && ...
        audit.ExcludedRecords.SourceCDFRecord(1) == 16 && ...
        isequal(audit.KeptSourceRows, find(expected)) && ...
        audit.MissingDeltaTRetained == 1, ...
        audit.MissingDeltaTRetained, 1, 'Original file/row mapping; NaN DeltaT retained.');
    report.Synthetic = struct('Input', data, 'Keep', keep, 'Audit', audit);

    %% nine daily and nine hourly annual CDFs
    catalog = Case1_Event_Catalog;
    catalog = catalog(catalog.Spacecraft == 1, :);
    allStart = datetime(2013, 1, 1, 'TimeZone', 'UTC');
    allEnd = datetime(2022, 1, 1, 'TimeZone', 'UTC');
    files = [cfg.LECPNativeDailyCDFs(:); cfg.LECPNativeHourlyCDFs(:)];
    cadence = [repmat("day", numel(cfg.LECPNativeDailyCDFs), 1); ...
        repmat("hour", numel(cfg.LECPNativeHourlyCDFs), 1)];
    addCheck('native_file_count', numel(cfg.LECPNativeDailyCDFs) == 9 && ...
        numel(cfg.LECPNativeHourlyCDFs) == 9, numel(files), 18, ...
        '2013 through 2021 original annual files for each official cadence.');
    allMaskExact = true; allEpochExact = true; allPayloadExact = true;
    allDeltaTExact = true; allAuditExact = true;
    retainedLong = 0; sourceLong = 0;
    for ii = 1:numel(files)
        obj = dataobj(char(files(ii)));
        rawEpoch = getv(obj, 'Epoch');
        rawDeltaT = getv(obj, 'DeltaT');
        rawFlux = getv(obj, 'FHDU_SectoredFluxes');
        rawSigma = getv(obj, 'FHDU_SectoredFluxUncertainties');
        product = Voyager_Read_CDF_Product(char(files(ii)), 'lecp_sector_daily');
        originalEpoch = datetime(EpochTT(int64(rawEpoch.data(:))).epochUnix, ...
            'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
        originalDeltaT = double(rawDeltaT.data(:));
        product.SourceFileIndex = repmat(ii, numel(originalEpoch), 1);
        product.SourceRecordNumber = (1:numel(originalEpoch)).';
        [selected, selectedAudit] = Case1_Select_LECP_Epoch(product, allStart, allEnd);
        inWindow = originalEpoch >= allStart & originalEpoch < allEnd;
        rejected = inWindow & originalDeltaT < 0;
        expectedSelected = inWindow & ~(originalDeltaT < 0);
        allMaskExact = allMaskExact && isequal(selected(:), expectedSelected);
        allEpochExact = allEpochExact && ...
            isequal(product.Epoch(selected), originalEpoch(expectedSelected));
        allDeltaTExact = allDeltaTExact && isequaln(product.DeltaT(:), originalDeltaT);
        allPayloadExact = allPayloadExact && ...
            validPayloadEqual(product.FHDU_SectoredFluxes, rawFlux, selected) && ...
            validPayloadEqual(product.FHDU_SectoredFluxUncertainties, rawSigma, selected);
        allAuditExact = allAuditExact && ...
            selectedAudit.InputRecords == numel(originalEpoch) && ...
            selectedAudit.RetainedRecords == sum(expectedSelected) && ...
            selectedAudit.NegativeDeltaTRejected == sum(rejected) && ...
            isequal(selectedAudit.ExcludedRecords.SourceRow(:), find(rejected)) && ...
            isequal(selectedAudit.ExcludedRecords.EpochUTC, originalEpoch(rejected));
        nominalSeconds = 86400;
        if cadence(ii) == "hour", nominalSeconds = 3600; end
        long = inWindow & originalDeltaT > nominalSeconds;
        sourceLong = sourceLong + sum(long);
        retainedLong = retainedLong + sum(long & selected);
        inEvents = false(size(originalEpoch));
        for jj = 1:height(catalog)
            inEvents = inEvents | (originalEpoch >= catalog.StartUTC(jj)-days(7) & ...
                originalEpoch < catalog.EndUTCExclusive(jj)+days(7));
        end
        report.SourceFiles(end+1, :) = {cadence(ii), files(ii), ...
            string(Case1_File_SHA256(char(files(ii)))), numel(originalEpoch), ...
            sum(rejected), sum(selected), sum(rejected & inEvents)};
        report.SelectionAudits{end+1, 1} = selectedAudit;
        rows = find(rejected);
        for jj = 1:numel(rows)
            kk = rows(jj);
            report.ExcludedRecords(end+1, :) = {files(ii), kk, ...
                originalEpoch(kk), product.DeltaT(kk), originalDeltaT(kk)};
        end
    end
    daily = report.SourceFiles.Cadence == "day";
    hourly = report.SourceFiles.Cadence == "hour";
    addCheck('native_daily_input_count', sum(report.SourceFiles.InputRecords(daily)) == 1870, ...
        sum(report.SourceFiles.InputRecords(daily)), 1870, 'Nine annual daily files.');
    addCheck('native_daily_retained_count', sum(report.SourceFiles.RetainedRecords(daily)) == 1867, ...
        sum(report.SourceFiles.RetainedRecords(daily)), 1867, 'Three negative DeltaT removed.');
    addCheck('native_hourly_retained_count', sum(report.SourceFiles.RetainedRecords(hourly)) == 14874, ...
        sum(report.SourceFiles.RetainedRecords(hourly)), 14874, 'No hourly negative DeltaT.');
    addCheck('native_selection_mask_exact', allMaskExact, double(allMaskExact), 1, ...
        'Original Epoch interval membership and signed DeltaT only.');
    addCheck('native_original_epoch_exact', allEpochExact, double(allEpochExact), 1, ...
        'Decoded TT2000 Epoch unchanged by selection.');
    addCheck('native_valid_flux_sigma_exact', allPayloadExact, double(allPayloadExact), 1, ...
        'All retained metadata-valid flux/sigma equal raw CDF values.');
    addCheck('native_signed_DeltaT_preserved', allDeltaTExact, double(allDeltaTExact), 1, ...
        'Real negatives retain their sign before selection; no DeltaT fill in these CDFs.');
    addCheck('native_audit_original_rows_exact', allAuditExact, double(allAuditExact), 1, ...
        'Exact original row numbers and Epoch in the exclusion audit.');
    addCheck('native_all_long_positive_retained', retainedLong == sourceLong && sourceLong > 0, ...
        retainedLong, sourceLong, 'Longer-than-cadence durations never cause rejection.');
    actualDates = report.ExcludedRecords.EpochUTC;
    actualDates.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    expectedDates = ["2015-05-13 10:40:19.005"; ...
        "2020-02-05 06:09:04.130"; "2020-03-04 04:10:39.030"];
    addCheck('native_three_negative_dates', isequal(string(actualDates), expectedDates), ...
        height(report.ExcludedRecords), 3, 'Verified original TT2000 dates.');
    addCheck('native_no_negative_rows_in_event_windows', ...
        sum(report.SourceFiles.NegativeDeltaTInEventWindows) == 0, ...
        sum(report.SourceFiles.NegativeDeltaTInEventWindows), 0, ...
        'None of the 45 current V1 event windows overlaps those three records.');

    %% merged-file source mapping and exact payload conservation
    mergedExact = true; mergedAuditExact = true; duplicateCount = 0;
    groups = {cfg.LECPNativeDailyCDFs(:), cfg.LECPNativeHourlyCDFs(:)};
    for ii = 1:numel(groups)
        merged = Case1_Read_LECP_CDFs(groups{ii});
        [mergedKeep, mergedAudit] = Case1_Select_LECP_Epoch(merged, allStart, allEnd);
        duplicateCount = duplicateCount + merged.DuplicateIdenticalRecordsRemoved;
        for jj = 1:numel(groups{ii})
            source = Voyager_Read_CDF_Product(char(groups{ii}(jj)), 'lecp_sector_daily');
            selected = mergedKeep & merged.SourceFileIndex == jj;
            originalRows = merged.SourceRecordNumber(selected);
            mergedExact = mergedExact && ...
                isequal(merged.Epoch(selected), source.Epoch(originalRows)) && ...
                isequaln(merged.FHDU_SectoredFluxes(selected, :, :), ...
                source.FHDU_SectoredFluxes(originalRows, :, :)) && ...
                isequaln(merged.FHDU_SectoredFluxUncertainties(selected, :, :), ...
                source.FHDU_SectoredFluxUncertainties(originalRows, :, :));
        end
        excluded = mergedAudit.ExcludedRecords;
        for jj = 1:height(excluded)
            mergedRow = excluded.SourceRow(jj);
            mergedAuditExact = mergedAuditExact && ...
                excluded.SourceCDF(jj) == merged.SourceManifest.SourceFile( ...
                merged.SourceFileIndex(mergedRow)) && ...
                excluded.SourceCDFRecord(jj) == merged.SourceRecordNumber(mergedRow);
        end
    end
    addCheck('merged_source_epoch_flux_sigma_exact', mergedExact, double(mergedExact), 1, ...
        'Merged retained records map exactly to each original file and record.');
    addCheck('merged_exclusion_source_mapping_exact', mergedAuditExact, double(mergedAuditExact), 1, ...
        'Excluded negative records identify the original CDF and original row.');
    addCheck('merged_no_duplicate_records_removed', duplicateCount == 0, duplicateCount, 0, ...
        'The current 18 distinct annual sources contain no duplicate Epoch records.');

    %% non-rounded Epoch and unchanged signed three-dimensional geometry
    sampleEpoch = startTime + seconds([17.125; 13725.2; 43199.99; 86399.5]);
    sampleB = [1 2 3; -2 3 1; 3 -1 2; 1 -3 -2];
    cfg.PADCadence = 'day';
    dailyPointing = Case1_Predicted_LECP_Pointing(sampleEpoch, sampleB, cfg);
    cfg.PADCadence = 'hour';
    hourlyPointing = Case1_Predicted_LECP_Pointing(sampleEpoch, sampleB, cfg);
    addCheck('pointing_original_epoch_exact', ...
        isequal(dailyPointing.TimeUTC, sampleEpoch) && isequal(hourlyPointing.TimeUTC, sampleEpoch), ...
        numel(sampleEpoch), numel(sampleEpoch), 'Fractional-second times survive both cadence modes.');
    addCheck('pointing_cadence_does_not_round', ...
        isequaln(dailyPointing.ParticleRTN, hourlyPointing.ParticleRTN) && ...
        isequaln(dailyPointing.PitchAngle_deg, hourlyPointing.PitchAngle_deg), ...
        double(isequaln(dailyPointing.ParticleRTN, hourlyPointing.ParticleRTN)), 1, ...
        'Same input Epoch and B yield the same attitude and angle.');
    attitude = Case1_Read_Predicted_Attitude(sampleEpoch, cfg);
    geometry = Case1_LECP_Geometry(cfg);
    expectedParticle = nan(size(hourlyPointing.ParticleRTN));
    expectedMu = nan(size(hourlyPointing.Mu));
    for ii = 1:numel(sampleEpoch)
        u = attitude.C_SC_to_RTN(:, :, ii)*geometry.ParticleSC;
        expectedParticle(ii, :, :) = reshape(u.', 1, 8, 3);
        expectedMu(ii, :) = (sampleB(ii, :)/norm(sampleB(ii, :)))*u;
    end
    matrixTolerance = 1e-12;
    angleTolerance = 1e-10;
    vectorError = max(abs(hourlyPointing.ParticleRTN-expectedParticle), [], 'all');
    unitError = max(abs(sqrt(sum(hourlyPointing.ParticleRTN.^2, 3))-1), [], 'all');
    signError = max(abs(hourlyPointing.ParticleRTN+hourlyPointing.LookRTN), [], 'all');
    muError = max(abs(hourlyPointing.Mu-expectedMu), [], 'all');
    angleError = max(abs(hourlyPointing.PitchAngle_deg-acosd(expectedMu)), [], 'all');
    addCheck('pointing_original_epoch_rotation', isfinite(vectorError) && vectorError < matrixTolerance, ...
        vectorError, matrixTolerance, 'Independent direct-reader evaluation at identical Epoch.');
    addCheck('pointing_unit_directions', isfinite(unitError) && unitError < matrixTolerance, ...
        unitError, matrixTolerance, 'All eight directions remain unit vectors.');
    addCheck('pointing_particle_negative_look', isfinite(signError) && signError < matrixTolerance, ...
        signError, matrixTolerance, 'Approved inward particle/outward aperture sign.');
    addCheck('pointing_full_3D_dot_product', isfinite(muError) && muError < matrixTolerance, ...
        muError, matrixTolerance, 'All R/T/N components contribute to mu.');
    addCheck('pointing_angle_formula', isfinite(angleError) && angleError < angleTolerance, ...
        angleError, angleTolerance, 'Alpha=acosd(mu); no angle binning or averaging.');
    normalComponent = max(abs(hourlyPointing.ParticleRTN(:, :, 3)), [], 'all');
    addCheck('pointing_normal_component_not_zeroed', normalComponent > matrixTolerance, ...
        normalComponent, matrixTolerance, 'Numerical regression only, never a science cut.');
    report.Pointing = hourlyPointing;
    report.PointingInput = struct('EpochUTC', sampleEpoch, 'B_RTN', sampleB);
    report.NumericalTolerances = struct('Matrix', matrixTolerance, 'Angle_deg', angleTolerance);
catch ME
    report.FatalError = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    addCheck('execution_completed', false, NaN, NaN, string(ME.message));
end

%% classified reproducible audit, no scientific figure or derived PAD file
report.CompletedUTC = datetime('now', 'TimeZone', 'UTC');
report.Passed = all(report.Checks.Passed) && strlength(report.FatalError) == 0;
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'epoch_anchor');
if ~isfolder(folder), mkdir(folder); end
report.AuditFile = fullfile(folder, ['epoch_anchor_regression_', ...
    char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS')), '.mat']);
save(report.AuditFile, 'report');
fprintf('Epoch-anchor regression: passed=%d, checks=%d/%d\n', ...
    report.Passed, sum(report.Checks.Passed), height(report.Checks));
fprintf('Test audit: %s\n', report.AuditFile);
disp(report.Checks);

    function addCheck(name, passed, measured, limit, detail)
        report.Checks(end+1, :) = {string(name), logical(passed), ...
            double(measured), double(limit), string(detail)};
    end
end

function same = validPayloadEqual(value, raw, keep)
% Compare valid physical values exactly; preserve original missing values.
original = double(raw.data(keep, :, :));
selected = value(keep, :, :);
valid = isfinite(original) & abs(original) < 1e30;
if isfield(raw, 'FILLVAL'), valid = valid & original ~= double(raw.FILLVAL); end
if isfield(raw, 'VALIDMIN'), valid = valid & original >= double(raw.VALIDMIN); end
if isfield(raw, 'VALIDMAX'), valid = valid & original <= double(raw.VALIDMAX); end
same = isequal(selected(valid), original(valid)) && all(isnan(selected(~valid)));
end
