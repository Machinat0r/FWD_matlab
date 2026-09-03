function report = Case1_Test_Peak_Plot_Audits()
%Case1_Test_Peak_Plot_Audits Independently check all 45 five-panel PADs.
%   Run after the final hourly overview and peak PNG/MAT files are complete.
%   Peak and neighbors are reconstructed from raw official hourly CDF rows.
%   No renderer or its selection routine is called to create the reference.
%   The test reads sources once and writes only a classified validation MAT.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% sources and audit containers
cfg = Case1_Config;
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
Case1_Add_IRFU_Path(cfg.IRFURoot);
catalog = Case1_Event_Catalog;
catalog = catalog(catalog.Spacecraft == 1, :);
overviewFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1h', ...
    'derived', 'pitch_angle', '2013-2021', 'predicted_ck');
currentEvent = "";
report = struct;
report.StartedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.TestSHA256 = Case1_File_SHA256(report.TestFile);
report.MATLABVersion = version;
report.Scope = ['Independent raw-CDF selection, normalization, provenance ', ...
    'and artifact regression; no new scientific processing or figures.'];
report.Checks = table('Size', [0 6], ...
    'VariableTypes', {'string', 'string', 'logical', 'double', 'double', 'string'}, ...
    'VariableNames', {'EventID', 'Check', 'Passed', 'Measured', 'Limit', 'Detail'});
report.Events = table('Size', [0 12], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'string', ...
    'string', 'double', 'double', 'logical', 'double', 'logical', 'string'}, ...
    'VariableNames', {'EventID', 'PeakMAT', 'PeakSHA256', 'FigurePNG', ...
    'FigureSHA256', 'OverviewSHA256', 'InputRows', 'SelectedCount', ...
    'PeakPADUsable', 'PlottedCount', 'Passed', 'Error'});
report.References = cell(height(catalog), 1);
report.FatalError = "";

try
    source = Case1_Read_LECP_CDFs(cfg.LECPNativeHourlyCDFs);
    report.SourceManifest = source.SourceManifest;
    report.SourceMetadata = source.SourceMetadata;
    allFlux = reshape(source.FHDU_SectoredFluxes(:, 10, 1:7), numel(source.Epoch), 7);
    allSigma = reshape(source.FHDU_SectoredFluxUncertainties(:, 10, 1:7), numel(source.Epoch), 7);
    allBand = reshape(source.FHDU_EnergyRange(:, :, 10), numel(source.Epoch), 2);
    sourceCodeHash = containers.Map('KeyType', 'char', 'ValueType', 'char');
    mats = dir(fullfile(cfg.PeakPADDataFolder, 'V1_*_P1_peak5_hourly.mat'));
    pngs = dir(fullfile(cfg.PeakPADFolder, 'V1_*_P1_peak5_hourly.png'));
    addCheck('45_v1_events', height(catalog) == 45, height(catalog), 45, 'Fixed event catalog.');
    addCheck('45_peak_mat_files', numel(mats) == 45, numel(mats), 45, 'One MAT per event.');
    addCheck('45_peak_png_files', numel(pngs) == 45, numel(pngs), 45, 'One PNG per event.');

    for ie = 1:height(catalog)
        currentEvent = string(catalog.EventID(ie));
        firstCheck = height(report.Checks) + 1;
        peakFile = ""; peakSHA = ""; figureFile = ""; figureSHA = "";
        overviewSHA = ""; eventError = "";
        nInput = NaN; nSelected = NaN; nPlotted = NaN; peakUsable = false;
        try
            stamp = char(datetime(catalog.StartUTC(ie), 'Format', 'yyyyMMdd'));
            name = ['V1_', char(currentEvent), '_', stamp, '_P1_peak5_hourly'];
            peakFile = string(fullfile(cfg.PeakPADDataFolder, [name, '.mat']));
            figureFile = string(fullfile(cfg.PeakPADFolder, [name, '.png']));
            addCheck('event_artifacts_exist', isfile(peakFile) && isfile(figureFile), ...
                double(isfile(peakFile))+double(isfile(figureFile)), 2, name);
            saved = load(char(peakFile));
            A = saved.peakAudit;
            peakSHA = string(Case1_File_SHA256(char(peakFile)));
            figureSHA = string(Case1_File_SHA256(char(figureFile)));
            info = imfinfo(char(figureFile));
            addCheck('valid_landscape_png', strcmpi(info.Format, 'png') && ...
                info.Width > info.Height && info.Height > 0 && A.FigureCreated && ...
                string(A.FigureFile) == figureFile, info.Width, info.Height, ...
                'PNG decodes, is landscape, and is linked to the audit.');
            overviewName = ['V1_', char(currentEvent), '_', stamp, '_', stamp, ...
                '_COHO1h_raw_LECP_P1_pitch_angle_predictedCK_1h_nativeCDF_Epoch.mat'];
            overviewFile = fullfile(overviewFolder, overviewName);
            O = load(overviewFile);
            T = O.pitchAngleTable;
            overviewSHA = string(Case1_File_SHA256(overviewFile));
            nInput = height(T);
            startTime = catalog.StartUTC(ie)-days(cfg.ContextDays);
            endTime = catalog.EndUTCExclusive(ie)+days(cfg.ContextDays);
            sourceInWindow = source.Epoch >= startTime & source.Epoch < endTime;
            sourceRows = find(sourceInWindow & ~(source.DeltaT(:) < 0));
            negativeRows = find(sourceInWindow & source.DeltaT(:) < 0);
            addCheck('window_and_input_identity', string(A.EventID) == currentEvent && ...
                A.WindowStartUTC == startTime && A.WindowEndUTCExclusive == endTime && ...
                A.Spacecraft == 1 && string(A.PADCadence) == "hour" && ...
                A.InputTableRowCount == numel(sourceRows) && height(T) == numel(sourceRows), ...
                A.InputTableRowCount, numel(sourceRows), 'Exact approved hourly event window.');
            addCheck('original_epoch_selection_provenance', ...
                isequaln(A.RecordSelection, O.recordSelection) && ...
                isequal(A.RecordSelection.KeptSourceRows(:), sourceRows) && ...
                isequal(A.RecordSelection.ExcludedRecords.SourceRow(:), negativeRows) && ...
                string(A.RecordSelection.Policy) == "epoch_drop_negative_deltat", ...
                numel(A.RecordSelection.KeptSourceRows), numel(sourceRows), ...
                'Original Epoch selection; negative DeltaT only.');
            addCheck('source_hashes_and_context_retained', ...
                isequaln(A.SourceLECP, source.SourceManifest) && ...
                isequaln(A.SourceLECP, O.sourceLECP) && ...
                isequaln(A.SourceMAG, O.sourceMAG) && ...
                isequaln(A.PointingAudit, T.Properties.UserData) && ...
                isequaln(A.SourceTableUserData, T.Properties.UserData) && ...
                string(A.PitchAngleAuditFile) == string(overviewFile), ...
                height(A.SourceLECP), height(source.SourceManifest), ...
                'Native CDF SHA256, MAG paths, and full attitude/overview context.');
            codeMatches = numel(A.CodeFiles) == numel(A.CodeSHA256) && ~isempty(A.CodeFiles);
            for ii = 1:numel(A.CodeFiles)
                file = char(A.CodeFiles(ii));
                if ~isKey(sourceCodeHash, file)
                    sourceCodeHash(file) = Case1_File_SHA256(file);
                end
                codeMatches = codeMatches && string(A.CodeSHA256(ii)) == string(sourceCodeHash(file));
            end
            addCheck('production_code_hashes_current', codeMatches, ...
                numel(A.CodeSHA256), numel(A.CodeFiles), 'Recorded production code matches final files.');
            addCheck('approved_display_energy_only', ...
                isequal(A.DisplayEnergyMeV, [0.57 1.78]) && ...
                isequal(saved.opts.P1DisplayEnergyMeV, [0.57 1.78]) && ...
                isfield(A, 'EnergyLabelPolicy') && strlength(string(A.EnergyLabelPolicy)) > 0, ...
                A.DisplayEnergyMeV(2), 1.78, 'Display convention retained separately from raw CDF bounds.');

            rawFlux = allFlux(sourceRows, :);
            rawSigma = allSigma(sourceRows, :);
            rawFlux(~isfinite(rawFlux) | rawFlux < 0) = NaN;
            rawSigma(~isfinite(rawSigma) | rawSigma < 0) = NaN;
            time = source.Epoch(sourceRows);
            pa = nan(nInput, 7);
            complete = false(nInput, 1);
            if ~isempty(T)
                actualRaw = sectorArray(T, 'RawFlux_S%d_1h');
                actualSigma = sectorArray(T, 'FluxUncertainty_S%d_1h');
                pa = sectorArray(T, 'PA_S%d_deg');
                complete = all(isfinite(rawFlux) & rawFlux > 0 & isfinite(pa), 2);
                addCheck('overview_raw_source_payload_exact', ...
                    isequal(T.EpochUTC, time) && ...
                    isequal(T.SourceRow(:), sourceRows) && ...
                    isequal(T.SourceCDFRecord(:), source.SourceRecordNumber(sourceRows)) && ...
                    isequal(string(T.SourceCDF), source.SourceManifest.SourceFile( ...
                    source.SourceFileIndex(sourceRows))) && ...
                    isequaln(actualRaw, rawFlux) && isequaln(actualSigma, rawSigma), ...
                    nInput, numel(sourceRows), 'Overview source rows, raw flux and sigma unchanged.');
                addCheck('overview_original_energy_and_display_distinct', ...
                    isequaln([T.P1EnergyLower_MeV T.P1EnergyUpper_MeV], allBand(sourceRows, :)) && ...
                    all(T.P1DisplayEnergyLower_MeV == 0.57) && ...
                    all(T.P1DisplayEnergyUpper_MeV == 1.78), ...
                    max(T.P1EnergyUpper_MeV), 0.89, 'Source energy remains the original CDF value.');
                addCheck('overview_no_new_validity_gate_or_background', ...
                    isequal(T.PADUsable, complete) && ~any(T.BackgroundCorrectionApplied) && ...
                    all(T.SourceToDifferentialFluxFactor == 1), nnz(T.PADUsable), nnz(complete), ...
                    'S1-S7 finite positive flux and finite PA; no recalibration.');
            end

            %% independent global maximum and chronological neighbor search
            candidate = any(isfinite(rawFlux) & rawFlux > 0, 2);
            candidateRows = find(candidate);
            expectedPeak = NaN; expectedPeakFlux = NaN;
            expectedPeakTime = NaT(1, 1, 'TimeZone', 'UTC');
            selected = nan(5, 1);
            rowMaximum = nan(nInput, 1);
            for ii = 1:nInput
                values = rawFlux(ii, isfinite(rawFlux(ii, :)) & rawFlux(ii, :) > 0);
                if ~isempty(values), rowMaximum(ii) = max(values); end
            end
            if any(candidate)
                positiveValues = rawFlux(isfinite(rawFlux) & rawFlux > 0);
                expectedPeakFlux = max(positiveValues);
                [tieRows, ~] = find(rawFlux == expectedPeakFlux);
                tieRows = unique(tieRows);
                [~, tieOrder] = sort(time(tieRows));
                expectedPeak = tieRows(tieOrder(1));
                expectedPeakTime = time(expectedPeak);
                selected(3) = expectedPeak;
                eligible = find(complete);
                [~, chronological] = sort(time(eligible));
                eligible = eligible(chronological);
                before = eligible(time(eligible) < expectedPeakTime);
                after = eligible(time(eligible) > expectedPeakTime);
                if ~isempty(before)
                    chosenBefore = before(max(1, numel(before)-1):end);
                    selected(3-numel(chosenBefore):2) = chosenBefore;
                end
                if ~isempty(after)
                    chosenAfter = after(1:min(2, numel(after)));
                    selected(4:3+numel(chosenAfter)) = chosenAfter;
                end
            end
            addCheck('raw_peak_ignores_S8_and_PAD_gate', ...
                isequaln(A.PeakTableRow, expectedPeak) && ...
                isequaln(A.PeakFlux, expectedPeakFlux) && ...
                isequaln(A.PeakEpochUTC, expectedPeakTime) && ...
                isequal(A.CandidateTableRows(:), candidateRows(:)) && ...
                isequaln(A.RowMaximumSectorFlux(:), rowMaximum), ...
                A.PeakTableRow, expectedPeak, 'Direct maximum of positive S1-S7 CDF values; earliest tie.');
            addCheck('nearest_two_neighbors_with_missing_skipped', ...
                isequaln(A.SelectedTableRows(:), selected) && ...
                isequal(A.PanelOffsets(:), (-2:2).'), ...
                nnz(isfinite(A.SelectedTableRows)), nnz(isfinite(selected)), ...
                'Chronological complete-PAD search on each side, within the window.');

            %% direct source values and fixed-denominator normalization
            expectedTime = NaT(5, 1, 'TimeZone', 'UTC');
            expectedComplete = false(5, 1);
            expectedFlux = nan(5, 7); expectedSigma = expectedFlux; expectedPA = expectedFlux;
            expectedNormalizedFlux = expectedFlux; expectedNormalizedSigma = expectedFlux;
            expectedScale = nan(5, 1); expectedBand = nan(5, 2);
            normalizers = cell(5, 1); rowsMatch = true;
            for ip = 1:5
                row = selected(ip);
                if ~isfinite(row)
                    rowsMatch = rowsMatch && isempty(A.SelectedRows{ip});
                    continue
                end
                expectedTime(ip) = time(row);
                expectedComplete(ip) = complete(row);
                expectedFlux(ip, :) = rawFlux(row, :);
                expectedSigma(ip, :) = rawSigma(row, :);
                expectedPA(ip, :) = pa(row, :);
                expectedScale(ip) = rowMaximum(row);
                normalizers{ip} = find(rawFlux(row, :) == rowMaximum(row));
                expectedNormalizedFlux(ip, :) = rawFlux(row, :)/rowMaximum(row);
                expectedNormalizedSigma(ip, :) = rawSigma(row, :)/rowMaximum(row);
                expectedBand(ip, :) = allBand(sourceRows(row), :);
                originalRow = T(row, :);
                originalRow.Properties.UserData = [];
                rowsMatch = rowsMatch && isequaln(A.SelectedRows{ip}, originalRow);
            end
            peakUsable = isfinite(expectedPeak) && complete(expectedPeak);
            nSelected = nnz(isfinite(selected)); nPlotted = nnz(expectedComplete);
            addCheck('selected_rows_and_original_epochs_exact', rowsMatch && ...
                isequaln(A.SelectedEpochUTC, expectedTime) && ...
                isequal(A.SelectedPADUsable, expectedComplete), nSelected, A.SelectedCount, ...
                'Selected row copies retain every overview field and original time.');
            addCheck('seven_unmerged_source_points_exact', isequal(A.Sectors, 1:7) && ...
                isequaln(A.RawFlux, expectedFlux) && isequaln(A.RawSigma, expectedSigma) && ...
                isequaln(A.PA_deg, expectedPA) && isequal(size(A.RawFlux), [5 7]), ...
                size(A.RawFlux, 2), 7, 'All seven original sectors; no angle merge or fit.');
            addCheck('own_max_fixed_scale_normalization_exact', ...
                isequaln(A.NormalizationFlux, expectedScale) && ...
                isequaln(A.NormalizingSectors, normalizers) && ...
                isequaln(A.NormalizedFlux, expectedNormalizedFlux) && ...
                isequaln(A.NormalizedSigma, expectedNormalizedSigma), ...
                nSelected, A.SelectedCount, 'J/Jmax and sigma/Jmax use each row own fixed denominator.');
            addCheck('selected_original_CDF_energy_retained', ...
                isequaln(A.SourceEnergyMeV, expectedBand), ...
                nnz(isfinite(A.SourceEnergyMeV)), nnz(isfinite(expectedBand)), ...
                'Original 0.57-0.89 metadata retained, independently of the display label.');
            expectedStatus = "no_retained_records";
            if nInput > 0, expectedStatus = "no_positive_sector_flux"; end
            if isfinite(expectedPeak)
                if nPlotted == 5
                    expectedStatus = "complete";
                elseif ~peakUsable
                    expectedStatus = "peak_pad_unavailable";
                else
                    expectedStatus = "partial_neighbors";
                end
            end
            addCheck('peak_and_blank_panel_status_correct', A.SelectedCount == nSelected && ...
                A.PlottedCount == nPlotted && A.PeakPADUsable == peakUsable && ...
                string(A.Status) == expectedStatus, A.PlottedCount, nPlotted, ...
                'An incomplete peak remains central; unavailable neighbors stay blank.');
            report.References{ie} = struct('EventID', currentEvent, 'SourceRows', sourceRows, ...
                'PeakTableRow', expectedPeak, 'SelectedTableRows', selected, ...
                'SelectedEpochUTC', expectedTime, 'SelectedPADUsable', expectedComplete, ...
                'RawFlux', expectedFlux, 'RawSigma', expectedSigma, 'PA_deg', expectedPA, ...
                'NormalizedFlux', expectedNormalizedFlux, 'NormalizedSigma', expectedNormalizedSigma);
        catch ME
            eventError = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
            addCheck('event_test_completed', false, NaN, NaN, string(ME.message));
        end
        eventPassed = all(report.Checks.Passed(firstCheck:end)) && strlength(eventError) == 0;
        report.Events(end+1, :) = {currentEvent, peakFile, peakSHA, figureFile, figureSHA, ...
            overviewSHA, nInput, nSelected, peakUsable, nPlotted, eventPassed, eventError};
    end
catch ME
    report.FatalError = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    addCheck('test_completed', false, NaN, NaN, string(ME.message));
end

%% classified verification output only
report.CompletedUTC = datetime('now', 'TimeZone', 'UTC');
report.Passed = all(report.Checks.Passed) && all(report.Events.Passed) && ...
    height(report.Events) == 45 && strlength(report.FatalError) == 0;
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'peak_pad');
if ~isfolder(folder), mkdir(folder); end
report.AuditFile = fullfile(folder, ['peak_plot_audits_regression_', ...
    char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS')), '.mat']);
save(report.AuditFile, 'report');
fprintf('Peak PAD integration regression: passed=%d; events=%d/45; checks=%d/%d\n', ...
    report.Passed, nnz(report.Events.Passed), nnz(report.Checks.Passed), height(report.Checks));
fprintf('Test audit: %s\n', report.AuditFile);
if any(~report.Checks.Passed), disp(report.Checks(~report.Checks.Passed, :)); end
if any(~report.Events.Passed), disp(report.Events(~report.Events.Passed, :)); end

    function addCheck(name, passed, measured, limit, detail)
        report.Checks(end+1, :) = {currentEvent, string(name), logical(passed), ...
            double(measured), double(limit), string(detail)};
    end
end

function value = sectorArray(T, pattern)
value = nan(height(T), 7);
for ii = 1:7, value(:, ii) = T.(sprintf(pattern, ii)); end
end
