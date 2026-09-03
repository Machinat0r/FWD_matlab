function report = Case1_QA_Context3D_Artifacts(dayRunFile, hourRunFile, maxEvents)
%Case1_QA_Context3D_Artifacts Read-only artifact and missing-PAD audit.
%   Supply the two completed Run_Case1_Context3Days_L1_First result MATs.
%   Production programs, plots and source data are never changed here.
%   PNGs are decoded and each event window is checked against this batch.
%   Missing peak causes use saved flux, B, attitude and PA independently.

%% the two explicit completed batches
cfg = Case1_Config;
if nargin < 3, maxEvents = Inf; end
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.Scope = "Read-only exact-report artifact inventory, +/-3-day windows, L1-first settings, and missing-center attribution; no replot or science-data mutation.";
report.ExpectedEventsPerCadence = min(45,maxEvents);
report.FullInventoryRequested = maxEvents >= 45;
report.TimestampPolicy = "Modification times are recorded exactly. PNG freshness is checked; MAT freshness is diagnostic because an overwritten HDF5 MAT can retain an earlier observed last-write timestamp. Saved options, data, window and hashes establish its current content.";
report.QAFile = string([mfilename('fullpath'), '.m']);
report.QASHA256 = string(Case1_File_SHA256(char(report.QAFile)));
report.RunFiles = string({dayRunFile; hourRunFile});
report.RunSHA256 = strings(2, 1);
report.Checks = table('Size', [0 5], ...
    'VariableTypes', {'string', 'string', 'string', 'logical', 'string'}, ...
    'VariableNames', {'Cadence', 'EventID', 'Check', 'Passed', 'Detail'});
report.Artifacts = table;
report.Events = table;
report.Peaks = table;
report.PeakCentralRows = cell(0, 1);
report.Errors = strings(0, 1);
currentCadence = ""; currentEvent = "";
expectedCadence = ["day"; "hour"];

for iRun = 1:2
    currentCadence = expectedCadence(iRun);
    resultFile = char(report.RunFiles(iRun));
    runSaved = load(resultFile, 'result');
    result = runSaved.result;
    R = result.Run.ReportV1;
    report.RunSHA256(iRun) = string(Case1_File_SHA256(resultFile));
    addCheck('45_report_events', height(R) == 45 && numel(unique(string(R.EventID))) == 45, ...
        sprintf('%d rows, %d unique EventID values',height(R),numel(unique(string(R.EventID)))));
    addCheck('run_scope', result.Cadence == currentCadence && result.ContextDays == 3 && ...
        string(result.SourcePriority) == "l1_first", 'Explicit current-batch cadence, three-day context and L1-first.');
    addCheck('no_report_errors', ~any(contains(string(R.Status), 'error')), ...
        char(join(unique(string(R.Status)), ', ')));
    if iRun == 1
        dayIDs = sort(string(R.EventID));
    else
        addCheck('same_45_event_ids', isequal(dayIDs, sort(string(R.EventID))), ...
            'Daily and hourly reports target the same 45 events.');
    end

    for iEvent = 1:min(height(R),maxEvents)
        currentEvent = string(R.EventID(iEvent));
        fprintf('QA %s %d/%d %s\n',currentCadence,iEvent,min(height(R),maxEvents),currentEvent);
        try
            pngFile = string(R.FigureFile(iEvent));
            matFile = string(R.PitchAngleAuditFile(iEvent));
            appendArtifact('overview_png', pngFile, result.CreatedUTC);
            appendArtifact('overview_mat', matFile, result.CreatedUTC);
            O = load(char(matFile), 'pitchAngleTable', 'opts');
            T = O.pitchAngleTable;
            expectedStart = dateshift(R.EventStartUTC(iEvent), 'start', 'day')-days(3);
            expectedEnd = expectedStart+days(7);
            addCheck('exact_seven_day_window', R.PlotStartUTC(iEvent) == expectedStart && ...
                R.PlotEndUTCExclusive(iEvent) == expectedEnd, ...
                sprintf('%s through %s, end exclusive', char(expectedStart), char(expectedEnd)));
            addCheck('saved_options_current', O.opts.ContextDays == 3 && ...
                string(O.opts.LECPSourcePriority) == "l1_first" && ...
                string(O.opts.PADCadence) == currentCadence, ...
                'Saved production options match this batch, independently of old folder names.');
            inWindow = true(height(T), 1);
            if ~isempty(T), inWindow = T.EpochUTC >= expectedStart & T.EpochUTC < expectedEnd; end
            addCheck('all_saved_rows_in_window', all(inWindow), sprintf('%d saved PAD input rows', height(T)));
            [fluxOK, paOK, bOK, attitudeFound, directionOK] = validity(T, currentCadence);
            gate = fluxOK & paOK;
            savedUsable = false(height(T), 1);
            sourceProduct = strings(height(T), 1);
            sourcePriority = strings(height(T), 1);
            if ~isempty(T)
                savedUsable = logical(T.PADUsable);
                sourceProduct = string(T.SourceProduct);
                sourcePriority = string(T.SourcePriority);
            end
            addCheck('only_seven_sector_gate', isequal(savedUsable, gate), ...
                'Finite positive S1-S7 flux and finite PA only; uncertainty is not a gate.');
            addCheck('source_priority_per_row', all(sourcePriority == "l1_first"), ...
                'Each saved input row identifies the current product-priority policy.');

            E = struct;
            E.Cadence = currentCadence; E.EventID = currentEvent;
            E.WindowStartUTC = expectedStart; E.WindowEndUTCExclusive = expectedEnd;
            E.WindowDays = days(expectedEnd-expectedStart);
            E.InputRows = height(T); E.PADUsableRows = sum(gate);
            E.L1Rows = sum(sourceProduct == "L1_UTC_mean");
            E.L2Rows = sum(sourceProduct == "L2");
            E.MissingFluxRows = sum(~fluxOK); E.MissingPARows = sum(~paOK);
            E.InvalidMagneticRows = sum(~bOK); E.MissingAttitudeRows = sum(~attitudeFound);
            E.MissingDirectionRows = sum(~directionOK);
            E.ReportStatus = string(R.Status(iEvent));
            E.FigurePNG = pngFile; E.OverviewMAT = matFile;
            report.Events = appendRow(report.Events, E);

            if currentCadence ~= "hour", continue; end
            peakPNG = string(R.PeakPADFigureFile(iEvent));
            peakMAT = string(R.PeakPADAuditFile(iEvent));
            appendArtifact('peak_png', peakPNG, result.CreatedUTC);
            appendArtifact('peak_mat', peakMAT, result.CreatedUTC);
            P = load(char(peakMAT), 'peakAudit', 'opts');
            A = P.peakAudit;
            addCheck('peak_window_current', A.WindowStartUTC == expectedStart && ...
                A.WindowEndUTCExclusive == expectedEnd && P.opts.ContextDays == 3 && ...
                string(P.opts.LECPSourcePriority) == "l1_first", ...
                'Peak search and all neighbors belong to this seven-day L1-first run.');
            addCheck('peak_manifest_matches', string(R.PeakPADStatus(iEvent)) == string(A.Status) && ...
                R.PeakPADSelectedCount(iEvent) == A.SelectedCount && ...
                R.PeakPADPlottedCount(iEvent) == A.PlottedCount, ...
                'Run report and actual peak MAT agree on status and counts.');
            selected = A.SelectedTableRows;
            use = isfinite(selected);
            selectedInside = all(selected(use) >= 1 & selected(use) <= height(T));
            if selectedInside
                selectedInside = all(inWindow(selected(use)));
            end
            addCheck('selected_neighbors_in_window', selectedInside, 'No selected index lies outside the saved new-window input.');
            expectedPlotted = false(5, 1);
            expectedPlotted(use) = gate(selected(use));
            addCheck('selected_actual_PAD_availability', ...
                isequal(logical(A.SelectedPADUsable), expectedPlotted) && ...
                A.PlottedCount == sum(expectedPlotted), 'Plotted count is reconstructed from saved seven-sector values.');

            F = fluxValues(T, currentCadence);
            positive = F; positive(~isfinite(positive) | positive <= 0) = NaN;
            metric = max(positive, [], 2, 'omitnan');
            candidates = find(isfinite(metric));
            expectedPeak = NaN;
            if ~isempty(candidates)
                ties = candidates(metric(candidates) == max(metric(candidates)));
                [~, earliest] = min(T.EpochUTC(ties));
                expectedPeak = ties(earliest);
            end
            addCheck('peak_reselected_in_new_window', isequaln(A.PeakTableRow, expectedPeak), ...
                'Peak independently reconstructed from this saved current-window L1-first sequence.');

            Q = struct;
            Q.EventID = currentEvent; Q.Status = string(A.Status);
            Q.PeakEpochUTC = A.PeakEpochUTC; Q.PeakTableRow = A.PeakTableRow;
            Q.SelectedCount = A.SelectedCount; Q.PlottedCount = A.PlottedCount;
            Q.CentralPADUsable = false; Q.CentralFluxComplete = false;
            Q.CentralPAComplete = false; Q.CentralBValid = false;
            Q.CentralAttitudeFound = false; Q.CentralDirectionComplete = false;
            Q.CentralMAGSampleCount = NaN; Q.CentralBMeanMagnitude_nT = NaN;
            Q.CentralSourceProduct = "none"; Q.CentralPADQuality = "no peak";
            Q.CentralMissingFluxSectorCount = NaN; Q.CentralMissingPASectorCount = NaN;
            Q.CentralCause = "no_peak"; Q.Classification = "no_peak";
            Q.SkippedMissingBeforeCount = A.SkippedMissingBeforeCount;
            Q.SkippedMissingAfterCount = A.SkippedMissingAfterCount;
            Q.PeakPNG = peakPNG; Q.PeakMAT = peakMAT;
            center = table;
            if isfinite(A.PeakTableRow)
                row = A.PeakTableRow; center = T(row, :);
                Q.CentralPADUsable = gate(row); Q.CentralFluxComplete = fluxOK(row);
                Q.CentralPAComplete = paOK(row); Q.CentralBValid = bOK(row);
                Q.CentralAttitudeFound = attitudeFound(row);
                Q.CentralDirectionComplete = directionOK(row);
                Q.CentralMAGSampleCount = T.MAGVectorSampleCount(row);
                Q.CentralBMeanMagnitude_nT = T.BMeanMagnitude_nT(row);
                Q.CentralSourceProduct = string(T.SourceProduct(row));
                Q.CentralPADQuality = string(T.PADQuality(row));
                Q.CentralMissingFluxSectorCount = sum(~isfinite(F(row, :)) | F(row, :) <= 0);
                pa = T{row, cellstr(compose('PA_S%d_deg', 1:7))};
                Q.CentralMissingPASectorCount = sum(~isfinite(pa));
                causes = strings(0, 1);
                if ~fluxOK(row), causes(end+1) = "flux_missing_or_nonpositive"; end
                if ~bOK(row), causes(end+1) = "magnetic_vector_missing_or_zero"; end
                if ~attitudeFound(row), causes(end+1) = "predicted_attitude_unavailable"; end
                if attitudeFound(row) && ~directionOK(row), causes(end+1) = "sector_direction_unavailable"; end
                if ~paOK(row) && bOK(row) && directionOK(row)
                    causes(end+1) = "PA_missing_despite_finite_field_and_direction";
                end
                if gate(row)
                    Q.CentralCause = "usable";
                    if A.PlottedCount == 5, Q.Classification = "complete_five";
                    else, Q.Classification = "partial_neighbors"; end
                else
                    Q.CentralCause = join(causes, '+');
                    Q.Classification = "center_missing";
                    addCheck('missing_center_cause_identified', ~isempty(causes), ...
                        char(Q.CentralCause));
                end
            end
            center.Properties.UserData = [];
            report.PeakCentralRows{end+1, 1} = center;
            report.Peaks = appendRow(report.Peaks, Q);
        catch ME
            report.Errors(end+1, 1) = currentCadence+" "+currentEvent+": "+string(ME.message);
            fprintf(2,'QA EVENT ERROR: %s\n',report.Errors(end));
            addCheck('event_QA_completed', false, ME.message);
        end
    end
end

%% complete current-run counts and categorized missing-center results
currentCadence = "all"; currentEvent = "all";
nExpected = report.ExpectedEventsPerCadence;
addCheck('requested_overview_records', height(report.Events) == 2*nExpected, sprintf('%d', height(report.Events)));
addCheck('requested_peak_records', height(report.Peaks) == nExpected, sprintf('%d', height(report.Peaks)));
addCheck('requested_pngs', sum(endsWith(report.Artifacts.Kind, '_png')) == 3*nExpected, ...
    sprintf('%d', sum(endsWith(report.Artifacts.Kind, '_png'))));
addCheck('requested_mat_artifacts', sum(endsWith(report.Artifacts.Kind, '_mat')) == 3*nExpected, ...
    sprintf('%d', sum(endsWith(report.Artifacts.Kind, '_mat'))));
addCheck('unique_artifact_paths', numel(unique(report.Artifacts.Path)) == height(report.Artifacts), ...
    'Each production PNG/MAT appears once in the exact report-backed inventory.');
report.Summary = struct;
report.Summary.OverviewEvents = height(report.Events);
report.Summary.PeakEvents = height(report.Peaks);
report.Summary.CompleteFive = NaN; report.Summary.PartialNeighbors = NaN;
report.Summary.CenterMissing = NaN; report.Summary.NoPeak = NaN;
report.Summary.CenterMagneticCause = NaN; report.Summary.CenterFluxCause = NaN;
report.Summary.CenterAttitudeCause = NaN;
if ~isempty(report.Peaks)
    report.Summary.CompleteFive = sum(report.Peaks.Classification == "complete_five");
    report.Summary.PartialNeighbors = sum(report.Peaks.Classification == "partial_neighbors");
    report.Summary.CenterMissing = sum(report.Peaks.Classification == "center_missing");
    report.Summary.NoPeak = sum(report.Peaks.Classification == "no_peak");
    report.Summary.CenterMagneticCause = sum(contains(report.Peaks.CentralCause, 'magnetic_vector'));
    report.Summary.CenterFluxCause = sum(contains(report.Peaks.CentralCause, 'flux_missing'));
    report.Summary.CenterAttitudeCause = sum(contains(report.Peaks.CentralCause, 'attitude_unavailable'));
end
report.Summary.AllArtifactsReadable = all(report.Artifacts.Readable);
report.Summary.AllArtifactsFresh = all(report.Artifacts.ModifiedWithinRun);
report.Summary.AllPNGsFresh = all(report.Artifacts.ModifiedWithinRun(endsWith(report.Artifacts.Kind,'_png')));
report.Summary.MATsWithObservedOlderTimestamp = sum(~report.Artifacts.ModifiedWithinRun & endsWith(report.Artifacts.Kind,'_mat'));
report.PassedCount = sum(report.Checks.Passed);
report.CheckCount = height(report.Checks);
report.Passed = all(report.Checks.Passed);
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'context3d_l1_first');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(folder, ['artifact_missing_PAD_QA_', stamp, '.mat']));
save(report.AuditFile, 'report', '-v7.3');
disp(report.Summary);
fprintf('Artifact/missing-PAD QA: %d/%d checks passed.\n%s\n', ...
    report.PassedCount, report.CheckCount, report.AuditFile);
if ~report.Passed, disp(report.Checks(~report.Checks.Passed, :)); end

    function addCheck(name, passed, detail)
        report.Checks(end+1, :) = {currentCadence, currentEvent, string(name), ...
            logical(passed), string(detail)};
    end

    function appendArtifact(kind, file, startedUTC)
        item = inspectArtifact(kind, file, startedUTC);
        item.Cadence = currentCadence; item.EventID = currentEvent;
        report.Artifacts = appendRow(report.Artifacts, item);
        current = true;
        if endsWith(kind,'_png'), current = item.ModifiedWithinRun; end
        addCheck([kind, '_readable'], item.Exists && item.Readable && ...
            current, string(item.Path)+" "+item.Error);
    end
end

function [fluxOK, paOK, bOK, attitudeFound, directionOK] = validity(T, cadence)
if isempty(T)
    fluxOK = false(0, 1); paOK = fluxOK; bOK = fluxOK;
    attitudeFound = fluxOK; directionOK = fluxOK;
    return
end
F = fluxValues(T, cadence);
fluxOK = all(isfinite(F) & F > 0, 2);
paOK = all(isfinite(T{:, cellstr(compose('PA_S%d_deg', 1:7))}), 2);
if cadence == "day", fieldNames = {'BR_daily_nT', 'BT_daily_nT', 'BN_daily_nT'};
else, fieldNames = {'BR_hourly_nT', 'BT_hourly_nT', 'BN_hourly_nT'}; end
B = T{:, fieldNames};
bOK = all(isfinite(B), 2) & vecnorm(B, 2, 2) > 0;
directionOK = true(height(T), 1);
for component = ["R", "T", "N"]
    names = cellstr(compose('ParticleU%s_S%d', component, 1:7));
    directionOK = directionOK & all(isfinite(T{:, names}), 2);
end
data = T.Properties.UserData;
attitudeFound = logical(data.Attitude.Found(:));
end

function F = fluxValues(T, cadence)
if isempty(T), F = zeros(0, 7); return; end
if cadence == "day", suffix = '1d'; else, suffix = '1h'; end
% Use scalar sprintf inputs; compose treats a char suffix as format data.
names = cell(1,7);
for sector = 1:7
    names{sector} = sprintf('RawFlux_S%d_%s',sector,suffix);
end
F = T{:, names};
end

function A = inspectArtifact(kind, file, startedUTC)
A = struct('Kind', string(kind), 'Path', string(file), 'Exists', isfile(file), ...
    'Readable', false, 'ModifiedUTC', NaT(1, 1, 'TimeZone', 'UTC'), ...
    'ModifiedWithinRun', false, 'Bytes', NaN, 'Width', NaN, 'Height', NaN, ...
    'Format', "", 'SHA256', "", 'Error', "");
if ~A.Exists, A.Error = "file_missing"; return; end
try
    jFile = java.io.File(char(file));
    A.ModifiedUTC = datetime(double(jFile.lastModified())/1000, ...
        'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    A.Bytes = double(jFile.length());
    A.ModifiedWithinRun = A.ModifiedUTC >= startedUTC && ...
        A.ModifiedUTC <= datetime('now', 'TimeZone', 'UTC');
    A.SHA256 = string(Case1_File_SHA256(char(file)));
    if endsWith(kind, '_png')
        info = imfinfo(char(file));
        pixels = imread(char(file));
        A.Format = string(info.Format); A.Width = info.Width; A.Height = info.Height;
        A.Readable = strcmpi(info.Format, 'png') && size(pixels, 1) == info.Height && ...
            size(pixels, 2) == info.Width && info.Width > 0 && info.Height > 0;
    else
        info = whos('-file', char(file));
        A.Format = "MAT"; A.Readable = ~isempty(info);
    end
catch ME
    A.Error = string(ME.message);
end
end

function T = appendRow(T, row)
if isempty(T), T = struct2table(row); else, T = [T; struct2table(row)]; end
end
