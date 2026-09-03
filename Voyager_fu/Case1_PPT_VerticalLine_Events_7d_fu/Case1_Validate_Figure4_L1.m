function report = Case1_Validate_Figure4_L1(runAuditFile)
%Case1_Validate_Figure4_L1 Independently verify the four regenerated days.
%   Read-only for source CDFs, production code and figures. This validator
%   does not call the L1 reader, fallback/average function or PAD renderer.
%   It directly reads original variables with IRFU dataobj, then computes
%   each selected UTC bin, uncertainty, B vector and 3-D dot product anew.

%% fixed requested audit and hashes
if nargin < 1
    runAuditFile = ['Z:/SPART-WORK/Data/Voyager/voyager1/lecp/validation/' ...
        'florinski_figure4/Florinski2008_Figure4_run_20260903_095215_226.mat'];
end
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
loaded = load(runAuditFile, 'result'); result = loaded.result;
startTime = datetime(2004, 8, 8, 'TimeZone', 'UTC');
endTime = startTime+days(4);
checkName = strings(0, 1); checkPassed = false(0, 1); checkDifference = zeros(0, 1);
code = result.CodeManifest;
code.CurrentSHA256 = strings(height(code), 1);
for ii = 1:height(code)
    code.CurrentSHA256(ii) = Case1_File_SHA256(char(code.SourceFile(ii)));
end
verify('Every recorded production-code hash is unchanged', ...
    all(code.SHA256 == code.CurrentSHA256), 0);
sources = [result.SourceManifests{1}; result.SourceManifests{2}];
[~, uniqueSource] = unique(sources.SourceFile, 'stable');
sources = sources(uniqueSource, :);
sources.CurrentSHA256 = strings(height(sources), 1);
for ii = 1:height(sources)
    sources.CurrentSHA256(ii) = Case1_File_SHA256(char(sources.SourceFile(ii)));
end
verify('All source CDF and attitude/support hashes are unchanged', ...
    all(sources.SHA256 == sources.CurrentSHA256), 0);
rateFile = sources.SourceFile(contains(sources.DatasetID, 'LEV-1-RATES'));
cohoFile = sources.SourceFile(contains(sources.Instrument, 'COHO'));
assert(numel(rateFile) == 1 && numel(cohoFile) == 1);

%% direct CDF variables, with explicit source fill/range handling
rateCDF = dataobj(char(rateFile));
epoch = readTime(getv(rateCDF, 'Epoch'));
delta = cleanValue(getv(rateCDF, 'DeltaT'));
rate = cleanValue(getv(rateCDF, 'FHDU_SectoredRates'));
sigma = cleanValue(getv(rateCDF, 'FHDU_SectoredRateUncertainties'));
inWindow = epoch >= startTime & epoch < endTime;
rawRecords = find(inWindow);
negativeRecords = find(inWindow & delta(:) < 0);
rows = find(inWindow & ~(delta(:) < 0));
rate = reshape(rate(rows, 10, :), numel(rows), 8);
sigma = reshape(sigma(rows, 10, :), numel(rows), 8);
epoch = epoch(rows); delta = delta(rows);
cohoCDF = dataobj(char(cohoFile));
magEpoch = readTime(getv(cohoCDF, 'Epoch'));
mag = [cleanValue(getv(cohoCDF, 'BR')), cleanValue(getv(cohoCDF, 'BT')), ...
    cleanValue(getv(cohoCDF, 'BN'))];
magKeep = magEpoch >= startTime & magEpoch < endTime;
magEpoch = magEpoch(magKeep); mag = mag(magKeep, :);
factor = 1/(0.44*(1.78-0.57));
tables = {result.DailyPAD, result.HourlyPAD};
cadences = ["day", "hour"]; suffixes = ["1d", "1h"];
derived = cell(2, 1);
for mode = 1:2
    tab = tables{mode}; cadence = cadences(mode); suffix = suffixes(mode);
    if cadence == "day", width = days(1); expectedCount = 4; bSuffix = 'daily';
    else, width = hours(1); expectedCount = 55; bSuffix = 'hourly'; end
    verify(cadence+" complete PAD count", ...
        height(tab) == expectedCount && nnz(tab.PADUsable) == expectedCount, 0);
    verify(cadence+" source identity and window", ...
        all(tab.SourceProduct == "L1_UTC_mean") && ...
        all(tab.EpochUTC >= startTime & tab.EpochUTC < endTime) && ...
        all(isnan(tab.SourceCDFRecord)) && all(isnan(tab.DeltaT_s)), 0);
    n = height(tab);
    means = nan(n, 8); errors = nan(n, 8); counts = zeros(n, 8);
    meanB = nan(n, 3); nB = zeros(n, 1); pa = nan(n, 8);
    normError = zeros(n, 8); nativeRecordIDs = cell(n, 1);
    un = nan(n, 8);
    for ii = 1:n
        first = dateshift(tab.EpochUTC(ii), 'start', char(cadence));
        selected = find(epoch >= first & epoch < first+width);
        nativeRecordIDs{ii} = rows(selected);
        verify(cadence+" source record mapping "+ii, ...
            isequal(tab.L1SourceRecords{ii}.SourceCDFRecord, rows(selected)), 0);
        verify(cadence+" source original DeltaT "+ii, ...
            isequaln(tab.L1SourceRecords{ii}.DeltaT_s, delta(selected)), 0);
        verify(cadence+" representative bin center "+ii, ...
            tab.EpochUTC(ii) == first+width/2, 0);
        for sector = 1:8
            available = selected(isfinite(rate(selected, sector)) & rate(selected, sector) >= 0);
            counts(ii, sector) = numel(available);
            if isempty(available), continue; end
            means(ii, sector) = sum(rate(available, sector))/numel(available);
            if all(isfinite(sigma(available, sector)) & sigma(available, sector) >= 0)
                errors(ii, sector) = sqrt(sum(sigma(available, sector).^2))/numel(available);
            end
        end
        magSelected = magEpoch >= first & magEpoch < first+width & all(isfinite(mag), 2);
        nB(ii) = nnz(magSelected);
        meanB(ii, :) = sum(mag(magSelected, :), 1)/nB(ii);
        bhat = meanB(ii, :)/sqrt(sum(meanB(ii, :).^2));
        for sector = 1:8
            u = [tab.(sprintf('ParticleUR_S%d', sector))(ii), ...
                tab.(sprintf('ParticleUT_S%d', sector))(ii), ...
                tab.(sprintf('ParticleUN_S%d', sector))(ii)];
            normError(ii, sector) = abs(sum(u.^2)-1);
            un(ii, sector) = u(3);
            pa(ii, sector) = acosd(max(-1, min(1, sum(u.*bhat))));
        end
    end
    fluxNames = compose("RawFlux_S%d_%s", (1:8).', suffix);
    sigmaNames = compose("FluxUncertainty_S%d_%s", (1:8).', suffix);
    countNames = compose("Samples_S%d_%s", (1:8).', suffix);
    rateNames = compose("RawRate_S%d_%s", (1:8).', suffix);
    bNames = {['BR_', bSuffix, '_nT'], ['BT_', bSuffix, '_nT'], ['BN_', bSuffix, '_nT']};
    compare(cadence+" independently computed sector means", means, tab{:, cellstr(rateNames)}, 1e-13);
    compare(cadence+" converted flux", means*factor, tab{:, cellstr(fluxNames)}, 1e-13);
    compare(cadence+" propagated sigma", errors*factor, tab{:, cellstr(sigmaNames)}, 1e-13);
    compare(cadence+" sector sample counts", counts, tab{:, cellstr(countNames)}, 0);
    compare(cadence+" COHO complete-vector mean", meanB, tab{:, bNames}, 1e-14);
    compare(cadence+" COHO vector count", nB, tab.MAGVectorSampleCount, 0);
    compare(cadence+" full three-dimensional pitch-angle dot product", ...
        pa, tab{:, cellstr(compose("PA_S%d_deg", (1:8).'))}, 1e-10);
    verify(cadence+" unit direction norms", max(normError, [], 'all') < 1e-12, max(normError, [], 'all'));
    verify(cadence+" nonzero normal components retained", all(abs(un) > 1e-6, 'all'), min(abs(un), [], 'all'));
    derived{mode} = struct('RateMean', means, 'RateSigma', errors, ...
        'Flux', means*factor, 'FluxSigma', errors*factor, 'SectorCount', counts, ...
        'MagneticVectorMean', meanB, 'MagneticVectorCount', nB, ...
        'PitchAngle_deg', pa, 'ParticleUN', un, 'OriginalRecordIDs', {nativeRecordIDs});
end

%% fixed four panels and optional old-rate-mean consistency
four = result.FourDayPAD;
expectedFlux = derived{1}.Flux(:, 1:7);
expectedSigma = derived{1}.FluxSigma(:, 1:7);
maximum = max(expectedFlux, [], 2);
compare('Four-panel raw flux', expectedFlux, four.RawFlux, 1e-13);
compare('Four-panel raw sigma', expectedSigma, four.RawSigma, 1e-13);
compare('Four-panel Jmax', maximum, four.NormalizationFlux, 1e-13);
compare('Four-panel J/Jmax', expectedFlux./maximum, four.NormalizedFlux, 1e-13);
compare('Four-panel sigma/Jmax', expectedSigma./maximum, four.NormalizedSigma, 1e-13);
compare('Four-panel 3-D angle', derived{1}.PitchAngle_deg(:, 1:7), four.PA_deg, 1e-10);
verify('Four fixed daily panels complete with 0-180 degree axis', ...
    four.PlottedCount == 4 && all(four.SelectedPADUsable) && ...
    isequal(four.DisplayXLimits_deg, [0, 180]), 0);
oldFile = ['C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/' ...
    'Voyager_Interstellar_Monthly/Voyager1_Selected_Event_Data/derived/' ...
    'lecp_pitch_angle_florinski2008_daily/' ...
    'V1_2004-Fig4-DOY221-224_20040808_20040811_LECP_P1_Florinski2008_pitch_angle_1d.csv'];
oldComparison = struct('File', string(oldFile), 'Available', isfile(oldFile));
if oldComparison.Available
    old = readtable(oldFile, 'VariableNamingRule', 'preserve');
    old.EpochUTC = datetime(old.EpochUTC, 'TimeZone', 'UTC');
    old = old(old.EpochUTC >= startTime & old.EpochUTC < endTime, :);
    oldFlux = old{:, cellstr(compose("Flux_S%d_1d", (1:8).'))};
    oldComparison.FluxMaxAbsDifference = max(abs(oldFlux-derived{1}.Flux), [], 'all');
    oldComparison.FactorMaxAbsDifference = max(abs(old.SourceToDifferentialFluxFactor-factor));
    compare('Old daily rate-to-flux results agree; old angles unused', oldFlux, derived{1}.Flux, 1e-12);
end

%% standalone audit; no prior output is modified
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.InputRunAudit = string(runAuditFile);
report.InputRunAuditSHA256 = string(Case1_File_SHA256(char(runAuditFile)));
report.CodeManifest = code; report.SourceManifest = sources;
report.Checks = table(checkName, checkPassed, checkDifference, ...
    'VariableNames', {'Check', 'Passed', 'MaxAbsDifferenceOrDiagnostic'});
report.Passed = all(checkPassed); report.CheckCount = numel(checkPassed);
report.DirectCDFRecordCount = numel(rawRecords);
report.NegativeDeltaTRejected = negativeRecords;
report.Daily = derived{1}; report.Hourly = derived{2}; report.OldDailyComparison = oldComparison;
report.Scope = 'Only 2004-08-08 through 2004-08-11 Figure 4; no 45-event audit read or source changes.';
report.ValidatorFile = string([mfilename('fullpath'), '.m']);
report.ValidatorSHA256 = string(Case1_File_SHA256(char(report.ValidatorFile)));
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(fileparts(runAuditFile), ...
    ['Figure4_L1_independent_validation_', stamp, '.mat']));
save(report.AuditFile, 'report', '-v7.3');
disp(report.Checks(~report.Checks.Passed, :));
fprintf('Figure4 independent validation: %d/%d checks passed.\n', nnz(checkPassed), numel(checkPassed));
fprintf('Direct source records: %d; negative DeltaT: %d.\n', numel(rawRecords), numel(negativeRecords));
fprintf('Old daily flux maximum absolute difference: %.16g.\n', oldComparison.FluxMaxAbsDifference);
disp(report.Daily.SectorCount);
disp(report.Daily.MagneticVectorCount);
fprintf('Daily |uN| range: %.8g to %.8g.\n', min(abs(report.Daily.ParticleUN), [], 'all'), max(abs(report.Daily.ParticleUN), [], 'all'));
disp(report.AuditFile);
assert(report.Passed, 'Figure4 independent validation failed; see the saved audit.');

    function verify(name, passed, difference)
        checkName(end+1, 1) = name;
        checkPassed(end+1, 1) = passed;
        checkDifference(end+1, 1) = difference;
    end

    function compare(name, expected, actual, tolerance)
        sameMissing = isequal(isnan(expected), isnan(actual));
        difference = abs(expected-actual);
        difference = difference(isfinite(difference));
        if isempty(difference), maximumDifference = 0; else, maximumDifference = max(difference); end
        verify(name, isequal(size(expected), size(actual)) && sameMissing && ...
            maximumDifference <= tolerance, maximumDifference);
    end
end

function value = cleanValue(raw)
value = double(raw.data);
if isfield(raw, 'FILLVAL'), value(value == double(raw.FILLVAL)) = NaN; end
if isfield(raw, 'VALIDMIN') && isscalar(raw.VALIDMIN)
    value(value < double(raw.VALIDMIN)) = NaN;
end
if isfield(raw, 'VALIDMAX') && isscalar(raw.VALIDMAX)
    value(value > double(raw.VALIDMAX)) = NaN;
end
value(~isfinite(value) | abs(value) >= 1e30) = NaN;
end

function epoch = readTime(raw)
if contains(lower(string(raw.type)), 'tt2000')
    secondsUnix = EpochTT(int64(raw.data(:))).epochUnix;
else
    secondsUnix = double(raw.data(:));
end
epoch = datetime(secondsUnix, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
end
