function audit = Run_Case1_Selected6_20200923_PAD()
%Run_Case1_Selected6_20200923_PAD Six hours boxed by the user on 2026-09-04.
%   Reuse the exact hourly overview values after checking original L1 CDF,
%   same-hour magnetic vectors and the existing full 3D pointing routine.
%   Validation recalculations do not replace the overview payload.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-04

%% paths and six explicit UTC bin centers
cfg = Case1_Config;
cfg.PADCadence = 'hour';
Case1_Add_IRFU_Path(cfg.IRFURoot);
eventID = 'Case1-S01-L02';
parentName = ['V1_Case1-S01-L02_20200923_20200923_', ...
    'COHO1h_raw_LECP_P1_pitch_angle_predictedCK_1h_nativeCDF_Epoch'];
parentFile = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1h', 'derived', ...
    'pitch_angle', '2013-2021', 'predicted_ck', [parentName, '.mat']);
parentFigure = fullfile(cfg.OutputRoot, 'Epoch_hourly', [parentName, '.png']);
rateFile = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'native', 'l1', ...
    'sectored_rates', '2020', 'voyager-1_lecp_lev-1-rates_20200101_v1.1.1-01.cdf');
magFile = fullfile(cfg.DataRoot, 'voyager1', 'coho', '1hr', 'l2', ...
    'merged_mag_plasma', '2020', '09', 'voyager1_coho1hr_merged_mag_plasma_20200901_v01.cdf');
name = 'V1_Case1-S01-L02_20200923_PAD_selected6_20200922T2230_20200923T0330';
figureFile = fullfile(cfg.OutputRoot, 'Selected_PAD_times_hourly', [name, '.png']);
dataFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1h', 'derived', ...
    'selected_pitch_angle', '2020');
auditFile = fullfile(dataFolder, [name, '.mat']);
selectedEpoch = datetime(2020, 9, 22, 22, 30, 0, 'TimeZone', 'UTC') + hours((0:5).');
sourceFiles = string({parentFile; parentFigure; rateFile; magFile});
sourceHashes = strings(size(sourceFiles));
for ii = 1:numel(sourceFiles)
    sourceHashes(ii) = Case1_File_SHA256(sourceFiles(ii));
end
S = load(parentFile);
T = S.pitchAngleTable;
rows = zeros(6, 1);
for ii = 1:6
    row = find(T.EpochUTC == selectedEpoch(ii));
    assert(isscalar(row), 'Requested epoch has no unique overview row.');
    rows(ii) = row;
end
selected = T(rows, :);
assert(all(selected.PADUsable) && all(selected.SourceProduct == "L1_UTC_mean"), ...
    'Expected six complete L1 UTC-hour PAD records.');
assert(S.opts.ContextDays == 3 && strcmpi(S.opts.LECPSourcePriority, 'l1_first'), ...
    'Expected the current three-day context and L1-first overview.');

%% validate existing means and uncertainties directly against original CDF
fprintf('Checking six selected hours against original L1 and COHO CDFs...\n');
rates = Voyager_Read_CDF_Product(rateFile, 'lecp_sector_daily');
mag = Voyager_Read_CDF_Product(magFile, 'coho');
rateValue = squeeze(rates.FHDU_SectoredRates(:, 10, :));
rateSigma = squeeze(rates.FHDU_SectoredRateUncertainties(:, 10, :));
factor = 1/(0.44*(1.78-0.57));
Braw = [mag.BR(:), mag.BT(:), mag.BN(:)];
B = nan(6, 3);
sampleCount = zeros(6, 7);
rawRecordIndices = cell(6, 1);
for ii = 1:6
    begin = dateshift(selectedEpoch(ii), 'start', 'hour');
    stop = begin + hours(1);
    keep = rates.Epoch >= begin & rates.Epoch < stop & ~(rates.DeltaT(:) < 0);
    rr = find(keep);
    assert(~isempty(rr), 'No L1 source records in a selected hour.');
    trace = selected.L1SourceRecords{ii};
    assert(isequal(rr(:), trace.SourceCDFRecord(:)), 'Original record identities changed.');
    rawRecordIndices{ii} = rr;
    for iSector = 1:7
        r = rateValue(rr, iSector);
        s = rateSigma(rr, iSector);
        good = isfinite(r) & r >= 0;
        n = nnz(good);
        sampleCount(ii, iSector) = n;
        assert(n > 0, 'A selected sector has no valid rate sample.');
        meanRate = mean(r(good));
        meanSigma = NaN;
        if all(isfinite(s(good)) & s(good) >= 0)
            meanSigma = sqrt(sum(s(good).^2))/n;
        end
        compareValue(selected.(sprintf('RawRate_S%d_1h', iSector))(ii), meanRate);
        compareValue(selected.(sprintf('RawFlux_S%d_1h', iSector))(ii), meanRate*factor);
        compareValue(selected.(sprintf('FluxUncertainty_S%d_1h', iSector))(ii), meanSigma*factor);
        assert(selected.(sprintf('Samples_S%d_1h', iSector))(ii) == n, 'Sample count changed.');
    end
    assert(selected.SourceToDifferentialFluxFactor(ii) == factor, 'Rate conversion changed.');
    validB = mag.Epoch >= begin & mag.Epoch < stop & all(isfinite(Braw), 2);
    assert(nnz(validB) == selected.MAGVectorSampleCount(ii), 'MAG sample count changed.');
    B(ii, :) = mean(Braw(validB, :), 1);
end
compareValue(B, [selected.BR_hourly_nT, selected.BT_hourly_nT, selected.BN_hourly_nT]);

%% verify approved three-component pointing without changing the six payloads
pointing = Case1_Predicted_LECP_Pointing(selectedEpoch, B, cfg);
for iSector = 1:7
    compareValue(selected.(sprintf('PA_S%d_deg', iSector)), pointing.PitchAngle_deg(:, iSector));
    compareValue([selected.(sprintf('ParticleUR_S%d', iSector)), ...
        selected.(sprintf('ParticleUT_S%d', iSector)), selected.(sprintf('ParticleUN_S%d', iSector))], ...
        squeeze(pointing.ParticleRTN(:, iSector, :)));
end
fprintf('CDF rates, propagated sigma, sample counts, magnetic vectors and 3D PA all match.\n');

%% render exact selected rows; preserve the old overview and five-time figures
opts = struct('Visible', false, 'ExportDPI', cfg.ExportDPI, ...
    'Spacecraft', 1, 'P1DisplayEnergyMeV', cfg.P1DisplayEnergyMeV);
audit = Case1_Plot_Selected_PAD(T, selectedEpoch, eventID, figureFile, opts);
assert(isequal(max(audit.NormalizedFlux, [], 2), ones(6, 1)), 'Normalization failed.');
audit.ParentOverviewMAT = string(parentFile);
audit.ParentOverviewPNG = string(parentFigure);
audit.ParentProductionOptions = S.opts;
audit.ParentCodeManifest = S.codeManifest;
audit.ParentSourceLECP = S.sourceLECP;
audit.ParentSourceMAG = S.sourceMAG;
audit.Validation = struct('Passed', true, 'RawCDFRecords', {rawRecordIndices}, ...
    'SectorSampleCount', sampleCount, 'B_RTN_nT', B, 'Pointing', pointing, ...
    'NumericComparisonRelativeTolerance', 1e-11, ...
    'Description', "Original CDF rates, errors and hourly means; same-hour complete B vectors; current approved full 3D pointing; exact selected rows and per-panel maximum normalization.");
audit.SourceManifest = table(sourceFiles, sourceHashes, ...
    'VariableNames', {'SourceFile', 'SHA256'});
for ii = 1:numel(sourceFiles)
    assert(Case1_File_SHA256(sourceFiles(ii)) == sourceHashes(ii), 'An input file changed during rendering.');
end
codeNames = ["Run_Case1_Selected6_20200923_PAD.m"; "Case1_Plot_Selected_PAD.m"; ...
    "Voyager_Read_CDF_Product.m"; "Case1_Predicted_LECP_Pointing.m"; ...
    "Case1_Read_Predicted_Attitude.m"; "Case1_LECP_Geometry.m"; ...
    "Case1_Config.m"; "Case1_Add_IRFU_Path.m"; "Case1_File_SHA256.m"];
codeFiles = fullfile(string(cfg.CodeRoot), codeNames);
codeHashes = strings(size(codeFiles));
for ii = 1:numel(codeFiles), codeHashes(ii) = Case1_File_SHA256(codeFiles(ii)); end
audit.CurrentCodeManifest = table(codeFiles, codeHashes, ...
    'VariableNames', {'SourceFile', 'SHA256'});
audit.FigureSHA256 = Case1_File_SHA256(figureFile);
info = imfinfo(figureFile);
audit.FigurePixelSize = [info.Width info.Height];
audit.AuditFile = string(auditFile);
if ~isfolder(dataFolder), mkdir(dataFolder); end
save(auditFile, 'audit', '-v7.3');
fprintf('Saved six PAD panels: %s\n', figureFile);
fprintf('Saved reproducibility audit: %s\n', auditFile);
disp(table(selectedEpoch, audit.NormalizationFlux, ...
    'VariableNames', {'EpochUTC', 'Jmax'}));
end

function compareValue(actual, expected)
% Numerical agreement check only; this tolerance is not a data quality cut.
assert(isequal(size(actual), size(expected)), 'Validation shape mismatch.');
assert(isequal(isnan(actual), isnan(expected)), 'Validation missing-value mismatch.');
good = ~isnan(expected);
assert(all(isfinite(actual(good))) && all(isfinite(expected(good))), 'Unexpected infinity.');
assert(all(abs(actual(good)-expected(good)) <= 1e-11*max(1, abs(expected(good)))), ...
    'The overview payload differs from the original CDF/current geometry.');
end

