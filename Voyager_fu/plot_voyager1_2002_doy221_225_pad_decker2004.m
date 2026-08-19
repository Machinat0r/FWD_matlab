% Reproduce the Voyager 1 LECP P1 pitch-angle distributions for the
% Decker et al. (2004) Figure 2 samples (1) and (3).
%
% Processing follows the paper: 24-hour LECP sector averages, a daily
% average MAG vector, seven active LECP sectors, normalization by the
% largest sector intensity of each day, and a second-order weighted fit in
% pitch cosine.  No interpolation or additional temporal smoothing is used.

programFolder = fileparts(mfilename('fullpath'));
monthlyRoot = fullfile(programFolder, 'Voyager_Interstellar_Monthly');
addpath(monthlyRoot);

dataRoot = fullfile(monthlyRoot, 'Voyager1_Selected_Event_Data');
lecpCDF = fullfile(dataRoot, 'voyager1', 'lecp', '1d', 'l2', ...
    'sectored_flux', '2002', '08', ...
    'voyager1_lecp_lev2_daily_sectored_20020808_20020814.cdf');
magCDF = fullfile(dataRoot, 'voyager1', 'coho', '1hr', 'l2', ...
    'merged_mag_plasma', '2002', '08', ...
    'voyager1_coho1hr_merged_mag_plasma_20020808_20020815.cdf');

if ~isfile(lecpCDF) || ~isfile(magCDF)
    downloader = fullfile(programFolder, ...
        'download_voyager1_2002_pad_data.py');
    pythonExe = fullfile(getenv('USERPROFILE'), '.cache', ...
        'codex-runtimes', 'codex-primary-runtime', 'dependencies', ...
        'python', 'python.exe');
    if ~isfile(pythonExe), pythonExe = 'python'; end
    command = sprintf('"%s" "%s"', pythonExe, downloader);
    [status, output] = system(command);
    fprintf('%s\n', output);
    assert(status == 0, 'VoyagerPAD:DownloadFailed', ...
        'The CDAWeb programmatic download failed.');
end

outputFolder = fullfile(monthlyRoot, ...
    'Voyager1_2002_PAD_Decker2004');
derivedFolder = fullfile(dataRoot, 'derived', 'lecp_pitch_angle_daily');
cacheRoot = fullfile(dataRoot, 'derived', 'cdf_cache');
if ~isfolder(outputFolder), mkdir(outputFolder); end
if ~isfolder(derivedFolder), mkdir(derivedFolder); end
if ~isfolder(cacheRoot), mkdir(cacheRoot); end

lecp = Voyager_Read_CDF_Product(lecpCDF, 'lecp_sector_daily', ...
    'CacheRoot', cacheRoot);
mag = Voyager_Read_CDF_Product(magCDF, 'coho', ...
    'CacheRoot', cacheRoot);

targetDOY = [221 225];
targetDate = datetime(2002, 1, 1, 'TimeZone', 'UTC') + ...
    days(targetDOY - 1);

energy = lecp.FHDU_Energy;
if size(energy, 1) == numel(lecp.Epoch)
    energy = energy(find(any(isfinite(energy), 2), 1), :);
end
energy = energy(:);
[~, p1Index] = min(abs(energy - 0.73));

% Published disk: +T right, +R down; S1 is upper-left.  These vectors
% point out of each telescope.  A detected particle entering the telescope
% has the opposite velocity direction.  The paper's pitch-cosine ordering
% uses the field azimuth in the LECP scan plane.
sector = (1:8).';
sectorPlotAngleDeg = (112.5 + (sector - 1) * 45);
telescopeUT = cosd(sectorPlotAngleDeg);
telescopeUR = -sind(sectorPlotAngleDeg);
particleUT = -telescopeUT;
particleUR = -telescopeUR;

allValues = table;
summary = table;
pad = repmat(struct, numel(targetDOY), 1);

for ii = 1:numel(targetDOY)
    dayStart = targetDate(ii);
    dayEnd = dayStart + days(1);

    lecpUse = lecp.Epoch >= dayStart & lecp.Epoch < dayEnd;
    assert(nnz(lecpUse) == 1, 'VoyagerPAD:LECPDailyRecord', ...
        'Expected one LECP daily record for DOY %d; found %d.', ...
        targetDOY(ii), nnz(lecpUse));
    lecpRecord = find(lecpUse, 1);
    sectorFlux = reshape(lecp.FHDU_SectoredFluxes( ...
        lecpRecord, p1Index, 1:8), 8, 1);
    sectorSigma = reshape(lecp.FHDU_SectoredFluxUncertainties( ...
        lecpRecord, p1Index, 1:8), 8, 1);

    magUse = mag.Epoch >= dayStart & mag.Epoch < dayEnd & ...
        isfinite(mag.BR) & isfinite(mag.BT) & isfinite(mag.BN);
    assert(any(magUse), 'VoyagerPAD:MAGDailyRecord', ...
        'No valid hourly MAG vector records for DOY %d.', targetDOY(ii));
    br = mean(mag.BR(magUse), 'omitnan');
    bt = mean(mag.BT(magUse), 'omitnan');
    bn = mean(mag.BN(magUse), 'omitnan');
    if isfield(mag, 'ABS_B')
        meanAbsB = mean(mag.ABS_B(magUse), 'omitnan');
    else
        meanAbsB = sqrt(br.^2 + bt.^2 + bn.^2);
    end
    vectorMagnitude = sqrt(br.^2 + bt.^2 + bn.^2);
    scanMagnitude = hypot(br, bt);
    elevationDeg = atan2d(bn, scanMagnitude);
    azimuthDeg = mod(atan2d(bt, br), 360);

    pitchCosine = (br * particleUR + bt * particleUT) / scanMagnitude;
    pitchCosine = max(-1, min(1, pitchCosine));
    pitchAngleDeg = acosd(pitchCosine);
    pitchCosineFullB = (br * particleUR + bt * particleUT) / ...
        vectorMagnitude;

    active = sector <= 7;
    reliable = active & isfinite(sectorFlux) & sectorFlux > 0 & ...
        isfinite(sectorSigma) & sectorSigma > 0;
    assert(nnz(reliable) == 7, 'VoyagerPAD:InvalidSector', ...
        'One or more active-sector samples are invalid on DOY %d.', ...
        targetDOY(ii));

    normalizer = max(sectorFlux(reliable));
    normalizedIntensity = sectorFlux / normalizer;
    normalizedSigma = sectorSigma / normalizer;

    muFit = pitchCosine(reliable);
    yFit = normalizedIntensity(reliable);
    sigmaFit = normalizedSigma(reliable);
    design = [ones(size(muFit)), muFit, muFit.^2];
    weight = 1 ./ sigmaFit.^2;
    coefficients = (design.' * (weight .* design)) \ ...
        (design.' * (weight .* yFit));
    residual = yFit - design * coefficients;
    chiSquare = sum((residual ./ sigmaFit).^2);
    reducedChiSquare = chiSquare / max(1, nnz(reliable) - 3);

    pad(ii).doy = targetDOY(ii);
    pad(ii).date = dayStart;
    pad(ii).mu = pitchCosine;
    pad(ii).normalizedIntensity = normalizedIntensity;
    pad(ii).normalizedSigma = normalizedSigma;
    pad(ii).coefficients = coefficients;
    pad(ii).meanAbsB = meanAbsB;
    pad(ii).elevationDeg = elevationDeg;
    pad(ii).azimuthDeg = azimuthDeg;
    pad(ii).nMAG = nnz(magUse);
    pad(ii).streamingTSign = sign(coefficients(2) * bt);

    quality = repmat("valid active-sector daily average", 8, 1);
    quality(8) = "excluded: S8 blocked by sun shield";
    current = table(repmat(dayStart, 8, 1), ...
        repmat(targetDOY(ii), 8, 1), sector, active, reliable, ...
        sectorPlotAngleDeg, particleUR, particleUT, pitchCosine, ...
        pitchAngleDeg, pitchCosineFullB, sectorFlux, sectorSigma, ...
        normalizedIntensity, normalizedSigma, quality, ...
        'VariableNames', {'DateUTC', 'DOY', 'Sector', 'ActiveSector', ...
        'UsedInPAD', 'DiskAngleDeg', 'ParticleUR', 'ParticleUT', ...
        'PitchCosine_ScanPlaneB', 'PitchAngleDeg_ScanPlaneB', ...
        'PitchCosine_FullB_Diagnostic', ...
        'Flux_counts_per_MeV_cm2_s_sr', 'FluxUncertainty_1sigma', ...
        'NormalizedIntensity', 'NormalizedUncertainty_1sigma', ...
        'Quality'});
    allValues = [allValues; current]; %#ok<AGROW>

    currentSummary = table(dayStart, targetDOY(ii), nnz(magUse), ...
        br, bt, bn, meanAbsB, vectorMagnitude, elevationDeg, azimuthDeg, ...
        energy(p1Index), normalizer, coefficients(1), coefficients(2), ...
        coefficients(3), reducedChiSquare, ...
        "LECP supplied 24-h daily average", ...
        "arithmetic mean of valid COHO 1-h vectors in UTC day", ...
        "weighted second-order polynomial in pitch cosine", ...
        'VariableNames', {'DateUTC', 'DOY', 'MAGHourlyRecords', ...
        'BR_daily_nT', 'BT_daily_nT', 'BN_daily_nT', ...
        'MeanAbsB_daily_nT', 'MagnitudeOfMeanVector_nT', ...
        'BElevationDeg', 'BAzimuthDeg', 'P1NominalEnergy_MeV', ...
        'FluxNormalizer', 'FitC0', 'FitC1', 'FitC2', ...
        'ReducedChiSquare', 'LECPAveraging', 'MAGAveraging', 'FitMethod'});
    summary = [summary; currentSummary]; %#ok<AGROW>
end

valuesCSV = fullfile(derivedFolder, ...
    'V1_2002_DOY221_225_LECP_P1_PAD_Decker2004_values.csv');
summaryCSV = fullfile(derivedFolder, ...
    'V1_2002_DOY221_225_LECP_P1_PAD_Decker2004_summary.csv');
matFile = fullfile(derivedFolder, ...
    'V1_2002_DOY221_225_LECP_P1_PAD_Decker2004.mat');
writetable(allValues, valuesCSV);
writetable(summary, summaryCSV);
save(matFile, 'allValues', 'summary', 'pad', 'lecpCDF', 'magCDF');

fig = figure('Color', 'white', 'Visible', 'on', ...
    'Name', 'Voyager 1 LECP P1 PADs: 2002 DOY 221 and 225', ...
    'Position', [120 120 960 455]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
panelLabel = {'(1)', '(3)'};

for ii = 1:numel(pad)
    ax = nexttile(layout);
    hold(ax, 'on');
    active = 1:7;
    errorbar(ax, pad(ii).mu(active), ...
        pad(ii).normalizedIntensity(active), ...
        pad(ii).normalizedSigma(active), 'o', ...
        'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
        'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5.5, ...
        'LineWidth', 1.0, 'CapSize', 6);
    muCurve = linspace(-1, 1, 401).';
    fitCurve = [ones(size(muCurve)), muCurve, muCurve.^2] * ...
        pad(ii).coefficients;
    plot(ax, muCurve, fitCurve, '--', 'Color', [0.15 0.15 0.15], ...
        'LineWidth', 1.25);

    xlim(ax, [-1.05 1.05]);
    ylim(ax, [0 1.12]);
    xticks(ax, [-1 -0.5 0 0.5 1]);
    yticks(ax, 0:0.2:1.0);
    xlabel(ax, 'pitch cosine');
    if ii == 1, ylabel(ax, 'normalized intensity'); end
    title(ax, sprintf('2002_%03d', pad(ii).doy), ...
        'Interpreter', 'none', 'FontWeight', 'normal');
    text(ax, 0.02, 0.98, panelLabel{ii}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontWeight', 'bold');
    text(ax, 0.50, 0.91, sprintf( ...
        'B=(%.3f, %.1f, %.1f)', pad(ii).meanAbsB, ...
        pad(ii).elevationDeg, pad(ii).azimuthDeg), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 9, ...
        'Interpreter', 'none');
    if pad(ii).streamingTSign < 0
        marker = '<'; markerX = -0.72;
    else
        marker = '>'; markerX = 0.72;
    end
    plot(ax, markerX, 0.79, marker, 'Color', [0 0 0], ...
        'MarkerFaceColor', [0 0 0], 'MarkerSize', 10, ...
        'HandleVisibility', 'off');
    set(ax, 'Box', 'on', 'TickDir', 'in', 'FontSize', 10, ...
        'LineWidth', 0.9, 'XMinorTick', 'on', 'YMinorTick', 'on');
    grid(ax, 'off');
end

sgtitle(layout, ['Voyager 1 LECP 0.57-1.78 MeV protons: ' ...
    '24-h pitch-angle distributions'], 'FontWeight', 'normal');

pngFile = fullfile(outputFolder, ...
    'Voyager1_LECP_P1_PAD_2002_DOY221_225_Decker2004.png');
pdfFile = fullfile(outputFolder, ...
    'Voyager1_LECP_P1_PAD_2002_DOY221_225_Decker2004.pdf');
exportgraphics(fig, pngFile, 'Resolution', 240);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

fprintf('\nVoyager 1 PAD reproduction completed.\n');
fprintf('Figure: %s\n', pngFile);
fprintf('Values: %s\n', valuesCSV);
fprintf('Summary: %s\n\n', summaryCSV);
disp(summary(:, {'DOY', 'MAGHourlyRecords', 'MeanAbsB_daily_nT', ...
    'BElevationDeg', 'BAzimuthDeg', 'ReducedChiSquare'}));
