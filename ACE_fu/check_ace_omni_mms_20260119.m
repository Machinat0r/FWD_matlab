%% MMS1 burst data for the 2026-01-19 strong-field event
% Directly executable script. The figure contains MMS1 FGM, FPI moments,
% FPI omnidirectional ion/electron energy flux, and FEEPS energetic
% electron omnidirectional intensity. No local functions are required.

clearvars
close all
clc

%% Paths and IRF-MATLAB
baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir)
    baseDir = 'C:\Users\Administrator\Documents\FWD_matlab\ACE_fu';
end
cacheDir = fullfile(baseDir, 'omni_mms_cache');
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

irfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';
if exist(irfDir, 'dir')
    addpath(genpath(irfDir));
end
if exist('irf_plot', 'file') ~= 2 || exist('dataobj', 'file') ~= 2
    error('IRF-MATLAB/dataobj was not found on the MATLAB path.');
end

burst0 = iso2epoch('2026-01-19T19:16:24Z');
burst1 = iso2epoch('2026-01-19T19:16:29Z');
burstTint = [burst0 burst1];

blue = [0.10 0.28 0.82];
green = [0.20 0.68 0.02];
red = [0.93 0.12 0.08];
black = [0.03 0.03 0.03];

%% Required local MMS files
fgmFile = fullfile(cacheDir, ...
    'mms1_fgm_brst_20260119T191620_191645.cdf');
disFile = fullfile(cacheDir, ...
    'mms1_fpi_brst_dis_20260119T191620_191645.cdf');
desFile = fullfile(cacheDir, ...
    'mms1_fpi_brst_des_20260119T191620_191645.cdf');
feepsFile = fullfile(cacheDir, ...
    'mms1_feeps_brst_l2_electron_20260119T191620_191645.cdf');

for requiredFile = {fgmFile, disFile, desFile}
    if ~isfile(requiredFile{1})
        error('Required MMS file is missing: %s', requiredFile{1});
    end
end

% FEEPS is a large daily file. If the local 25 s subset is absent, ask
% CDAWeb to generate a small subset and cache that CDF instead.
if ~isfile(feepsFile)
    feepsRequest = [ ...
        'https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/' ...
        'sp_phys/datasets/MMS1_FEEPS_BRST_L2_ELECTRON/data/' ...
        '20260119T191620Z,20260119T191645Z/ALL-VARIABLES?format=cdf'];
    fprintf('Preparing the 25 s MMS1 FEEPS subset from CDAWeb...\n');
    try
        requestOptions = weboptions('ContentType', 'json', 'Timeout', 180);
        feepsResult = webread(feepsRequest, requestOptions);
        if ~isfield(feepsResult, 'FileDescription') || ...
                isempty(feepsResult.FileDescription)
            error('CDAWeb returned no FEEPS subset file.');
        end
        downloadOptions = weboptions('Timeout', 180);
        websave(feepsFile, feepsResult.FileDescription(1).Name, ...
            downloadOptions);
    catch ME
        error('Unable to obtain the MMS1 FEEPS subset: %s', ME.message);
    end
end

%% MMS1 burst magnetic field
fgmObj = dataobj(fgmFile);
bVar = get_variable(fgmObj, 'mms1_fgm_b_gse_brst_l2_clean');
tB = EpochTT(bVar.DEPEND_0.data).epochUnix;
B = double(bVar.data);
B(abs(B) > 1e20) = NaN;

%% MMS1 FPI burst moments and omnidirectional spectra
disObj = dataobj(disFile);
desObj = dataobj(desFile);

niVar = get_variable(disObj, 'mms1_dis_numberdensity_brst');
viVar = get_variable(disObj, 'mms1_dis_bulkv_gse_brst');
tiParVar = get_variable(disObj, 'mms1_dis_temppara_brst');
tiPerpVar = get_variable(disObj, 'mms1_dis_tempperp_brst');
neVar = get_variable(desObj, 'mms1_des_numberdensity_brst');
teParVar = get_variable(desObj, 'mms1_des_temppara_brst');
tePerpVar = get_variable(desObj, 'mms1_des_tempperp_brst');

tNi = EpochTT(niVar.DEPEND_0.data).epochUnix;
tVi = EpochTT(viVar.DEPEND_0.data).epochUnix;
tTi = EpochTT(tiParVar.DEPEND_0.data).epochUnix;
tNe = EpochTT(neVar.DEPEND_0.data).epochUnix;
tTe = EpochTT(teParVar.DEPEND_0.data).epochUnix;

Ni = double(niVar.data);
Vi = double(viVar.data);
TiPar = double(tiParVar.data);
TiPerp = double(tiPerpVar.data);
Ne = double(neVar.data);
TePar = double(teParVar.data);
TePerp = double(tePerpVar.data);

Ni(abs(Ni) > 1e20) = NaN;
Vi(abs(Vi) > 1e20) = NaN;
TiPar(abs(TiPar) > 1e20) = NaN;
TiPerp(abs(TiPerp) > 1e20) = NaN;
Ne(abs(Ne) > 1e20) = NaN;
TePar(abs(TePar) > 1e20) = NaN;
TePerp(abs(TePerp) > 1e20) = NaN;

iFluxVar = get_variable(disObj, ...
    'mms1_dis_energyspectr_omni_brst');
iEnergyVar = get_variable(disObj, 'mms1_dis_energy_brst');
eFluxVar = get_variable(desObj, ...
    'mms1_des_energyspectr_omni_brst');
eEnergyVar = get_variable(desObj, 'mms1_des_energy_brst');

tIFpi = EpochTT(iFluxVar.DEPEND_0.data).epochUnix;
tEFpi = EpochTT(eFluxVar.DEPEND_0.data).epochUnix;
iFlux = double(iFluxVar.data);
iEnergy = double(iEnergyVar.data);
eFlux = double(eFluxVar.data);
eEnergy = double(eEnergyVar.data);

iFlux(iFlux <= 0 | iFlux > 1e30) = NaN;
eFlux(eFlux <= 0 | eFlux > 1e30) = NaN;
iEnergy(iEnergy <= 0 | iEnergy > 1e8) = NaN;
eEnergy(eEnergy <= 0 | eEnergy > 1e8) = NaN;

% Ensure every energy sweep is ordered from low to high energy.
for it = 1:size(iEnergy, 1)
    [iEnergy(it, :), order] = sort(iEnergy(it, :));
    iFlux(it, :) = iFlux(it, order);
end
for it = 1:size(eEnergy, 1)
    [eEnergy(it, :), order] = sort(eEnergy(it, :));
    eFlux(it, :) = eFlux(it, order);
end

iSpec = struct('t', tIFpi, 'p', iFlux, 'f', iEnergy, ...
    'p_label', {{'log_{10} DEF_i'}}, ...
    'f_label', {{'E_i', '[eV]'}});
eSpec = struct('t', tEFpi, 'p', eFlux, 'f', eEnergy, ...
    'p_label', {{'log_{10} DEF_e'}}, ...
    'f_label', {{'E_e', '[eV]'}});

%% MMS1 FEEPS burst omnidirectional energetic-electron intensity
feepsObj = dataobj(feepsFile);
sensorIds = [1:5 9:12];
sides = {'top', 'bottom'};

probeName = sprintf([ ...
    'mms1_epd_feeps_brst_l2_electron_top_' ...
    'intensity_sensorid_%d'], sensorIds(1));
probeVar = get_variable(feepsObj, probeName);
nFeepsTime = size(probeVar.data, 1);
nFeepsEnergy = size(probeVar.data, 2);
nFeepsEyes = numel(sensorIds) .* numel(sides);
feepsEyes = nan(nFeepsTime, nFeepsEnergy, nFeepsEyes);
feepsEnergyEyes = nan(nFeepsEyes, nFeepsEnergy);
tFeeps = EpochTT(probeVar.DEPEND_0.data).epochUnix;

eyeIndex = 0;
for iSide = 1:numel(sides)
    for sensorId = sensorIds
        eyeIndex = eyeIndex + 1;
        prefix = sprintf( ...
            'mms1_epd_feeps_brst_l2_electron_%s_', sides{iSide});

        intensityName = sprintf('%sintensity_sensorid_%d', ...
            prefix, sensorId);
        qualityName = sprintf('%squality_indicator_sensorid_%d', ...
            prefix, sensorId);
        energyName = sprintf('%senergy_centroid_sensorid_%d', ...
            prefix, sensorId);

        intensityVar = get_variable(feepsObj, intensityName);
        qualityVar = get_variable(feepsObj, qualityName);
        energyVar = get_variable(feepsObj, energyName);

        intensity = double(intensityVar.data);
        intensity(intensity < 0 | intensity > 1e30) = NaN;

        % Quality 0 is valid and 1 is caution. Values 2--4 are excluded.
        if ~isempty(qualityVar)
            quality = double(qualityVar.data);
            if isvector(quality)
                bad = repmat(quality(:) >= 2, 1, size(intensity, 2));
            else
                bad = quality >= 2;
            end
            intensity(bad) = NaN;
        end

        feepsEyes(:, :, eyeIndex) = intensity;
        feepsEnergyEyes(eyeIndex, :) = double(energyVar.data(:))';
    end
end

feepsFlux = mean(feepsEyes, 3, 'omitnan');
feepsFlux(feepsFlux <= 0 | feepsFlux > 1e30) = NaN;
feepsEnergy = mean(feepsEnergyEyes, 1, 'omitnan');
[feepsEnergy, order] = sort(feepsEnergy);
feepsFlux = feepsFlux(:, order);

feepsSpec = struct('t', tFeeps, 'p', feepsFlux, 'f', feepsEnergy, ...
    'p_label', {{'log_{10} J_e'}}, ...
    'f_label', {{'E_e', '[keV]'}});

%% Robust color ranges over the displayed five-second interval
iWindow = tIFpi >= burst0 & tIFpi <= burst1;
iLog = log10(iFlux(iWindow, :));
iLog = iLog(isfinite(iLog));
if isempty(iLog)
    iColor = [6 9];
else
    iColor = prctile(iLog, [5 99]);
end
iColor = [floor(iColor(1) .* 10) ceil(iColor(2) .* 10)] ./ 10;

eWindow = tEFpi >= burst0 & tEFpi <= burst1;
eLog = log10(eFlux(eWindow, :));
eLog = eLog(isfinite(eLog));
if isempty(eLog)
    eColor = [7 9];
else
    eColor = prctile(eLog, [5 99]);
end
eColor = [floor(eColor(1) .* 10) ceil(eColor(2) .* 10)] ./ 10;

feepsWindow = tFeeps >= burst0 & tFeeps <= burst1;
feepsLog = log10(feepsFlux(feepsWindow, :));
feepsLog = feepsLog(isfinite(feepsLog));
if isempty(feepsLog)
    feepsColor = [2 5];
else
    feepsColor = prctile(feepsLog, [5 99]);
end
feepsColor = [floor(feepsColor(1) .* 10) ...
    ceil(feepsColor(2) .* 10)] ./ 10;

%% MMS1 burst figure: moments plus low- and high-energy spectra
h = irf_plot(8, 'newfigure');
fig = gcf;
set(fig, 'Color', 'w', ...
    'Name', 'MMS1 FGM FPI FEEPS 2026-01-19', ...
    'Units', 'pixels', 'Position', [80 20 1100 1450], ...
    'Renderer', 'painters');
figData = get(fig, 'UserData');
if ~isstruct(figData), figData = struct; end
figData.t_start_epoch = burst0;
set(fig, 'UserData', figData);

% (a) MMS1 magnetic-field components and magnitude
set(h(1), 'ColorOrder', [blue; green; red; black], ...
    'NextPlot', 'replacechildren');
irf_plot(h(1), [tB B], 'LineWidth', 0.8);
ylabel(h(1), 'B [nT]');
ylim(h(1), [-55 55]);
lineHandles = flipud(findobj(h(1), 'Type', 'line'));
lgd = legend(h(1), lineHandles, {'B_x','B_y','B_z','|B|'}, ...
    'Location', 'northeast', 'NumColumns', 4, 'Box', 'off', ...
    'FontSize', 8, 'Interpreter', 'tex', 'AutoUpdate', 'off');
lgd.ItemTokenSize = [10 7];

% (b) FPI ion and electron number density
set(h(2), 'ColorOrder', [green; blue], ...
    'NextPlot', 'replacechildren');
irf_plot(h(2), [tNi Ni], 'LineWidth', 0.85);
hold(h(2), 'on');
irf_plot(h(2), [tNe Ne], 'LineWidth', 0.85, 'Color', blue);
hold(h(2), 'off');
ylabel(h(2), 'N [cm^{-3}]');
ylim(h(2), [5 45]);
lineHandles = flipud(findobj(h(2), 'Type', 'line'));
lgd = legend(h(2), lineHandles, {'N_i','N_e'}, ...
    'Location', 'northeast', 'NumColumns', 2, 'Box', 'off', ...
    'FontSize', 8, 'Interpreter', 'tex', 'AutoUpdate', 'off');
lgd.ItemTokenSize = [10 7];

% (c) FPI ion bulk velocity in GSE
set(h(3), 'ColorOrder', [blue; green; red], ...
    'NextPlot', 'replacechildren');
irf_plot(h(3), [tVi Vi], 'LineWidth', 0.8);
ylabel(h(3), 'V_i [km/s]');
ylim(h(3), [-650 150]);
lineHandles = flipud(findobj(h(3), 'Type', 'line'));
lgd = legend(h(3), lineHandles, {'V_x','V_y','V_z'}, ...
    'Location', 'northeast', 'NumColumns', 3, 'Box', 'off', ...
    'FontSize', 8, 'Interpreter', 'tex', 'AutoUpdate', 'off');
lgd.ItemTokenSize = [10 7];

% (d) FPI ion parallel and perpendicular temperatures
set(h(4), 'ColorOrder', [red; blue], ...
    'NextPlot', 'replacechildren');
irf_plot(h(4), [tTi TiPar], 'LineWidth', 0.85);
hold(h(4), 'on');
irf_plot(h(4), [tTi TiPerp], 'LineWidth', 0.85, 'Color', blue);
hold(h(4), 'off');
ylabel(h(4), 'T_i [eV]');
ylim(h(4), [150 1250]);
lineHandles = flipud(findobj(h(4), 'Type', 'line'));
lgd = legend(h(4), lineHandles, {'T_{i\parallel}','T_{i\perp}'}, ...
    'Location', 'northeast', 'NumColumns', 2, 'Box', 'off', ...
    'FontSize', 8, 'Interpreter', 'tex', 'AutoUpdate', 'off');
lgd.ItemTokenSize = [10 7];

% (e) FPI electron parallel and perpendicular temperatures
set(h(5), 'ColorOrder', [red; blue], ...
    'NextPlot', 'replacechildren');
irf_plot(h(5), [tTe TePar], 'LineWidth', 0.85);
hold(h(5), 'on');
irf_plot(h(5), [tTe TePerp], 'LineWidth', 0.85, 'Color', blue);
hold(h(5), 'off');
ylabel(h(5), 'T_e [eV]');
ylim(h(5), [140 290]);
lineHandles = flipud(findobj(h(5), 'Type', 'line'));
lgd = legend(h(5), lineHandles, {'T_{e\parallel}','T_{e\perp}'}, ...
    'Location', 'northeast', 'NumColumns', 2, 'Box', 'off', ...
    'FontSize', 8, 'Interpreter', 'tex', 'AutoUpdate', 'off');
lgd.ItemTokenSize = [10 7];

% (f) FPI ion omnidirectional differential energy flux
[~, cbI] = irf_spectrogram(h(6), iSpec, 'log');
set(h(6), 'YScale', 'log', 'YLim', [3 2e4], ...
    'YTick', 10 .^ (1:4));
ylabel(h(6), {'E_i', '[eV]'});
caxis(h(6), iColor);
cbI.Label.String = 'log_{10} DEF_i';
cbI.Label.Interpreter = 'tex';

% (g) FPI electron omnidirectional differential energy flux
[~, cbE] = irf_spectrogram(h(7), eSpec, 'log');
set(h(7), 'YScale', 'log', 'YLim', [5 3e4], ...
    'YTick', 10 .^ (1:4));
ylabel(h(7), {'E_e', '[eV]'});
caxis(h(7), eColor);
cbE.Label.String = 'log_{10} DEF_e';
cbE.Label.Interpreter = 'tex';

% (h) FEEPS energetic-electron omnidirectional intensity
[~, cbFeeps] = irf_spectrogram(h(8), feepsSpec, 'log');
set(h(8), 'YScale', 'log', 'YLim', [45 600], ...
    'YTick', [50 100 200 500]);
ylabel(h(8), {'E_e', '[keV]'});
caxis(h(8), feepsColor);
cbFeeps.Label.String = 'log_{10} J_e';
cbFeeps.Label.Interpreter = 'tex';

spectralMap = irf_colormap('waterfall');
colormap(h(6), spectralMap);
colormap(h(7), spectralMap);
colormap(h(8), spectralMap);

text(h(6), 0.985, 0.82, 'FPI ions', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'Color', 'w', 'FontName', 'Arial', 'FontSize', 9, ...
    'FontWeight', 'bold');
text(h(7), 0.985, 0.82, 'FPI electrons', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'Color', 'w', 'FontName', 'Arial', 'FontSize', 9, ...
    'FontWeight', 'bold');
text(h(8), 0.985, 0.82, 'FEEPS electrons', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'Color', 'w', 'FontName', 'Arial', 'FontSize', 9, ...
    'FontWeight', 'bold');

%% Shared axes, panel labels, and gap-free layout
panelLetters = 'abcdefgh';
for ip = 1:numel(h)
    if ip <= 5
        letterColor = 'k';
    else
        letterColor = 'w';
    end
    text(h(ip), 0.015, 0.82, panelLetters(ip), ...
        'Units', 'normalized', 'FontName', 'Arial', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'Color', letterColor);
    grid(h(ip), 'off');
    box(h(ip), 'on');
    set(h(ip), 'FontName', 'Arial', 'FontSize', 10, ...
        'LineWidth', 0.8, 'TickDir', 'in', 'Layer', 'top', ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
        'XMinorTick', 'off', 'XTickLabelRotation', 0);
end

irf_zoom(h, 'x', burstTint);
irf_timeaxis(h, 'date');
burstTicks = 0:1:(burst1 - burst0);
burstTickLabels = cellstr(datestr( ...
    datetime(burst0 + burstTicks, 'ConvertFrom', 'posixtime', ...
    'TimeZone', 'UTC'), 'HH:MM:SS'));
for ip = 1:numel(h)
    set(h(ip), 'XLim', [0 burst1 - burst0], ...
        'XTick', burstTicks, 'XMinorTick', 'off');
end
for ip = 1:(numel(h) - 1)
    set(h(ip), 'XTickLabel', []);
end
set(h(end), 'XTickLabel', burstTickLabels);
xlabel(h(end), '2026-01-19 UTC');

drawnow
panelLeft = 0.105;
panelRight = 0.875;
panelTop = 0.975;
panelBottom = 0.055;
panelGap = 0.0012;
panelWeights = [1 1 1 1 1 1.25 1.25 1.25];
availableHeight = panelTop - panelBottom - ...
    panelGap .* (numel(h) - 1);
heightUnit = availableHeight ./ sum(panelWeights);
panelCursor = panelTop;
panelPositions = nan(numel(h), 4);
for ip = 1:numel(h)
    panelHeight = panelWeights(ip) .* heightUnit;
    panelY = panelCursor - panelHeight;
    panelPositions(ip, :) = [panelLeft panelY ...
        panelRight - panelLeft panelHeight];
    set(h(ip), 'Units', 'normalized', ...
        'ActivePositionProperty', 'position', ...
        'Position', panelPositions(ip, :));
    panelCursor = panelY - panelGap;
end

colorbarLeft = 0.892;
colorbarWidth = 0.016;
set(cbI, 'Units', 'normalized', ...
    'Position', [colorbarLeft panelPositions(6, 2) ...
    colorbarWidth panelPositions(6, 4)], ...
    'FontName', 'Arial', 'FontSize', 8);
set(cbE, 'Units', 'normalized', ...
    'Position', [colorbarLeft panelPositions(7, 2) ...
    colorbarWidth panelPositions(7, 4)], ...
    'FontName', 'Arial', 'FontSize', 8);
set(cbFeeps, 'Units', 'normalized', ...
    'Position', [colorbarLeft panelPositions(8, 2) ...
    colorbarWidth panelPositions(8, 4)], ...
    'FontName', 'Arial', 'FontSize', 8);
drawnow

%% Save MATLAB, PNG, and editable vector PDF outputs
set(fig, 'PaperPositionMode', 'auto');
burstFig = fullfile(baseDir, ...
    'mms_burst_20260119_191624_191629.fig');
burstPng = fullfile(baseDir, ...
    'mms_burst_20260119_191624_191629.png');
burstPdf = fullfile(baseDir, ...
    'mms_burst_20260119_191624_191629.pdf');

savefig(fig, burstFig);
exportgraphics(fig, burstPng, 'Resolution', 300, ...
    'BackgroundColor', 'white');
exportgraphics(fig, burstPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');

fprintf('\nSaved MMS FGM/FPI/FEEPS figure: %s\n', burstPdf);
fprintf('Color ranges: FPI i %s, FPI e %s, FEEPS e %s\n', ...
    mat2str(iColor), mat2str(eColor), mat2str(feepsColor));
