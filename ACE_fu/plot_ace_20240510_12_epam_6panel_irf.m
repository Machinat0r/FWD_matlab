%% ACE 2024 May superstorm magnetic field and EPAM overview
% Run this file directly as a MATLAB script. It intentionally contains no
% local functions so that all loaded and derived variables remain visible
% in the base workspace.

clearvars
close all
clc

%% Paths and IRF-MATLAB
baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir)
    baseDir = pwd;
end

irfCandidates = {
    'C:\Users\Administrator\Documents\irfu-matlab-master'
    'C:\Users\Administrator\Documents\irf_matlab'
    'C:\Users\Administrator\Documents\irfu-matlab'
    };
for ii = 1:numel(irfCandidates)
    if exist(irfCandidates{ii}, 'dir')
        addpath(genpath(irfCandidates{ii}));
        break
    end
end
if exist('irf_plot', 'file') ~= 2 || exist('irf_spectrogram', 'file') ~= 2
    error('IRF-MATLAB was not found. Add the irfu-matlab folder to the MATLAB path.');
end

cacheDir = fullfile(baseDir, 'hapi_cache');
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

helper = fullfile(baseDir, 'hapi_csv_to_mat_v4.py');
if ~exist(helper, 'file')
    error('Missing HAPI cache helper: %s', helper);
end

pyExe = 'python';
bundledPython = fullfile(getenv('USERPROFILE'), '.cache', 'codex-runtimes', ...
    'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
if exist(bundledPython, 'file')
    pyExe = bundledPython;
elseif ~isempty(getenv('ACE_PYTHON')) && exist(getenv('ACE_PYTHON'), 'file')
    pyExe = getenv('ACE_PYTHON');
end

%% Time interval and HAPI cache files
plot0 = datenum(2024, 5, 10, 0, 0, 0);
plot1 = datenum(2024, 5, 12, 0, 0, 0);
epoch0 = datenum(1970, 1, 1, 0, 0, 0);
tint = ([plot0 plot1] - epoch0) .* 86400;

mfiCsv1 = fullfile(cacheDir, 'AC_H3_MFI_20240510_Magnitude_BGSEc.csv');
mfiMat1 = fullfile(cacheDir, 'AC_H3_MFI_20240510_Magnitude_BGSEc.mat');
mfiCsv2 = fullfile(cacheDir, 'AC_H3_MFI_20240511_Magnitude_BGSEc.csv');
mfiMat2 = fullfile(cacheDir, 'AC_H3_MFI_20240511_Magnitude_BGSEc.mat');
mfiUrl1 = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_MFI', ...
    '&parameters=Magnitude,BGSEc', ...
    '&time.min=2024-05-10T00:00:00Z', ...
    '&time.max=2024-05-11T00:00:00Z&format=csv'];
mfiUrl2 = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_MFI', ...
    '&parameters=Magnitude,BGSEc', ...
    '&time.min=2024-05-11T00:00:00Z', ...
    '&time.max=2024-05-12T00:00:00Z&format=csv'];

epmCsv = fullfile(cacheDir, ...
    'AC_H3_EPM_20240510T0000_20240512T0000_all_telescopes.csv');
epmMat = fullfile(cacheDir, ...
    'AC_H3_EPM_20240510T0000_20240512T0000_all_telescopes.mat');
% HAPI requires requested parameters to follow their dataset order.
epmParams = ['P1,P2,P3,P4,P5,P6,P7,P8,', ...
    'DE1,DE2,DE3,DE4,', ...
    'E1p,E2p,E3p,E4p,', ...
    'P1p,P2p,P3p,P4p,P5p,P6p,P7p,P8p,', ...
    'E1,E2,E3,E4'];
epmUrl = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_EPM', ...
    '&parameters=', epmParams, ...
    '&time.min=2024-05-10T00:00:00Z', ...
    '&time.max=2024-05-12T00:00:00Z&format=csv'];

for mfiJob = {mfiCsv1,mfiMat1,mfiUrl1; mfiCsv2,mfiMat2,mfiUrl2}'
    if ~exist(mfiJob{2}, 'file')
        fprintf('Preparing daily MFI cache from CDAWeb/HAPI...\n');
        cmd = sprintf('"%s" "%s" "%s" "%s" "%s"', ...
            pyExe, helper, mfiJob{1}, mfiJob{2}, mfiJob{3});
        [status, cmdout] = system(cmd);
        if status ~= 0
            error('MFI download/cache preparation failed:\n%s', cmdout);
        end
    end
end

if ~exist(epmMat, 'file')
    fprintf('Preparing EPAM cache from CDAWeb/HAPI...\n');
    cmd = sprintf('"%s" "%s" "%s" "%s" "%s"', ...
        pyExe, helper, epmCsv, epmMat, epmUrl);
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('EPAM download/cache preparation failed:\n%s', cmdout);
    end
end

%% Load and validate data
mfiPart1 = load(mfiMat1);
mfiPart2 = load(mfiMat2);
mfi.t = [mfiPart1.t(:); mfiPart2.t(:)];
mfi.y = [mfiPart1.y; mfiPart2.y];  % [|B| Bx_GSE By_GSE Bz_GSE]
epm = load(epmMat);  % y column order follows epmParams above

tMfi = double(mfi.t(:));
yMfi = double(mfi.y);
tEpm = double(epm.t(:));
yEpm = double(epm.y);
yMfi(abs(yMfi) > 1e30) = NaN;
yEpm(abs(yEpm) > 1e30) = NaN;

if size(yMfi, 2) < 4
    error('MFI cache must contain [Magnitude BGSEc_x BGSEc_y BGSEc_z].');
end
if size(yEpm, 2) < 28
    error('EPAM cache has %d columns; 28 columns are required.', size(yEpm, 2));
end

Bmag = yMfi(:, 1);
BGSE = yMfi(:, 2:4);
J_LEMS30 = yEpm(:, 1:8);
J_DE = yEpm(:, 9:12);
J_LEFS60 = yEpm(:, 13:16);
J_LEMS120 = yEpm(:, 17:24);
J_LEFS150 = yEpm(:, 25:28);

J_LEMS30(J_LEMS30 <= 0) = NaN;
J_LEMS120(J_LEMS120 <= 0) = NaN;
J_DE(J_DE <= 0) = NaN;
J_LEFS60(J_LEFS60 <= 0) = NaN;
J_LEFS150(J_LEFS150 <= 0) = NaN;

tMfiEpoch = (tMfi - epoch0) .* 86400;
tEpmEpoch = (tEpm - epoch0) .* 86400;

fprintf('\nValid EPAM samples in the 2024 May storm interval:\n');
fprintf('LEMS30  P1-P8  :');
fprintf(' %d', sum(isfinite(J_LEMS30), 1));
fprintf('\n');
fprintf('LEMS120 P1p-P8p:');
fprintf(' %d', sum(isfinite(J_LEMS120), 1));
fprintf('\n');
fprintf('DE30    DE1-DE4:');
fprintf(' %d', sum(isfinite(J_DE), 1));
fprintf('\n');
fprintf('LEFS60  E1p-E4p:');
fprintf(' %d', sum(isfinite(J_LEFS60), 1));
fprintf('\n');
fprintf('LEFS150 E1-E4  :');
fprintf(' %d', sum(isfinite(J_LEFS150), 1));
fprintf('\n\n');

%% Energy bins in keV
edges_LEMS30 = [46 67 115 193 315 580 1060 1880 4700];
edges_LEMS120 = [47 68 115 195 321 580 1060 1900 4800];
edges_DE = [38 53 103 175 315];
edges_LEFS60 = [45 62 102 175 290];
edges_LEFS150 = [44 58 104 180 295];

energy_LEMS30 = sqrt(edges_LEMS30(1:end-1) .* edges_LEMS30(2:end));
energy_LEMS120 = sqrt(edges_LEMS120(1:end-1) .* edges_LEMS120(2:end));
energy_DE = sqrt(edges_DE(1:end-1) .* edges_DE(2:end));
energy_LEFS60 = sqrt(edges_LEFS60(1:end-1) .* edges_LEFS60(2:end));
energy_LEFS150 = sqrt(edges_LEFS150(1:end-1) .* edges_LEFS150(2:end));

dtEpm = median(diff(tEpmEpoch), 'omitnan') ./ 2;
if ~isfinite(dtEpm) || dtEpm <= 0
    dtEpm = 6;
end

spec_LEMS30.t = tEpmEpoch;
spec_LEMS30.f = energy_LEMS30;
spec_LEMS30.p = J_LEMS30;
spec_LEMS30.dt = dtEpm;
spec_LEMS30.df = diff(edges_LEMS30) ./ 2;
spec_LEMS30.f_label = 'E_i [keV]';
spec_LEMS30.p_label = 'log_{10} J_i';

spec_LEMS120.t = tEpmEpoch;
spec_LEMS120.f = energy_LEMS120;
spec_LEMS120.p = J_LEMS120;
spec_LEMS120.dt = dtEpm;
spec_LEMS120.df = diff(edges_LEMS120) ./ 2;
spec_LEMS120.f_label = 'E_i [keV]';
spec_LEMS120.p_label = 'log_{10} J_i';

spec_DE.t = tEpmEpoch;
spec_DE.f = energy_DE;
spec_DE.p = J_DE;
spec_DE.dt = dtEpm;
spec_DE.df = diff(edges_DE) ./ 2;
spec_DE.f_label = 'E_e [keV]';
spec_DE.p_label = 'log_{10} J_e';

spec_LEFS60.t = tEpmEpoch;
spec_LEFS60.f = energy_LEFS60;
spec_LEFS60.p = J_LEFS60;
spec_LEFS60.dt = dtEpm;
spec_LEFS60.df = diff(edges_LEFS60) ./ 2;
spec_LEFS60.f_label = 'E_e [keV]';
spec_LEFS60.p_label = 'log_{10} J_e';

spec_LEFS150.t = tEpmEpoch;
spec_LEFS150.f = energy_LEFS150;
spec_LEFS150.p = J_LEFS150;
spec_LEFS150.dt = dtEpm;
spec_LEFS150.df = diff(edges_LEFS150) ./ 2;
spec_LEFS150.f_label = 'E_e [keV]';
spec_LEFS150.p_label = 'log_{10} J_e';

%% Shared robust color limits
ionLog = [log10(J_LEMS30(:)); log10(J_LEMS120(:))];
ionLog = sort(ionLog(isfinite(ionLog)));
if isempty(ionLog)
    ionClim = [0 1];
else
    nIon = numel(ionLog);
    ionClim = [ionLog(max(1, round(0.02 .* nIon))), ...
        ionLog(min(nIon, round(0.98 .* nIon)))];
end
if ionClim(2) <= ionClim(1)
    ionClim = ionClim(1) + [-0.5 0.5];
end

electronLog = [log10(J_DE(:)); log10(J_LEFS60(:)); log10(J_LEFS150(:))];
electronLog = sort(electronLog(isfinite(electronLog)));
if isempty(electronLog)
    electronClim = [0 1];
else
    nElectron = numel(electronLog);
    electronClim = [electronLog(max(1, round(0.02 .* nElectron))), ...
        electronLog(min(nElectron, round(0.98 .* nElectron)))];
end
if electronClim(2) <= electronClim(1)
    electronClim = electronClim(1) + [-0.5 0.5];
end

%% Six IRF panels
h = irf_plot(6, 'newfigure');
fig = gcf;
set(fig, 'Color', 'w', 'Name', 'ACE 2024 May storm EPAM six-panel', ...
    'Units', 'pixels', 'Position', [80 25 1050 1280], 'Renderer', 'painters');

figUserData = get(fig, 'UserData');
if ~isstruct(figUserData)
    figUserData = struct;
end
figUserData.t_start_epoch = tint(1);
set(fig, 'UserData', figUserData);

blue = [0.12 0.26 0.78];
green = [0.32 0.72 0.00];
red = [0.93 0.12 0.08];

set(h(1), 'ColorOrder', [blue; green; red; 0 0 0], ...
    'NextPlot', 'replacechildren');
irf_plot(h(1), [tMfiEpoch BGSE Bmag], 'LineWidth', 0.75);
ylabel(h(1), 'B [nT]');
set(h(1), 'ColorOrder', [blue; green; red; 0 0 0]);
irf_legend(h(1), {'B_x','B_y','B_z','|B|'}, [0.72 0.90]);

irf_spectrogram(h(2), spec_LEMS30, 'log');
irf_spectrogram(h(3), spec_LEMS120, 'log');
irf_spectrogram(h(4), spec_DE, 'log');
irf_spectrogram(h(5), spec_LEFS60, 'log');
irf_spectrogram(h(6), spec_LEFS150, 'log');

set(h(2:6), 'YScale', 'log', 'Layer', 'top');
ylabel(h(2), 'E_i [keV]');
ylabel(h(3), 'E_i [keV]');
ylabel(h(4), 'E_e [keV]');
ylabel(h(5), 'E_e [keV]');
ylabel(h(6), 'E_e [keV]');

ylim(h(2), [edges_LEMS30(1) edges_LEMS30(end)]);
ylim(h(3), [edges_LEMS120(1) edges_LEMS120(end)]);
ylim(h(4), [edges_DE(1) edges_DE(end)]);
ylim(h(5), [edges_LEFS60(1) edges_LEFS60(end)]);
ylim(h(6), [edges_LEFS150(1) edges_LEFS150(end)]);
set(h(2), 'YTick', [50 100 300 1000 3000]);
set(h(3), 'YTick', [50 100 300 1000 3000]);
set(h(4:6), 'YTick', [50 100 300]);

clim(h(2), ionClim);
clim(h(3), ionClim);
clim(h(4), electronClim);
clim(h(5), electronClim);
clim(h(6), electronClim);
colormap(fig, turbo(256));

%% Panel labels and physical channel labels
panelLetters = 'abcdef';
for ii = 1:numel(h)
    text(h(ii), 0.018, 0.88, panelLetters(ii), 'Units', 'normalized', ...
        'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', 'k', 'Interpreter', 'tex');
end
text(h(2), 0.985, 0.89, 'LEMS30 ions: P_{1-8}', 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'w');
text(h(2), 0.985, 0.73, 'P_{1-6}: fill (LEMS30 detector noise)', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 8, 'Color', 'w', 'Interpreter', 'tex');
text(h(3), 0.985, 0.89, 'LEMS120 ions: P^{\prime}_{1-8}', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'w', 'Interpreter', 'tex');
text(h(4), 0.985, 0.89, 'DE30 electrons: DE_{1-4}', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'w', 'Interpreter', 'tex');
text(h(5), 0.985, 0.89, 'LEFS60: e^- (+ions), E^{\prime}_{1-4}', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'w', 'Interpreter', 'tex');
text(h(6), 0.985, 0.89, 'LEFS150: e^- (+ions), E_{1-4}', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'w', 'Interpreter', 'tex');
text(h(6), 0.985, 0.73, 'E_{1-3}: fill', 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 8, ...
    'Color', 'w', 'Interpreter', 'tex');

%% Time range and axis styling
irf_zoom(h, 'x', tint);
for ii = 1:numel(h)
    grid(h(ii), 'off');
    box(h(ii), 'on');
    set(h(ii), 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 0.8, ...
        'TickDir', 'in', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
        'XMinorTick', 'off', 'XTickLabelRotation', 0);
end
irf_timeaxis(h, 'date');
hourTicks = 0:6*3600:(tint(2) - tint(1));
hourTickLabels = cellstr(datestr(plot0 + hourTicks(:) ./ 86400, 'mm-dd HH:MM'));
for ii = 1:numel(h)
    set(h(ii), 'XTick', hourTicks, 'XMinorTick', 'off');
end
for ii = 1:5
    set(h(ii), 'XTickLabel', []);
end
set(h(6), 'XTickLabel', hourTickLabels);
xlabel(h(6), '2024-05-10 to 12 UTC');

%% Final gap-free layout and aligned colorbars
drawnow
panelLeft = 0.105;
panelRight = 0.805;
panelTop = 0.955;
panelBottom = 0.070;
panelGap = 0.0015;
panelHeight = (panelTop - panelBottom - panelGap .* (numel(h) - 1)) ./ numel(h);
for ii = 1:numel(h)
    panelY = panelTop - ii .* panelHeight - (ii - 1) .* panelGap;
    set(h(ii), 'Units', 'normalized', 'ActivePositionProperty', 'position', ...
        'Position', [panelLeft panelY panelRight - panelLeft panelHeight], ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
end

cb = findall(fig, 'Type', 'colorbar');
if ~isempty(cb)
    cbY = zeros(numel(cb), 1);
    for ii = 1:numel(cb)
        cbPos = get(cb(ii), 'Position');
        cbY(ii) = cbPos(2);
    end
    [~, cbOrder] = sort(cbY, 'descend');
    cb = cb(cbOrder);
    nColorbars = min(5, numel(cb));
    for ii = 1:nColorbars
        panelPos = get(h(ii + 1), 'Position');
        set(cb(ii), 'Units', 'normalized', ...
            'Position', [0.830, panelPos(2) + 0.07 .* panelPos(4), ...
            0.018, 0.86 .* panelPos(4)], ...
            'FontName', 'Arial', 'FontSize', 8, 'TickDirection', 'in');
        cb(ii).Label.String = 'log_{10} J';
        cb(ii).Label.Interpreter = 'tex';
        cb(ii).Label.FontSize = 9;
    end
end

drawnow

%% Save editable outputs
outFig = fullfile(baseDir, 'ace_20240510_12_epam_spectrogram.fig');
outPng = fullfile(baseDir, 'ace_20240510_12_epam_spectrogram.png');
outPdf = fullfile(baseDir, 'ace_20240510_12_epam_spectrogram.pdf');
savefig(fig, outFig);
exportgraphics(fig, outPng, 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Saved EPAM spectrogram: %s\n', outPdf);

set(gcf,'paperpositionmode','auto')
%% Second figure: all particle channels as line plots
cyan = [0.00 0.53 0.62];
purple = [0.70 0.24 0.70];
orange = [0.94 0.53 0.05];
gray = [0.28 0.28 0.28];
brown = [0.68 0.40 0.08];
fluxColors = [blue; green; red; cyan; purple; orange; gray; brown];

hLine = irf_plot(6, 'newfigure');
figLine = gcf;
set(figLine, 'Color', 'w', 'Name', 'ACE 2024 May storm EPAM line plots', ...
    'Units', 'pixels', 'Position', [110 25 1050 1280], 'Renderer', 'painters');

lineUserData = get(figLine, 'UserData');
if ~isstruct(lineUserData)
    lineUserData = struct;
end
lineUserData.t_start_epoch = tint(1);
set(figLine, 'UserData', lineUserData);

set(hLine(1), 'ColorOrder', [blue; green; red; 0 0 0], ...
    'NextPlot', 'replacechildren');
irf_plot(hLine(1), [tMfiEpoch BGSE Bmag], 'LineWidth', 0.75);
ylabel(hLine(1), 'B [nT]');
set(hLine(1), 'ColorOrder', [blue; green; red; 0 0 0]);
irf_legend(hLine(1), {'B_x','B_y','B_z','|B|'}, [0.74 0.90]);

set(hLine(2), 'ColorOrder', fluxColors, 'NextPlot', 'replacechildren');
irf_plot(hLine(2), [tEpmEpoch J_LEMS30], 'LineWidth', 0.70);
ylabel(hLine(2), 'J_i');

set(hLine(3), 'ColorOrder', fluxColors, 'NextPlot', 'replacechildren');
irf_plot(hLine(3), [tEpmEpoch J_LEMS120], 'LineWidth', 0.70);
ylabel(hLine(3), 'J_i');

set(hLine(4), 'ColorOrder', fluxColors(1:4, :), 'NextPlot', 'replacechildren');
irf_plot(hLine(4), [tEpmEpoch J_DE], 'LineWidth', 0.75);
ylabel(hLine(4), 'J_e');

set(hLine(5), 'ColorOrder', fluxColors(1:4, :), 'NextPlot', 'replacechildren');
irf_plot(hLine(5), [tEpmEpoch J_LEFS60], 'LineWidth', 0.75);
ylabel(hLine(5), 'J_e');

set(hLine(6), 'ColorOrder', fluxColors(1:4, :), 'NextPlot', 'replacechildren');
irf_plot(hLine(6), [tEpmEpoch J_LEFS150], 'LineWidth', 0.75);
ylabel(hLine(6), 'J_e');

set(hLine(2:6), 'YScale', 'log');

% Set clean decade limits independently for each particle panel.
lineGroups = {J_LEMS30, J_LEMS120, J_DE, J_LEFS60, J_LEFS150};
for ii = 1:numel(lineGroups)
    values = lineGroups{ii};
    values = values(isfinite(values) & values > 0);
    if ~isempty(values)
        lowerLimit = 10 .^ floor(log10(min(values)));
        upperLimit = 10 .^ ceil(log10(max(values)));
        if upperLimit <= lowerLimit
            upperLimit = lowerLimit .* 10;
        end
        ylim(hLine(ii + 1), [lowerLimit upperLimit]);
    end
end

%% Line-panel labels and physical energy legends
lineLetters = 'abcdef';
for ii = 1:numel(hLine)
    text(hLine(ii), 0.018, 0.86, lineLetters(ii), 'Units', 'normalized', ...
        'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', 'k', 'Interpreter', 'tex');
end

text(hLine(2), 0.10, 0.88, 'LEMS30 ions', 'Units', 'normalized', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
lems30Labels = {'P_1: 46-67 keV','P_2: 67-115 keV', ...
    'P_3: 115-193 keV','P_4: 193-315 keV','P_5: 315-580 keV', ...
    'P_6: 0.58-1.06 MeV','P_7: 1.06-1.88 MeV','P_8: 1.88-4.70 MeV'};
legendX = [0.27 0.45 0.63 0.81];
for jj = 1:8
    legendRow = floor((jj - 1) ./ 4);
    legendColumn = mod(jj - 1, 4) + 1;
    text(hLine(2), legendX(legendColumn), 0.91 - 0.15 .* legendRow, ...
        lems30Labels{jj}, 'Units', 'normalized', 'FontSize', 6.5, ...
        'Color', fluxColors(jj, :), 'Interpreter', 'tex');
end
text(hLine(2), 0.10, 0.60, 'P_{1-6}: fill', ...
    'Units', 'normalized', 'FontSize', 7, ...
    'Color', [0.35 0.35 0.35], 'Interpreter', 'tex');

text(hLine(3), 0.10, 0.90, 'LEMS120 ions', 'Units', 'normalized', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
lems120Labels = {
    'P^{\prime}_1: 47-68 keV', 'P^{\prime}_2: 68-115 keV', ...
    'P^{\prime}_3: 115-195 keV', 'P^{\prime}_4: 195-321 keV', ...
    'P^{\prime}_5: 321-580 keV', 'P^{\prime}_6: 0.58-1.06 MeV', ...
    'P^{\prime}_7: 1.06-1.90 MeV', 'P^{\prime}_8: 1.90-4.80 MeV'};
for jj = 1:8
    legendRow = floor((jj - 1) ./ 4);
    legendColumn = mod(jj - 1, 4) + 1;
    text(hLine(3), legendX(legendColumn), 0.91 - 0.15 .* legendRow, ...
        lems120Labels{jj}, 'Units', 'normalized', 'FontSize', 6.5, ...
        'Color', fluxColors(jj, :), 'Interpreter', 'tex');
end

text(hLine(4), 0.10, 0.88, 'DE30 electrons', 'Units', 'normalized', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
deLabels = {'DE_1: 38-53 keV', 'DE_2: 53-103 keV', ...
    'DE_3: 103-175 keV', 'DE_4: 175-315 keV'};
legendX4 = [0.31 0.49 0.67 0.84];
for jj = 1:4
    text(hLine(4), legendX4(jj), 0.88, deLabels{jj}, ...
        'Units', 'normalized', 'FontSize', 7, ...
        'Color', fluxColors(jj, :), 'Interpreter', 'tex');
end
text(hLine(6), 0.78, 0.70, 'E_{1-3}: fill', 'Units', 'normalized', ...
    'FontSize', 7, 'Color', [0.35 0.35 0.35], 'Interpreter', 'tex');

text(hLine(5), 0.10, 0.88, 'LEFS60: e^- (+ions)', 'Units', 'normalized', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'tex');
lefs60Labels = {'E^{\prime}_1: 45-62 keV', 'E^{\prime}_2: 62-102 keV', ...
    'E^{\prime}_3: 102-175 keV', 'E^{\prime}_4: 175-290 keV'};
for jj = 1:4
    text(hLine(5), legendX4(jj), 0.88, lefs60Labels{jj}, ...
        'Units', 'normalized', 'FontSize', 7, ...
        'Color', fluxColors(jj, :), 'Interpreter', 'tex');
end

text(hLine(6), 0.10, 0.88, 'LEFS150: e^- (+ions)', 'Units', 'normalized', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'tex');
lefs150Labels = {'E_1: 44-58 keV','E_2: 58-104 keV', ...
    'E_3: 104-180 keV','E_4: 180-295 keV'};
for jj = 1:4
    text(hLine(6), legendX4(jj), 0.88, lefs150Labels{jj}, ...
        'Units', 'normalized', 'FontSize', 7, ...
        'Color', fluxColors(jj, :), 'Interpreter', 'tex');
end

%% Line-figure time axis and gap-free layout
irf_zoom(hLine, 'x', tint);
for ii = 1:numel(hLine)
    grid(hLine(ii), 'off');
    box(hLine(ii), 'on');
    set(hLine(ii), 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 0.8, ...
        'TickDir', 'in', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
        'XMinorTick', 'off', 'XTickLabelRotation', 0);
end
irf_timeaxis(hLine, 'date');
for ii = 1:numel(hLine)
    set(hLine(ii), 'XTick', hourTicks, 'XMinorTick', 'off');
end
for ii = 1:5
    set(hLine(ii), 'XTickLabel', []);
end
set(hLine(6), 'XTickLabel', hourTickLabels);
xlabel(hLine(6), '2024-05-10 to 12 UTC');

drawnow
linePanelLeft = 0.11;
linePanelRight = 0.965;
linePanelTop = 0.955;
linePanelBottom = 0.070;
linePanelGap = 0.0015;
linePanelHeight = (linePanelTop - linePanelBottom - ...
    linePanelGap .* (numel(hLine) - 1)) ./ numel(hLine);
for ii = 1:numel(hLine)
    linePanelY = linePanelTop - ii .* linePanelHeight - (ii - 1) .* linePanelGap;
    set(hLine(ii), 'Units', 'normalized', 'ActivePositionProperty', 'position', ...
        'Position', [linePanelLeft linePanelY ...
        linePanelRight - linePanelLeft linePanelHeight], ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off', 'XMinorTick', 'off');
end
drawnow

%% Save editable line-plot outputs
set(gcf,'paperpositionmode','auto')
lineFigFile = fullfile(baseDir, 'ace_20240510_12_epam_lines.fig');
linePngFile = fullfile(baseDir, 'ace_20240510_12_epam_lines.png');
linePdfFile = fullfile(baseDir, 'ace_20240510_12_epam_lines.pdf');
savefig(figLine, lineFigFile);
exportgraphics(figLine, linePngFile, 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(figLine, linePdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Saved EPAM line figure: %s\n', linePdfFile);
