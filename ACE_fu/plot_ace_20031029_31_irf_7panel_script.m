%% ACE 2003 Halloween storm, IRF-MATLAB 7-panel script
% Run this script directly in MATLAB.  It is intentionally a script, not a
% function, so that all variables remain in the base workspace for inspection.

clearvars
close all
clc

%% Paths
baseDir = 'C:\Users\Administrator\Documents\FWD_matlab\ACE_fu';
irfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';
if exist(irfDir, 'dir')
    addpath(genpath(irfDir));
end

if exist('irf_plot', 'file') ~= 2
    error('irf_matlab is not on the MATLAB path. Please add irfu-matlab-master first.');
end

cacheDir = fullfile(baseDir, 'hapi_cache');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end

helper = fullfile(baseDir, 'hapi_csv_to_mat_v4.py');
bundledPython = fullfile(getenv('USERPROFILE'), '.cache', 'codex-runtimes', ...
    'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
if exist(bundledPython, 'file'), pyExe = bundledPython; else, pyExe = 'python'; end

plot0 = datenum(2003, 10, 29, 0, 0, 0);
plot1 = datenum(2003, 10, 31, 0, 0, 0);
eventTag = '20031029T0000_20031031T0000';

mfiCsv = fullfile(cacheDir, ['AC_H3_MFI_' eventTag '_Magnitude_BGSEc.csv']);
mfiFile = fullfile(cacheDir, ['AC_H3_MFI_' eventTag '_Magnitude_BGSEc.mat']);
sweCsv = fullfile(cacheDir, ['AC_H0_SWE_' eventTag '_Vp_Tpr.csv']);
sweFile = fullfile(cacheDir, ['AC_H0_SWE_' eventTag '_Vp_Tpr.mat']);
epmCsv = fullfile(cacheDir, ['AC_H3_EPM_' eventTag '_all_telescopes.csv']);
epmFile = fullfile(cacheDir, ['AC_H3_EPM_' eventTag '_all_telescopes.mat']);

mfiUrl = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_MFI', ...
    '&parameters=Magnitude,BGSEc&time.min=2003-10-29T00:00:00Z', ...
    '&time.max=2003-10-31T00:00:00Z&format=csv'];
sweUrl = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H0_SWE', ...
    '&parameters=Vp,Tpr&time.min=2003-10-29T00:00:00Z', ...
    '&time.max=2003-10-31T00:00:00Z&format=csv'];
epmParams = ['P1,P2,P3,P4,P5,P6,P7,P8,DE1,DE2,DE3,DE4,', ...
    'E1p,E2p,E3p,E4p,P1p,P2p,P3p,P4p,P5p,P6p,P7p,P8p,E1,E2,E3,E4'];
epmUrl = ['https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_EPM', ...
    '&parameters=' epmParams '&time.min=2003-10-29T00:00:00Z', ...
    '&time.max=2003-10-31T00:00:00Z&format=csv'];

cacheJobs = {mfiCsv,mfiFile,mfiUrl; sweCsv,sweFile,sweUrl; ...
    epmCsv,epmFile,epmUrl};
for iJob = 1:size(cacheJobs, 1)
    if ~exist(cacheJobs{iJob, 2}, 'file')
        cmd = sprintf('"%s" "%s" "%s" "%s" "%s"', pyExe, helper, ...
            cacheJobs{iJob, 1}, cacheJobs{iJob, 2}, cacheJobs{iJob, 3});
        [status, cmdout] = system(cmd);
        if status ~= 0, error('HAPI cache preparation failed:\n%s', cmdout); end
    end
end

%% Load cached HAPI data
mfi = load(mfiFile);  % t: MATLAB datenum, y: [|B| Bx By Bz]
swe = load(sweFile);  % AC_H0_SWE Level-2 64-s data: [Vp Tpr(K)]
epm = load(epmFile);  % 28 EPAM channels in dataset order

tMfi = double(mfi.t(:));
yMfi = double(mfi.y);
tSwe = double(swe.t(:));
ySwe = double(swe.y);
tEpm = double(epm.t(:));
yEpm = double(epm.y);

yMfi(abs(yMfi) > 1e30) = NaN;
ySwe(abs(ySwe) > 1e30) = NaN;
yEpm(abs(yEpm) > 1e30) = NaN;

epoch0 = datenum(1970, 1, 1, 0, 0, 0);
tMfiEpoch = (tMfi - epoch0) .* 86400;
tSweEpoch = (tSwe - epoch0) .* 86400;
tEpmEpoch = (tEpm - epoch0) .* 86400;

tint = ([plot0 plot1] - epoch0) .* 86400;

%% Variables
Bmag = yMfi(:, 1);
BGSE = yMfi(:, 2:4);
Vp = ySwe(:, 1);
Ti_eV = ySwe(:, 2) ./ 11604.51812;  % AC_H0_SWE:Tpr unit is Kelvin

%% SWE availability report for panels b/c
sweInPlot = tSwe >= plot0 & tSwe <= plot1;
sweGood = sweInPlot & isfinite(Vp) & isfinite(Ti_eV);
fprintf('SWE source : %s\n', sweFile);
fprintf('SWE points : %d samples in event interval, %d valid Vp/Ti samples\n', ...
    sum(sweInPlot), sum(sweGood));
if any(sweGood)
    tSweGood = tSwe(sweGood);
    sweGapsMin = diff(tSweGood) .* 1440;
    sweGapIndex = find(sweGapsMin > 15);
    for iiGap = 1:numel(sweGapIndex)
        fprintf('SWE gap    : %s - %s (%.1f min)\n', ...
            datestr(tSweGood(sweGapIndex(iiGap)), 'HH:MM:SS'), ...
            datestr(tSweGood(sweGapIndex(iiGap) + 1), 'HH:MM:SS'), ...
            sweGapsMin(sweGapIndex(iiGap)));
    end
end

% Use LEMS120 for ions. LEMS30 P1-P6 became unreliable around DOY 302,
% exactly when the Halloween storm reached ACE.
J_i = yEpm(:, 17:24);
J_e = yEpm(:, 9:12);
J_i(J_i <= 0 | abs(J_i) > 1e30) = NaN;
J_e(J_e <= 0 | abs(J_e) > 1e30) = NaN;

inPlot = tMfi >= plot0 & tMfi <= plot1;
[~, iMaxRel] = max(Bmag(inPlot), [], 'omitnan');
iAll = find(inPlot);
if isempty(iMaxRel) || isempty(iAll)
    tBmax = plot0;
else
    tBmax = tMfi(iAll(iMaxRel));
end
tBmaxEpoch = (tBmax - epoch0) .* 86400;

%% Prepare spectrogram records for irf_spectrogram
E_e_edges = [38 53 103 175 315] .* 1e3;
E_i_edges = [47 68 115 195 321 580 1060 1900 4800] .* 1e3;
E_e = sqrt(E_e_edges(1:end-1) .* E_e_edges(2:end));
E_i = sqrt(E_i_edges(1:end-1) .* E_i_edges(2:end));
dE_e = (E_e_edges(2:end) - E_e_edges(1:end-1)) ./ 2;
dE_i = (E_i_edges(2:end) - E_i_edges(1:end-1)) ./ 2;

dtEpm = median(diff(tEpmEpoch), 'omitnan') ./ 2;
if ~isfinite(dtEpm) || dtEpm <= 0
    dtEpm = 6;
end

specrec_e.t = tEpmEpoch;
specrec_e.f = E_e;
specrec_e.p = J_e;
specrec_e.dt = dtEpm;
specrec_e.df = dE_e;
specrec_e.f_label = 'E_e [eV]';
specrec_e.p_label = 'log_{10} J_e';

specrec_i.t = tEpmEpoch;
specrec_i.f = E_i;
specrec_i.p = J_i;
specrec_i.dt = dtEpm;
specrec_i.df = dE_i;
specrec_i.f_label = 'E_i [eV]';
specrec_i.p_label = 'log_{10} J_i';

%% Figure and IRF panels
h = irf_plot(7, 'newfigure');
fig = gcf;
set(fig, 'Color', 'w', 'Name', 'ACE 2003 Halloween storm IRF 7-panel', ...
    'Units', 'pixels', 'Position', [80 30 980 1320]);
figUserData = get(fig, 'UserData');
if ~isstruct(figUserData)
    figUserData = struct;
end
figUserData.t_start_epoch = tint(1);
set(fig, 'UserData', figUserData);

blue = [0.12 0.26 0.78];
green = [0.32 0.72 0.00];
red = [0.93 0.12 0.08];
cyan = [0.00 0.53 0.62];
purple = [0.70 0.24 0.70];
brown = [0.65 0.43 0.00];
gray = [0.25 0.25 0.25];
orange = [0.92 0.53 0.10];
fluxColors = [blue; green; red; cyan; purple; brown; gray; orange];

%% a) Magnetic field
set(h(1), 'ColorOrder', [blue; green; red; 0 0 0], 'NextPlot', 'replacechildren');
irf_plot(h(1), [tMfiEpoch BGSE Bmag], 'LineWidth', 0.75);
ylabel(h(1), 'B [nT]');
set(h(1), 'ColorOrder', [blue; green; red; 0 0 0]);
irf_legend(h(1), {'B_x','B_y','B_z','|B|'}, [0.74 0.90]);
irf_legend(h(1), 'a', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

%% b) Proton speed
set(h(2), 'ColorOrder', blue, 'NextPlot', 'replacechildren');
irf_plot(h(2), [tSweEpoch Vp], 'LineWidth', 0.9);
ylabel(h(2), 'V_p [km/s]');
irf_legend(h(2), {'V_p'}, [0.94 0.88], 'color', blue);
irf_legend(h(2), 'b', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
if ~any(sweInPlot & isfinite(Vp))
    text(h(2), 0.50, 0.50, ...
        {'ACE/SWEPAM L2 unavailable', 'SEP contamination during Halloween storm'}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.35 0.35 0.35], 'FontSize', 10);
end

%% c) Proton temperature in eV
set(h(3), 'ColorOrder', red, 'NextPlot', 'replacechildren');
irf_plot(h(3), [tSweEpoch Ti_eV], 'LineWidth', 0.9);
ylabel(h(3), 'T_i [eV]');
irf_legend(h(3), {'T_i'}, [0.94 0.88], 'color', red);
irf_legend(h(3), 'c', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
if ~any(sweInPlot & isfinite(Ti_eV))
    text(h(3), 0.50, 0.50, ...
        'No standard SWEPAM proton temperature', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.35 0.35 0.35], 'FontSize', 10);
end

%% d/e) Particle flux spectrograms
irf_spectrogram(h(4), specrec_e, 'log');
set(h(4), 'YScale', 'log', 'Layer', 'top');
ylabel(h(4), 'E_e [eV]');
irf_legend(h(4), 'd', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

irf_spectrogram(h(5), specrec_i, 'log');
set(h(5), 'YScale', 'log', 'Layer', 'top');
ylabel(h(5), 'E_i [eV]');
irf_legend(h(5), 'e', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

%% f/g) Line-style particle flux panels, plotted as log10(J) with IRF_PLOT
logJe = log10(J_e);
logJi = log10(J_i);

set(h(6), 'ColorOrder', fluxColors(1:4, :), 'NextPlot', 'replacechildren');
irf_plot(h(6), [tEpmEpoch logJe], 'LineWidth', 0.75);
ylabel(h(6), 'J_e');
set(h(6), 'ColorOrder', fluxColors(1:4, :));
irf_legend(h(6), {'e^- 38-53 keV','e^- 53-103 keV','e^- 103-175 keV','e^- 175-315 keV'}, ...
    [0.58 0.90], 'FontSize', 8);
irf_legend(h(6), 'f', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

set(h(7), 'ColorOrder', fluxColors(1:8, :), 'NextPlot', 'replacechildren');
irf_plot(h(7), [tEpmEpoch logJi], 'LineWidth', 0.75);
ylabel(h(7), 'J_i');
set(h(7), 'ColorOrder', fluxColors(1:8, :));
irf_legend(h(7), {'i^+ 47-68 keV','i^+ 68-115 keV','i^+ 115-195 keV','i^+ 195-321 keV', ...
    'i^+ 321-580 keV','i^+ 0.58-1.06 MeV','i^+ 1.06-1.90 MeV','i^+ 1.90-4.80 MeV'}, ...
    [0.58 0.90], 'FontSize', 7);
irf_legend(h(7), 'g', [0.02 0.88], 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

%% Axes styling, Bmax marker, time axis
for ii = 1:numel(h)
    grid(h(ii), 'off');
    box(h(ii), 'on');
    set(h(ii), 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 0.8, ...
        'TickDir', 'in', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    % irf_pl_mark(h(ii), tBmaxEpoch, 'k', 'LineStyle', '--', 'LineWidth', 0.8);
end

irf_zoom(h, 'x', tint);
irf_timeaxis(h, 'date');
eventTicks = 0:6*3600:(tint(2) - tint(1));
eventTickLabels = cellstr(datestr(plot0 + eventTicks(:) ./ 86400, 'mm-dd HH:MM'));
for ii = 1:numel(h), set(h(ii), 'XTick', eventTicks, 'XMinorTick', 'off'); end
for ii = 1:numel(h)-1, set(h(ii), 'XTickLabel', []); end
set(h(end), 'XTickLabel', eventTickLabels);
xlabel(h(end), '2003-10-29 to 31 UTC');

% Flux-line panels are plotted in log10(J).  Relabel the y ticks as 10^n.
for ii = [6 7]
    yl = ylim(h(ii));
    yt = ceil(yl(1)):floor(yl(2));
    if numel(yt) < 2
        yt = round(linspace(yl(1), yl(2), 3));
    end
    set(h(ii), 'YTick', yt);
    ytxt = cell(size(yt));
    for jj = 1:numel(yt)
        ytxt{jj} = sprintf('10^{%d}', yt(jj));
    end
    set(h(ii), 'YTickLabel', ytxt);
end

% A compact color map close to the previous generated figure.
fluxMapStops = [46 45 95; 30 115 190; 98 190 185; 240 230 70; 250 125 25; 185 25 25] ./ 255;
fluxMapX = linspace(0, 1, size(fluxMapStops, 1));
fluxMapXi = linspace(0, 1, 256);
colormap(fig, interp1(fluxMapX, fluxMapStops, fluxMapXi, 'linear'));

% Color limits from robust percentiles of log10 flux values.
tmp = logJe(isfinite(logJe));
if ~isempty(tmp)
    tmp = sort(tmp(:));
    n = numel(tmp);
    lo = tmp(max(1, round(0.05 .* n)));
    hi = tmp(min(n, round(0.95 .* n)));
    caxis(h(4), [lo hi]);
end
tmp = logJi(isfinite(logJi));
if ~isempty(tmp)
    tmp = sort(tmp(:));
    n = numel(tmp);
    lo = tmp(max(1, round(0.05 .* n)));
    hi = tmp(min(n, round(0.95 .* n)));
    caxis(h(5), [lo hi]);
end

% Lock the final layout after IRF spectrogram/colorbar commands.  This keeps
% the two spectrogram panels aligned with the line panels while reserving a
% consistent strip on the right for colorbars.
drawnow
panelLeft = 0.12;
panelRight = 0.80;
panelTop = 0.93;
panelBottom = 0.09;
panelGap = 0.002;
panelHeight = (panelTop - panelBottom - panelGap .* (numel(h) - 1)) ./ numel(h);
for ii = 1:numel(h)
    panelY = panelTop - ii .* panelHeight - (ii - 1) .* panelGap;
    set(h(ii), 'Units', 'normalized', 'ActivePositionProperty', 'position', ...
        'Position', [panelLeft panelY panelRight - panelLeft panelHeight], ...
        'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
end

cb = findall(fig, 'Type', 'colorbar');
if numel(cb) >= 2
    cbY = zeros(numel(cb), 1);
    for ii = 1:numel(cb)
        cbPos = get(cb(ii), 'Position');
        cbY(ii) = cbPos(2);
    end
    [~, cbOrder] = sort(cbY, 'descend');
    cb = cb(cbOrder);
    cbLeft = panelRight + 0.025;
    cbWidth = 0.018;
    p4 = get(h(4), 'Position');
    p5 = get(h(5), 'Position');
    set(cb(1), 'Units', 'normalized', ...
        'Position', [cbLeft p4(2) + 0.08 .* p4(4) cbWidth 0.84 .* p4(4)]);
    set(cb(2), 'Units', 'normalized', ...
        'Position', [cbLeft p5(2) + 0.08 .* p5(4) cbWidth 0.84 .* p5(4)]);
end

drawnow

%% Save editable outputs
outFig = fullfile(baseDir, 'ace_20031029_31_solarwind_overview.fig');
outPng = fullfile(baseDir, 'ace_20031029_31_solarwind_overview.png');
outPdf = fullfile(baseDir, 'ace_20031029_31_solarwind_overview.pdf');

savefig(fig, outFig);
exportgraphics(fig, outPng, 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Saved overview: %s\n', outPdf);
set(gca,"XTickLabelRotation",0)
set(gcf,'Renderer','painters');
set(gcf,'paperpositionmode','auto')
