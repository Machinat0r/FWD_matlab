function out = plot_ace_20260119_matlab_7panel(baseDir)
%PLOT_ACE_20260119_MATLAB_7PANEL Draw the ACE 2026-01-19 overview in MATLAB.
%
% This standalone MATLAB renderer reads the local HAPI MAT caches generated
% in ACE_fu/hapi_cache and produces a 7-panel editable figure:
%   a) Bx, By, Bz, |B| in GSE
%   b) scalar proton speed Vp
%   c) proton temperature converted from Kelvin to eV
%   d) electron flux spectrogram
%   e) ion flux spectrogram
%   f) electron flux line plot
%   g) ion flux line plot
%
% Usage:
%   out = plot_ace_20260119_matlab_7panel
%   out = plot_ace_20260119_matlab_7panel('C:\Users\Administrator\Documents\FWD_matlab\ACE_fu')

if nargin < 1 || isempty(baseDir)
    baseDir = fileparts(mfilename('fullpath'));
end

cacheDir = fullfile(baseDir, 'hapi_cache');
mfiFile = fullfile(cacheDir, 'AC_H3_MFI_20260119_Magnitude_BGSEc.mat');
sweFile = fullfile(cacheDir, 'AC_K0_SWE_20260119_Vp_Tpr.mat');
epmFile = fullfile(cacheDir, 'AC_H3_EPM_20260119T1200_20260120T0000_P1_P2_P3_P4_P5_P6_P7_P8_DE1_DE2_DE3_DE4.mat');

plot0 = datenum(2026, 1, 19, 12, 0, 0);
plot1 = datenum(2026, 1, 20, 0, 0, 0);

mfi = loadNumericCache(mfiFile);
swe = loadNumericCache(sweFile);
epm = loadNumericCache(epmFile);

tm = mfi.t;
ym = mfi.y;
ts = swe.t;
ys = swe.y;
te = epm.t;
ye = epm.y;

bmag = ym(:, 1);
bgse = ym(:, 2:4);
vp = ys(:, 1);
tiEv = ys(:, 2) ./ 11604.51812;

ions = ye(:, 1:8);
electrons = ye(:, 9:12);
ions(ions <= 0 | abs(ions) > 1e30) = NaN;
electrons(electrons <= 0 | abs(electrons) > 1e30) = NaN;

inPlot = tm >= plot0 & tm <= plot1;
[~, imaxRel] = max(bmag(inPlot), [], 'omitnan');
idx = find(inPlot);
if isempty(imaxRel) || isempty(idx)
    tBmax = datenum(2026, 1, 19, 18, 59, 38.687);
else
    tBmax = tm(idx(imaxRel));
end

fig = figure('Color', 'w', 'Name', 'ACE 2026-01-19 overview 7-panel', ...
    'Units', 'pixels', 'Position', [80 40 980 1320]);
tl = tiledlayout(fig, 7, 1, 'TileSpacing', 'none', 'Padding', 'compact');
ax = gobjects(7, 1);

cc = componentColors();
fluxC = fluxColors();

ax(1) = nexttile(tl);
plot(ax(1), tm, bgse(:, 1), 'Color', cc(1, :), 'LineWidth', 0.75); hold(ax(1), 'on');
plot(ax(1), tm, bgse(:, 2), 'Color', cc(2, :), 'LineWidth', 0.75);
plot(ax(1), tm, bgse(:, 3), 'Color', cc(3, :), 'LineWidth', 0.75);
plot(ax(1), tm, bmag, 'k', 'LineWidth', 0.85);
ylabel(ax(1), 'B [nT]');
legend(ax(1), {'B_x', 'B_y', 'B_z', '|B|'}, 'Location', 'northeast', ...
    'Orientation', 'horizontal', 'Box', 'off', 'Interpreter', 'tex');
panelLabel(ax(1), 'a');

ax(2) = nexttile(tl);
plot(ax(2), ts, vp, 'Color', cc(1, :), 'LineWidth', 0.9);
ylabel(ax(2), 'V_p [km/s]');
legend(ax(2), {'V_p'}, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');
panelLabel(ax(2), 'b');

ax(3) = nexttile(tl);
plot(ax(3), ts, tiEv, 'Color', cc(3, :), 'LineWidth', 0.9);
ylabel(ax(3), 'T_i [eV]');
legend(ax(3), {'T_i'}, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');
panelLabel(ax(3), 'c');

ax(4) = nexttile(tl);
plotFluxSpectrogram(ax(4), te, electrons, [38 53 103 175 315] .* 1e3, 'E_e [eV]', 'log_{10} J_e');
panelLabel(ax(4), 'd');

ax(5) = nexttile(tl);
plotFluxSpectrogram(ax(5), te, ions, [46 68 115 195 321 580 1060 1900 4800] .* 1e3, 'E_i [eV]', 'log_{10} J_i');
panelLabel(ax(5), 'e');

ax(6) = nexttile(tl);
semilogy(ax(6), te, electrons, 'LineWidth', 0.75); hold(ax(6), 'on');
applyLineColors(ax(6), fluxC(1:4, :));
ylabel(ax(6), 'J_e');
legend(ax(6), {'e^- 38-53 keV', 'e^- 53-103 keV', 'e^- 103-175 keV', 'e^- 175-315 keV'}, ...
    'Location', 'northeast', 'Box', 'on', 'NumColumns', 2, 'Interpreter', 'tex', 'FontSize', 8);
panelLabel(ax(6), 'f');

ax(7) = nexttile(tl);
semilogy(ax(7), te, ions, 'LineWidth', 0.75); hold(ax(7), 'on');
applyLineColors(ax(7), fluxC(1:8, :));
ylabel(ax(7), 'J_i');
legend(ax(7), {'H^+ 46-68 keV', 'H^+ 67-115 keV', 'H^+ 115-195 keV', 'H^+ 193-321 keV', ...
    'H^+ 315-580 keV', 'H^+ 0.58-1.06 MeV', 'H^+ 1.06-1.90 MeV', 'H^+ 1.88-4.80 MeV'}, ...
    'Location', 'northeast', 'Box', 'on', 'NumColumns', 2, 'Interpreter', 'tex', 'FontSize', 7);
ylabel(ax(7), 'J_i');
xlabel(ax(7), '2026-01-19/20 UTC');
panelLabel(ax(7), 'g');

for ii = 1:numel(ax)
    styleAxis(ax(ii), plot0, plot1);
    addBmaxLine(ax(ii), tBmax);
end
linkaxes(ax, 'x');
formatTimeAxes(ax, plot0, plot1);

out.fig = fullfile(baseDir, 'ace_20260119_overview_matlab.fig');
out.png = fullfile(baseDir, 'ace_20260119_overview_matlab.png');
out.pdf = fullfile(baseDir, 'ace_20260119_overview_matlab.pdf');
out.bmax_time = tBmax;

drawnow;
savefig(fig, out.fig);
exportgraphics(fig, out.png, 'Resolution', 300);
exportgraphics(fig, out.pdf, 'ContentType', 'vector');

fprintf('Saved FIG: %s\n', out.fig);
fprintf('Saved PNG: %s\n', out.png);
fprintf('Saved PDF: %s\n', out.pdf);
fprintf('Bmax UTC : %s\n', datestr(tBmax, 'yyyy-mm-dd HH:MM:SS.FFF'));
end

function s = loadNumericCache(file)
if ~exist(file, 'file')
    error('Required cache file was not found: %s', file);
end
s = load(file);
s.t = double(s.t(:));
s.y = double(s.y);
s.y(abs(s.y) > 1e30) = NaN;
end

function plotFluxSpectrogram(ax, t, flux, energyEdgesEv, ylab, cblab)
flux(flux <= 0 | abs(flux) > 1e30) = NaN;
logFlux = log10(flux);
tEdges = timeEdges(t);
c = NaN(size(logFlux, 2) + 1, size(logFlux, 1) + 1);
c(1:size(logFlux, 2), 1:size(logFlux, 1)) = logFlux.';
pcolor(ax, tEdges, energyEdgesEv, c);
shading(ax, 'flat');
set(ax, 'YScale', 'log', 'YDir', 'normal', 'Layer', 'top');
ylabel(ax, ylab);
yt = 10 .^ (ceil(log10(min(energyEdgesEv))) : floor(log10(max(energyEdgesEv))));
set(ax, 'YTick', yt);
colormap(ax, aceFluxColormap(256));
vals = logFlux(isfinite(logFlux));
if ~isempty(vals)
    clim(ax, percentileLocal(vals, [5 95]));
end
cb = colorbar(ax, 'east');
cb.Label.String = cblab;
cb.Label.Interpreter = 'tex';
cb.Color = [1 1 1];
end

function e = timeEdges(t)
t = t(:);
if numel(t) < 2
    e = [t - 1/1440; t + 1/1440];
    return
end
mid = (t(1:end-1) + t(2:end)) ./ 2;
e = [t(1) - (mid(1) - t(1)); mid; t(end) + (t(end) - mid(end))];
end

function styleAxis(ax, plot0, plot1)
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [plot0 plot1]);
set(ax, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 0.8, ...
    'XTick', plot0:3/24:plot1, 'TickDir', 'in');
end

function formatTimeAxes(ax, plot0, plot1)
ticks = plot0:3/24:plot1;
for ii = 1:numel(ax)
    set(ax(ii), 'XLim', [plot0 plot1], 'XTick', ticks);
    if ii < numel(ax)
        set(ax(ii), 'XTickLabel', []);
    else
        labels = cell(size(ticks));
        for jj = 1:numel(ticks)
            labels{jj} = datestr(ticks(jj), 'HH:MM');
        end
        set(ax(ii), 'XTickLabel', labels);
    end
end
end

function addBmaxLine(ax, tBmax)
yl = ylim(ax);
plot(ax, [tBmax tBmax], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylim(ax, yl);
end

function panelLabel(ax, label)
text(ax, 0.018, 0.88, label, 'Units', 'normalized', 'FontSize', 18, ...
    'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
end

function c = componentColors()
c = [0.12 0.26 0.78; 0.32 0.72 0.00; 0.93 0.12 0.08];
end

function c = fluxColors()
c = [0.12 0.26 0.78; 0.32 0.72 0.00; 0.93 0.12 0.08; ...
    0.00 0.53 0.62; 0.70 0.24 0.70; 0.65 0.43 0.00; ...
    0.25 0.25 0.25; 0.92 0.53 0.10];
end

function applyLineColors(ax, colors)
h = findobj(ax, 'Type', 'line');
h = flipud(h(:));
for ii = 1:min(numel(h), size(colors, 1))
    set(h(ii), 'Color', colors(ii, :));
end
end

function cmap = aceFluxColormap(n)
stops = [46 45 95; 30 115 190; 98 190 185; 240 230 70; 250 125 25; 185 25 25] ./ 255;
x = linspace(0, 1, size(stops, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, stops, xi, 'linear');
end

function q = percentileLocal(x, p)
x = sort(x(:));
x = x(isfinite(x));
if isempty(x)
    q = [0 1];
    return
end
n = numel(x);
q = zeros(size(p));
for ii = 1:numel(p)
    pos = 1 + (n - 1) .* p(ii) ./ 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(ii) = x(lo);
    else
        q(ii) = x(lo) + (x(hi) - x(lo)) .* (pos - lo);
    end
end
end
