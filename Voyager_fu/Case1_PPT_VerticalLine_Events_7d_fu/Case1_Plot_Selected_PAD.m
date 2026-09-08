function audit = Case1_Plot_Selected_PAD(padTable, selectedEpoch, eventID, figureFile, opts)
%Case1_Plot_Selected_PAD Draw explicitly selected hourly PAD records.
%   No peak search, neighboring-record replacement, averaging or fitting.
%   Each S1--S7 PAD uses its own Jmax; sigma uses the same fixed scale.
%   The input remains unchanged and the exact selected rows are audited.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-05

%% exact time selection
if nargin < 5, opts = struct; end
if ~isfield(opts, 'Visible'), opts.Visible = false; end
if ~isfield(opts, 'ExportDPI'), opts.ExportDPI = 220; end
if ~isfield(opts, 'Spacecraft'), opts.Spacecraft = 1; end
if ~isfield(opts, 'P1DisplayEnergyMeV'), opts.P1DisplayEnergyMeV = [0.57 1.78]; end
selectedEpoch = selectedEpoch(:);
assert(~isempty(selectedEpoch) && ~any(isnat(selectedEpoch)), 'Supply valid UTC epochs.');
assert(~isempty(selectedEpoch.TimeZone), 'Supply timezone-aware epochs.');
selectedEpoch.TimeZone = 'UTC';
assert(numel(unique(selectedEpoch)) == numel(selectedEpoch), 'Duplicate requested epochs.');
nPanel = numel(selectedEpoch);
rows = zeros(nPanel, 1);
for ii = 1:nPanel
    row = find(padTable.EpochUTC == selectedEpoch(ii));
    assert(isscalar(row), 'Each requested epoch must identify exactly one source row.');
    rows(ii) = row;
end
selected = padTable(rows, :);
assert(all(selected.PADUsable), 'A requested PAD is unavailable; do not substitute another time.');
flux = nan(nPanel, 7); sigma = flux; pa = flux;
for iSector = 1:7
    flux(:, iSector) = selected.(sprintf('RawFlux_S%d_1h', iSector));
    sigma(:, iSector) = selected.(sprintf('FluxUncertainty_S%d_1h', iSector));
    pa(:, iSector) = selected.(sprintf('PA_S%d_deg', iSector));
end
assert(all(isfinite(flux) & flux > 0, 'all') && all(isfinite(pa), 'all'), ...
    'All seven sectors must have positive flux and finite pitch angle.');
jmax = max(flux, [], 2);
normalizedFlux = flux ./ jmax;
normalizedSigma = sigma ./ jmax;

%% retain original payload and normalization
audit = struct;
audit.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
audit.EventID = string(eventID);
audit.SelectedEpochUTC = selectedEpoch;
audit.SelectedTableRows = rows;
audit.SourceTableUserData = padTable.Properties.UserData;
selected.Properties.UserData = [];
audit.SelectedRows = selected;
audit.RawFlux = flux;
audit.RawSigma = sigma;
audit.PA_deg = pa;
audit.NormalizationFlux = jmax;
audit.NormalizedFlux = normalizedFlux;
audit.NormalizedSigma = normalizedSigma;
audit.Sectors = 1:7;
audit.DisplayEnergyMeV = opts.P1DisplayEnergyMeV;
audit.SelectionPolicy = "Exact user-selected epochs in input order; no peak search or time substitution.";
audit.NormalizationPolicy = "Each row J/Jmax and sigma/Jmax; denominator treated as fixed; no additional uncertainty propagation.";
audit.ProcessingPolicy = "Existing hourly payload and full 3D PA retained. No new average, interpolation, angle merging, fitting, background subtraction or S8 points.";
audit.UncertaintyPolicy = "Finite nonnegative sigma shown; unavailable sigma retains its flux point without an error bar.";
audit.FigureFile = string(figureFile);

%% common axes; keep all finite error bars visible
yLow = 0; yHigh = 1;
for ii = 1:nPanel
    y = normalizedFlux(ii, :); dy = normalizedSigma(ii, :);
    good = isfinite(dy) & dy >= 0;
    yLow = min([yLow, y(good)-dy(good)]);
    yHigh = max([yHigh, y(good)+dy(good)]);
end
yLow = floor(yLow*5)/5;
% Leave enough room for the date, time and Jmax above the error bars.
yHigh = max(1.6, ceil((yHigh+0.30*(yHigh-yLow))*5)/5);
audit.DisplayYLimits = [yLow yHigh];
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
figureWidth = 340*nPanel;
left = 0.050; right = 0.015; gap = 0.017;
headerFontSize = 17; footnoteFontSize = 14;
if nPanel <= 2
    % Two-panel figures need room for the y label and end-point tick labels.
    figureWidth = max(800, figureWidth);
    left = 0.11; right = 0.04; gap = 0.055;
    headerFontSize = 15; footnoteFontSize = 11;
end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [80 80 figureWidth 610]);
figureCleanup = onCleanup(@() closeHiddenFigure(fig, opts.Visible));
panelAxes = gobjects(nPanel, 1);
for ii = 1:nPanel
    panelAxes(ii) = irf_subplot(1, nPanel, -ii);
end
panelWidth = (1-left-right-(nPanel-1)*gap)/nPanel;
for ii = 1:nPanel
    ax = panelAxes(ii);
    ax.Position = [left+(ii-1)*(panelWidth+gap), 0.16, panelWidth, 0.70];
    hold(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 16, ...
        'LineWidth', 1.3, 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.025], ...
        'XLim', [0 180], 'XTick', 0:45:180, 'YLim', [yLow yHigh], ...
        'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
    if yLow == 0 && yHigh <= 1.8, ax.YTick = 0:0.2:yHigh; end
    if ii > 1, ax.YTickLabel = []; end
    xlabel(ax, 'Pitch angle (\circ)', 'FontSize', 18, 'Interpreter', 'tex');
    if ii == 1, ylabel(ax, 'Normalized intensity', 'FontSize', 19); end
    title(ax, sprintf('(%c)', 'a'+ii-1), 'FontSize', 16, 'FontWeight', 'normal');
    dateText = char(datetime(selectedEpoch(ii), 'Format', 'yyyy-MM-dd'));
    timeText = [char(datetime(selectedEpoch(ii), 'Format', 'HH:mm')), ' UTC'];
    panelText(ax, {dateText; timeText}, [0.96 0.97], ...
        'FontSize', 15, 'Interpreter', 'none');
    panelText(ax, sprintf('J_{max} = %.4g', jmax(ii)), [0.96 0.835], ...
        'FontSize', 14, 'Interpreter', 'tex');
    x = pa(ii, :); y = normalizedFlux(ii, :); dy = normalizedSigma(ii, :);
    good = isfinite(dy) & dy >= 0;
    errorbar(ax, x(good), y(good), dy(good), ...
        'LineStyle', 'none', 'Marker', 'none', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 7);
    plot(ax, x, y, 'ko', 'LineStyle', 'none', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 6.5, 'LineWidth', 1);
end
header = sprintf('Voyager %d  |  %s  |  LECP P1 %.2f-%.2f MeV', ...
    opts.Spacecraft, char(eventID), opts.P1DisplayEnergyMeV);
annotation(fig, 'textbox', [0.05 0.945 0.93 0.047], 'String', header, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', headerFontSize, 'Interpreter', 'none');
annotation(fig, 'textbox', [0.08 0.025 0.89 0.040], ...
    'String', 'J_{max} in cm^{-2} s^{-1} sr^{-1} MeV^{-1}; each panel normalized by its own J_{max}', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', footnoteFontSize, 'Interpreter', 'tex');
folder = fileparts(figureFile);
if ~isfolder(folder), mkdir(folder); end
exportgraphics(fig, char(figureFile), 'Resolution', opts.ExportDPI, 'BackgroundColor', 'white');
audit.FigureCreated = true;
end

function panelText(ax, labels, position, varargin)
ht = irf_legend(ax, labels, position, varargin{:});
set(ht, 'Color', 'k', 'FontName', 'Times New Roman');
end

function closeHiddenFigure(fig, visible)
if ~logical(visible) && isgraphics(fig), close(fig); end
end
