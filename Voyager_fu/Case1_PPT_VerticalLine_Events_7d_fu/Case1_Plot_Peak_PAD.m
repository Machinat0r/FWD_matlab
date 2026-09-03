function peakAudit = Case1_Plot_Peak_PAD(padTable, eventRow, figureFile, opts)
%Case1_Plot_Peak_PAD Five source-identified PADs around the flux maximum.
%   The input is the in-memory PAD table produced from the original CDFs.
%   Peak = largest positive S1--S7 sector flux in the event plotting window.
%   Ties use the earliest Epoch. The peak remains selected if its PAD is
%   unavailable. On either side, select the nearest two complete PAD records;
%   missing records are skipped without interpolation or time relabelling.
%   Each panel uses J/Jmax and sigma/Jmax, treating Jmax as a fixed scale.
%   S8, angle merging, fitting and new quality thresholds are not used.
%   opts.RenderFigure=false runs selection/auditing without creating a figure.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% options and audit structure
if nargin < 4, opts = struct; end
if ~isfield(opts, 'PADCadence'), opts.PADCadence = 'hour'; end
if ~isfield(opts, 'P1DisplayEnergyMeV'), opts.P1DisplayEnergyMeV = [0.57 1.78]; end
if ~isfield(opts, 'Visible'), opts.Visible = false; end
if ~isfield(opts, 'ExportDPI'), opts.ExportDPI = 220; end
if ~isfield(opts, 'RenderFigure'), opts.RenderFigure = true; end
if ~isfield(opts, 'Spacecraft'), opts.Spacecraft = 1; end
if strcmpi(opts.PADCadence, 'hour')
    suffix = '1h';
elseif strcmpi(opts.PADCadence, 'day')
    suffix = '1d';
else
    error('Case1:PeakPADCadence', 'PADCadence must be hour or day.');
end
if ~istable(padTable)
    error('Case1:PeakPADTable', 'The PAD input must be a table.');
end
eventID = string(eventRow.EventID);
if numel(eventID) ~= 1
    error('Case1:PeakPADEvent', 'Supply one event row.');
end
startTime = eventRow.PlotStartUTC;
endTime = eventRow.PlotEndUTCExclusive;
peakAudit = struct;
peakAudit.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
peakAudit.EventID = eventID;
peakAudit.Spacecraft = opts.Spacecraft;
peakAudit.PADCadence = string(opts.PADCadence);
peakAudit.WindowStartUTC = startTime;
peakAudit.WindowEndUTCExclusive = endTime;
peakAudit.Status = "no_retained_records";
peakAudit.PeakTableRow = NaN;
peakAudit.PeakEpochUTC = NaT(1, 1, 'TimeZone', 'UTC');
peakAudit.PeakFlux = NaN;
peakAudit.PeakPADUsable = false;
peakAudit.SelectedTableRows = nan(5, 1);
peakAudit.SelectedCount = 0;
peakAudit.PlottedCount = 0;
peakAudit.PanelOffsets = (-2:2).';
peakAudit.SelectedEpochUTC = NaT(5, 1, 'TimeZone', 'UTC');
peakAudit.SelectedPADUsable = false(5, 1);
peakAudit.SelectedRows = cell(5, 1);
peakAudit.RawFlux = nan(5, 7);
peakAudit.RawSigma = nan(5, 7);
peakAudit.PA_deg = nan(5, 7);
peakAudit.NormalizedFlux = nan(5, 7);
peakAudit.NormalizedSigma = nan(5, 7);
peakAudit.NormalizationFlux = nan(5, 1);
peakAudit.NormalizingSectors = cell(5, 1);
peakAudit.SourceEnergyMeV = nan(5, 2);
peakAudit.DisplayEnergyMeV = reshape(opts.P1DisplayEnergyMeV, 1, 2);
peakAudit.InputTableRowCount = height(padTable);
peakAudit.CandidateTableRows = zeros(0, 1);
peakAudit.RowMaximumSectorFlux = nan(height(padTable), 1);
peakAudit.SkippedMissingBeforeCount = 0;
peakAudit.SkippedMissingAfterCount = 0;
peakAudit.SkippedMissingBeforeTableRows = zeros(0, 1);
peakAudit.SkippedMissingAfterTableRows = zeros(0, 1);
peakAudit.Sectors = 1:7;
peakAudit.FluxUnits = "cm^-2 s^-1 sr^-1 MeV^-1";
peakAudit.PeakPolicy = "Maximum positive finite raw S1-S7 sector flux within [window start,window end); earliest Epoch on ties; peak retained even if PAD is incomplete.";
peakAudit.NeighborPolicy = "Nearest two complete S1-S7 PAD records strictly before and after the peak, within the same window; missing records skipped; absent neighbors left blank.";
peakAudit.NormalizationPolicy = "Each selected record divided by its own maximum positive S1-S7 flux; uncertainty divided by the same fixed scale; denominator uncertainty not additionally propagated.";
peakAudit.FluxPolicy = "Original source CDF flux and sigma, original Epoch, after the approved DeltaT<0 exclusion; no averaging, interpolation, angle merging, fitting or S8 subtraction here.";
peakAudit.EnergyLabelPolicy = "User-approved 0.57-1.78 MeV display convention; source CDF energy bounds and flux values are preserved without recalibration.";
peakAudit.UncertaintyDisplayPolicy = "Finite nonnegative uncertainties are displayed; a point with unavailable uncertainty is retained without an error bar. Common y limits contain all displayed bars.";
peakAudit.FigureFile = string(figureFile);
peakAudit.FigureCreated = false;
peakAudit.DisplayYLimits = [0 1.2];
peakAudit.SourceTableUserData = padTable.Properties.UserData;
if ismember('SourceProduct',padTable.Properties.VariableNames) && ...
        any(padTable.SourceProduct == "L1_UTC_mean")
    peakAudit.FluxPolicy = "Input uses unchanged L2 records where complete, plus approved L1 UTC-bin rate means where L2 is incomplete/absent; see each SelectedRows.SourceProduct and SourceTableUserData.L1FallbackAudit. This renderer adds no averaging, interpolation or fitting.";
    peakAudit.EnergyLabelPolicy = "Display 0.57-1.78 MeV; source CDF bounds unchanged. L1 separately uses the approved historical 0.44*(1.78-0.57) conversion; no L2 rescaling.";
    peakAudit.PeakPolicy = "Maximum positive S1-S7 flux in the L2-first/L1-fallback table; earliest time on ties; L1 means use UTC bin center. The selected peak stays selected when its PA is unavailable.";
    userData = padTable.Properties.UserData;
    if isstruct(userData) && isfield(userData,'L1FallbackAudit') && ...
            isfield(userData.L1FallbackAudit,'SourcePriority')
        peakAudit.SourcePriority = string(userData.L1FallbackAudit.SourcePriority);
        peakAudit.FluxPolicy = "Input source priority=" + peakAudit.SourcePriority + ...
            "; complete UTC bins from the selected level; L1 rate means at bin centers, retained L2 payload unchanged. Source records and replaced L2 are in L1FallbackAudit. No additional average/interpolation in this renderer.";
        peakAudit.PeakPolicy = "Maximum positive S1-S7 flux in the selected " + ...
            peakAudit.SourcePriority + " table within the event window; earliest tie; L1 uses bin center. An unavailable peak PAD stays selected.";
    end
end

%% peak and nearest available PAD records
if ~isempty(padTable)
    needed = {'EpochUTC', 'PADUsable'};
    for iSector = 1:7
        needed = [needed, {sprintf('RawFlux_S%d_%s', iSector, suffix), ...
            sprintf('FluxUncertainty_S%d_%s', iSector, suffix), ...
            sprintf('PA_S%d_deg', iSector)}]; %#ok<AGROW>
    end
    if ~all(ismember(needed, padTable.Properties.VariableNames))
        error('Case1:PeakPADVariables', 'The PAD table lacks required %s variables.', suffix);
    end
    time = padTable.EpochUTC;
    if isempty(time.TimeZone), time.TimeZone = 'UTC'; end
    time.TimeZone = 'UTC';
    inWindow = ~isnat(time) & time >= startTime & time < endTime;
    flux = nan(height(padTable), 7);
    sigma = flux; pa = flux;
    for iSector = 1:7
        flux(:, iSector) = padTable.(sprintf('RawFlux_S%d_%s', iSector, suffix));
        sigma(:, iSector) = padTable.(sprintf('FluxUncertainty_S%d_%s', iSector, suffix));
        pa(:, iSector) = padTable.(sprintf('PA_S%d_deg', iSector));
    end
    positiveFlux = flux;
    positiveFlux(~isfinite(positiveFlux) | positiveFlux <= 0) = NaN;
    rowMaximum = max(positiveFlux, [], 2, 'omitnan');
    completePAD = padTable.PADUsable == 1 & ...
        all(isfinite(flux) & flux > 0, 2) & all(isfinite(pa), 2);
    peakAudit.RowMaximumSectorFlux = rowMaximum;
    candidates = find(inWindow & isfinite(rowMaximum));
    peakAudit.CandidateTableRows = candidates;
    if isempty(candidates)
        peakAudit.Status = "no_positive_sector_flux";
    else
        peakFlux = max(rowMaximum(candidates));
        ties = candidates(rowMaximum(candidates) == peakFlux);
        [~, firstTime] = min(time(ties));
        peakRow = ties(firstTime);
        peakAudit.PeakTableRow = peakRow;
        peakAudit.PeakEpochUTC = time(peakRow);
        peakAudit.PeakFlux = peakFlux;
        peakAudit.PeakPADUsable = completePAD(peakRow);
        peakAudit.SelectedTableRows(3) = peakRow;

        before = find(inWindow & time < time(peakRow));
        after = find(inWindow & time > time(peakRow));
        [~, order] = sort(time(before), 'descend'); before = before(order);
        [~, order] = sort(time(after), 'ascend'); after = after(order);
        goodBefore = before(completePAD(before));
        goodAfter = after(completePAD(after));
        nBefore = min(2, numel(goodBefore));
        nAfter = min(2, numel(goodAfter));
        if nBefore > 0
            peakAudit.SelectedTableRows(3-(1:nBefore)) = goodBefore(1:nBefore);
        end
        if nBefore == 2
            stopBefore = find(before == goodBefore(nBefore), 1);
        else
            stopBefore = numel(before);
        end
        if nAfter > 0
            peakAudit.SelectedTableRows(3+(1:nAfter)) = goodAfter(1:nAfter);
        end
        if nAfter == 2
            stopAfter = find(after == goodAfter(nAfter), 1);
        else
            stopAfter = numel(after);
        end
        searchedBefore = before(1:stopBefore);
        searchedAfter = after(1:stopAfter);
        peakAudit.SkippedMissingBeforeTableRows = searchedBefore(~completePAD(searchedBefore));
        peakAudit.SkippedMissingAfterTableRows = searchedAfter(~completePAD(searchedAfter));
        peakAudit.SkippedMissingBeforeCount = numel(peakAudit.SkippedMissingBeforeTableRows);
        peakAudit.SkippedMissingAfterCount = numel(peakAudit.SkippedMissingAfterTableRows);

        for iPanel = 1:5
            row = peakAudit.SelectedTableRows(iPanel);
            if ~isfinite(row), continue; end
            peakAudit.SelectedEpochUTC(iPanel) = time(row);
            peakAudit.SelectedPADUsable(iPanel) = completePAD(row);
            selected = padTable(row, :);
            selected.Properties.UserData = [];
            peakAudit.SelectedRows{iPanel} = selected;
            peakAudit.RawFlux(iPanel, :) = flux(row, :);
            peakAudit.RawSigma(iPanel, :) = sigma(row, :);
            peakAudit.PA_deg(iPanel, :) = pa(row, :);
            peakAudit.NormalizationFlux(iPanel) = rowMaximum(row);
            peakAudit.NormalizingSectors{iPanel} = find(flux(row, :) == rowMaximum(row));
            peakAudit.NormalizedFlux(iPanel, :) = flux(row, :)/rowMaximum(row);
            peakAudit.NormalizedSigma(iPanel, :) = sigma(row, :)/rowMaximum(row);
            if all(ismember({'P1EnergyLower_MeV', 'P1EnergyUpper_MeV'}, ...
                    padTable.Properties.VariableNames))
                peakAudit.SourceEnergyMeV(iPanel, :) = ...
                    [padTable.P1EnergyLower_MeV(row), padTable.P1EnergyUpper_MeV(row)];
            end
        end
        peakAudit.SelectedCount = sum(isfinite(peakAudit.SelectedTableRows));
        peakAudit.PlottedCount = sum(peakAudit.SelectedPADUsable);
        if peakAudit.PlottedCount == 5
            peakAudit.Status = "complete";
        elseif ~peakAudit.PeakPADUsable
            peakAudit.Status = "peak_pad_unavailable";
        else
            peakAudit.Status = "partial_neighbors";
        end
    end
end

%% shared plotting range; no uncertainty-based exclusion
yLow = 0; yHigh = 1;
for iPanel = find(peakAudit.SelectedPADUsable).'
    y = peakAudit.NormalizedFlux(iPanel, :);
    dy = peakAudit.NormalizedSigma(iPanel, :);
    goodError = isfinite(dy) & dy >= 0;
    yLow = min([yLow, y(goodError)-dy(goodError)]);
    yHigh = max([yHigh, y(goodError)+dy(goodError)]);
end
yLow = floor(yLow*5)/5;
% Reserve a small header inside every panel, as in the supplied reference.
yHigh = max(1.2, ceil((yHigh+0.12*(yHigh-yLow))*5)/5);
peakAudit.DisplayYLimits = [yLow yHigh];
if ~logical(opts.RenderFigure), return; end

%% five portrait panels using the IRFU tight-panel and text helpers
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [80 80 1700 610], 'DefaultAxesPosition', [0.065 0.13 0.915 0.80]);
figureCleanup = onCleanup(@() closeHiddenFigure(fig, opts.Visible));
panelAxes = gobjects(5, 1);
for iPanel = 1:5
    panelAxes(iPanel) = irf_subplot(1, 5, -iPanel);
end
left = 0.060; right = 0.015; gap = 0.018;
panelWidth = (1-left-right-4*gap)/5;
for iPanel = 1:5
    ax = panelAxes(iPanel);
    ax.Position = [left+(iPanel-1)*(panelWidth+gap), 0.16, panelWidth, 0.70];
    hold(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 16, ...
        'LineWidth', 1.3, 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.025], ...
        'XLim', [0 180], 'XTick', 0:45:180, 'YLim', [yLow yHigh], ...
        'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
    if yLow == 0 && yHigh <= 1.4, ax.YTick = 0:0.2:yHigh; end
    if iPanel > 1, ax.YTickLabel = []; end
    xlabel(ax, 'Pitch angle (\circ)', 'FontSize', 18, 'Interpreter', 'tex');
    if iPanel == 1
        ylabel(ax, 'Normalized intensity', 'FontSize', 19);
    end
    titleText = sprintf('%+d', iPanel-3);
    if iPanel == 3, titleText = 'Peak'; end
    title(ax, titleText, 'FontSize', 16, 'FontWeight', 'normal');
    row = peakAudit.SelectedTableRows(iPanel);
    if ~isfinite(row)
        if isnat(peakAudit.PeakEpochUTC)
            message = {'No retained LECP'; 'record with positive flux'; 'in the event window'};
        else
            message = {'No additional valid PAD'; 'within the event window'};
        end
        panelText(ax, message, [0.5 0.50], ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'Interpreter', 'none');
        continue;
    end
    epoch = peakAudit.SelectedEpochUTC(iPanel);
    dateText = char(datetime(epoch, 'Format', 'yyyy-MM-dd'));
    timeText = [char(datetime(epoch, 'Format', 'HH:mm:ss')), ' UTC'];
    panelText(ax, {dateText; timeText}, [0.96 0.97], ...
        'FontSize', 15, 'Interpreter', 'none');
    intensityText = sprintf('J_{max} = %.3g', peakAudit.NormalizationFlux(iPanel));
    panelText(ax, intensityText, [0.96 0.835], 'FontSize', 14, 'Interpreter', 'tex');
    if ~peakAudit.SelectedPADUsable(iPanel)
        panelText(ax, {'Maximum-flux record'; 'complete S1-S7 PAD unavailable'}, ...
            [0.5 0.47], 'HorizontalAlignment', 'center', 'FontSize', 13, 'Interpreter', 'none');
        continue;
    end
    x = peakAudit.PA_deg(iPanel, :);
    y = peakAudit.NormalizedFlux(iPanel, :);
    dy = peakAudit.NormalizedSigma(iPanel, :);
    goodError = isfinite(dy) & dy >= 0;
    % All seven original points remain separate; even coincident angles stay.
    errorbar(ax, x(goodError), y(goodError), dy(goodError), ...
        'LineStyle', 'none', 'Marker', 'none', 'Color', 'k', ...
        'LineWidth', 1.2, 'CapSize', 7);
    plot(ax, x, y, 'ko', 'LineStyle', 'none', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 6.5, 'LineWidth', 1.0);
end
header = sprintf('Voyager %d  |  %s  |  LECP P1 %.2f-%.2f MeV', ...
    opts.Spacecraft, char(eventID), peakAudit.DisplayEnergyMeV);
annotation(fig, 'textbox', [0.06 0.945 0.92 0.047], 'String', header, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 17, 'Interpreter', 'none');
annotation(fig, 'textbox', [0.10 0.025 0.86 0.040], ...
    'String', 'J_{max} in cm^{-2} s^{-1} sr^{-1} MeV^{-1}; each panel normalized by its own J_{max}', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
folder = fileparts(figureFile);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
exportgraphics(fig, char(figureFile), 'Resolution', opts.ExportDPI, 'BackgroundColor', 'white');
peakAudit.FigureCreated = true;
end

function panelText(ax, labels, position, varargin)
% Apply one typography/color convention to every IRFU annotation line.
ht = irf_legend(ax, labels, position, varargin{:});
set(ht, 'Color', 'k', 'FontName', 'Times New Roman');
end

function closeHiddenFigure(fig, visible)
if ~logical(visible) && isgraphics(fig), close(fig); end
end
