function audit = Case1_Plot_Florinski_Figure4_PAD(padTable, figureFile, opts)
%Case1_Plot_Florinski_Figure4_PAD Four source-identified daily PAD records.
%   Use the existing CDF-to-PAD table for 2004/221--224 (August 8--11).
%   Original Epoch selects its UTC day; each day must contain at most one
%   retained daily-product record. Missing or incomplete PADs remain blank.
%   Each complete S1--S7 PAD uses its own Jmax and sigma/Jmax with a fixed
%   denominator. No flux averaging, angle merging, fitting or S8 is used.
%   opts.RenderFigure=false performs selection and auditing only.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% options and four-day audit
if nargin < 3, opts = struct; end
if ~isfield(opts, 'Visible'), opts.Visible = false; end
if ~isfield(opts, 'ExportDPI'), opts.ExportDPI = 220; end
if ~isfield(opts, 'RenderFigure'), opts.RenderFigure = true; end
if ~isfield(opts, 'P1DisplayEnergyMeV'), opts.P1DisplayEnergyMeV = [0.57 1.78]; end
if ~istable(padTable)
    error('Case1:FlorinskiPADTable', 'The PAD input must be a table.');
end
dayStart = datetime(2004, 8, 8, 'TimeZone', 'UTC') + days((0:3).');
audit = struct;
audit.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
audit.Reference = "Florinski et al. 2008, Figure 4";
audit.Spacecraft = 1;
audit.PADCadence = "day";
audit.WindowStartUTC = dayStart(1);
audit.WindowEndUTCExclusive = dayStart(end)+days(1);
audit.DayStartUTC = dayStart;
audit.DayOfYear = (221:224).';
audit.SelectedTableRows = nan(4, 1);
audit.SelectedEpochUTC = NaT(4, 1, 'TimeZone', 'UTC');
audit.SelectedRows = cell(4, 1);
audit.SelectedPADUsable = false(4, 1);
audit.RecordsPerDay = zeros(4, 1);
audit.StatusPerDay = repmat("no_retained_record", 4, 1);
audit.RawFlux = nan(4, 7);
audit.RawSigma = nan(4, 7);
audit.PA_deg = nan(4, 7);
audit.NormalizedFlux = nan(4, 7);
audit.NormalizedSigma = nan(4, 7);
audit.NormalizationFlux = nan(4, 1);
audit.NormalizingSectors = cell(4, 1);
audit.SourceCDF = strings(4, 1);
audit.SourceCDFRecord = nan(4, 1);
audit.SourceEnergyMeV = nan(4, 2);
audit.DisplayEnergyMeV = reshape(opts.P1DisplayEnergyMeV, 1, 2);
audit.SourceTableUserData = padTable.Properties.UserData;
audit.InputTableRowCount = height(padTable);
audit.SelectedCount = 0;
audit.PlottedCount = 0;
audit.Status = "no_retained_records";
audit.Sectors = 1:7;
audit.FluxUnits = "cm^-2 s^-1 sr^-1 MeV^-1";
audit.DaySelectionPolicy = "The original Epoch selects its UTC calendar day in 2004-08-08 through 2004-08-11 inclusive; more than one retained record in a day raises an error; no relabelling or reaveraging.";
audit.ValidityPolicy = "Use the existing PADUsable decision and require all seven original S1-S7 fluxes finite and positive and all seven PA values finite; missing or incomplete days remain blank.";
audit.NormalizationPolicy = "Each daily record divided by its own maximum positive finite S1-S7 flux; sigma divided by the same fixed scale; denominator uncertainty not additionally propagated.";
audit.FluxPolicy = "Use the current original-CDF PAD table after the approved DeltaT<0 exclusion; this renderer adds no averaging, interpolation, background subtraction, angle merging, fitting or S8 use.";
audit.EnergyLabelPolicy = "User-approved 0.57-1.78 MeV display convention; source CDF energy metadata and all flux values remain separate and unchanged.";
audit.UncertaintyDisplayPolicy = "Finite nonnegative sigma values receive error bars; unavailable sigma does not remove the point; common y limits contain all displayed bars.";
if ismember('SourceProduct',padTable.Properties.VariableNames) && ...
        any(padTable.SourceProduct == "L1_UTC_mean")
    audit.DaySelectionPolicy = "One L2-first or approved L1 UTC-day mean per fixed Figure4 day; L1 time is the UTC bin center. No additional averaging in this renderer.";
    audit.FluxPolicy = "Input L1 rate means use J=mean(R)/[0.44*(1.78-0.57)] without background subtraction; complete L2 values stay unchanged. Source records/counts/uncertainties are in SourceTableUserData.L1FallbackAudit.";
end
audit.FigureFile = string(figureFile);
audit.FigureCreated = false;
audit.DisplayXLimits_deg = [0 180];
audit.DisplayYLimits = [0 1.2];

%% retain one original daily record per reference day
if ~isempty(padTable)
    needed = {'EpochUTC', 'PADUsable'};
    for iSector = 1:7
        needed = [needed, {sprintf('RawFlux_S%d_1d', iSector), ...
            sprintf('FluxUncertainty_S%d_1d', iSector), ...
            sprintf('PA_S%d_deg', iSector)}]; %#ok<AGROW>
    end
    if ~all(ismember(needed, padTable.Properties.VariableNames))
        error('Case1:FlorinskiPADVariables', 'The PAD table lacks required daily-product variables.');
    end
    time = padTable.EpochUTC;
    if ~isdatetime(time)
        error('Case1:FlorinskiPADEpoch', 'EpochUTC must be a datetime vector.');
    end
    if isempty(time.TimeZone), time.TimeZone = 'UTC'; end
    time.TimeZone = 'UTC';
    for iPanel = 1:4
        rows = find(~isnat(time) & time >= dayStart(iPanel) & ...
            time < dayStart(iPanel)+days(1));
        audit.RecordsPerDay(iPanel) = numel(rows);
        if numel(rows) > 1
            error('Case1:FlorinskiMultipleDailyRecords', ...
                'UTC day %s contains %d retained records. No selection or reaveraging was applied.', ...
                char(datetime(dayStart(iPanel), 'Format', 'yyyy-MM-dd')), numel(rows));
        end
        if isempty(rows), continue; end
        row = rows(1);
        audit.SelectedTableRows(iPanel) = row;
        audit.SelectedEpochUTC(iPanel) = time(row);
        selected = padTable(row, :);
        selected.Properties.UserData = [];
        audit.SelectedRows{iPanel} = selected;
        for iSector = 1:7
            audit.RawFlux(iPanel, iSector) = padTable.(sprintf('RawFlux_S%d_1d', iSector))(row);
            audit.RawSigma(iPanel, iSector) = padTable.(sprintf('FluxUncertainty_S%d_1d', iSector))(row);
            audit.PA_deg(iPanel, iSector) = padTable.(sprintf('PA_S%d_deg', iSector))(row);
        end
        flux = audit.RawFlux(iPanel, :);
        positiveFlux = flux(isfinite(flux) & flux > 0);
        if ~isempty(positiveFlux)
            maximum = max(positiveFlux);
            audit.NormalizationFlux(iPanel) = maximum;
            audit.NormalizingSectors{iPanel} = find(flux == maximum);
            audit.NormalizedFlux(iPanel, :) = flux/maximum;
            audit.NormalizedSigma(iPanel, :) = audit.RawSigma(iPanel, :)/maximum;
        end
        audit.SelectedPADUsable(iPanel) = padTable.PADUsable(row) == 1 && ...
            all(isfinite(flux) & flux > 0) && all(isfinite(audit.PA_deg(iPanel, :)));
        audit.StatusPerDay(iPanel) = "pad_unavailable";
        if audit.SelectedPADUsable(iPanel), audit.StatusPerDay(iPanel) = "complete"; end
        if ismember('SourceCDF', padTable.Properties.VariableNames)
            audit.SourceCDF(iPanel) = string(padTable.SourceCDF(row));
        end
        if ismember('SourceCDFRecord', padTable.Properties.VariableNames)
            audit.SourceCDFRecord(iPanel) = padTable.SourceCDFRecord(row);
        end
        if all(ismember({'P1EnergyLower_MeV', 'P1EnergyUpper_MeV'}, padTable.Properties.VariableNames))
            audit.SourceEnergyMeV(iPanel, :) = ...
                [padTable.P1EnergyLower_MeV(row), padTable.P1EnergyUpper_MeV(row)];
        end
    end
end
audit.SelectedCount = sum(isfinite(audit.SelectedTableRows));
audit.PlottedCount = sum(audit.SelectedPADUsable);
if audit.PlottedCount == 4
    audit.Status = "complete";
elseif audit.SelectedCount > 0
    audit.Status = "partial_or_unavailable";
end

%% common limits retain large uncertainties without adding a quality gate
yLow = 0; yHigh = 1;
for iPanel = find(audit.SelectedPADUsable).'
    y = audit.NormalizedFlux(iPanel, :);
    dy = audit.NormalizedSigma(iPanel, :);
    goodError = isfinite(dy) & dy >= 0;
    yLow = min([yLow, y(goodError)-dy(goodError)]);
    yHigh = max([yHigh, y(goodError)+dy(goodError)]);
end
yLow = floor(yLow*5)/5;
% Leave the upper panel region free for original Epoch and normalization.
yHigh = max(1.2, ceil((yHigh+0.25*(yHigh-yLow))*5)/5);
audit.DisplayYLimits = [yLow yHigh];
if ~logical(opts.RenderFigure), return; end

%% four portrait panels with the existing IRFU layout and annotation helpers
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
visibility = 'off';
if logical(opts.Visible), visibility = 'on'; end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Position', [80 80 1450 660], 'DefaultAxesPosition', [0.065 0.13 0.915 0.80]);
figureCleanup = onCleanup(@() closeHiddenFigure(fig, opts.Visible));
panelAxes = gobjects(4, 1);
for iPanel = 1:4, panelAxes(iPanel) = irf_subplot(1, 4, -iPanel); end
left = 0.065; right = 0.015; gap = 0.022;
panelWidth = (1-left-right-3*gap)/4;
for iPanel = 1:4
    ax = panelAxes(iPanel);
    ax.Position = [left+(iPanel-1)*(panelWidth+gap), 0.155, panelWidth, 0.70];
    hold(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 17, ...
        'LineWidth', 1.3, 'Box', 'on', 'TickDir', 'in', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.025], ...
        'XLim', [0 180], 'XTick', 0:45:180, 'YLim', [yLow yHigh], ...
        'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
    if yLow == 0 && yHigh <= 1.4, ax.YTick = 0:0.2:yHigh; end
    if iPanel > 1, ax.YTickLabel = []; end
    xlabel(ax, 'Pitch angle (\circ)', 'FontSize', 19, 'Interpreter', 'tex');
    if iPanel == 1, ylabel(ax, 'Normalized intensity', 'FontSize', 20); end
    dateText = char(datetime(dayStart(iPanel), 'Format', 'yyyy-MM-dd'));
    title(ax, sprintf('2004/%d  |  %s', audit.DayOfYear(iPanel), dateText), ...
        'FontSize', 17, 'FontWeight', 'normal', 'Interpreter', 'none');
    if ~isfinite(audit.SelectedTableRows(iPanel))
        panelText(ax, {'No retained daily'; 'LECP record'}, [0.5 0.5], ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'none');
        continue;
    end
    timeText = [char(datetime(audit.SelectedEpochUTC(iPanel), 'Format', 'HH:mm:ss.SSS')), ' UTC'];
    panelText(ax, timeText, [0.96 0.97], 'FontSize', 16, 'Interpreter', 'none');
    panelText(ax, sprintf('J_{max} = %.3g', audit.NormalizationFlux(iPanel)), ...
        [0.96 0.90], 'FontSize', 16, 'Interpreter', 'tex');
    if ~audit.SelectedPADUsable(iPanel)
        panelText(ax, {'Complete S1-S7'; 'PAD unavailable'}, [0.5 0.47], ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'Interpreter', 'none');
        continue;
    end
    x = audit.PA_deg(iPanel, :);
    y = audit.NormalizedFlux(iPanel, :);
    dy = audit.NormalizedSigma(iPanel, :);
    goodError = isfinite(dy) & dy >= 0;
    errorbar(ax, x(goodError), y(goodError), dy(goodError), ...
        'LineStyle', 'none', 'Marker', 'none', 'Color', 'k', ...
        'LineWidth', 1.2, 'CapSize', 7);
    plot(ax, x, y, 'ko', 'LineStyle', 'none', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 7, 'LineWidth', 1.0);
end
header = sprintf('Voyager 1  |  Florinski 2008 Figure 4 interval  |  LECP P1 %.2f-%.2f MeV', ...
    audit.DisplayEnergyMeV);
annotation(fig, 'textbox', [0.06 0.945 0.92 0.045], 'String', header, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 18, 'Interpreter', 'none');
annotation(fig, 'textbox', [0.06 0.028 0.92 0.045], ...
    'String', 'J_{max} in cm^{-2} s^{-1} sr^{-1} MeV^{-1}; each daily record normalized by its own J_{max}', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 15, 'Interpreter', 'tex');
folder = fileparts(figureFile);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
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
