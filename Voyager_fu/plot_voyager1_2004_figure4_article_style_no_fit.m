% Voyager 1 LECP daily pitch-angle distributions for 2004 DOY 221--224.
% Layout follows the four-panel style of Figure 4 in Florinski et al. (2008).
%
% Each panel contains the seven usable LECP P1 (0.57--1.78 MeV proton)
% sector measurements after separate 24 h averaging.  The horizontal
% coordinate is mu = cos(pitch angle), calculated by the revised
% Florinski-style full-vector daily method.  Error bars are the propagated
% statistical uncertainties stored by that program.  No fit is drawn.

% Work around the R2026a recent-artifacts/Home-service crash on this PC.
try
    feature('homeSessionManagerImpl', false);
catch
end

programFolder = fileparts(mfilename('fullpath'));
monthlyRoot = fullfile(programFolder, 'Voyager_Interstellar_Monthly');
dataRoot = fullfile(monthlyRoot, 'Voyager1_Selected_Event_Data');

sourceFile = fullfile(dataRoot, 'derived', ...
    'lecp_pitch_angle_florinski2008_daily', ...
    ['V1_2004-Fig4-DOY221-224_20040808_20040811_' ...
     'LECP_P1_Florinski2008_pitch_angle_1d.csv']);

if ~isfile(sourceFile)
    run(fullfile(programFolder, ...
        'run_voyager1_2004_figure4_new_florinski_daily_pad_no_fit.m'));
end
assert(isfile(sourceFile), 'VoyagerPAD:MissingDerivedData', ...
    'Daily pitch-angle data were not generated: %s', sourceFile);

outputFolder = fullfile(monthlyRoot, ...
    'Voyager1_2004_Florinski_Figure4_ArticleStyle_NoFit');
if ~isfolder(outputFolder), mkdir(outputFolder); end

T = readtable(sourceFile, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
epoch = datetime(T.EpochUTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
    'TimeZone', 'UTC');

targetDays = datetime(2004, 8, 8:11, 'TimeZone', 'UTC');
targetDOY = 221:224;
nSector = 7;

mu = nan(4, nSector);
flux = nan(4, nSector);
sigma = nan(4, nSector);
pa = nan(4, nSector);
sampleCount = nan(4, nSector);
rowIndex = nan(4, 1);

for iDay = 1:4
    idx = find(dateshift(epoch, 'start', 'day') == targetDays(iDay), 1);
    assert(~isempty(idx), 'VoyagerPAD:MissingDay', ...
        'No daily row found for %s.', char(targetDays(iDay), 'yyyy-MM-dd'));
    assert(T.PADUsable(idx) == 1, 'VoyagerPAD:ExcludedDay', ...
        'Daily PAD was marked unusable for %s: %s', ...
        char(targetDays(iDay), 'yyyy-MM-dd'), T.PADQuality(idx));
    rowIndex(iDay) = idx;

    for iSector = 1:nSector
        mu(iDay, iSector) = T.(sprintf('Mu_S%d', iSector))(idx);
        pa(iDay, iSector) = T.(sprintf('PA_S%d_deg', iSector))(idx);
        flux(iDay, iSector) = ...
            T.(sprintf('Flux_S%d_1d', iSector))(idx);
        sigma(iDay, iSector) = ...
            T.(sprintf('FluxUncertainty_S%d_1d', iSector))(idx);
        sampleCount(iDay, iSector) = ...
            T.(sprintf('Samples_S%d_1d', iSector))(idx);
    end
end

% Keep the same vertical scale in all four panels, as in a comparison plot.
upperValues = flux + sigma;
yTop = ceil(1.12 * max(upperValues(:), [], 'omitnan') * 10) / 10;
if ~isfinite(yTop) || yTop <= 0, yTop = 1; end

fig = figure('Color', 'w', 'Visible', 'on', ...
    'Name', 'Voyager 1 2004 DOY 221-224 daily PAD', ...
    'Position', [120 80 960 780]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for iDay = 1:4
    ax = nexttile(tl, iDay);
    hold(ax, 'on');

    good = isfinite(mu(iDay, :)) & isfinite(flux(iDay, :)) & ...
        isfinite(sigma(iDay, :)) & sigma(iDay, :) >= 0;
    x = mu(iDay, good);
    y = flux(iDay, good);
    e = sigma(iDay, good);
    [x, order] = sort(x);
    y = y(order);
    e = e(order);

    % Observations and statistical uncertainty only; deliberately no fit.
    errorbar(ax, x, y, e, 'o', ...
        'Color', 'k', 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 5.2, ...
        'LineStyle', 'none', 'LineWidth', 1.05, ...
        'CapSize', 7);

    xlim(ax, [-1.05 1.05]);
    ylim(ax, [0 yTop]);
    xticks(ax, [-1 -0.5 0 0.5 1]);
    grid(ax, 'off');
    box(ax, 'on');
    ax.TickDir = 'in';
    ax.TickLength = [0.018 0.018];
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    xlabel(ax, '\mu', 'Interpreter', 'tex', 'FontSize', 14);
    ylabel(ax, {'j', '(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
        'Interpreter', 'tex', 'FontSize', 12);

    meanJ = mean(y, 'omitnan');
    text(ax, 0.04, 0.93, sprintf('(%d)  2004/%03d', ...
        iDay, targetDOY(iDay)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontName', 'Times New Roman', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', 'k');
    text(ax, 0.96, 0.93, sprintf('\\itj\\rm = %.3f', meanJ), ...
        'Units', 'normalized', 'Interpreter', 'tex', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontName', 'Times New Roman', 'FontSize', 11, 'Color', 'k');
end

% The paper-style panels carry their own dates and mean intensities; no
% super-title and no fit legend are added.

baseName = 'Voyager1_2004_DOY221_224_Florinski_Figure4_style_no_fit';
pngFile = fullfile(outputFolder, baseName + ".png");
pdfFile = fullfile(outputFolder, baseName + ".pdf");
figFile = fullfile(outputFolder, baseName + ".fig");
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
savefig(fig, figFile);

% Save the exact points displayed in the four panels in long-table form.
nRows = 4 * nSector;
DayUTC = NaT(nRows, 1, 'TimeZone', 'UTC');
DOY = nan(nRows, 1);
Sector = nan(nRows, 1);
Mu = nan(nRows, 1);
PitchAngle_deg = nan(nRows, 1);
Flux = nan(nRows, 1);
FluxUncertainty = nan(nRows, 1);
SectorSamples = nan(nRows, 1);

k = 0;
for iDay = 1:4
    for iSector = 1:nSector
        k = k + 1;
        DayUTC(k) = targetDays(iDay);
        DOY(k) = targetDOY(iDay);
        Sector(k) = iSector;
        Mu(k) = mu(iDay, iSector);
        PitchAngle_deg(k) = pa(iDay, iSector);
        Flux(k) = flux(iDay, iSector);
        FluxUncertainty(k) = sigma(iDay, iSector);
        SectorSamples(k) = sampleCount(iDay, iSector);
    end
end

plotData = table(DayUTC, DOY, Sector, Mu, PitchAngle_deg, Flux, ...
    FluxUncertainty, SectorSamples);
dataFile = fullfile(outputFolder, baseName + "_points.csv");
writetable(plotData, dataFile);

fprintf('Article-style Figure 4 plot written to:\n  %s\n', pngFile);
fprintf('Vector PDF written to:\n  %s\n', pdfFile);
fprintf('Displayed point table written to:\n  %s\n', dataFile);

