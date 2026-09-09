function result = Run_Case1_2020E01_BaselineComparison
%Run_Case1_2020E01_BaselineComparison Test alternate temporal baselines.
%   Current raw-CDF renderer, independent S1--S7, three-day centered means.
%   Fixed five-day intervals are sensitivity tests, not validated backgrounds.
%   Preserve the existing published figure; all new results have their own root.
%   Modified: 2026-09-09

%% paths and baseline definitions
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
root = fullfile(fileparts(cfg.OutputRoot),'2020E01_BaselineComparison');
if ~isfolder(root), mkdir(root); end
starts = datetime([2020 7 2;2020 7 22;2020 8 12;2020 8 24], ...
    'TimeZone','UTC');
ends = starts+days(5);
names = ["original_Jul02_06";"Jul22_26";"Aug12_16";"Aug24_28"];
allTables = cell(4,1);
allOptions = cell(4,1);
result = struct('Root',root,'GeneratedUTC',datetime('now','TimeZone','UTC'));
result.Runs = cell(4,1);
commonLimits = [];
referenceUpper = [];
baseline = nan(4,7);
negativePercent = nan(4,1);
s4HighestPercent = nan(4,1);
eventS4Median = nan(4,1);
eventS4Maximum = nan(4,1);

%% direct original CDF re-runs through the current formal program
for k = 1:4
    fprintf('\nBASELINE_COMPARISON %d/4: %s\n',k,names(k));
    [run, payload] = Run_Case1_2020E01_BottomVariants( ...
        'Modes','pad_difference','PADDisplayAverageDays',3, ...
        'BaselineUTC',[starts(k),ends(k)], ...
        'OutputRoot',fullfile(root,char(names(k))), ...
        'DifferenceColorLimits',commonLimits);
    T = payload{1}.Table;
    opts = payload{1}.Options;
    audit = T.Properties.UserData;
    if k == 1
        commonLimits = audit.ColorLimits;
        referenceUpper = audit.UpperPanelLines;
    else
        assert(isequaln(referenceUpper,audit.UpperPanelLines));
        assert(isequaln(commonLimits,audit.ColorLimits));
        assert(isequaln(T.PADUsable,allTables{1}.PADUsable));
        fields = [compose('Flux_S%d_1d',1:7),compose('PA_S%d_deg',1:7)];
        assert(isequaln(T{:,cellstr(fields)},allTables{1}{:,cellstr(fields)}));
    end
    opts.DifferenceColorLimits = commonLimits;
    allTables{k} = T;
    allOptions{k} = opts;
    result.Runs{k} = run;
    baseline(k,:) = audit.BaselineMeanBySector;
    D = T{:,cellstr(compose('DisplayValue_S%d',1:7))};
    usable = T.BottomPanelUsable;
    V = D(usable,:);
    negativePercent(k) = 100*nnz(V<0)/nnz(isfinite(V));
    % Diagnostic ranking only; original sector values remain independent.
    s4HighestPercent(k) = 100*mean(V(:,4)>=max(V(:,[1:3 5:7]),[],2));
    event = usable & T.EpochUTC>=datetime(2020,7,22,'TimeZone','UTC') & ...
        T.EpochUTC<datetime(2020,8,8,'TimeZone','UTC');
    eventS4Median(k) = median(D(event,4),'omitnan');
    eventS4Maximum(k) = max(D(event,4));
end

%% numerical sensitivity summary, retaining signed differences
summary = table(names,starts,ends,baseline,negativePercent,s4HighestPercent, ...
    eventS4Median,eventS4Maximum,'VariableNames', ...
    {'Name','StartUTC','EndUTCExclusive','BaselineS1toS7', ...
    'NegativeDisplayedCellPercent','S4HighestDayPercent', ...
    'EventS4Median3d','EventS4Max3d'});
T = allTables{1};
J = T{:,cellstr(compose('Flux_S%d_1d',1:7))};
otherMedian = median(J(:,[1:3 5:7]),2,'omitnan');
ratio = J(:,4)./otherMedian;
valid = isfinite(ratio) & isfinite(J(:,4)) & isfinite(otherMedian);
correlation = corrcoef(J(valid,4),otherMedian(valid));
result.S4ToOtherSectorMedianQuartiles = prctile(ratio(valid),[25 50 75]);
result.S4OtherMedianCorrelation = correlation(1,2);
result.DiagnosticPolicy = ['Median of the six other sectors is a diagnostic ', ...
    'reference only; no sector combination enters the PAD or baseline subtraction. ', ...
    'All finite days retained; no fitted correction, new threshold or outlier removal.'];
result.Summary = summary;
result.CommonColorLimits = commonLimits;
result.UpperPanelsPitchAnglesAndMasksIdentical = true;
result.BackgroundValidation = 'None; higher-flux intervals are sensitivity tests only.';
writetable(summary,fullfile(root,'baseline_summary.csv'));
daily = table(T.EpochUTC,J,otherMedian,ratio,T.PADUsable, ...
    'VariableNames',{'EpochUTC','S1toS7','OtherSixMedianDiagnostic', ...
    'S4RatioDiagnostic','PADUsable'});
writetable(daily,fullfile(root,'daily_sector_diagnostics.csv'));

%% matched PAD comparison, using the same current sector-cell renderer
fig = figure('Visible','off','Color','w','Position',[80 80 1550 1180]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
for k = 1:4
    ax = nexttile(layout); hold(ax,'on');
    plotted = Case1_Plot_LECP_Bottom_Variant(ax,allTables{k},allOptions{k});
    fields = cellstr(compose('DisplayValue_S%d',1:7));
    assert(isequaln(plotted{:,fields},allTables{k}{:,fields}));
    xlim(ax,datenum([datetime(2020,6,22),datetime(2020,9,8)]));
    datetick(ax,'x','dd-mmm','keeplimits');
    ylabel(ax,'PA (deg)'); box(ax,'on');
    title(ax,sprintf('%s to %s UTC | S4 baseline = %.3f', ...
        datestr(starts(k),'dd-mmm-yyyy'), ...
        datestr(ends(k)-days(1),'dd-mmm-yyyy'),baseline(k,4)));
end
xlabel(ax,'UTC');
title(layout,'Voyager 1 2020-E01 | Three-day PAD differences | Common color scale');
result.ComparisonFigure = fullfile(root,'PAD_baseline_comparison.png');
exportgraphics(fig,result.ComparisonFigure,'Resolution',180);
clear cleanup

%% signed S4 residuals expose any suppression by the zero-based color scale
fig = figure('Visible','off','Color','w','Position',[80 80 1450 880]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout); hold(ax,'on');
plot(ax,T.EpochUTC,J(:,4),'k-','LineWidth',1.2);
plot(ax,T.EpochUTC,otherMedian,'Color',[0.2 0.5 0.8],'LineWidth',1.2);
ylabel(ax,'Daily J'); legend(ax,{'S4','Median of other six (diagnostic)'},'Location','northwest');
title(ax,'Original daily flux: no baseline removed'); grid(ax,'on');
ax = nexttile(layout); plot(ax,T.EpochUTC,ratio,'k-');
ylabel(ax,'S4 / other-six median'); grid(ax,'on');
ax = nexttile(layout); hold(ax,'on');
for k = 1:4
    plot(ax,T.EpochUTC,allTables{k}.DisplayValue_S4,'LineWidth',1.2, ...
        'DisplayName',sprintf('%s--%s',datestr(starts(k),'dd-mmm'), ...
        datestr(ends(k)-days(1),'dd-mmm')));
end
yline(ax,0,'k--','HandleVisibility','off');
ylabel(ax,'Three-day S4 difference'); xlabel(ax,'UTC');
legend(ax,'Location','northwest'); grid(ax,'on');
title(layout,'S4 temporal-baseline test | Signed values retained');
result.DiagnosticFigure = fullfile(root,'S4_signed_diagnostics.png');
exportgraphics(fig,result.DiagnosticFigure,'Resolution',180);
clear cleanup
result.EntryPointSHA256 = Case1_File_SHA256([mfilename('fullpath'),'.m']);
save(fullfile(root,'comparison_result.mat'),'result','-v7.3');
disp(summary);
disp(result.S4ToOtherSectorMedianQuartiles);
fprintf('S4 / other-six median correlation: %.6f\n',result.S4OtherMedianCorrelation);
fprintf('BASELINE_COMPARISON_VERIFIED\n');
end
