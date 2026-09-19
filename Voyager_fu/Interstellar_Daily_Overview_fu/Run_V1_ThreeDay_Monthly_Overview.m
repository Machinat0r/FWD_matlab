function result = Run_V1_ThreeDay_Monthly_Overview(visible,requestedMode,outputRoot,makePlot)
% 六扇区版本：四个面板的非重叠三日平均、UTC自然月平均和Pearson相关。
% 原始CDF -> 既有日统计 -> 各面板独立按窗口算术平均。
if nargin<1
    visible = true;
end

if nargin<2
    requestedMode = 'both';
end
requestedMode = validatestring(requestedMode,{'both','three_day','monthly'});

%% 路径与分析时段
CodeDir = 'C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/';
ParentDir = 'Z:/SPART-WORK/Data/Voyager/';
OutputDir = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/averaged_no_S4/';
if nargin>=3 && ~isempty(outputRoot), OutputDir=char(outputRoot); end
if nargin<4, makePlot=true; end
addpath([CodeDir,'Interstellar_Daily_Overview_fu']);
if ~isfolder(OutputDir)
    mkdir(OutputDir);
end
prior = []; % 只用于保护未修改的输出和比较旧图，不提供科学输入。
priorFile = fullfile(OutputDir,'three_day_monthly_audit.mat');
if ~strcmp(requestedMode,'both') && isfile(priorFile)
    previous = load(priorFile,'result');
    prior = previous.result;
end
StartUTC = datetime(2012,8,25,'TimeZone','UTC');
EndUTC = datetime(2021,12,17,'TimeZone','UTC');

%% 读取原始CDF，沿用日统计及L1优先通路
[source,hourly] = Run_V1_Interstellar_Daily_Overview( ...
    'DataRoot',ParentDir,'OutputRoot',fullfile(OutputDir,'source_recompute'), ...
    'SyncArchive',false,'Visible',false,'MakePlot',false);
use = source.Daily.EpochUTC>=StartUTC & source.Daily.EpochUTC<EndUTC;
daily = source.Daily(use,:);
assert(all(diff(daily.EpochUTC)==days(1)),'Incomplete daily grid.');
values = [daily.BMean_nT daily.P1Mean daily.P1Median daily.P1SixSectorSum];
hourly = hourly(hourly.EpochUTC>=StartUTC & hourly.EpochUTC<EndUTC,:);

%% 两种分组：从起点开始每三日；UTC自然月
result = struct;
result.CreatedUTC = datetime('now','TimeZone','UTC');
result.Daily = daily;
result.HourlyP1 = hourly(:,{'EpochUTC','P1','FileIndex','CDFRecord'});
result.SourceMethod = source.Method;
result.SourceAudit = string(fullfile(source.OutputFolder,'V1_daily_overview_audit.mat'));
result.SourceCodeFiles = source.CodeFiles;
result.SourceCodeSHA256 = source.CodeSHA256;
result.COHOSources = source.COHOSources;
result.SectorAuditFiles = source.SectorAuditFiles;
result.Method = struct('PanelDSectors',[1 2 3 5 6 7], ...
    'StartUTC',StartUTC,'EndUTCExclusive',EndUTC, ...
    'Average','Unweighted arithmetic mean of finite daily panel values; each panel averaged independently.', ...
    'PanelC','Three-day: median of all finite original COHO hourly P1 values in the bin. Monthly: arithmetic mean of daily medians, unchanged.', ...
    'ThreeDay','Nonoverlapping [start+3*k days,start+3*(k+1) days), anchored at 2012-08-25 UTC.', ...
    'Monthly','UTC calendar months clipped to the requested interval at both ends.', ...
    'PartialWindows','Retained; actual day counts and finite contributing days per panel recorded.', ...
    'Missing','No filling or interpolation; empty windows NaN. Source-valid zero retained. No minimum coverage threshold.', ...
    'Correlation','Pearson on untransformed window statistics; same-window pairwise finite selection. No lag, detrending, weighting or outlier rejection.', ...
    'Plot','a linear; b/c/d logarithmic. Nonpositive values hidden only in plot. Thin lines break at NaN.', ...
    'IRFUReview','irf_waverage uses weighted fixed-point windows and replaces NaN with zero; explicit UTC grouping used here.');
mode = {'three_day','monthly'};
if ~strcmp(requestedMode,'both')
    mode = {requestedMode};
end
allCorrelations = table;
for m = 1:numel(mode)
    if strcmp(mode{m},'three_day')
        nominalStart = (StartUTC:caldays(3):EndUTC-seconds(1)).';
        nominalEnd = nominalStart+caldays(3);
        plotTitle = '3-day means (a,b,d) and 3-day P1 median (c)';
    else
        nominalStart = (dateshift(StartUTC,'start','month'):calmonths(1): ...
            dateshift(EndUTC-seconds(1),'start','month')).';
        nominalEnd = nominalStart+calmonths(1);
        plotTitle = 'Calendar-month means of daily values';
    end
    begin = nominalStart;
    finish = nominalEnd;
    begin(begin<StartUTC) = StartUTC;
    finish(finish>EndUTC) = EndUTC;
    [average,count,binIndex] = averageWindows(daily.EpochUTC,values,begin,finish);
    oldPanelCMean = average(:,3);
    medianCount = [];
    hourlyBinIndex = [];
    if strcmp(mode{m},'three_day')
        [average(:,3),medianCount,hourlyBinIndex] = medianWindows( ...
            hourly.EpochUTC,hourly.P1,begin,finish);
    end
    t = begin+(finish-begin)/2;
    actualDays = days(finish-begin);
    partial = begin~=nominalStart | finish~=nominalEnd;
    windows = table(nominalStart,nominalEnd,begin,finish,t,actualDays,partial, ...
        average(:,1),average(:,2),average(:,3),average(:,4), ...
        count(:,1),count(:,2),count(:,3),count(:,4), ...
        'VariableNames',{'NominalStartUTC','NominalEndUTC','StartUTC','EndUTCExclusive', ...
        'EpochUTC','CalendarDaysInWindow','PartialWindow','BMean_nT','MeanOfDailyP1Means', ...
        'MeanOfDailyP1Medians','MeanOfDailySixSectorSums', ...
        'ValidDays_A','ValidDays_B','ValidDays_C','ValidDays_D'});
    if strcmp(mode{m},'three_day')
        windows.Properties.VariableNames{'MeanOfDailyP1Medians'} = 'ThreeDayP1Median';
        windows.P1HourlySampleCount_C = medianCount;
        windows.PreviousMeanOfDailyP1Medians = oldPanelCMean;
    end
    [correlation,pairUsed] = correlateWindows(average,mode{m});
    allCorrelations = [allCorrelations;correlation]; %#ok<AGROW>
    files = strings(0,1);
    if makePlot
        files = drawPanels(t,average,StartUTC,EndUTC,plotTitle,mode{m},OutputDir,visible);
    end
    result.(mode{m}) = struct('Windows',windows,'DailyBinIndex',binIndex, ...
        'Correlation',correlation,'PairUsed',pairUsed,'OutputFiles',files, ...
        'HourlyP1BinIndex',hourlyBinIndex);
    writetable(windows,fullfile(OutputDir,[mode{m},'_averages.csv']),'Encoding','UTF-8');
    pairTable = table(t,pairUsed(:,1),pairUsed(:,2),pairUsed(:,3), ...
        'VariableNames',{'EpochUTC','Used_AB','Used_AC','Used_AD'});
    writetable(pairTable,fullfile(OutputDir,[mode{m},'_correlation_pairs.csv']));
end

if ~isempty(prior)
    preservedModes = setdiff({'three_day','monthly'},mode,'stable');
    for k = 1:numel(preservedModes)
        name = preservedModes{k};
        if isfield(prior,name)
            result.(name) = prior.(name); % 未修改结果的原样保留。
            allCorrelations = [allCorrelations;prior.(name).Correlation]; %#ok<AGROW>
        end
    end
    if strcmp(requestedMode,'three_day') && isfield(prior,'three_day')
        checkNames = {'BMean_nT','MeanOfDailyP1Means','MeanOfDailySixSectorSums', ...
            'ValidDays_A','ValidDays_B','ValidDays_D'};
        assert(isequaln(result.three_day.Windows(:,checkNames), ...
            prior.three_day.Windows(:,checkNames)),'Unrequested panels changed.');
        result.PanelsABDUnchanged = true;
    end
end
%% 保存独立结果与完整来源审计
result.Correlation = allCorrelations;
result.ProgramSHA256 = string(Case1_File_SHA256([mfilename('fullpath'),'.m']));
save(fullfile(OutputDir,'three_day_monthly_audit.mat'),'result','-v7.3');
writetable(allCorrelations,fullfile(OutputDir,'Pearson_correlations.csv'),'Encoding','UTF-8');
fid = fopen(fullfile(OutputDir,'三日与月平均_结果说明.md'),'w','n','UTF-8');
assert(fid>=0);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'# 六扇区版本：三日平均和月平均\n\n');
fprintf(fid,'使用原始CDF重算日统计；时段2012-08-25至2021-12-16（含）。d始终为S1、S2、S3、S5、S6、S7，排除S4和S8。\n\n');
fprintf(fid,'三日窗口从2012-08-25起连续且不重叠，最后不足三日的窗口保留。月平均按UTC自然月，首尾月仅使用时段内数据。a、b、d对有限日值等权算术平均。三日图c对窗口内所有原始COHO小时P1通量求中位数；月图c继续使用日中位数的算术平均。零值保留，不插值、不填缺测，不设置最低覆盖率。\n\n');
fprintf(fid,'Pearson相关基于未取对数的窗口统计量；每组按同一窗口双方有限配对，各面板在窗口内可有不同的有效日期。只报告描述性系数。\n\n');
fprintf(fid,'| 平均方式 | 配对 | Pearson r | 有效窗口数 |\n|---|---|---:|---:|\n');
for k = 1:height(allCorrelations)
    fprintf(fid,'| %s | %s | %.6f | %d |\n',char(allCorrelations.Mode(k)), ...
        char(allCorrelations.Pair(k)),allCorrelations.PearsonR(k),allCorrelations.ValidWindows(k));
end
fprintf(fid,'\n窗口起止、首尾不完整标志和每个面板实际贡献的天数见*_averages.csv；配对掩码见*_correlation_pairs.csv。均值与独立accumarray分组核对，三日中位数与排序后中间元素及分组结果核对，Pearson与中心化向量公式核对。\n');
disp(allCorrelations);
fprintf('Saved figures and correlation results: %s\n',OutputDir);
end

function [average,count,binIndex] = averageWindows(time,values,begin,finish)
%% 显式半开窗口，缺测保留，空窗口不赋零
n = numel(begin);
average = nan(n,4);
count = zeros(n,4);
binIndex = zeros(numel(time),1);
for k = 1:n
    rows = time>=begin(k) & time<finish(k);
    assert(all(binIndex(rows)==0),'Overlapping bins.');
    binIndex(rows) = k;
    for p = 1:4
        x = values(rows,p);
        x = x(isfinite(x));
        count(k,p) = numel(x);
        if ~isempty(x)
            average(k,p) = mean(x);
        end
    end
end
assert(all(binIndex>0),'Some requested days were omitted.');
for p = 1:4
    good = isfinite(values(:,p));
    independent = accumarray(binIndex(good),values(good,p),[n 1],@mean,NaN);
    assert(isequaln(average(:,p),independent),'Window means failed independent check.');
    assert(sum(count(:,p))==nnz(good),'Finite-day count mismatch.');
    assert(all(isnan(average(count(:,p)==0,p))),'Empty window must remain NaN.');
end
end

function [value,count,binIndex] = medianWindows(time,flux,begin,finish)
%% 对原始小时样本直接求三日中位数；保留源记录到窗口的映射。
n = numel(begin);
value = nan(n,1);
count = zeros(n,1);
binIndex = zeros(numel(time),1);
for k = 1:n
    rows = time>=begin(k) & time<finish(k);
    assert(all(binIndex(rows)==0),'Overlapping hourly bins.');
    binIndex(rows) = k;
    x = flux(rows & isfinite(flux));
    count(k) = numel(x);
    if isempty(x)
        continue
    end
    value(k) = median(x);
    ordered = sort(x);
    central = ordered([floor((numel(x)+1)/2),ceil((numel(x)+1)/2)]);
    assert(abs(value(k)-mean(central))<=8*eps(max(1,abs(value(k)))), ...
        'Direct median check failed.');
end
assert(all(binIndex>0),'Unassigned hourly records.');
good = isfinite(flux);
independent = accumarray(binIndex(good),flux(good),[n 1],@median,NaN);
assert(isequaln(value,independent),'Grouped median mismatch.');
assert(sum(count)==nnz(good),'Hourly sample count mismatch.');
end
function [summary,pairUsed] = correlateWindows(values,label)
%% 原始窗口均值上的同窗口Pearson相关
Mode = repmat(string(label),3,1);
Pair = ["a-b";"a-c";"a-d"];
PearsonR = nan(3,1);
ValidWindows = zeros(3,1);
FormulaDifference = nan(3,1);
pairUsed = false(size(values,1),3);
for k = 1:3
    use = isfinite(values(:,1)) & isfinite(values(:,k+1));
    pairUsed(:,k) = use;
    ValidWindows(k) = nnz(use);
    x = values(use,1);
    y = values(use,k+1);
    if numel(x)<2 || all(x==x(1)) || all(y==y(1))
        continue
    end
    r = corrcoef(x,y);
    PearsonR(k) = r(1,2);
    dx = x-mean(x);
    dy = y-mean(y);
    direct = sum(dx.*dy)/sqrt(sum(dx.^2)*sum(dy.^2));
    FormulaDifference(k) = abs(PearsonR(k)-direct);
    assert(FormulaDifference(k)<1e-12,'Pearson formula check failed.');
end
summary = table(Mode,Pair,PearsonR,ValidWindows,FormulaDifference);
end

function files = drawPanels(time,values,begin,finish,plotTitle,tag,out,visible)
%% 四面板，保留原配色、对数轴和细线
visibility = 'off';
if visible
    visibility = 'on';
end
fig = figure('Color','w','Position',[80 80 1600 1050],'Visible',visibility);
layout = tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
labels = {{'Mean of daily |B|','(nT)'}, ...
    {'Mean of daily P1 means','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'Mean of daily P1 medians','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'Mean of daily S1,S2,S3,S5,S6,S7 sums','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}};
if strcmp(tag,'three_day')
    labels{3}{1} = '3-day P1 median';
end
colors = [0.15 0.15 0.15;0.48 0.12 0.62;0.48 0.12 0.62;0.08 0.36 0.62];
ax = gobjects(4,1);
for p = 1:4
    ax(p) = nexttile(layout);
    y = values(:,p);
    if p>1
        y(y<=0) = NaN;
    end
    h = plot(ax(p),time,y,'.-','Color',colors(p,:),'MarkerSize',3,'LineWidth',0.4);
    assert(isequaln(h.YData(:),y));
    if p>1
        set(ax(p),'YScale','log');
    end
    ylabel(ax(p),labels{p});
    set(ax(p),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(p),[begin finish]);
    text(ax(p),0.008,0.89,sprintf('(%c)',96+p),'Units','normalized','FontWeight','bold');
    if p<4
        ax(p).XTickLabel = [];
    end
end
linkaxes(ax,'x');
xlabel(ax(4),'UTC');
xtickformat(ax(4),'yyyy');
title(layout,{'Voyager 1 | 2012-08-25 to 2021-12-16 | P1 0.57-1.78 MeV', ...
    [plotTitle,' | S4 and S8 excluded']},'FontWeight','bold','FontSize',14);
stem = ['V1_',tag,'_4panels_no_S4'];
files = string(fullfile(out,{[stem,'.png'],[stem,'.pdf'],[stem,'.fig']}));
exportgraphics(fig,files(1),'Resolution',220);
exportgraphics(fig,files(2),'ContentType','vector');
savefig(fig,files(3));
end

