function audit = Correlate_V1_Interstellar_Daily(includeS4)
% Voyager 1：日均磁场与另外三个面板的同日 Pearson 相关系数。
% 原始 CDF 读取、日统计和缺测处理沿用当前四面板程序。
% 原始日值计算，不取对数，不去趋势、不平滑、不移时。

% includeS4=true：panel d 包含 S4，S1-S7 求和；无参数保留原六扇区版本。
if nargin<1
    includeS4 = false;
end
assert(islogical(includeS4) && isscalar(includeS4),'includeS4 must be a logical scalar.');

%% 路径与时间范围
CodeDir = 'C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/';
ParentDir = 'Z:/SPART-WORK/Data/Voyager/';
OutputDir = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/';
addpath([CodeDir,'Interstellar_Daily_Overview_fu']);
ResultDir = fullfile(OutputDir,'correlation');
if includeS4
    ResultDir = fullfile(OutputDir,'with_S4','correlation');
end
if ~isfolder(ResultDir)
    mkdir(ResultDir);
end
StartUTC = datetime(2012,8,25,'TimeZone','UTC');
EndUTC = datetime(2021,12,17,'TimeZone','UTC');

%% 从原始 CDF 重新计算日统计
result = Run_V1_Interstellar_Daily_Overview( ...
    'DataRoot',ParentDir,'OutputRoot',fullfile(ResultDir,'source_recompute'), ...
    'SyncArchive',false,'Visible',false,'MakePlot',false);
useTime = result.Daily.EpochUTC>=StartUTC & result.Daily.EpochUTC<EndUTC;
daily = result.Daily(useTime,:);
assert(all(diff(daily.EpochUTC)==days(1)),'Daily UTC grid is incomplete.');
B = daily.BMean_nT;
if includeS4
    daily.P1SevenSectorSum = sum(daily.SectorDailyMean(:,1:7),2);
    daily.SevenSectorsComplete = all(isfinite(daily.SectorDailyMean(:,1:7)),2);
    daily.P1SevenSectorSum(~daily.SevenSectorsComplete) = NaN;
    result.Method.PreviousPanelD = result.Method.PanelD;
    result.Method.PanelD = 'Sum of finite UTC daily fluxes S1-S7, including S4 and excluding S8; all seven sectors required.';
end

%% 三组同日配对相关：每组分别排除缺测
fields = {'P1Mean','P1Median','P1SixSectorSum'};
Pair = ["a-b";"a-c";"a-d"];
Description = ["日均磁场 vs 质子日均通量"; ...
    "日均磁场 vs 质子日中位通量";"日均磁场 vs 六扇区日均通量之和"];
if includeS4
    fields{3} = 'P1SevenSectorSum';
    Description(3) = "日均磁场 vs 七扇区日均通量之和（含S4）";
end
PearsonR = nan(3,1);
ValidDays = zeros(3,1);
MissingPairDays = zeros(3,1);
IncludedZeroDays = zeros(3,1);
FirstPairedUTC = NaT(3,1,'TimeZone','UTC');
LastPairedUTC = FirstPairedUTC;
FormulaDifference = nan(3,1);
pairs = struct;
for k = 1:3
    J = daily.(fields{k});
    use = isfinite(B) & isfinite(J);
    x = B(use);
    y = J(use);
    ValidDays(k) = numel(x);
    MissingPairDays(k) = height(daily)-ValidDays(k);
    IncludedZeroDays(k) = nnz(x==0 | y==0);
    pairs.(fields{k}) = table(daily.EpochUTC,B,J,use, ...
        'VariableNames',{'EpochUTC','BMean_nT','ParticleDailyValue','Used'});
    if ~isempty(x)
        times = daily.EpochUTC(use);
        FirstPairedUTC(k) = times(1);
        LastPairedUTC(k) = times(end);
    end
    if numel(x)<2 || all(x==x(1)) || all(y==y(1))
        continue
    end
    R = corrcoef(x,y);
    PearsonR(k) = R(1,2);
    % 以中心化向量公式独立核对 corrcoef，阈值仅用于浮点校验。
    dx = x-mean(x);
    dy = y-mean(y);
    directR = sum(dx.*dy)/sqrt(sum(dx.^2)*sum(dy.^2));
    FormulaDifference(k) = abs(PearsonR(k)-directR);
    assert(FormulaDifference(k)<1e-12,'Pearson formula consistency failed.');
end
summary = table(Pair,Description,PearsonR,ValidDays,MissingPairDays, ...
    IncludedZeroDays,FirstPairedUTC,LastPairedUTC,FormulaDifference);

%% 保存结果、配对日和来源记录
method = struct;
method.PanelDSectors = [1 2 3 5 6 7];
if includeS4
    method.PanelDSectors = 1:7;
end
method.Coefficient = 'Pearson, zero lag, unweighted original daily values';
method.TimeWindow = '[2012-08-25 00:00 UTC,2021-12-17 00:00 UTC)';
method.Pairing = 'Pairwise finite same-UTC-day values; each pair may use different dates.';
method.Transform = 'No log transformation, detrending, smoothing, interpolation, lag search or outlier rejection; source-valid zero values included.';
method.Significance = 'Descriptive coefficients only. No independent-day p-values or causal interpretation.';
method.IRFUReview = 'irf_corr_deriv concerns derivatives/time markers; MATLAB corrcoef used for same-day Pearson correlation.';
audit = struct('CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'Summary',summary,'Method',method,'Daily',daily,'Pairs',pairs, ...
    'DailySourceMethod',result.Method,'COHOSources',result.COHOSources, ...
    'SectorAuditFiles',result.SectorAuditFiles, ...
    'DailyRunCodeFiles',result.CodeFiles,'DailyRunCodeSHA256',result.CodeSHA256, ...
    'CorrelationCodeSHA256',string(Case1_File_SHA256([mfilename('fullpath'),'.m'])));
save(fullfile(ResultDir,'V1_daily_Pearson_audit.mat'),'audit','-v7.3');
writetable(summary,fullfile(ResultDir,'V1_daily_Pearson_summary.csv'),'Encoding','UTF-8');
for k = 1:3
    writetable(pairs.(fields{k}),fullfile(ResultDir,['paired_days_',fields{k},'.csv']), ...
        'Encoding','UTF-8');
end
fid = fopen(fullfile(ResultDir,'磁场与三个面板的相关系数.md'),'w','n','UTF-8');
assert(fid>=0,'Cannot create correlation report.');
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 日均磁场与三个粒子面板的相关系数\n\n');
fprintf(fid,'时间：2012-08-25 至 2021-12-16（含）。Pearson 相关，原始日值，未取对数。每组按同日有限值分别配对，零值保留；缺测不插值。\n\n');
fprintf(fid,'| 配对 | Pearson r | 有效天数 N |\n|---|---:|---:|\n');
for k = 1:3
    fprintf(fid,'| %s | %.6f | %d |\n',char(Description(k)),PearsonR(k),ValidDays(k));
end
fprintf(fid,'\n未去趋势、未平滑、未移时、未删除大峰值。相关系数描述该时段整体线性关系；日序列的时间自相关没有用于显著性推断。\n');
fprintf(fid,'\n程序从原始 CDF 重新计算日统计。参与及未参与日期见 paired_days_*.csv；完整来源和方法见 V1_daily_Pearson_audit.mat。\n');
disp(summary(:,{'Pair','PearsonR','ValidDays','IncludedZeroDays'}));
fprintf('Correlation results: %s\n',ResultDir);
end


