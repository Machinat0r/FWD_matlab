function result = Run_V1_Interstellar_Daily_WithS4(visible)
% 另画一张四面板图：panel d 的日均通量求和包含 S4，仍排除 S8。
% 原始 CDF -> 既有日统计 -> S1-S7 求和 -> 独立结果目录。
if nargin<1
    visible = true;
end

%% 路径设置
CodeDir = 'C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/';
ParentDir = 'Z:/SPART-WORK/Data/Voyager/';
OutputDir = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/with_S4/';
addpath([CodeDir,'Interstellar_Daily_Overview_fu']);
if ~isfolder(OutputDir)
    mkdir(OutputDir);
end

%% 从原始 CDF 重算原来的八扇区日均值
source = Run_V1_Interstellar_Daily_Overview( ...
    'DataRoot',ParentDir,'OutputRoot',fullfile(OutputDir,'source_recompute'), ...
    'SyncArchive',false,'Visible',false,'MakePlot',false);
result = source;
result.OutputFolder = string(OutputDir);
result.SourceRecomputeAudit = string(fullfile(source.OutputFolder,'V1_daily_overview_audit.mat'));

%% panel d 加入 S4：七扇区有任一缺测则当日留空
sectorMean = source.Daily.SectorDailyMean;
result.Daily.P1SevenSectorSum = sum(sectorMean(:,1:7),2);
result.Daily.SevenSectorsComplete = all(isfinite(sectorMean(:,1:7)),2);
result.Daily.P1SevenSectorSum(~result.Daily.SevenSectorsComplete) = NaN;
result.Method.PanelD = ['Sum of UTC daily differential fluxes S1-S7, including S4 and excluding S8. ', ...
    'All seven sector daily means must be finite. Original L1-first source selection retained.'];
result.Method.PreviousPanelD = source.Method.PanelD;
result.Method.PreviousSixSectorSumRetained = true;

%% 检查：a/b/c 不变，七扇区求和等于六扇区加 S4
assert(isequaln(result.Daily.BMean_nT,source.Daily.BMean_nT));
assert(isequaln(result.Daily.P1Mean,source.Daily.P1Mean));
assert(isequaln(result.Daily.P1Median,source.Daily.P1Median));
valid = result.Daily.SevenSectorsComplete;
expected = source.Daily.P1SixSectorSum(valid)+sectorMean(valid,4);
actual = result.Daily.P1SevenSectorSum(valid);
roundoff = abs(actual-expected);
assert(all(roundoff<=32*eps(max(1,max(abs(actual),abs(expected))))), ...
    'Seven-sector sum check failed.');
result.Verification = struct('PanelsABCUnchanged',true, ...
    'SevenEqualsSixPlusS4',true,'MaxRoundoff',max(roundoff), ...
    'MissingDaysRetained',all(isnan(result.Daily.P1SevenSectorSum(~valid))));

%% 七扇区图与输出审计
result.DisplayAudit = Plot_V1_Interstellar_Daily_Overview(result,visible,true);
result.WithS4CodeSHA256 = string(Case1_File_SHA256([mfilename('fullpath'),'.m']));
inWindow = result.Daily.EpochUTC>=datetime(2012,8,25,'TimeZone','UTC') & ...
    result.Daily.EpochUTC<datetime(2021,12,17,'TimeZone','UTC');
x = result.Daily.P1SevenSectorSum(inWindow);
result.Verification.PlotWindowDays = nnz(inWindow);
result.Verification.SevenSectorValidDays = nnz(isfinite(x));
result.Verification.SevenSectorPositiveDays = nnz(isfinite(x)&x>0);
writetable(result.Daily,fullfile(OutputDir,'V1_daily_overview_with_S4.csv'),'Encoding','UTF-8');
save(fullfile(OutputDir,'V1_daily_overview_audit.mat'),'result','-v7.3');
disp(result.Verification);
fprintf('New figure including S4: %s\n',OutputDir);
end
