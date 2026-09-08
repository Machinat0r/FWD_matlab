clear;clc
% Voyager 1 日球层顶外四面板；参照 MMS_fu 手写程序的分节布局。
% 原始 CDF -> 日统计 -> 绘图。缺测、L1/L2 与能段处理见同目录 README。

%% 路径设置
CodeDir = 'C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/';
ParentDir = 'Z:/SPART-WORK/Data/Voyager/';
OutputDir = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/';
addpath([CodeDir,'Interstellar_Daily_Overview_fu']);

%% 运行设置
% true：先检查 NASA 目录并下载缺少的原始 CDF。
% false：直接读取本地原始 CDF；已核查的正式数据可以使用此项。
SyncArchive = false;
Visible = true;

%% 读取原始 CDF、计算日统计并画图
result = Run_V1_Interstellar_Daily_Overview( ...
    'DataRoot',ParentDir,'OutputRoot',OutputDir, ...
    'SyncArchive',SyncArchive,'Visible',Visible);

%% 结果
% 图中只显示 2012-08-25 至 2021-12-16；b/c/d 为对数轴。
% 细线连接相邻有效日，缺测处断开；不插值、不平滑。
Daily = result.Daily;
Coverage = result.Coverage;
disp(OutputDir);
