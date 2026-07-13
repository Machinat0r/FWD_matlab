%% ACE/MFI 全时段逐小时平均值与最大值
% 直接运行本脚本即可。原始 Z: 盘数据只读，不会被修改。

clear;
clc;

thisFolder = string(fileparts(mfilename('fullpath')));
addpath(thisFolder);

%% 用户设置
dataRoot = "Z:\SPART-WORK\Data\ACE";

% 结果与年度缓存默认放在本程序旁边。也可以改成本地 SSD 上的其他目录。
outputDir = fullfile(thisFolder, "ACE_MFI_hourly_output");

% NaT 表示自动采用数据库文件名中的第一天和最后一天。
startDate = NaT;
endDate = NaT;

% true：忽略已有年度缓存并重新计算；正常运行保持 false。
forceRebuild = false;

%% 执行
[Hourly, RunInfo, DailySelection, FileLog] = ace_mfi_hourly( ...
    dataRoot, outputDir, ...
    'StartDate', startDate, ...
    'EndDate', endDate, ...
    'ForceRebuild', forceRebuild, ...
    'WriteCSV', true, ...
    'MakeFigure', true, ...
    'FigureVisible', 'off', ...
    'CheckpointEveryDays', 30, ...
    'Verbose', true);

disp(RunInfo);

