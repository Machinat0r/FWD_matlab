%------written by Wending Fu, Nov.2023 in Beijing------------GSM坐标系
%------revised by Jiqiu Niu, Jan.2025 in Xi'an------------
%------modified by Wending Fu, Mar.2025 in Beijing------------
%% GSE
clear; %清除工作空间中的所有变量
clc; %清空 MATLAB 命令窗口
global ParentDir OutputDir %全局变量
ParentDir = 'Z:\Data\MMS\'; 
DownloadDir = 'C:\Users\pc\Downloads\MMS\';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir); %形成临时目录路径
% ParentDir = 'Z:/Data/MMS/'; 原版

Date = '2016-09-27/2016-09-28';
splitDate = regexp(Date,'/','split');%基于斜杠 /，分割 Date 字符串，返回开始日期和结束日期
ic = 1:4;iic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');%返回日期内所有有的数据，命令行输出"文件目录爬取完成"
filenames2= SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');

filenames = [filenames1,filenames2,filenames3]; %文件名列表
% filenames = filenames1;

% NameTags/FileGroups
expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
desiredLength2 = 11;
desiredLength = 17;
NameTags2 = NameTags(cellfun(@(x) strlength(x) == desiredLength2, NameTags));
NameTags = NameTags(cellfun(@(x) strlength(x) == desiredLength, NameTags));
%%
FileGroups2 = cell(1,length(NameTags2)); 
for j = 1:length(NameTags2)
    FileGroups2{j} = filenames(contains(filenames,NameTags2{j}));
end
FileGroups2 = cellfun(@cellstr,FileGroups2,'UniformOutput',false);%按时间分类整理后的文件名组
%%
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,'EHsearch/',splitDate{1},'To',splitDate{2},'/'];%Z盘里新建文件夹记录xxxx-xx-xxToxxxx-xx-xx[起始终止日期]【DFsearch记得改名】
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end
%分类和目录准备结束

units = irf_units; %%【！！】【常数统一】

%% 
% parpool;
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
%检查当前文件名是否包含特定的字符串（如'fgm'、'des-moms'或'dis-moms'）。如果包含，则执行接下来的代码块。
try
    SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 128, 'CheckSize', 0) %%% 不要直接把数据下到nas里，nas检索速度很慢
catch 
    writematrix([NameTags{TDT}(2:end-2),'的Z盘数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end
%下载文件到tempDir目录，并进行线程数设置、文件大小检查【后期可尝试删掉，影响运行速度】

SDCDataMove(TempDir,ParentDir)

%% get tint
mms.db_init('local_file_db',ParentDir);
formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
% 所得为字符串：两个NameTags元素的时间信息，并格式化为特定的日期时间格式
tempTint=irf.tint(tempDate);
% 字符串转换为时间间隔对象
try
    B1_ts=mms.get_data('B_gse_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
%     Pos = mms.get_data('R_gsm',tint);
%     Pos = Pos.gsmR1;
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    %将错误信息写入到OutputDir下的errorlog.txt文件中，使用追加模式。
    continue
end

%% 强变化电场判据
flag = 0;
mu0 = units.mu0;ep0=units.eps0;
ic=1:4;
try
    c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
    % [Exyz1, Exyz2, Exyz3, Exyz4] = Read_edp_data(FileGroups{TDT},TempDir, tint);
    c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
    c_eval('B_brst0?=irf.ts2mat(Bxyz?);',ic);
    c_eval('E_brst?=irf.ts2mat(Exyz?);',ic);
    c_eval('E_brst?=irf_resamp(E_brst?,E_brst1);',ic);
    
    c_eval('B_brst?=irf_resamp(B_brst0?,E_brst?);',ic);
    c_eval('B_brst?=irf_abs(B_brst?);',ic);
    c_eval('E_fac?=irf_convert_fac(E_brst?,B_brst?,[1,0,0]);',ic);
    c_eval('E_fac? = E_fac?(:, [1,4]);',ic);
    
 
    flagTime = []; tempB = 1; flag = 0;
    deltaTime = 33; % frequency 8192Hz * duration 4ms=4/1e3s，4ms窗口不断滑动
    
    c_eval('Emax? = movmax(E_fac?(:,2), [0, deltaTime]);',ic);
    c_eval('Emin? = movmin(E_fac?(:,2), [0, deltaTime]);',ic);
    c_eval('tempidx? = zeros(length(Emax?),2);',ic);
    c_eval('[~, tempidx?(:,1)] = ismember(Emax1, E_fac1(:, 2));',ic);
    c_eval('[~, tempidx?(:,2)] = ismember(Emin1, E_fac1(:, 2));',ic);
    tempmid=mean(tempidx1, 2);
    tempwindow = zeros(size(tempidx1));
    tempwindow_1=tempmid-17; tempwindow_2=tempmid+17;
    % exceed boundary question
    tempwindow_1(tempwindow_1 <=0) = 1; tempwindow_2(tempwindow_2 >= length(Emax1)) = length(Emax1);

    c_eval('condition1_? = Emax? - Emin? >= 100;',ic);
    condition1 = condition1_1 & condition1_2 & condition1_3 & condition1_4;
    c_eval('condition2_? = abs(tempidx?(:,1) - tempidx?(:,2)) <= 16;',ic);
    condition2 = condition2_1 & condition2_2 & condition2_3 & condition2_4;
    c_eval('condition3_? = all(bsxfun(@ge, tempidx?(:, 1:2), tempwindow_1) & bsxfun(@le, tempidx?(:, 1:2), tempwindow_2), 2);')
    condition3 = condition3_1 & condition3_2 & condition3_3 & condition3_4;
    validIdx = find(condition1 & condition2 & condition3);

    if ~isempty(validIdx)
        flag = 1;
        flagTime = [flagTime;E_fac1(validIdx, 1)];
        writematrix(['在',irf_time(E_fac1(validIdx(1), 1), 'epoch>utc'),'到',irf_time(E_fac1(validIdx(end), 1), 'epoch>utc'),...
            '期间找到了',num2str(length(validIdx)),'个EH'], [OutputDir, 'caselist.txt'], 'WriteMode', 'append', 'Encoding', 'UTF-8')
    end
    
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据检索出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
%% 符合判据的继续下载并出图 同1个时间段文件下载一次，逐个跑图
if flag == 1
try
    goal=[NameTags{TDT}(1:9),NameTags{TDT}(16:17)];
    mask = cellfun(@(x) ~isempty(strfind(x, goal)), NameTags2); %#ok<STREMP> 
    TDT2 = find(mask);
    SDCFilesDownload_NAS(FileGroups2{TDT2},TempDir, 'Threads', 32, 'CheckSize', 1)
    SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 32, 'CheckSize', 1)
    SDCDataMove(TempDir,ParentDir) %下载到tempDir的数据移动到ParentDir目录
    mms.db_init('local_file_db',ParentDir);
catch
     writematrix([NameTags{TDT}(2:end-2),'的画图-数据导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
    for ip = 1:length(flagTime)
        try
        PlotTint = irf_time([flagTime(ip)-0.15,flagTime(ip)+0.15],'epoch>epochTT'); %前后0.3s时间
        PlotTint2 = irf_time([flagTime(ip)-15,flagTime(ip)+15],'epoch>epochTT'); %前后30s时间
        % id_flagTime = SDCPlot(PlotTint,PlotTint2,ic,NameTags{TDT},flagTime(ip));需要返回值的话用这句
        SDCPlot(PlotTint,PlotTint2,ic,NameTags{TDT},flagTime(ip),ip);
        catch
         writematrix([irf_time(flagTime(ip),'epoch>utc'),'的画图调用出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
        continue
        end    
    end
end
end
%%
function [Exyz1, Exyz2, Exyz3, Exyz4] =Read_edp_data(FileGroups,TempDir, tint)
errortag = 1;
while errortag == 1
try
    c_eval(['Exyz?=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],1:4);
    errortag = 0;
catch 
    global ErrorFilePath_fwd_modified ParentDir
    delete(ErrorFilePath_fwd_modified);
    SDCFilesDownload_NAS(FileGroups,TempDir, 'Threads', 128, 'CheckSize', 0)
    SDCDataMove(TempDir,ParentDir)
    mms.db_init('local_file_db',ParentDir);
    clear ErrorFilePath_fwd_modified
end
end
end