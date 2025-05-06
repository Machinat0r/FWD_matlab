%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2019-08-01/2019-08-10';
splitDate = regexp(Date,'/','split');
ic = 1;
% filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');

% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy');
filenames_mec = SDCFilenames(Date,ic,'inst','mec','drm','srvy','dpt','epht89d'); 
filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','dis-moms');
% filenames = [filenames1,filenames2,filenames3,filenames4];
filenames = filenames_fast;
filenames_srvy = [filenames_srvy, filenames_mec];
%%
expr = '_\d+_v'; 
NameTags = regexp(filenames, expr, 'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1, length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = [filenames(contains(filenames, NameTags{j})),...
        filenames_srvy(contains(filenames_srvy, NameTags{j}(2:9)))];
end
FileGroups = cellfun(@cellstr, FileGroups, 'UniformOutput', false);

formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11),':00:00.000Z'];
splitDate = regexp(Date, '/', 'split');
%%
% 定义输出目录
global OutputDir
OutputDir = [ParentDir, 'BBF/', splitDate{1}, 'To', splitDate{2}, '/'];
mkdir(OutputDir);

mms.db_init('local_file_db', ParentDir);
for TDT = 1:length(NameTags) - 1 
clc;
fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2)...
    ,repmat('■',1,round(10*TDT/(length(NameTags)-1))),...
    repmat('□',1,10-round(10*TDT/(length(NameTags)-1))),'\n'])
flag = 0;
units = irf_units;
try
SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 32, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
catch
writematrix([NameTags{TDT}(2:end-2), '的Vi下载出现问题'], [OutputDir, 'errorlog.txt'],...
'WriteMode', 'append', 'Encoding', 'UTF-8')
    continue
end

tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tint = irf.tint(tempDate);

%%
try
    Pos = mms.get_data('R_gsm',tint, 1);
    Pos = Pos.data;
catch
    writematrix([NameTags{TDT}(2:end-2),'的mec数据出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end

%% 判据
if Pos(1,1) <= -10*units.RE/1e3 && abs(Pos(1,2)) <= 10*units.RE/1e3 && abs(Pos(1,3)) <= 10*units.RE/1e3
try
  
     % 获取第一颗卫星的速度数据

    V1_ts = mms.get_data('Vi_gse_fpi_fast_l2', tint, 1);
    V1 = irf.ts2mat(V1_ts);
    V1=irf_gse2gsm(V1);

    Vx = V1(:, 2);
    Vxx = Vx(1:end-1).*Vx(2:end);
    Vxmax = movmax(Vx, 11);
    id = find(Vxx <= 0);
    if ~isempty(id)
    for ii = 1:length(id)
    if Vxmax(id(ii)) >= 200
       flag = 1;
       writematrix([irf_time(V(id(ii),1), 'epoch>utc'), ' ', num2str(Vxmax(id(ii)))],...
            [OutputDir, 'caselist.txt'], 'WriteMode', 'append', 'Encoding', 'UTF-8')
    end
    end
    end
catch
    % 记录错误信息
    writematrix([NameTags{TDT}(2:end-2), '的数据导入2出现问题'], [OutputDir, 'errorlog.txt'],...
        'WriteMode', 'append', 'Encoding', 'UTF-8')
    continue
end
end
%% 符合判据的继续下载并出图
if flag == 1
try
    mms.db_init('local_file_db',ParentDir);
    for ii = 1:length(id)
        flagTime = V1(id(ii));
        PlotTint = irf_time([flagTime(end)-40,flagTime(end)+40],'epoch>epochTT');
        Plot_BV(PlotTint,ic,NameTags{TDT}(2:end-2),flagTime);
    end
catch
    writematrix([irf_time(flagTime,'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
end
end
end