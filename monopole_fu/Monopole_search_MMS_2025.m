%------written by Wending Fu, May.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = '/Volumes/SPART-NAS/Data/MMS/'; 
DownloadDir = '/Volumes/SPART-NAS/Data/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
% Date = '2019-09-24/2019-09-25';
Date = '2019-01-16/2019-01-17';
splitDate = regexp(Date,'/','split');
ic = 1:4;iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
% filenames = [filenames1,filenames2,filenames3,filenames4];
filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    % FileGroups{j} = [filenames(contains(filenames,NameTags{j})), filenames_srvy(contains(filenames_srvy,NameTags{j}(2:9)))];
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,'Monopole_Search/',splitDate{1},'To',splitDate{2},'/'];
mkdir(OutputDir)
% if ~isfolder([OutputDir,'OverviewFig/'])
%     mkdir([OutputDir,'OverviewFig/']);
% end

mms.db_init('local_file_db',ParentDir);

%%
units = irf_units;
% parfor_progress(length(NameTags)-1);
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;clear B1 B2 B3 B4 R1 R2 R3 R4 Pos;
PlotFlag = 0;
% parfor_progress;
% disp(['૮₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*TDT/length(NameTags))),repmat('□',1,10-round(10*TDT/length(NameTags))),'正在光速检索ing...'])
% fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])

if length(FileGroups{TDT})~=4, continue;end
% % % try
% % %     SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 32, 'CheckSize', 0)
% % % catch
% % %     writematrix([NameTags{TDT}(2:end-2),'的Z盘数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % %     continue
% % % end
% % % 
% % % SDCDataMove(TempDir,ParentDir)
% % % 
% % % mms.db_init('local_file_db',ParentDir);
formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
% formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T00:00:00.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tempTint=irf.tint(tempDate);

tryTimes = 0;
while tryTimes <= 10
try
    % B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    % tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));  
    mms.db_init('local_file_db','/Volumes/SPART-NAS/Data/MMS/');
    Pos = mms.get_data('R_gsm',tempTint,1);
    Pos = irf.ts2mat(Pos);
    tryTimes = 666;
catch
    tryTimes = tryTimes + 1;
    continue
end
end

if tryTimes ~= 666
writematrix([NameTags{TDT}(2:end-2),'的数据导入1出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
end

%% 判据
try
    tint = irf.tint(irf_time(Pos(1,1),'epoch>epochTT'),irf_time(Pos(end,1),'epoch>epochTT'));
    % c_eval can be used only in for not parfor
    B1_ts=mms.get_data('B_gsm_brst',tint,1);B2_ts=mms.get_data('B_gsm_brst',tint,2);
    B3_ts=mms.get_data('B_gsm_brst',tint,3);B4_ts=mms.get_data('B_gsm_brst',tint,4);
    B1 = irf.ts2mat(B1_ts);B2 = irf.ts2mat(B2_ts);
    B3 = irf.ts2mat(B3_ts);B4 = irf.ts2mat(B4_ts);
    c_eval('B? = irf_resamp(B?,B1);');
    
    Pos = mms.get_data('R_gsm',tint);
    R1 = Pos.gsmR1;R2 = Pos.gsmR2;R3 = Pos.gsmR3;R4 = Pos.gsmR4;
    R1 = [Pos.time.epochUnix R1(:,1:3)];R2 = [Pos.time.epochUnix R2(:,1:3)];
    R3 = [Pos.time.epochUnix R3(:,1:3)];R4 = [Pos.time.epochUnix R4(:,1:3)];

    B2 = irf_resamp(B2, B1); B3 = irf_resamp(B3, B1); B4 = irf_resamp(B4, B1);

    c_eval('R? = irf_resamp(R?,B1);')
    CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
    c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');
    
    LocPoint = zeros(length(B1),3)*nan;
    LocRes = cell(length(B1),1);
    Q = zeros(length(B1),1)*nan;
    resQ = cell(length(B1),1);
    
    Qerror = ones(length(B1),1)*1000;
    dLoc = ones(length(B1),15)*1000;
    
    % div
    gradB=c_4_grad(R1, R2, R3, R4, B1, B2, B3, B4, 'grad');
    divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

    flag_m = 0;
    time_flagm = 0;    
    
    RR12 = irf_abs(R1-R2); RR13 = irf_abs(R1-R3); RR14 = irf_abs(R1-R4); 
    RR23 = irf_abs(R2-R3); RR24 = irf_abs(R2-R4); RR34 = irf_abs(R3-R4); 
    RR_mean = (RR12(:,5) + RR13(:,5) + RR14(:,5) + RR23(:,5) + RR24(:,5) + RR34(:,5))/6;
    MultiPower = ceil(max([log10(RR_mean)]));
    if MultiPower > 3, continue; end
    
    id = nchoosek(1:6,2);
    % solve
    tic
    [Q_py, Loc_py] = py.cal_error_gpu.CalError(R1_list, R2_list, R3_list, R4_list, ...
                                                  B1_list, B2_list, B3_list, B4_list, ...
                                                  py.None, py.None, py.None);
    % % % parfor_progress(size(B1, 1)-1);
    % % % parfor i =1:size(B1, 1)
    % % % parfor_progress;
    % % % [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError(R1, R2, R3, R4,B1, B2, B3, B4,...
    % % %     i,i*sign(divB(i,2)),RR_mean(i),1);
    % % % 
    % % % tempd = irf_abs(LocRes{i}(id(:,1),:)-LocRes{i}(id(:,2),:));
    % % % tempd = tempd(:,4)/RR_mean(i);
    % % % tempd = tempd';
    % % % dLoc(i,:) = tempd * 100;
    % % % 
    % % % if ~isnan(resQ{i})
    % % %     Qerror(i) = abs(100*std(resQ{i})/Q(i));
    % % % else
    % % %     Qerror(i) = 1000;
    % % % end
    % % % end
    parfor_progress(0);
    toc
    meand = mean(dLoc,2);
    LocPoint = irf_abs(LocPoint);


    writematrix([irf_time(B1(1),'epoch>utc'),' ',num2str(min(meand)),' ', num2str(min(Qerror))...
        ,' ', num2str(min(LocPoint(:,4)))],[OutputDir,'DateList.txt'],'WriteMode','append','Encoding','UTF-8')
    TimeUTC = irf_time(B1(1),'epoch>utc');TimeUTC = erase(TimeUTC,{':','.'});
    save([OutputDir, TimeUTC, '.mat'], 'Q', 'resQ', 'LocPoint', 'LocRes')
    
    try
    tag = find(meand <= 66.1 & Qerror <= 100 & LocPoint(:,4) <= 3*mean(RR_mean));
    if ~isempty(tag)
    PlotFlag = 1;
    flagTime = B1(tag, 1);
    flagTimeUTC = irf_time(B1(1),'epoch>utc');flagTimeUTC = erase(flagTimeUTC,{':','.'});
    writematrix([flagTimeUTC,' ',num2str(meand(tag)),' ', num2str(Qerror(tag))...
    ,' ', num2str(LocPoint(tag,4))],[OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    catch
    writematrix([NameTags{TDT}(2:end-2),'的caselist导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入2出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end

%% 符合判据的继续下载并出图
% if PlotFlag == 1
% try
%     SDCFilesDownload(FileGroups{TDT},TempDir) 
%     SDCDataMove(TempDir,ParentDir)
%     mms.db_init('local_file_db',ParentDir);
%     PlotTint = irf_time([flagTime(end)-20,flagTime(end)+20],'epoch>epochTT');
%     id_flagTime = SDCPlot(PlotTint,desmoms,ic,NameTags{TDT},flagTime(end));
% catch
%     writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
%     'WriteMode','append','Encoding','UTF-8')
%     continue
% end
% end
%% 删除文件夹并生成记录文件
% % % try
% % %     cd(OutputDir)
% % %     rmdir(TempDir,'s');    
% % % catch
% % %     writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % % end
end
parfor_progress(0);