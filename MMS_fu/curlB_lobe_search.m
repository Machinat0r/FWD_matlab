%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2016-01-01/2023-01-01';
splitDate = regexp(Date,'/','split');
ic = 1:4;iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','srvy');
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
OutputDir = [ParentDir,'CurlB_Search/',splitDate{1},'To',splitDate{2},'/'];
mkdir(OutputDir)
% if ~isfolder([OutputDir,'OverviewFig/'])
%     mkdir([OutputDir,'OverviewFig/']);
% end

mms.db_init('local_file_db',ParentDir);

parfor_progress(length(NameTags)-1);
parfor TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;
parfor_progress;
units = irf_units;
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
% formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T00:00:00.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tempTint=irf.tint(tempDate);

tryTimes = 0;
while tryTimes <= 10
try
    % B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    % tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));  
    mms.db_init('local_file_db','D:/MMS/');
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

Id = find(Pos(:,2) <= -10*units.RE/1e3 & abs(Pos(:,3)) <= 10*units.RE/1e3...
    & abs(Pos(:,4)) >= 5*units.RE/1e3);
%% 判据
if ~isempty(Id)
try
    tint = irf.tint(irf_time(Pos(Id(1),1),'epoch>epochTT'),irf_time(Pos(Id(end),1),'epoch>epochTT'));
    % c_eval can be used only in for not parfor
    % % % c_eval("B?_ts=mms.get_data('B_gsm_brst',tint,?);",1:4);
    % % % c_eval('B? = irf.ts2mat(B?_ts);',ic);
    % % % c_eval('R? = Pos.gsmR?;',ic)
    % % % c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];',ic)

    B1_ts=mms.get_data('B_gsm_srvy',tint,1);B2_ts=mms.get_data('B_gsm_srvy',tint,2);
    B3_ts=mms.get_data('B_gsm_srvy',tint,3);B4_ts=mms.get_data('B_gsm_srvy',tint,4);
    B1 = irf.ts2mat(B1_ts);B2 = irf.ts2mat(B2_ts);
    B3 = irf.ts2mat(B3_ts);B4 = irf.ts2mat(B4_ts);

    % % % tStart = max([B1(1),B2(1),B3(1),B4(1)]);
    % % % B1(B1(:,1)<=tStart,:)=[];B2(B2(:,1)<=tStart,:)=[];
    % % % B3(B3(:,1)<=tStart,:)=[];B4(B4(:,1)<=tStart,:)=[];
    % % % 
    % % % tEnd = min([B1(end,1),B2(end,1),B3(end,1),B4(end,1)]);
    % % % B1(B1(:,1)>=tEnd,:)=[];B2(B2(:,1)>=tEnd,:)=[];
    % % % B3(B3(:,1)>=tEnd,:)=[];B4(B4(:,1)>=tEnd,:)=[];
    
    Pos = mms.get_data('R_gsm',tint);
    R1 = Pos.gsmR1;R2 = Pos.gsmR2;R3 = Pos.gsmR3;R4 = Pos.gsmR4;
    R1 = [Pos.time.epochUnix R1(:,1:3)];R2 = [Pos.time.epochUnix R2(:,1:3)];
    R3 = [Pos.time.epochUnix R3(:,1:3)];R4 = [Pos.time.epochUnix R4(:,1:3)];

    % [curlB,~]=c_4_grad('R?','B?','curl');
    [Jcurl,~] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);
    Jcurl = irf_abs(Jcurl)*1e9;

    IdSep = round(size(Jcurl,1)/5);
    [maxJ1,maxId1] = max(Jcurl(1:IdSep,5));         flagTime1 = B1(maxId1,1);
    [maxJ2,maxId2] = max(Jcurl(IdSep:2*IdSep,5));   flagTime2 = B1(IdSep + maxId2,1);
    [maxJ3,maxId3] = max(Jcurl(2*IdSep:3*IdSep,5)); flagTime3 = B1(2*IdSep + maxId3,1);
    [maxJ4,maxId4] = max(Jcurl(3*IdSep:4*IdSep,5)); flagTime4 = B1(3*IdSep + maxId4,1);
    [maxJ5,maxId5] = max(Jcurl(4*IdSep:end,5));     flagTime5 = B1(4*IdSep + maxId5,1);

    writematrix([irf_time(flagTime1,'epoch>utc'),' ',num2str(maxJ1)],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    writematrix([irf_time(flagTime2,'epoch>utc'),' ',num2str(maxJ2)],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    writematrix([irf_time(flagTime3,'epoch>utc'),' ',num2str(maxJ3)],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    writematrix([irf_time(flagTime4,'epoch>utc'),' ',num2str(maxJ4)],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    writematrix([irf_time(flagTime5,'epoch>utc'),' ',num2str(maxJ5)],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入2出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
end
%% 符合判据的继续下载并出图
% % % if flag == 1
% % % try
% % %     SDCFilesDownload(FileGroups{TDT},TempDir) 
% % %     SDCDataMove(TempDir,ParentDir)
% % %     mms.db_init('local_file_db',ParentDir);
% % %     PlotTint = irf_time([flagTime(end)-20,flagTime(end)+20],'epoch>epochTT');
% % %     id_flagTime = SDCPlot(PlotTint,desmoms,ic,NameTags{TDT},flagTime(end));
% % % catch
% % %     writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
% % %     'WriteMode','append','Encoding','UTF-8')
% % %     continue
% % % end
% % % end
%% 删除文件夹并生成记录文件
% % % try
% % %     cd(OutputDir)
% % %     rmdir(TempDir,'s');    
% % % catch
% % %     writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % % end
end
parfor_progress(0);