%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir OutputDir
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2016-01-01/2017-06-27';
splitDate = regexp(Date,'/','split');
ic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
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
OutputDir = [ParentDir,'TCS_Search/',splitDate{1},'To',splitDate{2},'/'];
mkdir(OutputDir)
mkdir([OutputDir,'Figs/']);

mms.db_init('local_file_db',ParentDir);

% parfor_progress(length(NameTags)-1);
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;
% parfor_progress;
units = irf_units;
% disp(['૮₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*TDT/length(NameTags))),repmat('□',1,10-round(10*TDT/length(NameTags))),'正在光速检索ing...'])
% fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])

flag = 0;
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
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tempTint=irf.tint(tempDate);

try
    B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
    Pos = mms.get_data('R_gsm',tint,1);
    Pos = Pos.data;
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入1出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end

%% 判据
if Pos(1,1) <= -10*units.RE/1e3 && abs(Pos(1,2)) <= 15*units.RE/1e3 && abs(Pos(1,3)) <= 5*units.RE/1e3
try
    B1_ts=mms.get_data('B_gsm_brst',tint,1);B2_ts=mms.get_data('B_gsm_brst',tint,2);
    B3_ts=mms.get_data('B_gsm_brst',tint,3);B4_ts=mms.get_data('B_gsm_brst',tint,4);
    B1 = irf.ts2mat(B1_ts);B2 = irf.ts2mat(B2_ts);
    B3 = irf.ts2mat(B3_ts);B4 = irf.ts2mat(B4_ts);

    tStart = max([B1(1),B2(1),B3(1),B4(1)]);
    B1(B1(:,1)<=tStart,:)=[];B2(B2(:,1)<=tStart,:)=[];
    B3(B3(:,1)<=tStart,:)=[];B4(B4(:,1)<=tStart,:)=[];

    tEnd = min([B1(end,1),B2(end,1),B3(end,1),B4(end,1)]);
    B1(B1(:,1)>=tEnd,:)=[];B2(B2(:,1)>=tEnd,:)=[];
    B3(B3(:,1)>=tEnd,:)=[];B4(B4(:,1)>=tEnd,:)=[];
    
    Pos = mms.get_data('R_gsm',tint);
    R1 = Pos.gsmR1;R2 = Pos.gsmR2;R3 = Pos.gsmR3;R4 = Pos.gsmR4;
    R1 = [Pos.time.epochUnix R1(:,1:3)];R2 = [Pos.time.epochUnix R2(:,1:3)];
    R3 = [Pos.time.epochUnix R3(:,1:3)];R4 = [Pos.time.epochUnix R4(:,1:3)];

    % [curlB,~]=c_4_grad('R?','B?','curl');
    [Jcurl,~] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);
    Jcurl = irf_abs(Jcurl);
    Jcurl(:,2:5) = Jcurl(:,2:5)*1e9;
    Jt = Jcurl(:,5);
    
    width = 257; % 128Hz, 2s
    J_std = movstd(Jt, width);
    J_mean = movmean(Jt, width);
    J_CV = J_std./J_mean;
    criter = 0.85; %%%%%%%%CRITERION OF CV%%%%%%%%
    Id_CV = find(Jt >= 10 & J_CV >= criter);
    
    Id_start = max([ones(size(Id_CV)), Id_CV - 128], [], 2);
    Id_end = min([size(J_CV,1).*ones(size(Id_CV)), Id_CV + 128], [], 2);
    Id_Duration = [];Time_Duration = [];
    for ii = 1:length(Id_CV)
    [J_max, Id_max] = max(Jcurl(Id_start(ii):Id_end(ii), 5));
    Id_max = Id_max + Id_start(ii) -1;
    J_max2 = J_max/sqrt(2);

    left_window = Jcurl(Id_start(ii):Id_max, 5);
    mask_left = left_window > J_max2;
    Id_left = find(~mask_left, 1, 'last');
    if isempty(Id_left), Id_left = Id_start(ii);
    else,Id_left = Id_left + Id_start(ii);end

    right_window = Jcurl(Id_max:Id_end(ii), 5);
    mask_right = right_window > J_max2;
    Id_right = find(~mask_right, 1, 'first');
    if isempty(Id_right), Id_right = Id_end(ii);
    else,Id_right = Id_right + Id_max - 3;end

    TimeLength = Jcurl(Id_right, 1) - Jcurl(Id_left, 1);


    if TimeLength <= 0.2 %%%%%%%%CRITERION OF DURATION%%%%%%%%
        Id_Duration = [Id_Duration, Id_CV(ii)];
        Time_Duration = [Time_Duration, TimeLength];
    end
    end  

    Time_Duration = Time_Duration(diff([0,Id_Duration]) ~= 1);
    Id_Duration = Id_Duration(diff([0,Id_Duration]) ~= 1); % only save the first one of the time series 
    if ~isempty(Id_Duration)
    flag = 1;
    flagTime = B1(Id_Duration, 1);
    for iii = 1:length(Id_Duration)
    writematrix([irf_time(flagTime(iii),'epoch>utc'),' ',num2str(Jt(Id_Duration(iii))),...
        ' ', num2str(Time_Duration(iii))],...
        [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入2出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
end
%% 符合判据的继续下载并出图
if flag == 1
try
    % mms.db_init('local_file_db',ParentDir);
    for ii = 1:length(Id_Duration)
        PlotTint = irf_time([flagTime(ii)-10,flagTime(ii)+10],'epoch>epochTT');
        TimeStamp = irf_time(flagTime(ii),'epoch>utc');
        TimeStamp = regexprep(TimeStamp, '[:.]', '-');
        Plot_TCS(PlotTint,ic,TimeStamp,flagTime(ii));
    end
catch
    writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
end
end
end
% parfor_progress(0);