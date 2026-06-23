%------written by Wending Fu, Oct.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
ic = 1:4;iic = 1:4;
load([pwd,'/NameTags.mat']);

Date = '2024-01-01T00:00:00.000Z/2025-01-01T00:00:00.000Z';
% Date = '2016-02-05T00:00:00.000Z/2016-03-01T00:00:00.000Z';
splitDate = regexp(Date,'/','split');
startDate = splitDate{1};
endDate   = splitDate{2};
startDT = datetime(startDate, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
endDT   = datetime(endDate,   'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
NameTagsDT = datetime(NameTags, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
mask = (NameTagsDT >= startDT) & (NameTagsDT <= endDT);
NameTags = NameTags(mask);
%% download data
% filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
% filenames = [filenames1,filenames2];
% filenames = filenames1;
% 
% expr = '_[0-9]+\_v';
% NameTags = regexp(filenames,expr,'match');
% 
% NameTags = unique(cellfun(@cellstr,NameTags));
% % % FileGroups = cell(1,length(NameTags)); 
% % % for j = 1:length(NameTags)
% % %     % FileGroups{j} = [filenames(contains(filenames,NameTags{j})), filenames_srvy(contains(filenames_srvy,NameTags{j}(2:9)))];
% % %     FileGroups{j} = filenames(contains(filenames,NameTags{j}));
% % % end
% % % FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组
%%
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
splitDate{1} = erase(splitDate{1},{':','.'});splitDate{2} = erase(splitDate{2},{':','.'});
OutputDir = [ParentDir,'Monopole_Search/maxB_Search/',splitDate{1},'To',splitDate{2},'/'];
clear mask NameTagsDT
mkdir(OutputDir)
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

% % % if length(FileGroups{TDT})~=4, continue;end
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
% % % formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
% % % tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
tempTint=irf.tint(tempDate);

tryTimes = 0;
while tryTimes <= 10
try
    % B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    % tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));  
    mms.db_init('local_file_db',ParentDir);
    Pos = mms.get_data('R_gsm',tempTint,1);
    Pos = irf.ts2mat(Pos);
    Pos = irf_abs(Pos);
    tryTimes = 666;
catch
    tryTimes = tryTimes + 1;
    continue
end
end

if tryTimes ~= 666
writematrix([NameTags{TDT}(1:end-1),'的数据导入1出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
end

%% 判据
try
    if Pos(1,2) < 0 || Pos(1,5)./units.RE*1e3 < 10, continue; end
    tint = irf.tint(irf_time(Pos(1,1),'epoch>epochTT'),irf_time(Pos(end,1),'epoch>epochTT'));

    % c_eval can be used only in for not parfor
    B1_ts=mms.get_data('B_gsm_brst',tint,1);B2_ts=mms.get_data('B_gsm_brst',tint,2);
    B3_ts=mms.get_data('B_gsm_brst',tint,3);B4_ts=mms.get_data('B_gsm_brst',tint,4);
    B1 = irf.ts2mat(B1_ts);B2 = irf.ts2mat(B2_ts);
    B3 = irf.ts2mat(B3_ts);B4 = irf.ts2mat(B4_ts);
    c_eval('Bt? = irf_abs(B?);',1:4)
    
    Pos = mms.get_data('R_gsm',tint);
    R1 = Pos.gsmR1;R2 = Pos.gsmR2; R3 = Pos.gsmR3;R4 = Pos.gsmR4;
    R1 = [Pos.time.epochUnix R1(:,1:3)]; R2 = [Pos.time.epochUnix R2(:,1:3)];
    R3 = [Pos.time.epochUnix R3(:,1:3)]; R4 = [Pos.time.epochUnix R4(:,1:3)];

    % c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
    % c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
    % c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
    % c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
    % 
    % c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
    % c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
    % c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
    % c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
    % c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);
    % 
    % 
    % 
    % c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
    % c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
    % c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
    % c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
    % c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);
    % 
    % c_eval('Ne? = irf_resamp(Ne?, B?);',ic);c_eval('Te? = irf_resamp(Te?, B?);',ic);
    % c_eval('Ni? = irf_resamp(Ni?, B?);',ic);c_eval('Ti? = irf_resamp(Ti?, B?);',ic);
    % 
    % c_eval('Pthe? = units.e*10^(15)*Ni?(:,2).*Ti?(:,2) + units.e*10^(15)*Ne?(:,2).*Te?(:,2);',ic);
    % c_eval('Pm? = 10^(-9)*Bt?(:,5).^2 / (2*units.mu0);',ic); %nPa
    % c_eval('Beta? = Pthe? ./ Pm?;',ic);
    % 
    % c_eval('Flag_msp? = find(Beta? <= 1);',ic);

    c_eval('Bmsh? = Bt?(:,[1,5]);',ic);
    % c_eval('Bmsh?(Flag_msp?,2) = 0;',ic);
    c_eval('Flag_StrongB? = find(Bmsh?(:,2) >= 100);',ic);

    if isempty(Flag_StrongB1) & isempty(Flag_StrongB2) & isempty(Flag_StrongB3) & isempty(Flag_StrongB4)
        writematrix([NameTags{TDT}(1:end-1) , '无maxB>100事件'],...
            [OutputDir,'Log.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
    end 

    if ~isempty(Flag_StrongB1)
        c_eval('[maxB?, id_maxB?] = max(Bmsh?(Flag_StrongB?,2));',1);
        writematrix([irf_time(B1(id_maxB1,1),'epoch>utc'),' maxB1 = ', num2str(maxB1)],...
            [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    if ~isempty(Flag_StrongB2)
        c_eval('[maxB?, id_maxB?] = max(Bmsh?(Flag_StrongB?,2));',2);
        writematrix([irf_time(B2(id_maxB2,1),'epoch>utc'),' maxB2 = ', num2str(maxB2)],...
            [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    if ~isempty(Flag_StrongB3)
        c_eval('[maxB?, id_maxB?] = max(Bmsh?(Flag_StrongB?,2));',3);
        writematrix([irf_time(B3(id_maxB3,1),'epoch>utc'),' maxB3 = ', num2str(maxB3)],...
            [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end

    if ~isempty(Flag_StrongB4)
        c_eval('[maxB?, id_maxB?] = max(Bmsh?(Flag_StrongB?,2));',3);
        writematrix([irf_time(B4(id_maxB4,1),'epoch>utc'),' maxB4 = ', num2str(maxB4)],...
            [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end

catch
    writematrix([NameTags{TDT}(1:end-1),'的数据导入2出现问题'],[OutputDir,'errorlog.txt'],...
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