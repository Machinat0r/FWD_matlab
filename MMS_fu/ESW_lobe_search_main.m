%------written by Wending Fu, July.2025 in Singapore------------
clear;clc;
global ParentDir OutputDir
ParentDir = 'D:\MMS/'; 
DownloadDir = 'C:\MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

Date = '2019-01-01/2025-01-01';
splitDate = regexp(Date,'/','split');
ic = 1;iic = 1;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2, filenames3, filenames4];
% filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [DownloadDir,'ESWSearch/',splitDate{1},'To',splitDate{2},'/'];
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end

units = irf_units;
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
% SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 64, 'CheckSize', 0)
% SDCDataMove(TempDir,ParentDir)

formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
% tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
tempTint=irf.tint(tempDate);

mms.db_init('local_file_db',ParentDir);
try
    B1_ts=mms.get_data('B_gsm_brst',tempTint,ic);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
    Pos = mms.get_data('R_gsm',tint,ic);
    Pos = Pos.data;
    if isempty(Pos)
    writematrix([NameTags{TDT}(2:end-2),'无位置数据'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
    end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end

%% lobe location
flag = 0;
if Pos(1,1) <= -10*units.RE/1e3 && abs(Pos(1,3)) >= 5*units.RE/1e3 % X<-10, |Z|>5
try
    % FLAG 1, discrete location
    Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
    Ne1=irf.ts2mat(Ne1_ts);
    
    if min(Ne1(:,2)) > 0.1, continue; else, flag = 1; end

    % FLAG 2, discrete E zero
    B1_ts=mms.get_data('B_gsm_brst',tint,ic);
    B1 = irf.ts2mat(B1_ts);

    c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic); 
    c_eval(['E?_ts = irf_gse2gsm(E?_ts);'],ic);
    c_eval('E? = irf.ts2mat(E?_ts);',ic);
    c_eval(['Efac=irf_convert_fac(E?,B?,[1,0,0]);'],ic);
    Efre = 1/median(diff(Efac(:,1)));

    Efac = irf_filt(Efac, 10, 0, Efre, 3);

    Ecri = Efac(1:end-1,4) .* Efac(2:end,4); 
    index_Ezero = find(Ecri <= 0);  % E过零点

    index_Elen = index_Ezero + [-round(0.2 * Efre), round(0.2 * Efre)];
    index_Elen(index_Elen(:,1) < 1 , 1) = 1;
    index_Elen(index_Elen(:,2) > size(Efac,1) , 2) = size(Efac,1);
    
    flagTime = []; i_indexE = 1;
    while i_indexE <= size(index_Elen, 1)
        tempE = Efac(index_Elen(i_indexE, 1):index_Elen(i_indexE, 2), 4);
        if max(tempE) >= 10
            flag = 2; 
            flagTime(end+1) = Efac(index_Ezero(i_indexE), 1);
            writematrix([irf_time(flagTime(end),'epoch>utc'),'找到了ESW, Emax=', num2str(max(tempE))],...
                [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
            i_indexE = find(index_Ezero>index_Ezero(i_indexE)+Efre, 1);
        else
            i_indexE = i_indexE+1;
            continue
        end
    end
  
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
else
    continue
end
%% 符合判据的继续下载并出图
if flag == 2
try
    mms.db_init('local_file_db',ParentDir);
    for i_plot = 1:length(flagTime)
    PlotTint = irf_time([flagTime(i_plot)-0.4,flagTime(i_plot)+0.4],'epoch>epochTT');
    ESWPlot(PlotTint,ic,[NameTags{TDT}(2:end-1),num2str(i_plot)]);
    end
catch
    writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
end
end
%% 删除文件夹并生成记录文件
% % % try
% % %     cd(OutputDir)
% % %     rmdir(tempDir,'s');    
% % % catch
% % %     writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % % end
end