%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir OutputDir
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2016-01-01/2023-01-01';
splitDate = regexp(Date,'/','split');
OutputDir = [ParentDir,'CurlB_Search/',splitDate{1},'To',splitDate{2},'/'];
CaseListPath = [OutputDir, 'caselist.txt'];

CaseList = split(importdata(CaseListPath),' ');
Jcurl = str2double(CaseList(:, 2));
[sortJ, sortId] = sort(Jcurl, 'descend'); % sortJ = Jcurl(sortId);
sortT = CaseList(sortId,1);

%%
Ne = zeros(size(sortJ));
% parfor_progress(length(sortT));
for tt = 1:length(sortT)
clc
% parfor_progress;
disp(['₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*tt/length(sortT))),repmat('□',1,10-round(10*tt/length(sortT))),'正在光速检索ing...'])
fprintf(['当前处理时间为:',sortT{tt},'\n'])

    TT = sortT{tt};
    TT = datetime(TT, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS''Z');
    TT1 = datestr(TT - seconds(1), 'yyyy-mm-ddTHH:MM:SS.FFFZ');
    TT2 = datestr(TT + seconds(1), 'yyyy-mm-ddTHH:MM:SS.FFFZ');
    TT = [TT1, '/', TT2];
    tint=irf.tint(TT);

    Datelist = regexp(TT,'\d+-\d+-\d+','match');
    Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
    Date = [Datelist{1},'/',Datelist{2}];

    ic = 1;
    filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms');
    filenames = [filenames2];

    try
    [filenames,~,~] = findFilenames(TT,filenames,'fast',ic);
    SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
    SDCDataMove(TempDir,ParentDir)
    catch
    writematrix([TT,'的Ne数据缺失'],[OutputDir,'errorlog2.txt'],'WriteMode','append','Encoding','UTF-8')
    end

    try
    Ne1_ts = mms.get_data('Ne_fpi_fast_l2',tint,1);
    Ne1=irf.ts2mat(Ne1_ts);
    catch
    global ErrorFilePath_fwd_modified
    delete(ErrorFilePath_fwd_modified)
    SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
    SDCDataMove(TempDir,ParentDir)
    clear ErrorFilePath_fwd_modified
    try
    Ne1_ts = mms.get_data('Ne_fpi_fast_l2',tint,1);
    Ne1=irf.ts2mat(Ne1_ts);
    catch
    writematrix([TT,'的Ne导入出现问题'],[OutputDir,'errorlog2.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    end

    if size(Ne1,1) > 1
        Ne(tt) = mean(Ne1(:,2));
    else
        Ne(tt) = 9999;
    end

end
% parfor_progress(0);

figure(1);
Nebinranges = 0:0.1:1;
Nebincounts = histc(Ne,Nebinranges);
Nebar = bar(Nebinranges+0.05,Nebincounts,'FaceColor',"#D95319");
Nebar.Labels = Nebar.YData;
Nebar.BarWidth = 1;
xlabel('Ne [cc]');ylabel('Counts');

figure(2);
Jlobe = sortJ(Ne<0.1);
Tlobe = sortT(Ne<0.1);
Jbinranges = 0:10:150;
Jbincounts = histc(Jlobe,Jbinranges);
Jbar = bar(Jbinranges+5,Jbincounts,'FaceColor',"#D95319");
Jbar.Labels = Jbar.YData;
Jbar.BarWidth = 1;
xlabel('J [nA/m^2]');ylabel('Counts');
xlim([0,150]);xticks([0:10:150])
set(gca,"XTickLabelRotation",0)
% plot(Jlobe)

save([OutputDir, 'Ne.mat'])