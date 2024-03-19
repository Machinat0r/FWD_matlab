%------written by Wending Fu, Dec.2023 in Beijing------------
clear; clc
close all
%% parent dir
parent_folder = '/Volumes/172.17.190.41/Data/MMS/DFSearch/';
global OutputDir ParentDir
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
TempDir = [ParentDir,'temp/'];
if ~isfolder(TempDir);mkdir(TempDir);end
OutputDir = [parent_folder, 'TailwardFlowPlot/'];
if ~isfolder(OutputDir);mkdir(OutputDir);end

caselist = importdata([parent_folder, 'combination/caselist.txt']);
caselist = sort(unique(caselist));

for caseid = 1:length(caselist)
try
flagTime = caselist{caseid};
TT = irf_time(flagTime, 'utc>epoch');
TT = irf_time(TT + 30 * [-1, 1], 'epoch>utc');
TT = [TT(1,1:23),'Z/',TT(2,1:23),'Z'];

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames_moms = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
filenames_all = [filenames1,filenames2,filenames3,filenames4];


[filenames_moms,~,~] = findFilenames(TT,filenames_moms,'brst',ic);
% [filenames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

SDCFilesDownload_NAS(filenames_moms,TempDir)
% SDCFilesDownload_NAS(filenames_fast,TempDir)
% SDCFilesDownload_NAS(filenames_srvy,TempDir)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
SDCDataMove(TempDir,ParentDir)
%%
c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],ic);
if gsmVi1(100,2) <= 0 || gsmVi1(200,2) <= 0 || gsmVi1(300,2) <= 0
[filenames_all,desmoms1,desmoms2] = findFilenames(TT,filenames_all,'brst',ic);
SDCFilesDownload_NAS(filenames_all,TempDir)
SDCDataMove(TempDir,ParentDir)
id_flagTime = SDCPlot(tint,desmoms1,desmoms2,ic,flagTime,irf_time(flagTime, 'epoch>utc'));
end
catch
writematrix([flagTime,'出现问题'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
end
end