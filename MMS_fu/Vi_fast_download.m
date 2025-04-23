%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2019-06-01/2019-06-02';
splitDate = regexp(Date,'/','split');
ic = 1;iic = 1;
% filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
% filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames_fast = SDCFilenames(Date,iic,'inst','fpi','drm','fast','dpt','dis-moms');
% filenames = [filenames_fast];
% [fileames_fast,~,~] = findFilenames(Date,filenames,'fast',ic);

SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 48, 'CheckSize', 0)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);