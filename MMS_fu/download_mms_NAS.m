close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

TT = '2024-01-01T00:00:00.000Z/2025-01-01T00:00:00.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
iic = 1:4;
try
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,iic,'inst','edp','drm','brst','dpt','dce');
filenames = [filenames1, filenames2, filenames3, filenames4];

filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 

[filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);
[filenames,~,~] = findFilenames(TT,filenames,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'CheckSize', 0, 'Threads', 16)
SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 64, 'CheckSize', 0)
catch
    warning('no files have been downloaded')
end
%% load data
SDCDataMove(TempDir,ParentDir)