%------written by Wending Fu, Apr.2025 in Beijing------------
clear;clc;
global ParentDir OutputDir
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
Date = '2015-01-01/2019-01-01';
splitDate = regexp(Date,'/','split');
OutputDir = [ParentDir,'CurlB_Search/',splitDate{1},'To',splitDate{2},'/'];
CaseListPath = [OutputDir, 'errorlog.txt'];
ErrorList = load(CaseListPath);

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms');
filenames = [filenames1, filenames2, filenames3, filenames4];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);