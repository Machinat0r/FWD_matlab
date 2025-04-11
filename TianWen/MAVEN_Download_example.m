%------written by Wending Fu, Apr.2024 in Beijing------------
% Example：SDCFilenames('2017-08-23/2017-08-24',[1:4],'inst','fgm','drm','brst')
% Input variable
% Date:'2017-08-23/2017-08-23';
% ic:[1:4]; %Satellite id
% inst:'mag'; %Instrument（optional）
% drm:'srvy'; %mode（optional）
% dpt:''; %description（optional)
% level:'l2'; %description（optional)
% -----------------------------------------------------------
% Notice!!!!
% If a file haen't been fully download and this script is forced to close
% for some reason, Please delete the recently file and let it download again (⊙﹏⊙)
% Some time i may update this script ヾ(≧▽≦*)o
% -----------------------------------------------------------
% see also SDCFilenames_modifed, SDCFilesDownload
%%
close all
clear;clc

DownloadDir = '/Users/fwd/Documents/MATLAB/MAVEN/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

TT = '2020-08-02T00:00:00.000Z/2020-08-03T00:00:00.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames = SDCFilenames_modified(Date,'MAVEN','inst','mag','level','l2');
% [filenames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
expression = 'pc1s_';
filenames = filenames(contains(filenames,expression));
fprintf([expression(1:end-1),'文件共',num2str(length(filenames)),'个(๑•̀ㅂ•́)و✧\n'])
SDCFilesDownload(filenames,TempDir)