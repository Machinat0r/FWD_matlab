%main
%gsm
clear;
clc;
%%
Date = '2015-09-01/2023-04-30';
% Date = '2022-07-21/2023-04-30';
% Date = '2017-04-30/2017-05-01';
% Date = '2017-01-01/2021-01-01';
% Date = '2022-06-09/2022-07-02';

splitDate = regexp(Date,'/','split');
ic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames3 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %To get loaction
filenames = [filenames1, filenames_srvy];

% expr = '_[0-9]+\_v';
expr = '[0-9]{8}';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

global ParentDir
ParentDir = 'Z:\Data\MMS\'; 

%%
NameTags{end+1} = ['_' strrep(splitDate{2},'-','') '235959_v'];
for TDT = 1:length(FileGroups)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [ParentDir,'temp\'];
clc
fprintf(['当前处理时间为:',NameTags{TDT},'\n'])


    SDCFilesDownload_NAS(FileGroups{TDT},tempDir);
    SDCDataMove(tempDir,ParentDir); mms.db_init('local_file_db',ParentDir);
end