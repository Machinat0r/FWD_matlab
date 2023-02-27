%main
%gsm
clear;
clc;
Date = '2019-08-16/2019-08-17';
splitDate = regexp(Date,'/','split');
ic = 1;
% filenames = SDCFilenames(Date,1,'inst','fgm,fpi,scm,edp','drm','brst');
filenames1 = SDCFilenames(Date,ic,'inst','fgm,edp','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3];

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags));
for j = 1:length(NameTags)
    temp = 0;
    for  i = 1:length(filenames)
        if strfind(filenames{i},NameTags{j}) > 0
            temp = temp + 1;
            FileGroups{j}{1,temp} = filenames{i};
        end
    end
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);

global OutputDir ParentDir
ParentDir = 'D:\MMS\'; 
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,splitDate{1},'To',splitDate{2},'\'];
if ~isfolder([OutputDir,'OverviewFig\'])
    mkdir([OutputDir,'OverviewFig\']);
end

for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'\'];
fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
for i = 1:length(FileGroups{TDT})
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'des-moms'))
        SDCFilesDownload(FileGroups{TDT}(i),tempDir)
    end
end
for j = 1:length(filenames_srvy)
    if ~isempty(strfind(filenames_srvy{j},NameTags{TDT}(2:9)))
        SDCFilesDownload(filenames_srvy(j),tempDir)
    end
end
end
SDCDataMove(tempDir,ParentDir)