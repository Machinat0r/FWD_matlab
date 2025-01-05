function [filenames,desmoms1,desmoms2] =  findFilenames(TT,filenames,datamode,ic)
%------modified by Wending Fu, Nov.2024 in Beijing------------
% Add a new mode to download multiple files at once.
% In this mode, the output of desmoms1 & desmoms2 is meaningless.
%------written by Wending Fu, Nov.2023 in Beijing------------
% This is the function to find the filenames which include the time
% interval, and need to be used with the SDCFilenames.m
% TT is the time interval need to use, and should be format as 'yyyy-mm-ddTHH:MM:SS.sssZ'
% filenames is the cell of whole-day filenames from the SDCFilenames
% datamode is to distinct 'srvy' or 'brst', becuse of its different naming type
%% check
if ~iscell(filenames)
    help findFilenames
    return
end
%% extract from TT
expr1 = '_\d+\_v';
NameTags = regexp(filenames,expr1,'match');
NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);

TT1 = strrep(TT, '-', '');TT2 = strrep(TT1, ':', '');
TT1 = regexp(TT1,'\d+T','match');TT2 = regexp(TT2,'\d+T\d+','match');

if datamode == 'srvy'
    Tm1 = str2double(TT1{1}(1:end-1));Tm2 = str2double(TT1{2}(1:end-1));
elseif datamode == 'brst' 
    Tm1 = str2double(strrep(TT2{1},'T',''));Tm2 = str2double(strrep(TT2{2},'T',''));
elseif datamode == 'fast' 
    Tm1 = str2double(strrep(TT2{1},'T',''));Tm2 = str2double(strrep(TT2{2},'T',''));
else
    help findFilenames
    return
end
%% 
NameTags_mat = cell2mat(NameTags);
NameTags_idx = NameTags_mat(Tm1<=NameTags_mat & NameTags_mat<=Tm2);
switch length(NameTags_idx)
    case 0, errordlg('当天无数据，请检查所选时间')
    case {1,2}, flag = length(NameTags_idx) - 1; i = find(NameTags_mat == NameTags_idx(1));
    otherwise, flag = 2; i = find(NameTags_mat == NameTags_idx);
end

%%
% '/' for MacOs, '\' for Windows
global ParentDir
if datamode == 'brst'
if flag == 0
    tempTag = num2str(NameTags{i});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag(1:4),'/',tempTag(5:6),'/',...
            tempTag(7:8),'/',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
elseif flag == 1
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NameTags{i});
    tempTag2 = num2str(NameTags{i+1});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag1(1:4),'/',tempTag1(5:6),'/',...
            tempTag1(7:8),'/',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag2(1:4),'/',tempTag2(5:6),'/',...
            tempTag2(7:8),'/',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
else
    filenames = [];
    for ii = 1:length(i)
    tempTag = num2str(NameTags{i});
    tempfilenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    filenames = [filenames, tempfilenames];
    end
    desmoms1 = nan; desmoms2 = desmoms1;
end
else
    if flag == 0
    tempTag = num2str(NameTags{i});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag(1:4),'/',tempTag(5:6),'/',...
            '/',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
    elseif flag == 1
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NameTags{i});
    tempTag2 = num2str(NameTags{i+1});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag1(1:4),'/',tempTag1(5:6),'/',...
            '/',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'/fpi/',datamode,'/l2/des-moms/',tempTag2(1:4),'/',tempTag2(5:6),'/',...
            '/',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
    else
    filenames = [];
    for ii = 1:length(i)
    tempTag = num2str(NameTags{i});
    tempfilenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    filenames = [filenames, tempfilenames];
    end 
    desmoms1 = nan; desmoms2 = desmoms1;
    end
end
end