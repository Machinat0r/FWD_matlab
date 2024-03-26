%------written by Wending Fu, Dec.2023 in Beijing------------
clear; clc
close all
%% parent dir
parent_folder = '/Volumes/172.17.190.41/Data/MMS/DFSearch';
combination_folder = [parent_folder, '/combination/'];
combination_caselist = [combination_folder, 'caselist.txt'];
combination_overview = [combination_folder,'overview/'];
if ~isfolder(combination_folder), mkdir(combination_folder); end
if ~isfolder(combination_overview), mkdir(combination_overview); end
dir_names = ls(parent_folder);
dir_names = strsplit(dir_names);
%% sub dir
for i = 1:length(dir_names)
if ~isempty(dir_names{i}) && ~contains(dir_names{i},'combination')
dir = dir_names{i};
sub_folder = [parent_folder, '/', dir];
sub_dir_names = ls(sub_folder);
%% combine the case list
if contains(sub_dir_names,'caselist.txt')
    filename_path = [sub_folder, '/', 'caselist.txt'];
    DF_list_temp = importdata(filename_path);
    DF_time_list = cellfun(@(x)([char(regexp(x,'\d+-\d+-\d+T\d+:\d+:\d+.\d{3}','match')), 'Z'])...
        , DF_list_temp,'UniformOutput',false);
    writecell(DF_time_list, combination_caselist, 'WriteMode', 'append', 'Encoding', 'UTF-8')
end

%% combine the overview figures
if contains(sub_dir_names, 'OverviewFig')
    overview_path = [sub_folder, '/', 'OverviewFig'];
    copyfile(overview_path,combination_overview)
end
end
end