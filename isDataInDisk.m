function [ToF, varargout] = isDataInDisk(filename,destDir)
%------modified by Wending Fu, Jan.2024 in Beijing------------
% this function can discriminate temporary folder now
%------modified by Wending Fu, Jan.2024 in Beijing------------
% this function can return the file path now
%------written by Wending Fu, Aug.2023 in Beijing------------
% script description:
% destDir: the path of the directory you restore the datas
%% filename
if isstring(filename)
    convertStringsToChars(filename);
elseif  ~ischar(filename)
    help isDataInDisk
    return
end
%% temp folder distinguis
temp_folder_flag = 1;
dir_list = regexp(ls(destDir), '\s+', 'split');
for file_i = 1:length(dir_list)
if ~isempty(dir_list{file_i})
if isfolder([destDir, dir_list{file_i}])
    temp_folder_flag = 0;break
end
end
end

if ~temp_folder_flag
%% folder path
filename_sep = regexp(filename,'_','split');

ii=1; tempDir = destDir;
while isempty(regexp(filename_sep{ii}(1),'[0-9]', 'once'))
    tempDir = [tempDir, filename_sep{ii}, '/'];
    ii = ii + 1;
end
%% in or not in, this is a question ╮(๑•́ ₃•̀๑)╭
ToF = false;
if isfolder(tempDir)
    if ismember('brst',filename_sep)
        filepath = [tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            ,'/', filename_sep{ii}(7:8), '/', filename];
        ToF = isfile(filepath);
    elseif ismember('srvy',filename_sep) || ismember('fast', filename_sep)
        filepath = [tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            , '/', filename];
        ToF = isfile(filepath);
    else
        filepath = [tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            , '/', filename];
        ToF = isfile(filepath);
        warning('Unrecongnized File Type (′?н?‘)')
    end
end

if ~ToF
    if isfolder([destDir, 'temp'])
        filepath = [destDir, 'temp','/', filename];
        ToF = isfile(filepath);
    end
end
else
    filepath = [destDir, filename];
    ToF = isfile(filepath);
end

fclose all;
varargout{1} = filepath;
Identification(mfilename('fullpath'));
end