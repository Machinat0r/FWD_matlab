function ToF = isDataInDisk(filename,destDir)
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
        ToF = isfile([tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            ,'/', filename_sep{ii}(7:8), '/', filename]);
    elseif ismember('srvy',filename_sep) || ismember('fast', filename_sep)
        ToF = isfile([tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            , '/', filename]);
    else
        ToF = isfile([tempDir, filename_sep{ii}(1:4),'/', filename_sep{ii}(5:6)...
            , '/', filename]);
        warning('Unrecongnized File Type (′?н?‘)')
    end
end

if ~ToF
    if isfolder([destDir, 'temp'])
        ToF = isfile([destDir, 'temp','/', filename]);
    end
end
fclose all;
Identification(mfilename('fullpath'));
end