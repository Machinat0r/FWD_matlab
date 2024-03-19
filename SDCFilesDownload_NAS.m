function SDCFilesDownload_NAS(filenames,OutputFiles_dir, varargin)
%------modified by Wending Fu, Jan.2024 in Beijing------------
% now it can use multi threads to download!
% NOTICE: multi threads may cause the progroam crash! 
% And it can NOT be STOPPED!
%------modified by Wending Fu, Jan.2024 in Beijing------------
% Solve the problem in first version
% This function can check the file size now!
% This check will slow the download speed, if you in a harry, please set
% the 'CheckSize' to 0, default is 1.
%------modified by Wending Fu, Aug.2023 in Beijing------------
% To adopt the low speed of SPART-NAS folder checking
% Remove the FULL-FOLDER-CHECK
% This version ONLY work after SDCDatamove running
%------written by Wending Fu, May.2021 in Beijing------------
% One problem may occur!
% If a file is not fully downloaded during the last download and the download breaks, 
% the file will not be re-downloaded when the download continues
% because matlab crawler cannot read the file size
% so to interrupt, delete the recently downloaded files and re-download them at the next run
% see also SDCFilenames,SDCPlot,SDCDataMove
%% input parser
parser = inputParser;
defaultCheckSize = 1;
defaultThreads = 0;
addParameter(parser, 'CheckSize', defaultCheckSize, @isnumeric);
addParameter(parser, 'Threads', defaultThreads, @isnumeric);
parse(parser, varargin{:});
CheckSize = parser.Results.CheckSize;
Threads = parser.Results.Threads;
%%
if isfolder(OutputFiles_dir) == 0
    mkdir(OutputFiles_dir);
end

global ParentDir
% '/' for MacOs, '\' for Windows
disp(['□□□□□□□□□□','开始下载✧(≖ ◡ ≖✿)'])
i = 1;
while i <= length(filenames)
%% 
    flag = 0;
    if isDataInDisk(filenames{i}, ParentDir) || isDataInDisk(filenames{i}, OutputFiles_dir)
        if CheckSize
            if isDataInDisk(filenames{i}, ParentDir), [~, filepath] = isDataInDisk(filenames{i}, ParentDir);
            else, [~, filepath] = isDataInDisk(filenames{i}, OutputFiles_dir); end
        [~, ContentLength] = FileContentLength(filenames{i});
        if dir(filepath).bytes == ContentLength
        s1 = ['本次任务文件夹或总文件夹中已有文件:',num2str(i),'/',num2str(length(filenames))];
%             waitbar(i/length(filenames),h,s1);
        disp([repmat('■',1,round(10*i/length(filenames))),repmat('□',1,10-round(10*i/length(filenames))),s1])
            if i < length(filenames)
                i = i + 1;
                flag = 1;
            else
                flag = 2;
            end
        else
            warning('The old file is imcompleted! It will be re-download now!')
            delete(filepath)
        end
        else
            s1 = ['本次任务文件夹或总文件夹中已有文件:',num2str(i),'/',num2str(length(filenames))];
%             waitbar(i/length(filenames),h,s1);
        disp([repmat('■',1,round(10*i/length(filenames))),repmat('□',1,10-round(10*i/length(filenames))),s1])
            if i < length(filenames)
                i = i + 1;
                flag = 1;
            else
                flag = 2;
            end
        end
    end
%%
    switch flag
        case 2
            break
        case 1
            continue
        case 0
    [file_url, ContentLength] = FileContentLength(filenames{i});
    output_filename = [OutputFiles_dir,filenames{i}];   
    
    %网站接口需挂vpn下载，如果matlab下载速度过慢则需先挂vpn再打开matlab  
    options = weboptions('Timeout',10);        
    ErrorTimes = 0;
    
    while ErrorTimes <= 20
        tic
        if Threads == 0
            websave(output_filename,file_url,options) ; 
        else
            command = sprintf('python3 /Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/download_files.py "%s" "%s" "%s" "%s"',...
                file_url, output_filename, filenames{i}, num2str(Threads)); % 若电脑发生卡顿，可以下调线程数（最后一个变量）
            system(command);
        end   
        TimeInterval = toc;   

        Dir = dir(OutputFiles_dir); 
        FileIndex = strcmp({Dir.name},filenames{i});
        SizeOfFile = Dir(FileIndex).bytes;
        AverageSpeed = SizeOfFile/(TimeInterval*1024*1024);
    
        if SizeOfFile == ContentLength
        %显示进度条
        s2 = ['已下载文件数: ',num2str(i),'/',num2str(length(filenames)), ...
            char(13,10)','最近文件下载速度:',num2str(AverageSpeed),'M/s'];  
        disp([repmat('■',1,round(10*i/length(filenames))),repmat('□',1,10-round(10*i/length(filenames))),s2])
        i = i + 1;
        break
        else
            fprintf('连接超时，请检查网络配置\n')
            ErrorTimes = ErrorTimes + 1;
        end
    end
    end
    if exist('ContentLength', 'var'), clear ContentLength; end
end

if isempty(filenames)
    disp([repmat('□',1,10),'无可下载项'])
else
    disp([repmat('■',1,10),'下载完毕ヽ(✿ﾟ▽ﾟ)ノ'])
end
Identification(mfilename('fullpath'));
end
%%
function [file_url, varargout] = FileContentLength(filename)
file_url = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?', ...
    'file=',filename]; 

try
info = py.requests.get(file_url, stream=true);
ContentLength = str2double(char(info.headers{'Content-Length'}));
varargout{1} = ContentLength;
catch
    warning('Get content length failed! Please check the Python version.')
    warning('This may cause the file download imcompleted!')
end
end