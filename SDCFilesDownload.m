function SDCFilesDownload(filenames,OutputFiles_dir)
%------written by Wending Fu, May.2021 in Beijing------------
% One problem may occur!
% If a file is not fully downloaded during the last download and the download breaks, 
% the file will not be re-downloaded when the download continues
% because matlab crawler cannot read the file size
% so to interrupt, delete the recently downloaded files and re-download them at the next run
% see also SDCFilenames,SDCPlot,SDCDataMove

if isfolder(OutputFiles_dir) == 0
    mkdir(OutputFiles_dir);
end

global ParentDir
spPath = cell(1);
for ii = 1:4
    Maindir = [ParentDir,'mms',num2str(ii)];
    tempPath =  regexp(genpath(Maindir),';','split');
    spPath = [spPath,tempPath];
end

% % h = waitbar(0,'正在下载✧(≖ ◡ ≖✿)');
disp(['□□□□□□□□□□','开始下载✧(≖ ◡ ≖✿)'])
i = 1;
while i <= length(filenames)
    %检查文件是否已下载过   
    flag = 0;
    for ii = 1:length(spPath) 
        if isfile([OutputFiles_dir,filenames{i}]) || isfile([spPath{ii},'\',filenames{i}])
            %这样判别会导致重复计算，当前任务文件夹被重复判别了很多次
            s1 = ['本次任务文件夹或总文件夹中已有文件:',num2str(i),'/',num2str(length(filenames))];
%             waitbar(i/length(filenames),h,s1);
            disp([repmat('■',1,round(10*i/length(filenames))),repmat('□',1,10-round(10*i/length(filenames))),s1])
            if i < length(filenames)
                i = i + 1;
                flag = 1;
            else
                flag = 2;
            end
            break
        end
    end
    
    switch flag
        case 2
            break
        case 1
            continue
        case 0
    url_file = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?', ...
        'file=',filenames{i}];
    output_filename = [OutputFiles_dir,filenames{i}];   
    
    %网站接口需挂vpn下载，如果matlab下载速度过慢则需先挂vpn再打开matlab  
    options = weboptions('Timeout',10);        
    ErrorTimes = 0;
    while ErrorTimes <= 20
    try
        tic
        websave(output_filename,url_file,options) ;    
        TimeInterval = toc;   
        ErrorTimes = 666;
    catch
        fprintf('连接超时，请挂vpn\n')
        ErrorTimes = ErrorTimes + 1;
    end
    end
    
  
    Dir = dir(OutputFiles_dir); 
    FileIndex = strcmp({Dir.name},filenames{i});
    SizeOfFile = Dir(FileIndex).bytes;
    AverageSpeed = SizeOfFile/(TimeInterval*1024*1024);
    
    %显示进度条
    s2 = ['已下载文件数: ',num2str(i),'/',num2str(length(filenames)), ...
        char(13,10)','最近文件下载速度:',num2str(AverageSpeed),'M/s'];
%     waitbar(i/length(filenames),h,s2);   
    disp([repmat('■',1,round(10*i/length(filenames))),repmat('□',1,10-round(10*i/length(filenames))),s2])
    i = i + 1;
    
    end
end

if isempty(filenames)
%     waitbar(0,h,'无可下载项');
    disp([repmat('□',1,10),'无可下载项'])
else
%     waitbar(1,h,'下载完毕ヽ(✿ﾟ▽ﾟ)ノ');
    disp([repmat('■',1,10),'下载完毕ヽ(✿ﾟ▽ﾟ)ノ'])
end
% pause(0.5)
% close(h)
Identification(mfilename('fullpath'));
end