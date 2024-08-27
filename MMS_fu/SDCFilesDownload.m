function SDCFilesDownload(filenames,OutputFiles_dir)
%------written by Wending Fu, May.2021 in Beijing------------
% One problem may occur!
% If a file is not fully downloaded during the last download and the download breaks, 
% the file will not be re-downloaded when the download continues
% because matlab crawler cannot read the file size
% so to interrupt, delete the recently downloaded files and re-download them at the next run
% see also SDCFilenames,SDCPlot,SDCDataMove
% ------modified by Wending Fu, Apr.2024 in Beijing------------
% A little modification for downloading the MAVEN data
% Without call the Python script, only Matlab, so it can't monitor the real-time speed! 
% May be sometime I will update it o(￣▽�?)�?

if isfolder(OutputFiles_dir) == 0
    mkdir(OutputFiles_dir);
end

% % % global ParentDir
% % % spPath = cell(1);
% % % for ii = 1:4
% % %     Maindir = [ParentDir,'mms',num2str(ii)];
% % %     tempPath =  regexp(genpath(Maindir),':','split'); % ':' for MacOs, ';' for Windows
% % %     spPath = [spPath,tempPath];
% % % end

% '/' for MacOs, '\' for Windows
% % h = waitbar(0,'正在下载�?(�? �? ≖✿)');
disp(['□□□□□□□□□□','�?始下载✧(�? �? ≖✿)'])
i = 1;
while i <= length(filenames)
    %�?查文件是否已下载�?   
    flag = 0;
    % for ii = 1:length(spPath) 
    if isfile([OutputFiles_dir,filenames{i}(1:end-3),'zip']) % || isfile([spPath{ii},'/',filenames{i}])
        %这样判别会导致重复计算，当前任务文件夹被重复判别了很多次
        s1 = ['本次任务文件夹或总文件夹中已有文�?:',num2str(i),'/',num2str(length(filenames))];
%             waitbar(i/length(filenames),h,s1);
        disp([repmat('�?',1,round(10*i/length(filenames))),repmat('�?',1,10-round(10*i/length(filenames))),s1])
        if i < length(filenames)
            i = i + 1;
            flag = 1;
        else
            flag = 2;
        end
        % break
    end
    % end
    
    switch flag
        case 2
            break
        case 1
            continue
        case 0
    if contains(filenames,'mms')
        url_file = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?', ...
            'file=',filenames{i}]; 
    elseif contains(filenames,'mvn')
        url_file = ['https://lasp.colorado.edu/maven/sdc/public/files/api/v1/search/science/fn_metadata/download_zip?', ...
            'file=',filenames{i}]; 
        filenames{i} = [filenames{i}(1:end-3),'zip'];
    end
        % % %     url_file = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?', ...
        % % % 'file=',filenames{i}];
        
    output_filename = [OutputFiles_dir,filenames{i}];   
    
    %网站接口�?挂vpn下载，如果matlab下载速度过慢则需先挂vpn再打�?matlab  
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
    
    %显示进度�?
    s2 = ['�?( ＾皿�?)っ\n','已下载文件数: ',num2str(i),'/',num2str(length(filenames)), ...
        char(13,10)','�?近文件下载�?�度:',num2str(AverageSpeed),'M/s'];
%     waitbar(i/length(filenames),h,s2);   
    disp([repmat('�?',1,round(10*i/length(filenames))),repmat('�?',1,10-round(10*i/length(filenames))),s2])
    i = i + 1;
    
    end
end

if isempty(filenames)
%     waitbar(0,h,'无可下载�?');
    disp([repmat('�?',1,10),'无可下载�?'])
else
%     waitbar(1,h,'下载完毕�?(✿ﾟ▽ﾟ)�?');
    disp([repmat('�?',1,10),'下载完毕�?(✿ﾟ▽ﾟ)�?'])
end
% pause(0.5)
% close(h)
Identification(mfilename('fullpath'));
end