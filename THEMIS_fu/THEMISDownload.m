function THEMISDownload(Date,Ids,Insts,Output_dir)
%%
% Created by fwd on 2021.12.27

% 已知可能会出现的问题：如果某文件在上一次下载时未下载完全而下载断开，在继续下载时不会重新下载该文件
% 原因是matlab自带的爬虫无法读取文件大小
% 因此若要中断请删除最近下载的文件，在下次运行时会重新下载
% 但若上一次不是ctrl C中断而是关掉进度框之后其报错自动退出，则文件下载应该是完整的

% 对于THEMIS的这个下载程序，由于没有文件名爬取过程，可能会出现文件名版本更新不匹配的情况
% example: THEMISDownload('20090301','tha','fgm','D:\THEMIS\temp\')
%         THEMISDownload('20090305',3,'ssc','D:\THEMIS\temp\')
%%
if ~isfolder(Output_dir)
    mkdir(Output_dir);
end
% transform date
Year = Date(1:4);
Mon = Date(5:6);
% transform satellites Ids
if isnumeric(Ids)
    Id = 'abcde';
    Ids = ['th',Id(Ids)];  
end

% Define Parameters
if ismember(Insts,{'fgm';'esa';'efi';'mom';'scm';'fbk'})
    filename = [Ids '_l2_' Insts '_' Date '_v01.cdf'];
    url = ['http://themis.ssl.berkeley.edu/data/themis/' Ids '/l2/' Insts '/' Year '/'];
elseif ismember(Insts,{'ssc'})
    filename = [Ids '_or_ssc_' Year Mon '01_v01.cdf'];
    url = ['https://cdaweb.gsfc.nasa.gov/pub/data/themis/' Ids '/ssc/' Year '/'];
else 
    disp('Input instrument has not been added')
end
url_file = [url filename];
output_filename = [Output_dir filename];
options = weboptions('Timeout',20);        
ErrorTimes = 0;
flag = 0;


global ParentDir
spPath = cell(1);
Maindir = [ParentDir,Ids,'\'];
tempPath =  regexp(genpath(Maindir),';','split');
spPath = [spPath,tempPath];

disp(['Begin Download (●ˇ∀ˇ●)'])


% Check file exist
for ii = 2:length(spPath) 
    if isfile([Output_dir,filename]) || isfile([spPath{ii},'\',filename])
        flag = 1;
        break
    end  
end

if flag == 0
%download
while ErrorTimes <= 20
    try
        tic
        websave(output_filename,url_file,options) ;    
        TimeInterval = toc;   
        ErrorTimes = 999;
    catch
        fprintf('Connection timed out. Please hang up vpn\n')
        ErrorTimes = ErrorTimes + 1;
    end
end

% calculate speed
if ErrorTimes == 21
    disp('Download failed (+_+)? Please check the input format and network connection')
else
    Dir = dir(Output_dir); 
    FileIndex = strcmp({Dir.name},filename);
    SizeOfFile = Dir(FileIndex).bytes;
    AverageSpeed = SizeOfFile/(TimeInterval*1024*1024);
    disp([ 'Download success (❤ ω ❤): ',filename...
        newline,'Latest file download speed:',num2str(AverageSpeed),'M/s'])
end
else
    disp(['The files are already in the temp folder or general folder:' filename])
end
end