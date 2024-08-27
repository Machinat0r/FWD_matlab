function filenames = SDCFilenames_modified(Date,ic,varargin)
%------written by Wending Fu, May.2021 in Beijing------------
% Example：SDCFilenames('2017-08-23/2017-08-24',[1:4],'inst','fgm','drm','brst')
% Input variable
% Date:'2017-08-23/2017-08-23';
% ic:[1:4]; %Satellite id
% inst:'fgm'; %Instrument（optional�?
% drm:'srvy'; %mode（optional�?
% dpt:''; %description（optional)
% Noticing!! 
% If an instrument does NOT have this descriptor, ALL data for that instrument will not be downloaded !!!
% see also SDCDataMove,SDCPlot,SDCFilesDownload
%------modified by Wending Fu, Apr.2024 in Beijing------------
% Modified for MAVEN download

splitDate = regexp(Date,'/','split');
% startDate = regexp(splitDate{1},'(.+?)\T','match');
% endDate = regexp(splitDate{2},'(.+?)\T','match');
if lower(ic) == "maven"
    sc_ids = 'maven';
else
sc_ids = '';
for i = ic(1:end) 
    sc_ids = ['mms',num2str(i),',',sc_ids];
end
end

flag = 1;inst='';drm='';dpt='';
while flag <= length(varargin)
    switch varargin{flag}
        case 'inst'
            flag = flag + 1;
            if sc_ids == 'maven'
            inst = ['&instrument=',varargin{flag}];
            else
            inst = ['&instrument_ids=',varargin{flag}];
            end
        case 'drm'
            flag = flag + 1;
            drm = ['&data_rate_modes=',varargin{flag}];
        case 'dpt'
            flag = flag + 1;
            dpt = ['&descriptors=',varargin{flag}];
        case 'level'
            flag = flag + 1;
            level = ['level=',varargin{flag}];
    end
    flag = flag + 1;
end
%% 
if sc_ids == "maven"
url = ['https://lasp.colorado.edu/maven/sdc/public/files/api/v1/search/science/fn_metadata/file_info?',...
    level,inst,'&start_date=',splitDate{1},'&end_date=',splitDate{2}];
else
sc_ids=['&sc_ids=',sc_ids(1:end-1)];
url = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_names/science?start_date=', ...
    splitDate{1},'&end_date=',splitDate{2},sc_ids,inst,drm,dpt];
end

options = weboptions('Timeout',5); 

ErrorTimes = 0;
while ErrorTimes <= 5
try
    sourcefile = webread(url,options);
    fprintf('文件目录爬取已完成\n')
    ErrorTimes = 666;
catch
    fprintf('连接超时，请挂vpn\n')
    ErrorTimes = ErrorTimes + 1;
end
end

% expression = '<a.+?href=\"(.+?)\">(.+?)</a>'; %识别网页中的链接
if sc_ids == "maven"
    sourcefile = sourcefile.files;
    expression = 'sts';
    filenames = {};
    for i_files = 1:length(sourcefile)
    if contains(sourcefile(i_files).file_name,expression)
        filenames{end+1} = sourcefile(i_files).file_name;
    end
    end
else
expression = 'mms[1234]_(.+?)\.cdf'; %识别cdf文件�?
filenames = regexp(sourcefile,expression,'match');
filenames = unique(filenames);
end

fprintf(['�?',num2str(length(filenames)),'个文件\n'])

% celldisp(filenames)
Identification(mfilename('fullpath'));
end