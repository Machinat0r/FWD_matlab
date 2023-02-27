function filenames = SDCFilenames(Date,ic,varargin)
%------written by Wending Fu, May.2021 in Beijing------------
% Example：SDCFilenames('2017-08-23/2017-08-24',[1:4],'inst','fgm','drm','brst')
% Input variable
% Date:'2017-08-23/2017-08-23';
% ic:[1:4]; %Satellite id
% inst:'fgm'; %Instrument（optional）
% drm:'srvy'; %mode（optional）
% dpt:''; %description（optional)
% Noticing!! 
% If an instrument does NOT have this descriptor, ALL data for that instrument will not be downloaded !!!
% see also SDCDataMove,SDCPlot,SDCFilesDownload


splitDate = regexp(Date,'/','split');
% startDate = regexp(splitDate{1},'(.+?)\T','match');
% endDate = regexp(splitDate{2},'(.+?)\T','match');
sc_ids = '';
for i = ic(1:end) 
    sc_ids = ['mms',num2str(i),',',sc_ids];
end

flag = 1;inst='';drm='';dpt='';
while flag <= length(varargin)
    switch varargin{flag}
        case 'inst'
            flag = flag + 1;
            inst = ['&instrument_ids=',varargin{flag}];
        case 'drm'
            flag = flag + 1;
            drm = ['&data_rate_modes=',varargin{flag}];
        case 'dpt'
            flag = flag + 1;
            dpt = ['&descriptors=',varargin{flag}];
    end
    flag = flag + 1;
end

sc_ids=['&sc_ids=',sc_ids(1:end-1)];
url = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_names/science?start_date=', ...
    splitDate{1},'&end_date=',splitDate{2},sc_ids,inst,drm,dpt];
options = weboptions('Timeout',5); 

ErrorTimes = 0;
while ErrorTimes <= 20
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
expression = 'mms[1234]_(.+?)\.cdf'; %识别cdf文件名
filenames = regexp(sourcefile,expression,'match');
filenames = unique(filenames);
% celldisp(filenames)
Identification(mfilename('fullpath'));
end