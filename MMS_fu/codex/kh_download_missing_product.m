function downloaded = kh_download_missing_product(plotStart,plotEnd,sc,instrument,descriptor)
%KH_DOWNLOAD_MISSING_PRODUCT Download a missing public MMS burst product.
%
% downloaded = kh_download_missing_product(t0,t1,3,'fpi','dis-moms')
%
% Uses the existing SDCFilenames, SDCFilesDownload and SDCDataMove tools in
% MMS_fu.  Only files whose start stamps lie near the requested plot window
% are fetched, and the latest version for each file start stamp is retained.

validateattributes(sc,{'numeric'},{'scalar','integer','>=',1,'<=',4});
if ~isdatetime(plotStart) || ~isdatetime(plotEnd)
    error('plotStart and plotEnd must be datetime values.');
end
if isempty(plotStart.TimeZone), plotStart.TimeZone = 'UTC'; end
if isempty(plotEnd.TimeZone), plotEnd.TimeZone = 'UTC'; end

dataRoot = 'Z:\SPART-WORK\Data\MMS';
codeDir = fileparts(mfilename('fullpath'));
tempDir = fullfile(codeDir,'download_temp');
if ~exist(tempDir,'dir'), mkdir(tempDir); end
tempDirWithSep = [tempDir filesep];

day0 = dateshift(plotStart,'start','day');
day1 = dateshift(plotEnd,'start','day') + days(1);
dateRange = sprintf('%s/%s',char(string(day0,'yyyy-MM-dd')),char(string(day1,'yyyy-MM-dd')));

allNames = SDCFilenames(dateRange,sc,'inst',char(instrument), ...
    'drm','brst','dpt',char(descriptor));
keep = false(size(allNames));
fileTimes = NaT(size(allNames),'TimeZone','UTC');
for k = 1:numel(allNames)
    token = regexp(allNames{k},'_(\d{14})_v','tokens','once');
    if isempty(token), continue; end
    fileTimes(k) = datetime(token{1},'InputFormat','yyyyMMddHHmmss','TimeZone','UTC');
    keep(k) = fileTimes(k) >= plotStart-minutes(10) && fileTimes(k) <= plotEnd;
end
names = allNames(keep);
fileTimes = fileTimes(keep);

if isempty(names)
    downloaded = table(strings(0,1),NaT(0,1,'TimeZone','UTC'),strings(0,1), ...
        'VariableNames',{'Filename','FileStartUTC','Status'});
    warning('No public %s/%s files were listed near %s to %s.', ...
        instrument,descriptor,char(plotStart),char(plotEnd));
    return
end

% Keep the lexicographically latest version for each timestamp/product.
[names,order] = sort(names);
fileTimes = fileTimes(order);
bases = regexprep(names,'_v[^_]+\.cdf$','');
[~,latest] = unique(bases,'last');
latest = sort(latest);
names = names(latest);
fileTimes = fileTimes(latest);

fprintf('KH_DOWNLOAD_START sc=%d %s/%s files=%d\n',sc,instrument,descriptor,numel(names));
SDCFilesDownload(names,tempDirWithSep);
SDCDataMove(tempDirWithSep,[dataRoot filesep]);
mms.db_init('local_file_db',[dataRoot filesep]);

downloaded = table(string(names(:)),fileTimes(:),repmat("downloaded",numel(names),1), ...
    'VariableNames',{'Filename','FileStartUTC','Status'});
fprintf('KH_DOWNLOAD_DONE files=%d\n',height(downloaded));
end
