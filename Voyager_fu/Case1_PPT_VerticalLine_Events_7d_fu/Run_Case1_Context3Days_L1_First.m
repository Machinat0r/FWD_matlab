function result = Run_Case1_Context3Days_L1_First(cadence)
%Run_Case1_Context3Days_L1_First Replot all requested V1 events at +/-3 days.
%   'day' overwrites 45 daily overview images; 'hour' overwrites 45 hourly
%   overviews and their five-time PADs. Complete L1 UTC rate means take
%   priority over L2. No argument runs both groups sequentially.
%   Filenames stay unchanged so previous event figures are overwritten.
if nargin < 1, cadence = 'both'; end
cadence = validatestring(cadence,{'day','hour','both'});
if strcmp(cadence,'both')
    result.Day = Run_Case1_Context3Days_L1_First('day');
    result.Hour = Run_Case1_Context3Days_L1_First('hour');
    return
end
catalog = Case1_Event_Catalog;
catalog = catalog(catalog.Spacecraft == 1,:);
assert(height(catalog)==45,'Unexpected requested V1 event catalog.');
result = struct('CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'Cadence',string(cadence),'ContextDays',3,'SourcePriority',"l1_first", ...
    'Catalog',catalog);
result.Run = Run_Case1_PPT_VerticalLine_Events_7d('RunPlots',true, ...
    'Overwrite',true,'Visible',false,'PADCadence',cadence, ...
    'EventIDs',string(catalog.EventID),'ContextDays',3, ...
    'LECPSourcePriority','l1_first','ExportPeakPAD',strcmp(cadence,'hour'));
cfg = Case1_Config;
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','context3d_l1_first');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
result.AuditFile = string(fullfile(folder,['replot_',cadence,'_',stamp,'.mat']));
result.CodeFiles = string(fullfile(cfg.CodeRoot, ...
    {'Run_Case1_Context3Days_L1_First.m','Run_Case1_PPT_VerticalLine_Events_7d.m', ...
    'Case1_Config.m','Case1_Check_Data.m','Case1_Write_Audit.m'})).';
result.CodeSHA256 = strings(size(result.CodeFiles));
for ii = 1:numel(result.CodeFiles)
    result.CodeSHA256(ii) = string(Case1_File_SHA256(char(result.CodeFiles(ii))));
end
save(result.AuditFile,'result','-v7.3');
fprintf('Three-day context, L1-first run audit: %s\n',result.AuditFile);
assert(~any(contains(string(result.Run.ReportV1.Status),'error')), ...
    'One or more production plots failed; inspect the saved run report.');
end
