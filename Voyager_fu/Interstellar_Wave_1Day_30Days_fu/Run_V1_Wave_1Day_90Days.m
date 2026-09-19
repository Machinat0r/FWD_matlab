function result=Run_V1_Wave_1Day_90Days(varargin)
% Daily Morlet figure, 2--90 days, from 2012-08-25 to latest valid public MAG.
% Direct original CDF input; inherited linear gap interpolation and [0,3] colors.
% UseLatestAvailable=false reproduces the former end date 2021-12-16.
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_90Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',true,@islogical);
addParameter(p,'UseLatestAvailable',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
stopUTC=datetime(2021,12,17,'TimeZone','UTC');
availability=struct;
verifiedMissing=strings(0,1);
if cfg.UseLatestAvailable
    codeRoot=fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
    Case1_Add_IRFU_Path(cfg.IRFURoot);
    availability=V1_Sync_Morlet_Archive(cfg.DataRoot,datetime(2012,8,25,'TimeZone','UTC'));
    stopUTC=availability.StopExclusiveUTC;
    verifiedMissing=availability.VerifiedMissingMonths;
end
result=Run_V1_Wave_1Day_30Days('DataRoot',cfg.DataRoot, ...
    'OutputRoot',cfg.OutputRoot,'IRFURoot',cfg.IRFURoot,'Visible',cfg.Visible, ...
    'MaxPeriodDays',90,'StopUTC',stopUTC, ...
    'VerifiedMissingMonths',verifiedMissing,'AvailabilityAudit',availability);
end

