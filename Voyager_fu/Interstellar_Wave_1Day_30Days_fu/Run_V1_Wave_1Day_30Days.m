function result = Run_V1_Wave_1Day_30Days(varargin)
% Voyager 1 daily |B| Morlet wavelet, following the supplied MMS_fu/Wave.m.
% User clarified on 2026-09-10: use the daily magnetic field in the figure.
% Daily sampling supports periods >=2 days. MaxPeriodDays controls the lower frequency.
% Read README_波动分析.md for source definitions and native IRFU edge treatment.

%% Paths and parameters
codeRoot=fileparts(mfilename('fullpath'));
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'MaxPeriodDays',30,@(x)isscalar(x)&&isfinite(x)&&x>2);
addParameter(p,'Visible',true,@islogical);
addParameter(p,'StopUTC',datetime(2021,12,17,'TimeZone','UTC'),@(x)isdatetime(x)&&isscalar(x)&&~isnat(x));
addParameter(p,'VerifiedMissingMonths',strings(0,1),@(x)isstring(x)||iscellstr(x));
addParameter(p,'AvailabilityAudit',struct,@isstruct);
parse(p,varargin{:}); cfg=p.Results;
cfg.StartUTC=datetime(2012,8,25,'TimeZone','UTC');
cfg.StopUTC.TimeZone='UTC';
assert(cfg.StopUTC>cfg.StartUTC && cfg.StopUTC==dateshift(cfg.StopUTC,'start','day'),'StopUTC must be an exclusive UTC midnight after the start.');
cfg.NumFrequencies=100;
cfg.WaveletWidth=5.36;
cfg.RequestedFrequencyHz=[1/(cfg.MaxPeriodDays*86400),1/86400];
cfg.ComputableFrequencyHz=[1/(cfg.MaxPeriodDays*86400),1/(2*86400)];
cfg.DisplayPeriodDays=[2 cfg.MaxPeriodDays];
cfg.PSDColorLimits=[0 3]; % User requested color range 0 to 3.
if ~isfolder(cfg.OutputRoot), mkdir(cfg.OutputRoot); end
diary(fullfile(cfg.OutputRoot,'run_daily_morlet.log'));
logCleanup=onCleanup(@()diary('off'));
addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path(cfg.IRFURoot);
fprintf('Daily Morlet run UTC: %s\n',string(datetime('now','TimeZone','UTC')));

%% Read original CDF and reproduce the figure's scalar daily mean
[raw,sources,duplicates]=V1_Wave_Read_Hourly(cfg);
day=(cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin]=ismember(dateshift(raw.EpochUTC,'start','day'),day);
good=matched & isfinite(raw.Bmag_nT);
dailyMean=accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count=accumarray(bin(good),1,[numel(day) 1],@sum,0);
daily=table(day+hours(12),dailyMean,count, ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
assert(sum(count)==nnz(good),'Input sample accounting failed.');
fprintf('%d UTC days: %d measured daily means, %d missing days.\n', ...
    height(daily),nnz(isfinite(daily.BMean_nT)),nnz(~isfinite(daily.BMean_nT)));

%% User-authorized linear interpolation of the missing daily means
observedDaily=daily;
[daily,interpolationAudit]=V1_Interpolate_Daily_Gaps(observedDaily);
fprintf('Linearly interpolated %d internal daily gaps. Observed values unchanged.\n',height(interpolationAudit));

%% Morlet transform: same IRFU function as Wave.m
[spec,segments]=V1_Daily_Morlet(daily,cfg);
spec.InterpolatedInputMask=daily.IsInterpolated;
assert(height(segments)==1,'Filled input must be analyzed as one continuous series.');

%% Reproducibility and figures
method=struct;
method.UserInstruction='User requested the supplied daily Morlet figure from 2012-08-25 through the last valid publicly available Voyager 1 magnetic day. Retain the established scalar daily mean, linear interior-gap interpolation, 2--90 day display and [0,3] colors for the 90-day entry.';
method.Source='VOYAGER1_COHO1HR_MERGED_MAG_PLASMA original CDF; Epoch and scalar ABS_B in nT.';
method.Daily='UTC [00:00,next 00:00), arithmetic mean of source-valid hourly ABS_B; point at 12:00 UTC. Identical to the supplied overview panel a definition. No count/coverage threshold.';
method.Quantity='Scalar daily mean |B| only. Its power is not the sum of BR/BT/BN wavelet powers.';
method.Wavelet=sprintf('Direct irf_wavelet; 100 logarithmic frequencies from 1/(%g days) to daily Nyquist 1/(2 days); Morlet width 5.36, returnpower=1, cutedge=1. No fixed-duration sliding window.',cfg.MaxPeriodDays);
method.Normalization='Unmodified IRFU convention: P=(2*pi)*abs(W)^2/f; nT^2/Hz. No new normalization, significance, noise subtraction, detrending, prefilter, FAC transform or temporal averaging of power.';
method.Gaps=sprintf('Fill only missing scalar daily means by explicit irf_resamp(...,''linear'') between nearest observed daily values. Observed values and original sample counts remain unchanged. Endpoint assertions prohibit extrapolation. All %d daily points form one continuous input; native cutedge applies only at the two full-series ends.',height(daily));
method.Interpolation=sprintf('The %d inserted values are estimates, recorded with both observed brackets, linear weights, gap lengths and IsInterpolated flags. Original NaNs remain in ObservedDaily and BMeanObserved_nT. The established linear interpolation is retained across the extended interval without a new gap-length cutoff. Long gaps have no independent observed wave information.',height(interpolationAudit));
method.Availability='End at the day of the last finite source ABS_B. Distinguish catalog/file endpoints from finite magnetic observations. Missing monthly CDFs are allowed only after direct official-directory verification. No exterior extrapolation.';
method.EdgeRule='Native cutedge c=floor(Fs/f): first c and last c+1 samples masked. Evaluate only scales with N_even>=2*c+2. Native odd final sample omission is retained and recorded. Authorized linear gap fill precedes this step; no extra physical threshold.';
method.Frequency=sprintf('Longest period %g days, minimum frequency %.12e Hz. Daily sampling Nyquist=5.787037037037e-6 Hz. Display 2--%g days, 2 days at the top. Values exactly at Nyquist have limiting sampling support.',cfg.MaxPeriodDays,cfg.ComputableFrequencyHz(1),cfg.MaxPeriodDays);
method.Electric='PWS sample verified: 10--56200 Hz; no electric-field observations in requested band. No electric spectrum or reconstructed E.';
method.Limit='Wavelet localization and endpoint effects depend on scale. Linear interpolation affects local spectral power; the interpolated days do not add independent observations. Spectral power alone does not establish propagating waves, modes or significance.';
method.Display=sprintf('Daily flat color cells; no spectral averaging or interpolation. Period axis in days: top 2, bottom %g. Fixed log10 color limits [0,3]. No method annotations inside the figure. Frequency-to-period conversion changes only the display coordinate; power values remain in nT^2/Hz.',cfg.MaxPeriodDays);
files=string(fullfile(codeRoot,{'Run_V1_Wave_1Day_30Days.m','V1_Wave_Read_Hourly.m', ...
    'V1_Daily_Morlet.m','V1_Interpolate_Daily_Gaps.m','Plot_V1_Daily_Morlet.m','Validate_V1_Daily_Morlet.m','V1_Sync_Morlet_Archive.m','Run_V1_Wave_1Day_60Days.m','Run_V1_Wave_1Day_90Days.m'})).';
files=[files;string({which('irf_wavelet');which('irf_resamp');which('Voyager_Read_CDF_Product')})];
hashes=strings(size(files));
for k=1:numel(files), hashes(k)=Case1_File_SHA256(files(k)); end
result=struct('Config',cfg,'Daily',daily,'ObservedDaily',observedDaily,'InterpolationAudit',interpolationAudit,'Spectrum',spec,'Segments',segments, ...
    'Method',method,'Sources',sources,'IdenticalDuplicateRecords',duplicates, ...
    'CodeAudit',table(files,hashes),'MATLABVersion',version,'Availability',cfg.AvailabilityAudit);
availability=cfg.AvailabilityAudit;
save(fullfile(cfg.OutputRoot,'official_availability_audit.mat'),'availability','-v7.3');
writetable(daily,fullfile(cfg.OutputRoot,'V1_daily_magnetic_input.csv'));
writetable(interpolationAudit,fullfile(cfg.OutputRoot,'V1_daily_interpolation_audit.csv'));
writetable(segments,fullfile(cfg.OutputRoot,'V1_daily_morlet_segments.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(cfg.OutputRoot,'source_manifest.csv'));
[result.FigureFiles,result.PlotAudit]=Plot_V1_Daily_Morlet(daily,spec,cfg);
result.OutputFiles=[string(fullfile(cfg.OutputRoot,{'V1_daily_morlet_analysis.mat', ...
    'official_availability_audit.mat','V1_daily_magnetic_input.csv','V1_daily_interpolation_audit.csv','V1_daily_morlet_segments.csv','source_manifest.csv','run_daily_morlet.log'})).';result.FigureFiles(:)];
save(fullfile(cfg.OutputRoot,'V1_daily_morlet_analysis.mat'),'result','raw','-v7.3');
fprintf('Finite coefficients: %d. Files saved in %s\n',nnz(isfinite(spec.Power)),cfg.OutputRoot);
end






