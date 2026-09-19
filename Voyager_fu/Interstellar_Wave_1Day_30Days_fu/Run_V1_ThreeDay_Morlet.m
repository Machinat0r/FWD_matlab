function result = Run_V1_ThreeDay_Morlet(varargin)
% Voyager 1: existing nonoverlapping three-day |B| means -> IRFU Morlet.
% Same original-CDF, scalar |B| and wavelet conventions as the daily entry.
% 2026-09-11: repeat on the earlier three-day dataset; carry forward linear
% gap interpolation and period-in-days axis. User revised color limits to [-1,3].

%% Paths and parameters
codeRoot=fileparts(mfilename('fullpath'));
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_ThreeDay_30Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'ReferenceAuditFile','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/averaged_no_S4/three_day_monthly_audit.mat',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
cfg.StartUTC=datetime(2012,8,25,'TimeZone','UTC');
cfg.StopUTC=datetime(2021,12,17,'TimeZone','UTC'); % Exclusive.
cfg.SampleDays=3;
cfg.SampleSeconds=3*86400;
cfg.SamplingFrequencyHz=1/cfg.SampleSeconds;
cfg.NumFrequencies=100;
cfg.WaveletWidth=5.36;
cfg.OriginalRequestedFrequencyHz=[1/(30*86400),1/86400];
cfg.ComputableFrequencyHz=[1/(30*86400),1/(6*86400)];
cfg.DisplayPeriodDays=[6 30];
cfg.PSDColorLimits=[-1 3];
if ~isfolder(cfg.OutputRoot), mkdir(cfg.OutputRoot); end
diary(fullfile(cfg.OutputRoot,'run_three_day_morlet.log'));
logCleanup=onCleanup(@()diary('off'));
addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path(cfg.IRFURoot);
fprintf('Three-day Morlet run UTC: %s\n',string(datetime('now','TimeZone','UTC')));

%% Original hourly CDF -> observed UTC daily means -> original three-day bins
[raw,sources,duplicates]=V1_Wave_Read_Hourly(cfg);
day=(cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin]=ismember(dateshift(raw.EpochUTC,'start','day'),day);
good=matched & isfinite(raw.Bmag_nT);
dailyMean=accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count=accumarray(bin(good),1,[numel(day) 1],@sum,0);
observedDaily=table(day+hours(12),dailyMean,count, ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
assert(sum(count)==nnz(good),'Hourly sample accounting failed.');
[observedWindows,dailyBinIndex]=threeDayMeans(observedDaily,cfg);
fprintf('%d daily rows; %d missing daily means. %d three-day bins; %d empty bins; %d partial bins.\n', ...
    height(observedDaily),nnz(~isfinite(dailyMean)),height(observedWindows), ...
    nnz(~isfinite(observedWindows.BMean_nT)),nnz(observedWindows.PartialWindow));

%% Interpolate empty three-day means, after averaging the observed daily values
windows=observedWindows;
fullRows=find(~windows.PartialWindow);
assert(isequal(fullRows,(1:numel(fullRows)).'),'Full windows must be a contiguous prefix.');
uniform=windows(fullRows,:);
[uniform,interpolationAudit]=interpolateThreeDayGaps(uniform,cfg);
windows.BMeanObserved_nT=windows.BMean_nT;
windows.IsInterpolated=false(height(windows),1);
windows.BMean_nT(fullRows)=uniform.BMean_nT;
windows.IsInterpolated(fullRows)=uniform.IsInterpolated;
windows.IsWaveletInput=~windows.PartialWindow;
fprintf('Linearly interpolated %d empty three-day bins; observed bin means unchanged.\n',height(interpolationAudit));

%% Direct Morlet transform, at the true three-day sample rate
assert(all(seconds(diff(uniform.EpochUTC))==cfg.SampleSeconds),'Nonuniform wavelet input.');
assert(all(isfinite(uniform.BMean_nT)),'No missing values may reach irf_wavelet.');
w=irf_wavelet([posixtime(uniform.EpochUTC),uniform.BMean_nT], ...
    'fs',cfg.SamplingFrequencyHz,'nf',cfg.NumFrequencies, ...
    'f',cfg.ComputableFrequencyHz,'wavelet_width',cfg.WaveletWidth,'cutedge',1);
[f,order]=sort(w.f);
nUsed=height(uniform)-mod(height(uniform),2);
assert(isequal(w.t,posixtime(uniform.EpochUTC(1:nUsed))),'Unexpected IRFU time handling.');
power=nan(height(uniform),numel(f));
power(1:nUsed,:)=w.p{1}(:,order);
windows.WaveletSampleReturned=false(height(windows),1);
windows.WaveletSampleReturned(fullRows(1:nUsed))=true;
spec=struct('TimeUTC',uniform.EpochUTC,'FrequencyHz',f(:),'Power',power, ...
    'SamplingFrequencyHz',cfg.SamplingFrequencyHz,'NyquistHz',cfg.SamplingFrequencyHz/2, ...
    'Units','nT^2/Hz','Estimator','irf_wavelet Morlet', ...
    'InterpolatedInputMask',uniform.IsInterpolated);
segments=table(uniform.EpochUTC(1),uniform.EpochUTC(end),height(uniform),nUsed, ...
    height(uniform)-nUsed,nnz(windows.PartialWindow),nnz(isfinite(power)), ...
    'VariableNames',{'FirstUTC','LastUTC','FullThreeDaySamples','EvenSamplesUsed', ...
    'OddFinalSampleOmitted','PartialWindowsExcluded','FiniteCoefficients'});

%% Verify the previous dataset and the sampling, interpolation and edge handling
validation=validateThreeDay(observedDaily,observedWindows,windows,dailyBinIndex, ...
    interpolationAudit,spec,cfg,nUsed);
disp(validation);
method=struct;
method.Source='Original VOYAGER1_COHO1HR_MERGED_MAG_PLASMA CDF: Epoch, scalar ABS_B (nT); fill and validity metadata handled by the existing reader.';
method.Daily='Arithmetic mean of finite hourly ABS_B in each UTC day. No filling of the daily series before three-day grouping.';
method.ThreeDay='Same panel-a values as Run_V1_ThreeDay_Monthly_Overview: nonoverlapping [2012-08-25+3*k days,2012-08-25+3*(k+1) days). Equal weight for each finite daily mean; no minimum coverage threshold. Preserve observed values and counts.';
method.Time='Epoch is the actual clipped window midpoint. The final two-day partial window is retained in the magnetic curve and table, excluded from the uniformly sampled wavelet input. Native irf_wavelet odd final-sample omission is also retained and recorded.';
method.Interpolation='Carry forward user-authorized linear gap interpolation, now applied only to empty complete three-day windows after the original bin averaging. Explicit irf_resamp(...,''linear'') with observed brackets; no extrapolation or changes to observed bins. All brackets, weights and gap lengths saved.';
method.Wavelet='Direct irf_wavelet as in MMS_fu/Wave.m; fs=1/(3*86400), 100 logarithmic frequencies, Morlet width 5.36, default returnpower=1 and cutedge=1. Whole continuous record; no fixed-duration sliding window.';
method.Normalization='Native IRFU P=2*pi*abs(W)^2/f, nT^2/Hz; scalar |B| power. No vector trace, prefilter, mean subtraction, detrending, extra spectral averaging, or renormalization.';
method.Frequency='One month=30 days; computed band [1/(30*86400),1/(6*86400)] Hz = [3.858024691358e-7,1.929012345679e-6] Hz. Period display 6--30 days; 6 days at top. Six days is the Nyquist boundary, with limiting sampling support.';
method.Display='Same two-panel style; period label in days, native three-day time cells and log10 color limits [-1,3]. Power remains nT^2/Hz. No added method annotations in the figure.';
method.Limit='Three-day averaging suppresses short-period variations; linear interpolation contributes estimates and can alter local wavelet power. It adds no independent observations. Native scale-dependent cutedge is retained at the full-record endpoints.';
method.Electric='Prior PWS availability check found no electric-field data in this low-frequency band; magnetic-only analysis is retained.';
files=string({[mfilename('fullpath'),'.m'];which('Plot_V1_ThreeDay_Morlet'); ...
    which('V1_Wave_Read_Hourly');which('Voyager_Read_CDF_Product'); ...
    which('irf_wavelet');which('irf_resamp')});
hashes=strings(size(files));
for k=1:numel(files), hashes(k)=Case1_File_SHA256(files(k)); end
result=struct('Config',cfg,'ObservedDaily',observedDaily,'ObservedWindows',observedWindows, ...
    'Windows',windows,'DailyBinIndex',dailyBinIndex,'InterpolationAudit',interpolationAudit, ...
    'Spectrum',spec,'Segments',segments,'Validation',validation,'Method',method, ...
    'Sources',sources,'IdenticalDuplicateRecords',duplicates, ...
    'CodeAudit',table(files,hashes),'MATLABVersion',version);

%% Export figures, input tables and reproducibility audit
writetable(observedDaily,fullfile(cfg.OutputRoot,'V1_observed_daily_magnetic_input.csv'));
writetable(windows,fullfile(cfg.OutputRoot,'V1_three_day_magnetic_input.csv'));
writetable(interpolationAudit,fullfile(cfg.OutputRoot,'V1_three_day_interpolation_audit.csv'));
writetable(segments,fullfile(cfg.OutputRoot,'V1_three_day_morlet_segments.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(cfg.OutputRoot,'source_manifest.csv'));
[result.FigureFiles,result.PlotAudit]=Plot_V1_ThreeDay_Morlet(windows,spec,cfg);
result.OutputFiles=[string(fullfile(cfg.OutputRoot,{'V1_three_day_morlet_analysis.mat', ...
    'V1_observed_daily_magnetic_input.csv','V1_three_day_magnetic_input.csv', ...
    'V1_three_day_interpolation_audit.csv','V1_three_day_morlet_segments.csv', ...
    'source_manifest.csv','run_three_day_morlet.log'})).';result.FigureFiles(:)];
save(fullfile(cfg.OutputRoot,'V1_three_day_morlet_analysis.mat'),'result','raw','-v7.3');
fprintf('Complete: %d finite wavelet coefficients. Outputs: %s\n',nnz(isfinite(power)),cfg.OutputRoot);
end

function [windows,binIndex]=threeDayMeans(daily,cfg)
%% Match the previously approved nonoverlapping three-day means exactly
begin=(cfg.StartUTC:caldays(3):cfg.StopUTC-seconds(1)).';
nominalEnd=begin+caldays(3);
finish=nominalEnd; finish(finish>cfg.StopUTC)=cfg.StopUTC;
value=nan(numel(begin),1); count=zeros(numel(begin),1);
hourlyCount=count; binIndex=zeros(height(daily),1);
for k=1:numel(begin)
    rows=daily.EpochUTC>=begin(k) & daily.EpochUTC<finish(k);
    assert(all(binIndex(rows)==0),'Overlapping bins.');
    binIndex(rows)=k;
    x=daily.BMean_nT(rows); x=x(isfinite(x));
    count(k)=numel(x); hourlyCount(k)=sum(daily.MAGSampleCount(rows));
    if ~isempty(x), value(k)=mean(x); end
end
assert(all(binIndex>0),'Unassigned daily rows.');
windows=table(begin,nominalEnd,begin,finish,begin+(finish-begin)/2, ...
    days(finish-begin),finish~=nominalEnd,value,count,hourlyCount, ...
    'VariableNames',{'NominalStartUTC','NominalEndUTC','StartUTC','EndUTCExclusive', ...
    'EpochUTC','CalendarDaysInWindow','PartialWindow','BMean_nT','ValidDays_A','MAGSampleCount'});
end

function [filled,audit]=interpolateThreeDayGaps(observed,cfg)
%% Fill internal empty complete bins only; use no extrapolation
filled=observed;
x=observed.BMean_nT; missing=~isfinite(x);
t=seconds(observed.EpochUTC-observed.EpochUTC(1));
assert(all(diff(t)==cfg.SampleSeconds),'Expected a complete three-day grid.');
assert(~missing(1) && ~missing(end),'Endpoint extrapolation is not authorized.');
filled.BMeanObserved_nT=x; filled.IsInterpolated=missing;
if any(missing)
    q=irf_resamp([t(~missing),x(~missing)],t(missing),'linear');
    filled.BMean_nT(missing)=q(:,2);
end
rows=find(missing); previous=zeros(size(rows)); following=previous;
valid=find(~missing);
for k=1:numel(rows)
    previous(k)=valid(find(valid<rows(k),1,'last'));
    following(k)=valid(find(valid>rows(k),1,'first'));
end
alpha=(t(rows)-t(previous))./(t(following)-t(previous));
gapBins=following-previous-1;
audit=table(observed.EpochUTC(rows),observed.EpochUTC(previous), ...
    observed.EpochUTC(following),x(previous),x(following),alpha, ...
    filled.BMean_nT(rows),gapBins,cfg.SampleDays*gapBins, ...
    'VariableNames',{'InterpolatedEpochUTC','PreviousObservedUTC','NextObservedUTC', ...
    'PreviousB_nT','NextB_nT','WeightOfNextValue','InterpolatedB_nT', ...
    'EmptyWindowsInGap','CalendarDaysInEmptyWindows'});
end

function report=validateThreeDay(daily,observed,filled,binIndex,audit,spec,cfg,nUsed)
%% Independent bin aggregation, original-product comparison and frequency check
good=isfinite(daily.BMean_nT);
independent=accumarray(binIndex(good),daily.BMean_nT(good),[height(observed) 1],@mean,NaN);
assert(isequaln(independent,observed.BMean_nT),'Independent bin average mismatch.');
assert(sum(observed.ValidDays_A)==nnz(good),'Finite daily count mismatch.');
assert(sum(observed.MAGSampleCount)==sum(daily.MAGSampleCount),'Hourly count mismatch.');
measured=isfinite(observed.BMean_nT);
assert(isequal(filled.BMean_nT(measured),observed.BMean_nT(measured)),'Observed bin changed.');
assert(isequal(filled.ValidDays_A,observed.ValidDays_A),'Observed counts changed.');
expected=(1-audit.WeightOfNextValue).*audit.PreviousB_nT+audit.WeightOfNextValue.*audit.NextB_nT;
assert(all(abs(expected-audit.InterpolatedB_nT)<1e-12),'Interpolation differs from bracket formula.');
assert(all(audit.WeightOfNextValue>0 & audit.WeightOfNextValue<1),'Extrapolation detected.');
assert(height(audit)==nnz(filled.IsInterpolated),'Interpolation count mismatch.');
assert(all(isfinite(filled.BMean_nT(filled.IsWaveletInput))),'Unfilled transform input.');
assert(~any(filled.PartialWindow & filled.IsWaveletInput),'Partial bin included in transform.');
assert(all(spec.FrequencyHz<=spec.NyquistHz*(1+1e-12)),'Frequency exceeds Nyquist.');
assert(all(diff(spec.FrequencyHz)>0),'Frequency ordering invalid.');
% Construct the native censor mask independently for each frequency.
scale=logspace(log10(0.5*cfg.SamplingFrequencyHz/cfg.ComputableFrequencyHz(2)), ...
    log10(0.5*cfg.SamplingFrequencyHz/cfg.ComputableFrequencyHz(1)),cfg.NumFrequencies);
edgeCounts=flipud(floor(2*scale(:)));
expectedFinite=false(size(spec.Power));
for j=1:numel(spec.FrequencyHz)
    expectedFinite(edgeCounts(j)+1:nUsed-edgeCounts(j)-1,j)=true;
end
assert(isequal(isfinite(spec.Power),expectedFinite),'Unexpected internal holes or incorrect edge mask.');
referenceMatched=false; referenceSHA="";
if isfile(cfg.ReferenceAuditFile)
    reference=load(cfg.ReferenceAuditFile,'result'); % Validation only; never a science input.
    columns={'NominalStartUTC','NominalEndUTC','StartUTC','EndUTCExclusive', ...
        'EpochUTC','CalendarDaysInWindow','PartialWindow','BMean_nT','ValidDays_A'};
    assert(isequaln(observed(:,columns),reference.result.three_day.Windows(:,columns)), ...
        'Three-day data do not match the previous overview.');
    referenceMatched=true; referenceSHA=string(Case1_File_SHA256(cfg.ReferenceAuditFile));
end
% Physical sampling test: a 12-day sinusoid sampled every 3 days peaks near 12 days.
sample=(0:599).'; testTime=posixtime(cfg.StartUTC)+sample*cfg.SampleSeconds;
testWave=irf_wavelet([testTime,sin(2*pi*sample*cfg.SampleDays/12)], ...
    'fs',cfg.SamplingFrequencyHz,'nf',cfg.NumFrequencies,'f',cfg.ComputableFrequencyHz, ...
    'wavelet_width',cfg.WaveletWidth,'cutedge',1);
testPower=mean(testWave.p{1}(100:500,:),1);
[~,peak]=max(testPower); peakDays=1/testWave.f(peak)/86400;
assert(abs(peakDays-12)<0.5,'Three-day sampling frequency conversion failed.');
report=struct('IndependentBinMeansExact',true,'ObservedValuesUnchanged',true, ...
    'InterpolationChecked',true,'NativeEdgeMaskExact',true, ...
    'PreviousThreeDayDatasetMatched',referenceMatched,'ReferenceAuditSHA256',referenceSHA, ...
    'TotalWindows',height(observed),'MeasuredWindows',nnz(measured), ...
    'InterpolatedWindows',height(audit),'WaveletInputSamples',nnz(filled.IsWaveletInput), ...
    'WaveletReturnedSamples',nUsed,'PartialWindowsExcluded',nnz(filled.PartialWindow), ...
    'Sinusoid12DayRecoveredPeriodDays',peakDays);
end


