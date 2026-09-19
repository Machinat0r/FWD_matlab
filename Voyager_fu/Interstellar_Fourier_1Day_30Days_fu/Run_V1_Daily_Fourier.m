function result=Run_V1_Daily_Fourier(varargin)
% Voyager 1 daily scalar |B|: Fourier time-period power spectrum.
% Raw CDF -> original UTC daily means -> previously authorized linear fill.
% Complete Hann windows, per-window mean removal, one-sided FFT PSD.
% See README_日均傅里叶时频图.md for sampling and window tradeoffs.

%% Paths and configurable parameters
codeRoot=fileparts(mfilename('fullpath'));
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_1Day_30Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'ReferenceAuditFile','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days/V1_daily_morlet_analysis.mat',@(x)ischar(x)||isstring(x));
addParameter(p,'WindowDays',90,@(x)isscalar(x)&&isfinite(x)&&x>=60&&x==round(x));
addParameter(p,'MaxPeriodDays',30,@(x)isscalar(x)&&isfinite(x)&&x>2);
addParameter(p,'PSDColorLimits',[-2 3.5],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)<x(2));
addParameter(p,'Visible',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
cfg.StartUTC=datetime(2012,8,25,'TimeZone','UTC');
cfg.StopUTC=datetime(2021,12,17,'TimeZone','UTC');
cfg.StepDays=1;
cfg.SampleSeconds=86400;
cfg.SamplingFrequencyHz=1/cfg.SampleSeconds;
cfg.OriginalRequestedFrequencyHz=[1/(30*86400),1/86400];
cfg.ComputableFrequencyHz=[1/(cfg.MaxPeriodDays*86400),1/(2*86400)];
cfg.DisplayPeriodDays=[2 cfg.MaxPeriodDays];
cfg.PSDColorLimits=reshape(cfg.PSDColorLimits,1,[]);
assert(cfg.MaxPeriodDays<=cfg.WindowDays,'The longest period exceeds the nonzero FFT frequency support; increase WindowDays.');
cfg.WindowFunction='Periodic Hann';
cfg.WindowMeanRemoval=true;
cfg.NFFT=cfg.WindowDays;
waveRoot=fullfile(fileparts(codeRoot),'Interstellar_Wave_1Day_30Days_fu');
addpath(waveRoot);
addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path(cfg.IRFURoot);
if ~isfolder(cfg.OutputRoot), mkdir(cfg.OutputRoot); end
diary(fullfile(cfg.OutputRoot,'run_daily_fourier.log'));
logCleanup=onCleanup(@()diary('off'));
fprintf('Daily Fourier run UTC: %s\n',string(datetime('now','TimeZone','UTC')));
fprintf('Window: %d days; step: 1 day; periodic Hann; window mean removal; NFFT=%d.\n', ...
    cfg.WindowDays,cfg.NFFT);

%% Read original CDF and reproduce the previous daily magnetic field
[raw,sources,duplicates]=V1_Wave_Read_Hourly(cfg);
day=(cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin]=ismember(dateshift(raw.EpochUTC,'start','day'),day);
good=matched & isfinite(raw.Bmag_nT);
dailyMean=accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count=accumarray(bin(good),1,[numel(day) 1],@sum,0);
observedDaily=table(day+hours(12),dailyMean,count, ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
assert(sum(count)==nnz(good),'Hourly sample accounting failed.');
[daily,interpolationAudit]=V1_Interpolate_Daily_Gaps(observedDaily);
fprintf('%d daily values, %d observed, %d linearly interpolated.\n', ...
    height(daily),nnz(~daily.IsInterpolated),height(interpolationAudit));
referenceMatched=false; referenceSHA="";
if isfile(cfg.ReferenceAuditFile)
    reference=load(cfg.ReferenceAuditFile,'result'); % Audit only, never science input.
    assert(isequaln(observedDaily,reference.result.ObservedDaily),'Original daily data changed.');
    assert(isequaln(daily,reference.result.Daily),'Interpolated daily data changed.');
    assert(isequaln(interpolationAudit,reference.result.InterpolationAudit),'Gap interpolation changed.');
    referenceMatched=true;
    referenceSHA=string(Case1_File_SHA256(cfg.ReferenceAuditFile));
    clear reference
end

%% Fourier spectrogram and independent normalization checks
[spec,windows,validation]=V1_Daily_Fourier(daily,cfg);
validation.PreviousDailyDatasetMatched=referenceMatched;
validation.PreviousDailyAuditSHA256=referenceSHA;
disp(validation);

%% Method and source audit
method=struct;
method.UserInstruction='2026-09-11: user requested drawing after the explicit proposal of 90-day windows, 1-day hop, Hann taper and per-window mean removal; execute that proposed configuration on the prior interpolated daily input.';
method.Source='Original VOYAGER1_COHO1HR_MERGED_MAG_PLASMA CDF, Epoch and scalar ABS_B in nT. Existing CDF fill and valid-range handling.';
method.Daily='Same UTC daily arithmetic means as the prior daily wavelet figure: finite hourly scalar ABS_B, equal sample weights, no minimum count threshold; daily timestamps at noon.';
method.Interpolation='Same previously user-authorized linear filling of the 79 internal missing daily values using irf_resamp(...,''linear''). Observed values, counts, original NaNs, brackets and weights retained. No extrapolation.';
method.Input='All 3401 daily samples retained. Scalar |B| only. Magnetic panel shows the same filled daily values used by the previous wavelet figure.';
method.Fourier=sprintf('Short-time Fourier transform: %d-day complete windows, 1-day hop; subtract each window arithmetic mean, apply periodic Hann taper, NFFT=%d with no zero padding. No linear detrend or additional temporal/frequency averaging.',cfg.WindowDays,cfg.NFFT);
method.IRFU='Reuse irf_psd for one complete segment per call, explicit window vector, internal overlap=0 and dflag=mean. The caller advances windows by one day. No averaging of distinct window spectra.';
method.Normalization='One-sided PSD: abs(FFT((x-mean(x)).*w)).^2/(Fs*sum(w.^2)); double strictly positive non-Nyquist bins. Correct irf_psd doubling at DC and even-NFFT Nyquist only. Units nT^2/Hz. No 2*pi factor; this is the standard Fourier periodogram normalization.';
method.Frequency=sprintf('Fs=1/86400 Hz. Period 2--%g days corresponds to %.12e--5.787037037037e-6 Hz. Nyquist=1/(2 days). FFT bin spacing Fs/NFFT. DC retained in audit but excluded from the displayed band.',cfg.MaxPeriodDays,cfg.ComputableFrequencyHz(1));
method.Time='Only full windows. Timestamp is midpoint of the actual window; no zero/reflection padding or endpoint extrapolation. Daily color-cell width denotes output step, while each estimate uses the whole window.';
method.Display=sprintf('Time on horizontal axis, logarithmic period in days, 2 days top and %g days bottom. Native FFT bins shown as flat cells with linear-frequency midpoint boundaries. Fixed log10 color limits [%g,%g]; no extra processing annotations inside figure.',cfg.MaxPeriodDays,cfg.PSDColorLimits(1),cfg.PSDColorLimits(2));
method.Limits='Hann taper and window length set time/frequency resolution; the one-day hop does not make spectra independent daily estimates. No data outside full-window coverage are inferred. Gap interpolation may modify local spectral power; inserted values add no independent measurements. Fourier and Morlet use different estimators and normalization conventions.';
method.Electric='Retain magnetic-only output: earlier PWS check found no electric-field observations in the requested low-frequency band.';
files=string({[mfilename('fullpath'),'.m'];which('V1_Daily_Fourier'); ...
    which('Plot_V1_Daily_Fourier');which('V1_Wave_Read_Hourly'); ...
    which('V1_Interpolate_Daily_Gaps');which('Voyager_Read_CDF_Product'); ...
    which('irf_psd');which('irf_resamp');which('Run_V1_Daily_Fourier_90Days')});
hashes=strings(size(files));
for k=1:numel(files), hashes(k)=Case1_File_SHA256(files(k)); end
result=struct('Config',cfg,'Daily',daily,'ObservedDaily',observedDaily, ...
    'InterpolationAudit',interpolationAudit,'Spectrum',spec,'Windows',windows, ...
    'Validation',validation,'Method',method,'Sources',sources, ...
    'IdenticalDuplicateRecords',duplicates,'CodeAudit',table(files,hashes),'MATLABVersion',version);

%% Export plots and reproducibility products
writetable(daily,fullfile(cfg.OutputRoot,'V1_daily_magnetic_input.csv'));
writetable(interpolationAudit,fullfile(cfg.OutputRoot,'V1_daily_interpolation_audit.csv'));
writetable(windows,fullfile(cfg.OutputRoot,'V1_fourier_window_audit.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(cfg.OutputRoot,'source_manifest.csv'));
frequency=table(spec.FrequencyHz,1./spec.FrequencyHz/86400, ...
    'VariableNames',{'FrequencyHz','PeriodDays'});
writetable(frequency,fullfile(cfg.OutputRoot,'V1_fourier_frequency_bins.csv'));
[result.FigureFiles,result.PlotAudit]=Plot_V1_Daily_Fourier(daily,spec,cfg);
result.OutputFiles=[string(fullfile(cfg.OutputRoot,{'V1_daily_fourier_analysis.mat', ...
    'V1_daily_magnetic_input.csv','V1_daily_interpolation_audit.csv', ...
    'V1_fourier_window_audit.csv','V1_fourier_frequency_bins.csv','source_manifest.csv', ...
    'run_daily_fourier.log'})).';result.FigureFiles(:)];
save(fullfile(cfg.OutputRoot,'V1_daily_fourier_analysis.mat'),'result','raw','-v7.3');
fprintf('Complete: %d windows, %d displayed frequencies. Outputs: %s\n', ...
    height(windows),numel(spec.FrequencyHz),cfg.OutputRoot);
end



