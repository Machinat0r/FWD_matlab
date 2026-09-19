function result=Run_V1_Fourier_PSD_2015(varargin)
% Whole-record Fourier PSD of the 2015 daily scalar magnetic magnitude.
% Original CDF -> UTC daily means -> established linear interpolation.
% One periodic Hann window over the entire one-year record, one mean removed.
% No sliding windows, zero padding, linear detrending or spectral smoothing.

%% Paths and analysis interval
codeRoot=fileparts(mfilename('fullpath'));
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_PSD_2015',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'ReferenceAuditFile','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_PSD_2014_2015/V1_Fourier_PSD_2014_2015_analysis.mat',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
cfg.StartUTC=datetime(2015,1,1,'TimeZone','UTC');
cfg.StopUTC=datetime(2016,1,1,'TimeZone','UTC'); % Includes all of 2015-12-31.
cfg.SampleSeconds=86400;
cfg.SamplingFrequencyHz=1/cfg.SampleSeconds;
cfg.DisplayFrequencyHz=[1/(90*86400),1/(2*86400)];
cfg.DisplayPeriodDays=[2 90];
cfg.ComputableFrequencyHz=cfg.DisplayFrequencyHz;
cfg.WindowFunction='Periodic Hann over the entire selected interval';
cfg.WindowMeanRemoval=true;
cfg.StepDays=1; % One output only because WindowDays equals the entire record.
addpath(fullfile(fileparts(codeRoot),'Interstellar_Wave_1Day_30Days_fu'));
addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path(cfg.IRFURoot);
if ~isfolder(cfg.OutputRoot), mkdir(cfg.OutputRoot); end
diary(fullfile(cfg.OutputRoot,'run_fourier_psd_2015.log'));
logCleanup=onCleanup(@()diary('off'));
fprintf('Whole-record PSD run UTC: %s\n',string(datetime('now','TimeZone','UTC')));

%% Read only the original CDF records in the selected year
[raw,sources,duplicates]=V1_Wave_Read_Hourly(cfg);
day=(cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin]=ismember(dateshift(raw.EpochUTC,'start','day'),day);
good=matched & isfinite(raw.Bmag_nT);
meanB=accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count=accumarray(bin(good),1,[numel(day) 1],@sum,0);
observedDaily=table(day+hours(12),meanB,count, ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
assert(height(observedDaily)==365,'Expected exactly 365 UTC days.');
assert(sum(count)==nnz(good),'Hourly sample accounting failed.');
[daily,interpolationAudit]=V1_Interpolate_Daily_Gaps(observedDaily);
fprintf('%d days: %d observed, %d previously authorized linear fills.\n', ...
    height(daily),nnz(~daily.IsInterpolated),height(interpolationAudit));
referenceMatched=false; referenceSHA="";
if isfile(cfg.ReferenceAuditFile)
    prior=load(cfg.ReferenceAuditFile,'result'); % Read only to verify the selected values.
    take=prior.result.Daily.EpochUTC>=cfg.StartUTC & prior.result.Daily.EpochUTC<cfg.StopUTC;
    assert(isequaln(daily,prior.result.Daily(take,:)),'Selected daily series changed.');
    takeObserved=prior.result.ObservedDaily.EpochUTC>=cfg.StartUTC & prior.result.ObservedDaily.EpochUTC<cfg.StopUTC;
    assert(isequaln(observedDaily,prior.result.ObservedDaily(takeObserved,:)),'Observed values changed.');
    oldAudit=prior.result.InterpolationAudit;
    takeAudit=oldAudit.InterpolatedDayUTC>=cfg.StartUTC & oldAudit.InterpolatedDayUTC<cfg.StopUTC;
    assert(isequaln(interpolationAudit,oldAudit(takeAudit,:)),'Interpolation changed.');
    referenceMatched=true; referenceSHA=string(Case1_File_SHA256(cfg.ReferenceAuditFile));
    clear prior oldAudit
end

%% One Fourier transform over the full selected interval
cfg.WindowDays=height(daily);
cfg.NFFT=height(daily);
[spec,windowAudit,validation]=V1_Daily_Fourier(daily,cfg);
assert(height(windowAudit)==1,'Whole-record PSD must contain exactly one spectrum.');
assert(windowAudit.StartUTC==cfg.StartUTC && windowAudit.EndUTCExclusive==cfg.StopUTC, ...
    'Fourier input interval differs from the requested dates.');
assert(spec.FFTLength==365,'Unexpected FFT length.');
validation.PreviousDailySubsetExact=referenceMatched;
validation.ReferenceAuditSHA256=referenceSHA;
validation.WholeRecordOnly=true;
disp(validation);

%% Full PSD table (DC and out-of-display-band values remain in the audit)
frequency=spec.AllFrequencyHz(:);
power=spec.AllPower(:);
period=nan(size(frequency));
period(frequency>0)=1./frequency(frequency>0)/86400;
psdTable=table(frequency,period,power,spec.BandMask(:), ...
    'VariableNames',{'FrequencyHz','PeriodDays','PSD_nT2_per_Hz','InDisplayedBand'});

%% Reproducibility record
method=struct;
method.UserInstruction='Use 2015-01-01 through 2015-12-31 only. Latest user request: display period in days with the horizontal axis reversed, 90 days on the left and 2 days on the right.';
method.Source='12 original monthly VOYAGER1_COHO1HR_MERGED_MAG_PLASMA CDF files; Epoch and scalar ABS_B (nT). Source fill/valid-range handling unchanged.';
method.Daily='Same finite hourly ABS_B arithmetic means in complete UTC days, timestamps at noon; no added count threshold. Only the selected 365 days are included.';
method.Gaps='Use the already authorized linear interpolation between original daily observations. No extrapolation. Original NaNs, observed sample counts, interpolation endpoints and weights are retained. Values match the same subset of the earlier daily analysis exactly.';
method.Fourier='One 365-sample transform over the entire selected record. Subtract the full-record arithmetic mean once and apply the previously used periodic Hann taper over this record. No 90-day sliding, no averaging of window PSDs, no linear detrending, no frequency smoothing, no zero padding.';
method.Normalization='Reuse V1_Daily_Fourier and irf_psd with one complete segment. PSD=abs(FFT((x-mean(x)).*w)).^2/(Fs*sum(w.^2)); double positive non-Nyquist bins only. Correct IRFU DC and even-length Nyquist weights. nT^2/Hz.';
method.Frequency='Fs=1/86400 Hz; NFFT=365; native spacing Fs/365. Display [1/(90*86400),1/(2*86400)] Hz, retaining the preceding 2--90 day band. No forced frequency point at the exact 90-day boundary; the full original FFT grid is saved, including DC.';
method.Plot='One black PSD line, logarithmic period (days) and PSD axes. T_days=1/f_Hz/86400; sort by increasing period and reorder matching PSD values. This changes only the horizontal coordinate; PSD remains nT^2/Hz. Display 90 days at the left and 2 days at the right (XDir=reverse), using the same native FFT bins. No spectral interpolation, colorbar or method footnotes.';
method.Limits='The whole-record spectrum describes the selected year together. Periodic Hann weighting and linear gap interpolation affect the estimate. Unsmooth periodogram values can fluctuate substantially between adjacent frequencies. No peak significance or fitted spectral index is inferred.';
files=string({[mfilename('fullpath'),'.m'];which('V1_Daily_Fourier'); ...
    which('V1_Wave_Read_Hourly');which('V1_Interpolate_Daily_Gaps'); ...
    which('Voyager_Read_CDF_Product');which('irf_psd');which('irf_resamp')});
hashes=strings(size(files));
for k=1:numel(files), hashes(k)=Case1_File_SHA256(files(k)); end
result=struct('Config',cfg,'Daily',daily,'ObservedDaily',observedDaily, ...
    'InterpolationAudit',interpolationAudit,'Spectrum',spec,'PSDTable',psdTable, ...
    'WindowAudit',windowAudit,'Validation',validation,'Sources',sources, ...
    'IdenticalDuplicateRecords',duplicates,'Method',method, ...
    'CodeAudit',table(files,hashes),'MATLABVersion',version);

%% Plot the PSD line and export
[result.FigureFiles,result.PlotAudit]=drawPSD(spec,cfg);
writetable(daily,fullfile(cfg.OutputRoot,'V1_daily_magnetic_input_2015.csv'));
writetable(interpolationAudit,fullfile(cfg.OutputRoot,'V1_daily_interpolation_audit_2015.csv'));
writetable(psdTable,fullfile(cfg.OutputRoot,'V1_Fourier_PSD_2015.csv'));
writetable(windowAudit,fullfile(cfg.OutputRoot,'V1_full_record_window_audit.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(cfg.OutputRoot,'source_manifest.csv'));
result.OutputFiles=[string(fullfile(cfg.OutputRoot,{'V1_Fourier_PSD_2015_analysis.mat', ...
    'V1_daily_magnetic_input_2015.csv','V1_daily_interpolation_audit_2015.csv', ...
    'V1_Fourier_PSD_2015.csv','V1_full_record_window_audit.csv', ...
    'source_manifest.csv','run_fourier_psd_2015.log'})).';result.FigureFiles(:)];
save(fullfile(cfg.OutputRoot,'V1_Fourier_PSD_2015_analysis.mat'),'result','raw','-v7.3');
fprintf('Saved one whole-record PSD, %d plotted bins; df=%.12e Hz.\n',numel(spec.FrequencyHz),spec.FrequencySpacingHz);
end

function [files,audit]=drawPSD(spec,cfg)
visibility='off'; if cfg.Visible, visibility='on'; end
fig=figure('Color','w','Visible',visibility,'Position',[70 70 1350 780]);
ax=axes(fig,'Position',[0.12 0.14 0.83 0.76]);
f=spec.FrequencyHz(:); p=spec.Power(:);
[periodDays,order]=sort(1./f/86400);
periodPower=p(order);
h=loglog(ax,periodDays,periodPower,'k-','LineWidth',0.9);
assert(isequal(h.XData(:),periodDays) && isequal(h.YData(:),p(order)), ...
    'Plotted coordinates or PSD values changed.');
set(ax,'XLim',cfg.DisplayPeriodDays,'XDir','reverse','FontSize',14,'TickDir','out','Box','on', ...
    'XGrid','on','YGrid','on','GridAlpha',0.13,'XMinorTick','on','YMinorTick','on');
periodTicks=[2 3 5 10 20 30 60 90];
set(ax,'XTick',periodTicks,'XTickLabel',string(periodTicks));
xlabel(ax,'周期（天）','FontName','Microsoft YaHei','FontSize',15);
ylabel(ax,'PSD_{|B|} (nT^2 Hz^{-1})','FontSize',15);
title(ax,'Voyager 1 | 2015-01-01 to 2015-12-31','FontSize',16,'FontWeight','bold');
base=fullfile(cfg.OutputRoot,'V1_Fourier_PSD_2015_2day_90days');
exportgraphics(fig,[base '.png'],'Resolution',200);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
savefig(fig,[base '.fig']);
files=string({[base '.png'];[base '.pdf'];[base '.fig']});
audit=struct('XScale','log','YScale','log','XDirection',get(ax,'XDir'),'XUnits','days','YUnits','nT^2/Hz', ...
    'PeriodLimitsDays',cfg.DisplayPeriodDays,'FrequencyLimitsHz',cfg.DisplayFrequencyHz, ...
    'PlottedBins',numel(f),'PlotOrderInSpectrum',order, ...
    'ActualMinimumPlottedPeriodDays',periodDays(1),'ActualMaximumPlottedPeriodDays',periodDays(end), ...
    'MethodAnnotationsAdded',false,'PlottedValuesExact',true);
if ~cfg.Visible, close(fig); end
end



