function result = Run_V1_PSD_ThreeMethods_20140601_20150601(varargin)
% Reproduce the supplied FFT, wavelet-mean and PMTM methods for Voyager 1.
% Raw CDF -> UTC daily scalar |B| means -> previously authorized linear fills.
% Scientific adaptations and original reference hashes are saved with results.

%% Paths and common settings
codeRoot = fileparts(mfilename('fullpath'));
p = inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager');
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_PSD_ThreeMethods_20140601_20150601');
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master');
addParameter(p,'ReferenceRoot','C:/Users/Administrator/.codex/attachments/70180ef1-3e39-4b3f-924d-a251573827cd');
addParameter(p,'ReferenceAuditFile','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_PSD_2014_2015/V1_Fourier_PSD_2014_2015_analysis.mat');
addParameter(p,'Visible',false,@islogical);
parse(p,varargin{:}); cfg = p.Results;
cfg.StartUTC = datetime(2014,6,1,'TimeZone','UTC');
cfg.StopUTC = datetime(2015,6,2,'TimeZone','UTC');
cfg.SampleSeconds = 86400;
cfg.SamplingFrequencyHz = 1/cfg.SampleSeconds;
cfg.DisplayPeriodDays = [2 90];
cfg.DisplayFrequencyHz = [1/(90*86400),1/(2*86400)];
cfg.NW = 3.5;
cfg.WaveletFrequencyCount = 200;
cfg.WaveletWidth = 5.36;
cfg.RemoveMean = false;
cfg.RemoveTrend = false;
addpath(fullfile(fileparts(codeRoot),'Interstellar_Wave_1Day_30Days_fu'));
addpath(fullfile(fileparts(codeRoot),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path(cfg.IRFURoot);
assert(exist('pmtm','file')==2 && exist('dpss','file')==2,'Signal Processing Toolbox is required.');
if ~isfolder(cfg.OutputRoot), mkdir(cfg.OutputRoot); end
diary(fullfile(cfg.OutputRoot,'run_three_methods.log'));
logCleanup = onCleanup(@()diary('off'));
fprintf('Three-method analysis UTC: %s\n',string(datetime('now','TimeZone','UTC')));

%% Original CDF input and the established daily interpolation
[raw,sources,duplicates] = V1_Wave_Read_Hourly(cfg);
day = (cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin] = ismember(dateshift(raw.EpochUTC,'start','day'),day);
good = matched & isfinite(raw.Bmag_nT);
meanB = accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count = accumarray(bin(good),1,[numel(day) 1],@sum,0);
observedDaily = table(day+hours(12),meanB,count, ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
[daily,interpolationAudit] = V1_Interpolate_Daily_Gaps(observedDaily);
assert(height(daily)==366 && sum(count)==nnz(good),'Daily sample accounting failed.');
assert(nnz(~daily.IsInterpolated)==360 && height(interpolationAudit)==6, ...
    'Observed/interpolated day counts differ from the previous figure; inspect source changes.');
assert(all(isfinite(daily.BMean_nT)) && all(seconds(diff(daily.EpochUTC))==86400));
prior = load(cfg.ReferenceAuditFile,'result'); % Validation only; never science input.
inPrior=prior.result.Daily.EpochUTC>=cfg.StartUTC & prior.result.Daily.EpochUTC<cfg.StopUTC;
inObserved=prior.result.ObservedDaily.EpochUTC>=cfg.StartUTC & prior.result.ObservedDaily.EpochUTC<cfg.StopUTC;
inGaps=prior.result.InterpolationAudit.InterpolatedDayUTC>=cfg.StartUTC & ...
    prior.result.InterpolationAudit.InterpolatedDayUTC<cfg.StopUTC;
assert(isequaln(daily,prior.result.Daily(inPrior,:)),'Daily input differs from the prior interval subset.');
assert(isequaln(observedDaily,prior.result.ObservedDaily(inObserved,:)),'Observed daily input changed.');
assert(isequaln(interpolationAudit,prior.result.InterpolationAudit(inGaps,:)),'Gap treatment changed.');
clear prior
fprintf('366 daily values: 360 observed + 6 established linear fills; prior subset matches exactly.\n');

%% Three attached methods on the same complete daily series
[spectra,details] = V1_Three_Spectral_Methods(posixtime(daily.EpochUTC),daily.BMean_nT,cfg);
validation = Validate_V1_Three_Spectral_Methods(daily,spectra,details,cfg);
validation.PreviousDailyInputExact = true;
validation.ReferenceAuditSHA256 = Case1_File_SHA256(cfg.ReferenceAuditFile);
disp(validation);

%% Processing and provenance audit
method = struct;
method.Interval = '2014-06-01 00:00 UTC <= Epoch < 2015-06-02 00:00 UTC; include all of the requested final day 2015-06-01.';
method.Input = 'Original monthly VOYAGER1_COHO1HR_MERGED_MAG_PLASMA CDF, Epoch and scalar ABS_B, nT. The quantity is |B|, not a component or the sum of vector PSDs.';
method.Daily = 'Arithmetic mean of finite hourly ABS_B in each complete UTC day, no new minimum count threshold; timestamps at noon. Original fill/valid-range handling is unchanged.';
method.Interpolation = 'Previously user-authorized linear interpolation of missing daily observations only: 6 of 366 days. No extrapolation. Original values, sample counts, neighbors and interpolation weights are retained.';
method.InputConditioning = 'The three supplied scripts do not explicitly remove mean or trend. Feed the same daily |B| directly to all three methods; no extra detrending, demeaning, filtering or spectral smoothing. This differs from the earlier periodic-Hann, mean-removed figure.';
method.FFT = 'Adapt the Hamming-window FFT core of turbulence_power_law_fft.m to one full 366-day record. The source selects an event frame; no Voyager event frame is requested here. Local irf_wavefft converts its frame parameters from milliseconds to samples despite the source-script point comment, and omits an exact whole-length frame and even Nyquist bin. Implement the same fft(x.*hamming(N)) core directly, retaining all 366 inputs and the Nyquist bin. No zero padding.';
method.FFTNormalization = 'Source script uses raw abs(FFT)^2. To retain the requested nT^2/Hz label, divide by Fs*sum(hamming(N).^2), then double positive non-Nyquist bins only. Hamming is MATLAB symmetric default. Full DC and out-of-band powers remain in the saved table.';
method.Wavelet = 'Call irf_wavelet with the supplied-script defaults nf=200, Morlet width=5.36, returnpower=1, cutedge=1; explicitly set frequency limits to 2--90 days. Follow nanmean(pB,1) by averaging finite power over the requested interval at each frequency. No B-to-E resampling is needed for a magnetic-only daily series.';
method.WaveletNormalization = 'Preserve native IRFU power=2*pi*abs(W).^2/f, nT^2/Hz. Preserve its default edge NaNs: scale-dependent valid counts and exact times are audited. Do not renormalize the curve to the FFT, fill edge powers, or interpret different frequency grids as equal resolving power. 366 is even, so no last sample is removed.';
method.PMTM = 'Exactly pmtm(x,3.5,[],Fs) as supplied: DPSS tapers, default adaptive weighting, one-sided PSD, default NFFT=512 for 366 inputs. The default drops the final of seven generated DPSS tapers, using six. Numerical equality to the explicit six-taper call is checked.';
method.PMTMBandwidth = 'NW=3.5 is dimensionless. W=NW/(N*dt) Hz is the half-bandwidth; nominal full smoothing width is 2W. Padding from 366 to 512 densifies the output frequency grid without adding observations or improving intrinsic resolution.';
method.Plot = 'Separate black PSD lines; log-log axes, period in days decreasing left-to-right from 90 to 2. PSD remains per Hz when the horizontal coordinate changes to period. Shared y limits cover every positive displayed value. Only titles and normal axis labels are added.';
method.Limits = 'Daily averaging and the authorized interpolation affect the spectra. Input mean is retained as in the attached scripts, so its spectral leakage can affect estimates. Wavelet edge removal changes time weights with scale; PMTM and FFT use different taper weights. No significance test, power-law fit, physical mode attribution or matched-frequency interpolation is applied. Daily sampling reaches Nyquist at two days.';
refNames = {'turbulence_power_law_fft.m';'turbulence_power_law_wavlet.m';'turbulence_power_law_pmtm.m'};
referenceFiles = string(fullfile(cfg.ReferenceRoot,refNames));
referenceHashes = strings(size(referenceFiles));
for k=1:numel(referenceFiles), referenceHashes(k)=Case1_File_SHA256(referenceFiles(k)); end
codeFiles = string({[mfilename('fullpath'),'.m'];which('V1_Three_Spectral_Methods'); ...
    which('Validate_V1_Three_Spectral_Methods');which('V1_Wave_Read_Hourly'); ...
    which('V1_Interpolate_Daily_Gaps');which('Voyager_Read_CDF_Product'); ...
    which('irf_wavelet');which('irf_wavefft');which('irf_resamp');which('pmtm');which('dpss');which('hamming')});
codeHashes = strings(size(codeFiles));
for k=1:numel(codeFiles), codeHashes(k)=Case1_File_SHA256(codeFiles(k)); end
result = struct('Config',cfg,'Daily',daily,'ObservedDaily',observedDaily, ...
    'InterpolationAudit',interpolationAudit,'Spectra',spectra,'Details',details, ...
    'Validation',validation,'Sources',sources,'IdenticalDuplicateRecords',duplicates, ...
    'Method',method,'CodeAudit',table(codeFiles,codeHashes), ...
    'ReferenceScripts',table(referenceFiles,referenceHashes),'MATLABVersion',version);

%% Three separate figures with common axes
allPower = [];
for k=1:numel(spectra)
    allPower = [allPower; spectra(k).Table.PSD_nT2_per_Hz(spectra(k).Table.InDisplayedBand)]; %#ok<AGROW>
end
positivePower = allPower(isfinite(allPower) & allPower>0);
ylimits = [min(1e-2,10^floor(log10(min(positivePower)))), ...
    max(1e4,10^ceil(log10(max(positivePower))))];
figureFiles = strings(0,1); plotAudit = struct([]);
for k=1:numel(spectra)
    [files,audit] = drawPSD(spectra(k),cfg,ylimits);
    figureFiles = [figureFiles;files]; %#ok<AGROW>
    if k==1, plotAudit=audit; else, plotAudit(k)=audit; end
    writetable(spectra(k).Table,fullfile(cfg.OutputRoot,['V1_PSD_' char(spectra(k).ID) '_20140601_20150601.csv']));
end
result.FigureFiles = figureFiles;
result.PlotAudit = plotAudit;
writetable(daily,fullfile(cfg.OutputRoot,'V1_daily_input_20140601_20150601.csv'));
writetable(interpolationAudit,fullfile(cfg.OutputRoot,'V1_daily_interpolation_audit.csv'));
writetable(details.WaveletAveragingAudit,fullfile(cfg.OutputRoot,'wavelet_averaging_audit.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(cfg.OutputRoot,'source_manifest.csv'));
writetable(result.CodeAudit,fullfile(cfg.OutputRoot,'code_manifest.csv'));
writetable(result.ReferenceScripts,fullfile(cfg.OutputRoot,'reference_scripts_manifest.csv'));
fid = fopen(fullfile(cfg.OutputRoot,'processing_audit.txt'),'w','n','UTF-8');
assert(fid>=0);
fields = fieldnames(method);
for k=1:numel(fields), fprintf(fid,'%s:\n%s\n\n',fields{k},method.(fields{k})); end
fprintf(fid,'Frequency range (Hz): %.15g to %.15g\n',cfg.DisplayFrequencyHz);
fprintf(fid,'PMTM: NFFT=%d, K=%d, half-bandwidth=%.15g Hz, full-width=%.15g Hz\n', ...
    details.PMTM.NFFT,details.PMTM.TaperCount,details.PMTM.HalfBandwidthHz,2*details.PMTM.HalfBandwidthHz);
fprintf(fid,'Wavelet contributing samples: %d to %d per frequency\n', ...
    min(details.WaveletAveragingAudit.ValidTimeCount),max(details.WaveletAveragingAudit.ValidTimeCount));
fclose(fid);
outMat = fullfile(cfg.OutputRoot,'V1_PSD_ThreeMethods_20140601_20150601_analysis.mat');
listing = dir(fullfile(cfg.OutputRoot,'*'));
result.OutputFiles = [string(fullfile({listing(~[listing.isdir]).folder},{listing(~[listing.isdir]).name})).';string(outMat)];
save(outMat,'result','raw','-v7.3');
fprintf('Saved three independent PSD figures. Shared y limits: %.3g to %.3g nT^2/Hz.\n',ylimits);
end

function [files,audit] = drawPSD(spec,cfg,ylimits)
visibility='off'; if cfg.Visible, visibility='on'; end
fig = figure('Color','w','Visible',visibility,'Position',[70 70 1350 780]);
ax = axes(fig,'Position',[0.12 0.14 0.83 0.76]);
use = spec.Table.InDisplayedBand;
period = spec.Table.PeriodDays(use); power = spec.Table.PSD_nT2_per_Hz(use);
[period,order] = sort(period); power = power(order);
h = loglog(ax,period,power,'k-','LineWidth',0.9);
set(ax,'XLim',cfg.DisplayPeriodDays,'XDir','reverse','YLim',ylimits, ...
    'FontSize',14,'TickDir','out','Box','on','XGrid','on','YGrid','on', ...
    'GridAlpha',0.13,'XMinorTick','on','YMinorTick','on');
ticks=[2 3 5 10 20 30 60 90]; set(ax,'XTick',ticks,'XTickLabel',string(ticks));
xlabel(ax,'周期（天）','FontName','Microsoft YaHei','FontSize',15);
ylabel(ax,'PSD_{|B|} (nT^2 Hz^{-1})','FontSize',15);
title(ax,['Voyager 1 | 2014-06-01 to 2015-06-01 | ' char(spec.Label)], ...
    'FontSize',16,'FontWeight','bold');
assert(isequal(h.XData(:),period) && isequal(h.YData(:),power));
base = fullfile(cfg.OutputRoot,['V1_PSD_' char(spec.ID) '_20140601_20150601_2day_90days']);
exportgraphics(fig,[base '.png'],'Resolution',200);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
savefig(fig,[base '.fig']);
files = string({[base '.png'];[base '.pdf'];[base '.fig']});
audit = struct('ID',spec.ID,'PlottedBins',numel(period),'XDirection',get(ax,'XDir'), ...
    'PeriodLimitsDays',cfg.DisplayPeriodDays,'PowerLimits',ylimits, ...
    'ActualPeriodBoundsDays',[min(period) max(period)],'PlottedValuesExact',true, ...
    'FootnotesAdded',false);
if ~cfg.Visible, close(fig); end
end


