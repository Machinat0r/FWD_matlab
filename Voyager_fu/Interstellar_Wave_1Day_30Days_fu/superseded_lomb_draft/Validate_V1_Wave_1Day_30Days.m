function validation = Validate_V1_Wave_1Day_30Days
% Scientific numerical checks with explicitly synthetic signals only.
cfg.StartUTC=datetime(2000,1,1,'TimeZone','UTC');
cfg.StopUTC=cfg.StartUTC+days(91); cfg.WindowDays=90; cfg.StepDays=45;
n=90*24; df=1/(n*3600); cfg.FrequencyHz=(1:n/2-1).'*df;
t=(cfg.StartUTC:hours(1):cfg.StopUTC-hours(1)).';
secondsFromStart=seconds(t-cfg.StartUTC);
amplitude=[0.02 0.01 0.015]; periods=[5 15 30]*86400;
x=sin(2*pi*secondsFromStart./periods).*amplitude;
raw=table(t,x+[0.1 0.3 -0.2],0.45+x(:,1), ...
    'VariableNames',{'EpochUTC','B_RTN_nT','Bmag_nT'});
s=V1_Wave_Lomb_Spectra(raw,cfg); center=2;
expected=sum(amplitude.^2)/2; recovered=sum(s.TracePSD(center,:))*df;
assert(abs(recovered/expected-1)<1e-9,'PSD normalization failed.');
[~,peak]=max(squeeze(s.ComponentPSD(center,:,:)),[],1);
recoveredDays=1./s.FrequencyHz(peak).'/86400;
assert(max(abs(recoveredDays-periods/86400))<1e-9,'Synthetic frequency recovery failed.');
assert(isequaln(s.TracePSD,sum(s.ComponentPSD,3)),'Trace identity failed.');

% Missing data must agree with direct plomb at retained timestamps.
cfg.FrequencyHz=logspace(log10(1/(30*86400)),log10(1/86400),100).';
gapRaw=raw; missing=mod((0:height(raw)-1).',24)>=8;
gapRaw.B_RTN_nT(missing,:)=NaN; gapRaw.Bmag_nT(missing)=NaN;
g=V1_Wave_Lomb_Spectra(gapRaw,cfg);
start=g.WindowAudit.WindowStartUTC(center); stop=g.WindowAudit.WindowStopExclusiveUTC(center);
good=gapRaw.EpochUTC>=start & gapRaw.EpochUTC<stop & all(isfinite(gapRaw.B_RTN_nT),2);
gx=gapRaw.B_RTN_nT(good,:); gx=gx-mean(gx,1);
direct=plomb(gx,seconds(gapRaw.EpochUTC(good)-start),cfg.FrequencyHz,'psd');
actual=squeeze(g.ComponentPSD(center,:,:));
assert(max(abs(actual-direct),[],'all')<1e-9,'Gap handling changed retained measurements.');

% An unobserved central day stays blank despite nearby observations.
dayStart=dateshift(g.TimeUTC(center),'start','day');
emptyDay=gapRaw.EpochUTC>=dayStart & gapRaw.EpochUTC<dayStart+days(1);
gapRaw.B_RTN_nT(emptyDay,:)=NaN; gapRaw.Bmag_nT(emptyDay)=NaN;
blank=V1_Wave_Lomb_Spectra(gapRaw,cfg);
assert(all(isnan(blank.TracePSD(center,:))),'Missing day was filled.');
assert(all(isnan(blank.MagnitudePSD(center,:))),'Missing magnitude day was filled.');
assert(all(isnan(blank.TracePSD(1,:))),'Incomplete boundary window was computed.');
validation=struct('Status','PASS','Signals','Synthetic only', ...
    'ExpectedVariance_nT2',expected,'RecoveredVariance_nT2',recovered, ...
    'RecoveredPeriodsDays',recoveredDays,'MissingSamples',nnz(missing), ...
    'PSDNormalization',true,'GapExclusion',true,'NoDataDayBlank',true,'BoundaryBlank',true);
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days';
if ~isfolder(out), mkdir(out); end
save(fullfile(out,'synthetic_validation.mat'),'validation');
disp(validation);
end
