function validation=Validate_V1_Daily_Morlet_60Days
% Verify the newly included long-period band against direct IRFU and a known signal.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(root),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
cfg=struct('NumFrequencies',100,'ComputableFrequencyHz',[1/(60*86400),1/(2*86400)]);
n=(0:719).';
t=datetime(2000,1,1,'TimeZone','UTC')+days(n);
x=0.45+0.02*sin(2*pi*n/45);
daily=table(t,x,'VariableNames',{'EpochUTC','BMean_nT'});
[spec,~]=V1_Daily_Morlet(daily,cfg);
reference=irf_wavelet([posixtime(t),x],'fs',1/86400,'nf',100, ...
    'f',cfg.ComputableFrequencyHz,'wavelet_width',5.36,'cutedge',1);
[f,order]=sort(reference.f);
assert(isequaln(spec.Power,reference.p{1}(:,order)),'60-day band differs from direct IRFU.');
assert(isequal(spec.FrequencyHz,f(:)),'Frequency grid mismatch.');
assert(abs(1/spec.FrequencyHz(1)/86400-60)<1e-10,'Missing 60-day lower-frequency endpoint.');
[~,j]=max(mean(spec.Power(180:540,:),1));
peakDays=1/spec.FrequencyHz(j)/86400;
assert(abs(peakDays/45-1)<0.07,'Known 45-day signal outside the old band was not recovered.');
nativeScales=logspace(log10(0.5/86400/cfg.ComputableFrequencyHz(2)),log10(0.5/86400/cfg.ComputableFrequencyHz(1)),100);
edgeCount=floor(2*nativeScales(end)); % Preserve IRFU floating-point floor exactly.
assert(all(isnan(spec.Power(1:edgeCount,1))),'Native start edge cut missing.');
assert(all(isfinite(spec.Power(edgeCount+1:end-edgeCount-1,1))),'Unexpected interior blanks at 60 days.');
assert(all(isnan(spec.Power(end-edgeCount:end,1))),'60-day native end edge cut missing.');
% The pre-existing 30-day configuration must still be available.
oldCfg=struct('NumFrequencies',100);
[oldSpec,~]=V1_Daily_Morlet(daily,oldCfg);
oldReference=irf_wavelet([posixtime(t),x],'fs',1/86400,'nf',100, ...
    'f',[1/(30*86400),1/(2*86400)],'wavelet_width',5.36,'cutedge',1);
[~,oldOrder]=sort(oldReference.f);
assert(isequaln(oldSpec.Power,oldReference.p{1}(:,oldOrder)),'Previous 30-day behavior changed.');
validation=struct('Status','PASS','DirectIRFUExactMatch',true, ...
    'LowestFrequencyHz',spec.FrequencyHz(1),'LongestPeriodDays',60, ...
    'KnownSignalPeriodDays',45,'RecoveredPeakPeriodDays',peakDays, ...
    'Native60DayEdgesCorrect',true,'LowFrequencyStartEdgeSamples',edgeCount, ...
    'LowFrequencyEndEdgeSamples',edgeCount+1,'Previous30DayBehaviorPreserved',true);
disp(validation);
end

