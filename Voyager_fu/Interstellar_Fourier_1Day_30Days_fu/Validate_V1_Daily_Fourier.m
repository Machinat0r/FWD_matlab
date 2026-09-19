function report=Validate_V1_Daily_Fourier
% Synthetic checks only; no Voyager data processing.
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Case1_PPT_VerticalLine_Events_7d_fu');
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
cfg=struct('WindowDays',90,'StepDays',1,'ComputableFrequencyHz',[1/(30*86400),1/(2*86400)]);
n=(0:299).';
t=datetime(2012,8,25,12,0,0,'TimeZone','UTC')+days(n);
A=0.02;
daily=table(t,0.5+A*sin(2*pi*n/10),false(size(n)), ...
    'VariableNames',{'EpochUTC','BMean_nT','IsInterpolated'});
[spec,windows,validation]=V1_Daily_Fourier(daily,cfg);
[~,peak]=max(mean(spec.Power,1));
peakDays=1/spec.FrequencyHz(peak)/86400;
assert(abs(peakDays-10)<1e-10,'Known 10-day signal not recovered.');
assert(max(abs(windows.IntegratedPower_nT2-A^2/2))<1e-12,'Sinusoid amplitude/power normalization failed.');
assert(windows.StartUTC(1)==datetime(2012,8,25,'TimeZone','UTC'));
assert(windows.CenterUTC(1)==windows.StartUTC(1)+days(45));
assert(all(diff(windows.CenterUTC)==days(1)));
daily.BMean_nT=0.5+A*(-1).^n;
[~,nyqWindows]=V1_Daily_Fourier(daily,cfg);
assert(max(abs(nyqWindows.IntegratedPower_nT2-A^2))<1e-12,'Nyquist PSD endpoint is doubled or halved incorrectly.');
daily.BMean_nT(:)=0.5;
[dc,~]=V1_Daily_Fourier(daily,cfg);
assert(max(dc.AllPower(:))<1e-20,'Constant background leaked through window mean removal.');
report=struct('TenDayPeakDays',peakDays,'SinusoidAmplitudePowerCorrect',true, ...
    'NyquistPowerCorrect',true,'ConstantBackgroundRemoved',true, ...
    'WindowTimesCorrect',true,'BaseValidation',validation);
disp(report);
end
