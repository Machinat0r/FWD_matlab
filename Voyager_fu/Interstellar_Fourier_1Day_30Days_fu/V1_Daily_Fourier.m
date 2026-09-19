function [spec,windows,validation] = V1_Daily_Fourier(daily,cfg)
% Explicit complete daily windows; IRFU FFT PSD per window, no PSD averaging.
% Periodic Hann taper; remove each window's arithmetic mean before tapering.
% Correct only IRFU's doubling of DC / even-length Nyquist in one-sided PSD.
dt=86400; Fs=1/dt; L=cfg.WindowDays; hop=cfg.StepDays;
assert(all(seconds(diff(daily.EpochUTC))==dt),'Daily grid is not uniform.');
assert(all(isfinite(daily.BMean_nT)),'Fill authorized internal daily gaps first.');
assert(L==round(L) && L>=60 && hop==round(hop) && hop>=1,'Invalid FFT window.');
starts=(1:hop:height(daily)-L+1).';
assert(~isempty(starts),'No complete FFT window.');
last=starts+L-1;
t=posixtime(daily.EpochUTC);
taper=0.5-0.5*cos(2*pi*(0:L-1).'/L); % Periodic Hann; no toolbox dependency.
allPower=nan(numel(starts),floor(L/2)+1);
meanB=nan(size(starts)); nInterpolated=zeros(size(starts));
integral=meanB; weightedMeanSquare=meanB; maxDirectDifference=0;
for k=1:numel(starts)
    rows=starts(k):last(k);
    x=daily.BMean_nT(rows);
    [p,f]=irf_psd([t(rows),x],L,Fs,taper,0,'mean');
    p(1)=p(1)/2; % DC has no separate negative-frequency partner.
    if mod(L,2)==0, p(end)=p(end)/2; end
    allPower(k,:)=p(:).';
    meanB(k)=mean(x);
    nInterpolated(k)=nnz(daily.IsInterpolated(rows));
    % Independent FFT and Parseval checks catch unit and endpoint errors.
    y=(x-mean(x)).*taper;
    z=fft(y,L);
    direct=abs(z(1:floor(L/2)+1)).^2/(Fs*sum(taper.^2));
    if mod(L,2)==0
        direct(2:end-1)=2*direct(2:end-1);
    else
        direct(2:end)=2*direct(2:end);
    end
    maxDirectDifference=max(maxDirectDifference,max(abs(direct-p)));
    integral(k)=sum(p)*(Fs/L);
    weightedMeanSquare(k)=sum(y.^2)/sum(taper.^2);
end
tol=1e-11*max(1,max(allPower(:)));
assert(maxDirectDifference<tol,'IRFU and independent FFT PSD disagree.');
assert(max(abs(integral-weightedMeanSquare))<1e-12,'PSD normalization fails Parseval.');
assert(all(isfinite(allPower(:))) && all(allPower(:)>=0),'Invalid Fourier power.');
bounds=cfg.ComputableFrequencyHz;
selected=f>=bounds(1)*(1-1e-12) & f<=bounds(2)*(1+1e-12);
assert(nnz(selected)>1,'Insufficient Fourier bins.');
center=daily.EpochUTC(starts)+days((L-1)/2);
windowStart=daily.EpochUTC(starts)-hours(12);
windowStop=daily.EpochUTC(last)+hours(12);
spec=struct('TimeUTC',center,'FrequencyHz',f(selected),'Power',allPower(:,selected), ...
    'AllFrequencyHz',f,'AllPower',allPower,'BandMask',selected, ...
    'SamplingFrequencyHz',Fs,'NyquistHz',Fs/2,'FFTLength',L, ...
    'FrequencySpacingHz',Fs/L,'Taper',taper,'TaperEnergy',sum(taper.^2), ...
    'Units','nT^2/Hz','Estimator','Short-time Fourier one-sided power spectral density');
windows=table((1:numel(starts)).',windowStart,windowStop,center,starts,last, ...
    L-nInterpolated,nInterpolated,meanB,integral,weightedMeanSquare, ...
    abs(integral-weightedMeanSquare), ...
    'VariableNames',{'WindowIndex','StartUTC','EndUTCExclusive','CenterUTC', ...
    'FirstDailyRow','LastDailyRow','ObservedDayCount','InterpolatedDayCount', ...
    'SubtractedMeanB_nT','IntegratedPower_nT2','HannWeightedMeanSquare_nT2','ParsevalAbsError'});
validation=struct('InputSamples',height(daily),'CompleteWindows',numel(starts), ...
    'RetainedFrequencyBins',nnz(selected),'FrequencySpacingHz',Fs/L, ...
    'IndependentFFTMaxAbsDifference',maxDirectDifference, ...
    'ParsevalMaxAbsError',max(abs(integral-weightedMeanSquare)), ...
    'NoInternalMissingSpectrum',true,'NyquistHz',Fs/2);
end
