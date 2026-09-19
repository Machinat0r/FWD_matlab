function [spectra,details] = V1_Three_Spectral_Methods(epochSeconds,x,cfg)
% Computational cores adapted from the three supplied turbulence scripts.
% Input nT, uniform daily time in seconds; PSD output nT^2/Hz.
x=x(:); epochSeconds=epochSeconds(:); N=numel(x); Fs=cfg.SamplingFrequencyHz;
assert(N==numel(epochSeconds) && N>=2);
assert(all(isfinite(x)) && all(abs(diff(epochSeconds)-1/Fs)<1e-6));
assert(~cfg.RemoveMean && ~cfg.RemoveTrend,'The supplied methods do not explicitly demean/detrend.');

%% Hamming FFT: full-record counterpart of the source event-frame FFT
w=hamming(N);
Y=fft(x.*w,N);
f=(0:floor(N/2))'*Fs/N;
P=abs(Y(1:numel(f))).^2/(Fs*sum(w.^2));
if mod(N,2)==0, P(2:end-1)=2*P(2:end-1); else, P(2:end)=2*P(2:end); end
spectra(1)=makeSpectrum("FFT","FFT (Hamming)",f,P,cfg);
details.FFT=struct('Window',w,'FFT',Y,'NFFT',N,'TaperEnergy',sum(w.^2), ...
    'SubtractedMean_nT',0,'InputMean_nT',mean(x),'FrequencySpacingHz',Fs/N, ...
    'IntegratedPower_nT2',sum(P)*(Fs/N), ...
    'HammingWeightedMeanSquare_nT2',sum((x.*w).^2)/sum(w.^2));

%% Morlet wavelet, then time average as nanmean(pB,1) in the source
wave=irf_wavelet([epochSeconds,x],'fs',Fs,'f',cfg.DisplayFrequencyHz, ...
    'nf',cfg.WaveletFrequencyCount,'wavelet_width',cfg.WaveletWidth, ...
    'returnpower',1,'cutedge',1);
% Preserve the supplied IRFU rule: only an odd final input sample is omitted.
Nwave=N-mod(N,2);
assert(isequal(wave.t(:),epochSeconds(1:Nwave)),'Unexpected wavelet sample selection.');
details.WaveletInputAudit=struct('RequestedSamples',N,'UsedSamples',Nwave, ...
    'OddFinalSampleOmitted',mod(N,2)==1, ...
    'OmittedEpochUTC',datetime(epochSeconds(Nwave+1:end),'ConvertFrom','posixtime','TimeZone','UTC'), ...
    'OmittedB_nT',x(Nwave+1:end),'FirstUsedEpochUTC', ...
    datetime(epochSeconds(1),'ConvertFrom','posixtime','TimeZone','UTC'), ...
    'LastUsedEpochUTC',datetime(epochSeconds(Nwave),'ConvertFrom','posixtime','TimeZone','UTC'));
power=wave.p{1};
P=mean(power,1,'omitnan').';
validCount=sum(isfinite(power),1).';
% A short record can have no interior times at long periods; retain NaN.
assert(isequal(isnan(P),validCount==0),'Unexpected missing wavelet power.');
spectra(2)=makeSpectrum("Wavelet","Morlet wavelet",wave.f,P,cfg);
firstUTC=NaT(size(validCount),'TimeZone','UTC'); lastUTC=firstUTC;
for j=find(validCount>0).'
    rows=find(isfinite(power(:,j)));
    firstUTC(j)=datetime(epochSeconds(rows(1)),'ConvertFrom','posixtime','TimeZone','UTC');
    lastUTC(j)=datetime(epochSeconds(rows(end)),'ConvertFrom','posixtime','TimeZone','UTC');
end
details.Wavelet=wave;
details.WaveletAveragingAudit=table(wave.f(:),1./wave.f(:)/86400,validCount, ...
    firstUTC,lastUTC, ...
    'VariableNames',{'FrequencyHz','PeriodDays','ValidTimeCount','FirstValidUTC','LastValidUTC'});

%% PMTM: preserve the supplied parameters and MATLAB defaults exactly
[P,f]=pmtm(x,cfg.NW,[],Fs);
spectra(3)=makeSpectrum("PMTM","PMTM",f,P,cfg);
[eigenvectors,eigenvalues]=dpss(N,cfg.NW);
details.PMTM=struct('NW',cfg.NW,'NFFT',max(256,2^nextpow2(N)), ...
    'TaperCount',size(eigenvectors,2)-1,'Tapers',eigenvectors(:,1:end-1), ...
    'Eigenvalues',eigenvalues(1:end-1),'Weighting','adapt','DropLastTaper',true, ...
    'HalfBandwidthHz',cfg.NW*Fs/N,'FrequencySpacingHz',Fs/max(256,2^nextpow2(N)));
end

function out=makeSpectrum(id,label,f,p,cfg)
f=f(:); p=p(:); period=nan(size(f)); period(f>0)=1./f(f>0)/86400;
band=f>=cfg.DisplayFrequencyHz(1)*(1-1e-12) & f<=cfg.DisplayFrequencyHz(2)*(1+1e-12);
if id=="Wavelet"
    assert(all(isfinite(p) | isnan(p)) && all(p(isfinite(p))>=0),'Invalid wavelet power.');
else
    assert(all(isfinite(p)) && all(p>=0),'Invalid spectral power.');
end
out=struct('ID',id,'Label',label,'Table',table(f,period,p,band, ...
    'VariableNames',{'FrequencyHz','PeriodDays','PSD_nT2_per_Hz','InDisplayedBand'}));
end

