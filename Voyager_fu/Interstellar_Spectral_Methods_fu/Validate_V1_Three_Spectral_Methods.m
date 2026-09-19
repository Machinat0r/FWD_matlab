function validation=Validate_V1_Three_Spectral_Methods(daily,spectra,details,cfg)
% Check physical units, reference equivalence, edge averaging and peak mapping.
x=daily.BMean_nT; N=numel(x); Fs=cfg.SamplingFrequencyHz;
[a,b]=periodogram(x,details.FFT.Window,N,Fs,'onesided');
fftTable=spectra(1).Table;
validation.FFTPeriodogramRelativeError=norm(a-fftTable.PSD_nT2_per_Hz)/norm(a);
assert(validation.FFTPeriodogramRelativeError<1e-12);
assert(max(abs(b-fftTable.FrequencyHz))<Fs*1e-12);
validation.FFTParsevalAbsError=abs(details.FFT.IntegratedPower_nT2-details.FFT.HammingWeightedMeanSquare_nT2);
assert(validation.FFTParsevalAbsError<1e-12);
mt=details.PMTM;
[a,b]=pmtm(x,mt.Tapers,mt.Eigenvalues,mt.NFFT,Fs,'adapt','DropLastTaper',false);
mtTable=spectra(3).Table;
validation.PMTMExplicitTapersRelativeError=norm(a-mtTable.PSD_nT2_per_Hz)/norm(a);
assert(validation.PMTMExplicitTapersRelativeError<1e-12);
assert(isequal(b,mtTable.FrequencyHz) && mt.TaperCount==6 && mt.NFFT==max(256,2^nextpow2(N)));
q=details.Wavelet.p{1}; counts=sum(isfinite(q),1).';
s=nan(size(counts));
supported=counts>0;
for j=find(supported).', values=q(isfinite(q(:,j)),j); s(j)=sum(values)/numel(values); end
wavePSD=spectra(2).Table.PSD_nT2_per_Hz;
assert(isequal(isnan(wavePSD),~supported),'Unsupported scales must remain NaN.');
validation.WaveletFiniteMeanRelativeError=norm(s(supported)-wavePSD(supported))/norm(s(supported));
validation.WaveletUnsupportedScales=nnz(~supported);
audit=details.WaveletAveragingAudit;
assert(isequal(audit.ValidTimeCount,counts));
assert(all(isnat(audit.FirstValidUTC(~supported))) && all(isnat(audit.LastValidUTC(~supported))));
assert(validation.WaveletFiniteMeanRelativeError<1e-12);
scale=(Fs/2)./details.Wavelet.f(:);
% Reconstruct the native logspace scales, preserving floating floor behavior.
nativeScale=logspace(log10(0.5*Fs/cfg.DisplayFrequencyHz(2)), ...
    log10(0.5*Fs/cfg.DisplayFrequencyHz(1)),cfg.WaveletFrequencyCount).';
assert(max(abs(scale-nativeScale))<1e-10);
Nwave=N-mod(N,2);
expectedCounts=max(0,Nwave-2*floor(2*nativeScale)-1);
assert(isequal(counts,expectedCounts));
validation.WaveletEdgeCountsExact=true;
assert(numel(details.Wavelet.t)==Nwave && details.WaveletInputAudit.UsedSamples==Nwave);
assert(isequal(details.Wavelet.t(:),posixtime(daily.EpochUTC(1:Nwave))));
validation.WaveletUsedSamples=Nwave;
validation.WaveletOddFinalSampleOmitted=mod(N,2)==1;
validation.WaveletValidTimeCountRange=[min(counts) max(counts)];
validation.PMTMTaperCount=mt.TaperCount;
validation.PMTMNFFT=mt.NFFT;

%% A known 20-day sinusoid checks time-unit and period conversion
t=(0:N-1)'/Fs; testPeriodDays=20;
signal=0.02*cos(2*pi*t/(testPeriodDays*86400));
[test,~]=V1_Three_Spectral_Methods(t,signal,cfg);
peakPeriods=zeros(1,3);
for j=1:3
    tab=test(j).Table; use=tab.InDisplayedBand & isfinite(tab.PSD_nT2_per_Hz);
    ff=tab.FrequencyHz(use); pp=tab.PSD_nT2_per_Hz(use);
    [~,idx]=max(pp); peakPeriods(j)=1/ff(idx)/86400;
end
assert(abs(1/(peakPeriods(1)*86400)-1/(testPeriodDays*86400))<=Fs/(2*N)+Fs*1e-12);
assert(abs(peakPeriods(2)-testPeriodDays)<1);
assert(abs(1/(peakPeriods(3)*86400)-1/(testPeriodDays*86400))<=mt.HalfBandwidthHz+mt.FrequencySpacingHz);
validation.Synthetic20DayPeakPeriods=peakPeriods;
validation.AllChecksPassed=true;
end


