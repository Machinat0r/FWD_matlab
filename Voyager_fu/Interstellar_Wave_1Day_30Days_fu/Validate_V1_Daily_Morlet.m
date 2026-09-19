function validation = Validate_V1_Daily_Morlet
% Compare exactly with the original IRFU Morlet, then test known periods/gaps.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(root),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
cfg.NumFrequencies=100;
t=(datetime(2000,1,1,'TimeZone','UTC'):days(1):datetime(2000,1,1,'TimeZone','UTC')+days(359)).';
x=0.45+0.02*sin(2*pi*(0:359).'/5);
daily=table(t,x,'VariableNames',{'EpochUTC','BMean_nT'});
[s,seg]=V1_Daily_Morlet(daily,cfg);
w=irf_wavelet([posixtime(t),x],'fs',1/86400,'nf',100, ...
    'f',[1/(30*86400),1/(2*86400)],'wavelet_width',5.36,'cutedge',1);
[~,order]=sort(w.f); expected=w.p{1}(:,order);
assert(isequaln(s.Power,expected),'Output differs from original irf_wavelet.');
[~,peak]=max(s.Power(180,:)); period=1/s.FrequencyHz(peak)/86400;
assert(abs(period/5-1)<0.06,'Known 5-day signal recovery failed.');
gapDaily=daily; gapDaily.BMean_nT(140:145)=NaN;
[g,gapSegments]=V1_Daily_Morlet(gapDaily,cfg);
assert(height(gapSegments)==2,'Gap was crossed.');
assert(all(isnan(g.Power(140:145,:)),'all'),'Missing day was filled.');
assert(all(isnan(g.Power(139,:))),'Native pre-gap edge cut absent.');
assert(all(isnan(g.Power(146,:))),'Native post-gap edge cut absent.');
assert(all(g.FrequencyHz<=g.NyquistHz*(1+1e-12)),'Above Nyquist computation.');
short=daily(1:3,:); [small,~]=V1_Daily_Morlet(short,cfg);
assert(all(isnan(small.Power),'all'),'Unsupported short segment was analyzed.');
% User-authorized interpolation: preserve observations and bridge the gap once.
[filled,interpolationAudit]=V1_Interpolate_Daily_Gaps(gapDaily);
assert(height(interpolationAudit)==6 && nnz(filled.IsInterpolated)==6);
assert(isequal(filled.BMean_nT(~filled.IsInterpolated),gapDaily.BMean_nT(~filled.IsInterpolated)));
weights=(1:6).'/7;
expectedGap=(1-weights)*gapDaily.BMean_nT(139)+weights*gapDaily.BMean_nT(146);
assert(max(abs(filled.BMean_nT(140:145)-expectedGap))<1e-12,'Linear bridge incorrect.');
assert(all(isnan(filled.BMeanObserved_nT(140:145))),'Original missing values not preserved.');
[filledSpec,filledSegments]=V1_Daily_Morlet(filled,cfg);
assert(height(filledSegments)==1,'Interpolated series was still split.');
assert(all(isfinite(filledSpec.Power(140:145,:)),'all'),'Internal gap still blank.');
filledReference=irf_wavelet([posixtime(t),filled.BMean_nT],'fs',1/86400,'nf',100, ...
    'f',[1/(30*86400),1/(2*86400)],'wavelet_width',5.36,'cutedge',1);
[~,filledOrder]=sort(filledReference.f);
assert(isequaln(filledSpec.Power,filledReference.p{1}(:,filledOrder)),'Full-series wavelet mismatch.');
validation=struct('Status','PASS','OriginalIRFUExactMatch',true, ...
    'KnownSignalPeriodDays',5,'RecoveredPeakPeriodDays',period, ...
    'NoGapFilling',true,'GapBoundaryBlank',true,'AboveNyquistExcluded',true, ...
    'ShortSegmentBlank',true,'FullSegmentSamples',seg.EvenSamplesUsed, ...
    'LinearGapInterpolation',true,'OriginalObservationsUnchanged',true, ...
    'FullSeriesAfterInterpolation',true,'NoInternalBlankAfterInterpolation',true);
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days';
save(fullfile(out,'morlet_synthetic_validation.mat'),'validation');
disp(validation);
end

