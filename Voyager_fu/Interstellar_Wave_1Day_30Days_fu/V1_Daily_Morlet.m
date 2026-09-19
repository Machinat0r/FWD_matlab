function [spec,segments] = V1_Daily_Morlet(daily,cfg)
% Direct irf_wavelet on continuous finite daily |B| segments.
% Same nf, Morlet width, normalization and default cutedge as MMS_fu/Wave.m.
% This helper never fills gaps. The formal entry first applies the user-authorized daily interpolation.
frequencyBounds=[1/(30*86400),1/(2*86400)]; Fs=1/86400;
if isfield(cfg,'ComputableFrequencyHz'), frequencyBounds=cfg.ComputableFrequencyHz; end
assert(frequencyBounds(1)>0 && frequencyBounds(2)<=Fs/2 && frequencyBounds(1)<frequencyBounds(2),'Invalid daily frequency band.');
a=logspace(log10(0.5*Fs/frequencyBounds(2)),log10(0.5*Fs/frequencyBounds(1)),cfg.NumFrequencies);
f=flipud((0.5*Fs./a).');
t=posixtime(daily.EpochUTC);
assert(all(abs(diff(t)-86400)<1e-6),'Daily time axis must have one row per UTC day.');
power=nan(height(daily),numel(f)); finite=isfinite(daily.BMean_nT);
edge=diff([false;finite;false]); starts=find(edge==1); stops=find(edge==-1)-1;
segments=table;
for k=1:numel(starts)
    rows=(starts(k):stops(k)).'; n=numel(rows);
    nUsed=n-mod(n,2); % irf_wavelet itself removes an odd final sample.
    % Exact native edge rule: first floor(Fs/f) samples and last floor(Fs/f)+1.
    nEdge=floor((1/86400)./f);
    allowed=nUsed>=2*nEdge+2;
    mapped=find(allowed); minFrequency=NaN; nFinite=0;
    if ~isempty(mapped)
        input=[t(rows),daily.BMean_nT(rows)];
        bounds=[f(mapped(1)),f(mapped(end))];
        if mapped(1)==1, bounds(1)=frequencyBounds(1); end
        if mapped(end)==numel(f), bounds(2)=frequencyBounds(2); end
        w=irf_wavelet(input,'fs',1/86400,'nf',numel(mapped), ...
            'f',bounds, ...
            'wavelet_width',5.36,'cutedge',1);
        [wf,order]=sort(w.f);
        assert(max(abs(wf(:)-f(mapped))./f(mapped))<1e-12,'Frequency grid mismatch.');
        assert(numel(w.t)==nUsed && isequal(w.t,t(rows(1:nUsed))),'Unexpected IRFU time handling.');
        power(rows(1:nUsed),mapped)=w.p{1}(:,order);
        minFrequency=f(mapped(1)); nFinite=nnz(isfinite(w.p{1}));
    end
    row=table(k,daily.EpochUTC(rows(1)),daily.EpochUTC(rows(end)),n,nUsed,n-nUsed, ...
        minFrequency,nFinite, ...
        'VariableNames',{'Segment','FirstUTC','LastUTC','OriginalDays','EvenSamplesUsed', ...
        'OddFinalSampleOmitted','LowestEvaluatedFrequencyHz','FiniteCoefficients'});
    segments=[segments;row]; %#ok<AGROW>
end
spec=struct('TimeUTC',daily.EpochUTC,'FrequencyHz',f,'Power',power, ...
    'RequestedFrequencyHz',[frequencyBounds(1),1/86400], ...
    'NyquistHz',1/(2*86400),'Units','nT^2/Hz','Estimator','irf_wavelet Morlet');
assert(all(all(isnan(power(~finite,:)))),'Missing day must remain missing.');
end



