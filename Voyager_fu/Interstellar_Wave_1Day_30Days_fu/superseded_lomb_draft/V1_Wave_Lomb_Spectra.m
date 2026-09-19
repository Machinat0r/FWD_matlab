function spec = V1_Wave_Lomb_Spectra(raw,cfg)
% Sliding Lomb--Scargle PSD; keep timestamps of all retained measurements.
% No call to fillmissing, interp1, resample, filtfilt, or uniform-grid FFT.
centers=(cfg.StartUTC+hours(12):days(cfg.StepDays):cfg.StopUTC-hours(12)).';
f=cfg.FrequencyHz(:); nt=numel(centers); nf=numel(f);
component=nan(nt,nf,3); scalar=nan(nt,nf); samplingWindow=nan(nt,nf);
vectorValid=all(isfinite(raw.B_RTN_nT),2); magValid=isfinite(raw.Bmag_nT);
nVector=zeros(nt,1); nMag=nVector; centralVector=nVector; centralMag=nVector;
spanVector=nan(nt,1); maxGap=spanVector; means=nan(nt,4);
status=repmat("outside_full_window",nt,1);
half=days(cfg.WindowDays/2);
starts=centers-half; stops=centers+half;
for k=1:nt
    dayStart=dateshift(centers(k),'start','day');
    dayMask=raw.EpochUTC>=dayStart & raw.EpochUTC<dayStart+days(1);
    centralVector(k)=nnz(dayMask & vectorValid); centralMag(k)=nnz(dayMask & magValid);
    if starts(k)<cfg.StartUTC || stops(k)>cfg.StopUTC, continue; end
    window=raw.EpochUTC>=starts(k) & raw.EpochUTC<stops(k);
    iv=window & vectorValid; im=window & magValid;
    nVector(k)=nnz(iv); nMag(k)=nnz(im);
    status(k)="no_central_day_vector";
    if nVector(k)>=4
        t=seconds(raw.EpochUTC(iv)-starts(k));
        spanVector(k)=(t(end)-t(1))/86400;
        maxGap(k)=max(diff(t))/3600;
        if centralVector(k)>0
            x=raw.B_RTN_nT(iv,:); means(k,1:3)=mean(x,1);
            px=columnPSD(x-means(k,1:3),t,f);
            component(k,:,:)=reshape(px,[1 nf 3]);
            % Dimensionless spectral window: diagnostic of actual observation times.
            samplingWindow(k,:)=abs(mean(exp(-2i*pi*t*f.'),1)).^2;
            status(k)="computed";
        end
    elseif centralVector(k)>0
        status(k)="fewer_than_4_vector_samples";
    end
    if nMag(k)>=4 && centralMag(k)>0
        t=seconds(raw.EpochUTC(im)-starts(k));
        x=raw.Bmag_nT(im); means(k,4)=mean(x);
        scalar(k,:)=columnPSD(x-means(k,4),t,f).';
    end
    if mod(k,365)==0, fprintf('Spectra: day %d / %d\n',k,nt); end
end
audit=table(centers,starts,stops,nVector,nMag,centralVector,centralMag, ...
    spanVector,maxGap,means,status, ...
    'VariableNames',{'CenterUTC','WindowStartUTC','WindowStopExclusiveUTC', ...
    'VectorSamples','MagnitudeSamples','CentralDayVectorSamples','CentralDayMagnitudeSamples', ...
    'VectorSampleSpanDays','LargestVectorSampleSeparationHours','SubtractedMeans_nT','VectorStatus'});
spec=struct('TimeUTC',centers,'FrequencyHz',f,'ComponentPSD',component, ...
    'TracePSD',sum(component,3),'MagnitudePSD',scalar, ...
    'SamplingWindow',samplingWindow,'WindowAudit',audit,'PSDUnits','nT^2/Hz');
end

function p=columnPSD(x,t,f)
p=zeros(numel(f),size(x,2));
active=any(x~=0,1);
if any(active), p(:,active)=plomb(x(:,active),t,f,'psd'); end
end
