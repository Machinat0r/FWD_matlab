function audit=Diagnose_V1_23Day_WindowChoice
% Exploratory interval sensitivity at the user's approximate 23-day feature.
% Raw CDF input and the existing daily interpolation/Hamming FFT are retained.
% Candidate intervals are fixed before inspecting their PSD; no window search.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(root),'Interstellar_Wave_1Day_30Days_fu'));
addpath(fullfile(fileparts(root),'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
cfg.DataRoot='Z:/SPART-WORK/Data/Voyager';
cfg.StartUTC=datetime(2014,1,1,'TimeZone','UTC');
cfg.StopUTC=datetime(2016,1,1,'TimeZone','UTC');
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_23Day_Window_Diagnosis';
if ~isfolder(out),mkdir(out);end
diary(fullfile(out,'run_diagnosis.log')); cleanup=onCleanup(@()diary('off'));
[raw,sources,duplicates]=V1_Wave_Read_Hourly(cfg);
day=(cfg.StartUTC:days(1):cfg.StopUTC-days(1)).';
[matched,bin]=ismember(dateshift(raw.EpochUTC,'start','day'),day);
good=matched & isfinite(raw.Bmag_nT);
meanB=accumarray(bin(good),raw.Bmag_nT(good),[numel(day) 1],@mean,NaN);
count=accumarray(bin(good),1,[numel(day) 1],@sum,0);
observed=table(day+hours(12),meanB,count,'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
[daily,gapAudit]=V1_Interpolate_Daily_Gaps(observed);
ref='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_PSD_ThreeMethods_2014_2015/V1_PSD_ThreeMethods_2014_2015_analysis.mat';
old=load(ref,'result'); assert(isequaln(daily,old.result.Daily)); clear old
names=["Previous_two_years";"Previous_one_year";"Current_box";"Extended_candidate";"Centered_candidate"];
datePairs=["2014-01-01","2015-12-31";"2014-06-01","2015-06-01"; ...
    "2014-11-01","2015-05-15";"2014-10-01","2015-06-30";"2014-12-01","2015-05-31"];
starts=datetime(datePairs(:,1),'TimeZone','UTC'); ends=datetime(datePairs(:,2),'TimeZone','UTC');
N=zeros(5,1); peakPeriod=N; peakPSD=N; resolution23=N; fills=N; parseval=N;
psdTables=cell(5,1); allPeriods=cell(5,1); allPowers=cell(5,1);
Fs=1/86400;
for k=1:5
    use=daily.EpochUTC>=starts(k) & daily.EpochUTC<ends(k)+days(1);
    x=daily.BMean_nT(use); N(k)=numel(x); fills(k)=nnz(daily.IsInterpolated(use));
    w=hamming(N(k)); Y=fft(x.*w);
    f=(0:floor(N(k)/2))'*Fs/N(k);
    P=abs(Y(1:numel(f))).^2/(Fs*sum(w.^2));
    if mod(N(k),2)==0, P(2:end-1)=2*P(2:end-1); else,P(2:end)=2*P(2:end);end
    [pcheck,fcheck]=periodogram(x,w,N(k),Fs,'onesided');
    assert(norm(P-pcheck)/norm(P)<1e-12 && max(abs(f-fcheck))<Fs*1e-12);
    parseval(k)=abs(sum(P)*Fs/N(k)-sum((x.*w).^2)/sum(w.^2));
    assert(parseval(k)<1e-12);
    T=nan(size(f)); T(f>0)=1./f(f>0)/86400;
    near=T>=20 & T<=26; idx=find(near);
    [peakPSD(k),j]=max(P(near)); peakPeriod(k)=T(idx(j));
    resolution23(k)=23^2/N(k); % Native grid spacing mapped locally to period; not FWHM.
    tab=table(f,T,P,'VariableNames',{'FrequencyHz','PeriodDays','PSD_nT2_per_Hz'});
    psdTables{k}=tab;
    writetable(tab,fullfile(out,names(k)+"_FFT.csv"));
    allPeriods{k}=T; allPowers{k}=P;
end
summary=table(names,starts,ends,N,fills,peakPeriod,peakPSD,resolution23, ...
 'VariableNames',{'Window','StartUTC','EndUTCInclusive','DailyCount','FilledDays', ...
 'StrongestNativePeriod_20to26Days','PSDAtThatBin_nT2_per_Hz','ApproxNativePeriodStepAt23Days'});
disp(summary);

% Inspect existing native 23-day wavelet coefficients only as result audit.
waveFile='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_90Days/V1_daily_morlet_analysis.mat';
wv=load(waveFile,'result'); q=wv.result.Spectrum;
[~,j]=min(abs(1./q.FrequencyHz/86400-23)); nativePeriod=1/q.FrequencyHz(j)/86400;
use=q.TimeUTC>=cfg.StartUTC & q.TimeUTC<cfg.StopUTC;
tt=q.TimeUTC(use); pp=q.Power(use,j);
peakRows=find(islocalmax(pp)); [~,order]=sort(pp(peakRows),'descend'); peakRows=peakRows(order);
wavePeaks=table(tt(peakRows),pp(peakRows),'VariableNames',{'TimeUTC','PSD_nT2_per_Hz'});
disp(wavePeaks(1:min(8,height(wavePeaks)),:));
writetable(wavePeaks,fullfile(out,'Existing_native_23day_wavelet_local_peaks.csv'));
writetable(table(tt,pp,'VariableNames',{'TimeUTC','NativeWaveletPSD_nT2_per_Hz'}),fullfile(out,'Existing_native_23day_wavelet_trace.csv'));

% The same plot scale avoids the previous PMTM-driven visual compression.
fig=figure('Visible','off','Color','w','Position',[70 70 1100 680]);
ax=axes(fig,'Position',[.12 .14 .84 .77]); hold(ax,'on');
show=[3 4 5]; styles={'k-','b-','r-'};
for k=1:3
    z=show(k); use=allPeriods{z}>=15 & allPeriods{z}<=40;
    plot(ax,allPeriods{z}(use),allPowers{z}(use),styles{k},'LineWidth',1.2);
end
set(ax,'XScale','log','YScale','log','XDir','reverse','XLim',[15 40], ...
    'YLim',[1 1e4],'XTick',[15 18 20 23 26 30 40],'XTickLabel',string([15 18 20 23 26 30 40]), ...
    'FontSize',13,'TickDir','out','Box','on','XGrid','on','YGrid','on','GridAlpha',.15);
xlabel(ax,'周期（天）','FontName','Microsoft YaHei','FontSize',15);
ylabel(ax,'PSD_{|B|} (nT^2 Hz^{-1})','FontSize',15);
title(ax,'Voyager 1 | FFT (Hamming)','FontSize',16);
labels=datePairs(show,1)+" to "+datePairs(show,2);
legend(ax,labels,'Location','southwest','FontSize',11,'Box','off');
base=fullfile(out,'V1_FFT_Window_Comparison_15day_40days');
exportgraphics(fig,[base '.png'],'Resolution',200);
exportgraphics(fig,[base '.pdf'],'ContentType','vector'); savefig(fig,[base '.fig']); close(fig);
writetable(summary,fullfile(out,'window_comparison.csv'));
writetable(sources(:,{'SourceFile','SHA256','RecordsInRange'}),fullfile(out,'source_manifest.csv'));
method=struct('UserRequest','Explain why the approximate 23-day feature becomes less clear, and suggest a time interval.', ...
    'Input','Original CDF -> established daily scalar means -> prior authorized linear interpolation; exact equality to previous daily input checked.', ...
    'Candidates','Three earlier intervals and two fixed candidates: 2014-10-01 through 2015-06-30, 2014-12-01 through 2015-05-31. Inclusive UTC days. No exhaustive window search.', ...
    'FFT','Same symmetric Hamming, no demeaning/detrending, no padding. Correct even/odd one-sided normalization. Methods checked against periodogram and Parseval.', ...
    'PeakSummary','Largest native FFT value within 20--26 days for descriptive comparison of the user approximate 23-day feature; no peak significance or claim of optimal interval. Native grids, record lengths and time weights differ. Peak PSD heights alone are not a significance comparison.', ...
    'WaveletAudit','Read existing full-record wavelet coefficients at nearest native period to 23 days; no recalculation from MAT daily input, no smoothing or time averaging, local maxima listed in descending power.', ...
    'Plot','Comparison zoom 15--40 days with fixed y limits 1--1e4, plotting only existing native bins. Display-only zoom does not filter the signal. No method footnotes.');
audit=struct('Config',cfg,'Method',method,'Summary',summary,'Spectra',{psdTables}, ...
    'Daily',daily,'ObservedDaily',observed,'InterpolationAudit',gapAudit,'Sources',sources, ...
    'Duplicates',duplicates,'ParsevalErrors',parseval,'NativeWaveletPeriodDays',nativePeriod, ...
    'ExistingWaveletLocalPeaks',wavePeaks,'ExistingWaveletFile',waveFile, ...
    'ExistingWaveletSHA256',Case1_File_SHA256(waveFile),'PriorDailyReference',ref, ...
    'PriorDailyReferenceSHA256',Case1_File_SHA256(ref), ...
    'CodeFile',[mfilename('fullpath') '.m'],'CodeSHA256',Case1_File_SHA256([mfilename('fullpath') '.m']), ...
    'MATLABVersion',version);
save(fullfile(out,'window_diagnosis_audit.mat'),'audit','raw','-v7.3');
fid=fopen(fullfile(out,'processing_audit.txt'),'w','n','UTF-8');
fields=fieldnames(method); for k=1:numel(fields),fprintf(fid,'%s:\n%s\n\n',fields{k},method.(fields{k}));end
fprintf(fid,'Existing native wavelet period: %.9f days\n',nativePeriod);fclose(fid);
end
