function result = Correlate_V1_Sunspots_B
% 原始CDF和SILSO原始CSV -> 日/三日/月统计 -> 同窗口Pearson。
%% 路径与原始数据
CodeDir='C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/sunspots_B_correlation/';
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
addpath(CodeDir); if ~isfolder(OutputDir), mkdir(OutputDir); end
source=Run_V1_ThreeDay_Monthly_Overview(false,'both',fullfile(OutputDir,'source_recompute'),false);
a=readmatrix(SolarFile,'Delimiter',';','FileType','text');
t=datetime(a(:,1),a(:,2),a(:,3),'TimeZone','UTC')+hours(12);
use=t>=source.Method.StartUTC & t<source.Method.EndUTCExclusive;
a=a(use,:); t=t(use); sourceRows=find(use);
sn=a(:,5); sn(sn==-1)=NaN;
assert(numel(unique(t))==numel(t) && all(isnan(sn)|sn>=0));
[matched,idx]=ismember(t,source.Daily.EpochUTC); assert(all(matched));
dailySun=NaN(height(source.Daily),1); dailySun(idx)=sn;
%% 与已有图完全相同的独立平均窗口
result=struct;
result.SourceAudit=fullfile(OutputDir,'source_recompute','three_day_monthly_audit.mat');
result.SolarSource=SolarFile; result.SolarSHA256=Case1_File_SHA256(SolarFile);
result.SolarSourceRows=sourceRows; result.SolarOriginalColumns=a;
result.Method='Pearson of untransformed values, same-date/window pairwise finite. Zero retained. Independent finite-day arithmetic means per variable; no common-day restriction within a bin. No lag, detrending, interpolation, weighting or outlier rejection. Partial endpoint bins retained.';
result.VoyagerMethod=source.Method;
modes={'daily','three_day','monthly'};
Mode=string(modes).'; PearsonR=NaN(3,1); ValidPairs=zeros(3,1); ZeroSunspotPairs=zeros(3,1); FormulaDifference=NaN(3,1);
for m=1:3
    tag=modes{m};
    if m==1
        epoch=source.Daily.EpochUTC; solar=dailySun; b=source.Daily.BMean_nT;
        begin=dateshift(epoch,'start','day'); finish=begin+days(1);
        solarCount=double(isfinite(solar)); bCount=double(isfinite(b));
    else
        w=source.(tag).Windows; epoch=w.EpochUTC; begin=w.StartUTC; finish=w.EndUTCExclusive;
        b=w.BMean_nT; bCount=w.ValidDays_A;
        solar=NaN(height(w),1); solarCount=zeros(height(w),1); bins=zeros(numel(t),1);
        for k=1:height(w)
            rows=t>=begin(k) & t<finish(k); assert(all(bins(rows)==0)); bins(rows)=k;
            v=sn(rows & isfinite(sn)); solarCount(k)=numel(v);
            if ~isempty(v), solar(k)=mean(v); end
        end
        assert(all(bins>0)); good=isfinite(sn);
        check=accumarray(bins(good),sn(good),[height(w) 1],@mean,NaN);
        assert(isequaln(solar,check) && sum(solarCount)==nnz(good));
    end
    paired=isfinite(solar)&isfinite(b); x=solar(paired); y=b(paired);
    ValidPairs(m)=nnz(paired); ZeroSunspotPairs(m)=nnz(x==0);
    assert(numel(x)>2 && std(x)>0 && std(y)>0);
    r=corrcoef(x,y); PearsonR(m)=r(1,2);
    dx=x-mean(x); dy=y-mean(y);
    independent=sum(dx.*dy)/sqrt(sum(dx.^2)*sum(dy.^2));
    FormulaDifference(m)=abs(PearsonR(m)-independent); assert(FormulaDifference(m)<1e-12);
    pairs=table(epoch,begin,finish,solar,b,solarCount,bCount,paired, ...
        'VariableNames',{'EpochUTC','StartUTC','EndUTCExclusive','SunspotNumber','B_nT','SunspotValidDays','BValidDays','UsedInCorrelation'});
    result.(tag)=pairs;
    writetable(pairs,fullfile(OutputDir,[tag,'_pairs.csv']));
end
result.Summary=table(Mode,PearsonR,ValidPairs,ZeroSunspotPairs,FormulaDifference);
result.CreatedUTC=datetime('now','TimeZone','UTC');
result.CodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
writetable(result.Summary,fullfile(OutputDir,'Pearson_sunspots_B.csv'));
save(fullfile(OutputDir,'sunspots_B_correlation_audit.mat'),'result','-v7.3');
format long g; disp(result.Summary);
end
