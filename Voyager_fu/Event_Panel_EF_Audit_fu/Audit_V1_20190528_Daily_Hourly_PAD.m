function report = Audit_V1_20190528_Daily_Hourly_PAD
%Audit_V1_20190528_Daily_Hourly_PAD Trace the existing daily/hourly difference.
% Read-only scientific audit. No changes to production figures or processing.
% Verify the already-approved daily mean using its cached original records
% and the equivalent sample-count-weighted identity over all hourly bins.
addpath(genpath('C:/Users/Administrator/Documents/irfu-matlab-master'));
root = 'Z:/SPART-WORK/Data/Voyager/voyager1';
stem = 'V1_Case1-S03-L02_20190527_20190527_COHO1h_raw_LECP_P1_pitch_angle_predictedCK_';
dailyFile = fullfile(root,'lecp/1d/derived/pitch_angle/2013-2021/predicted_ck',[stem '1d_nativeCDF_Epoch.mat']);
hourlyFile = fullfile(root,'lecp/1h/derived/pitch_angle/2013-2021/predicted_ck',[stem '1h_nativeCDF_Epoch.mat']);
out = fullfile(root,'lecp/derived/daily_hourly_audit/2019/20190528');
if ~isfolder(out), mkdir(out); end
d = load(dailyFile); h = load(hourlyFile);
day = datetime(2019,5,28,'TimeZone','UTC');
daily = d.pitchAngleTable(dateshift(d.pitchAngleTable.EpochUTC,'start','day')==day,:);
hourly = h.pitchAngleTable(dateshift(h.pitchAngleTable.EpochUTC,'start','day')==day,:);
candidates = h.l1FallbackAudit.Candidates;
candidates = candidates(dateshift(candidates.BinStartUTC,'start','day')==day,:);
assert(height(daily)==1 && height(hourly)==2 && height(candidates)==5);
assert(d.opts.ContextDays==15 && h.opts.ContextDays==3);
records = daily.L1SourceRecords{1};
assert(height(records)==20);
rate = records.Rate_S1_to_S8;
rows = records.SourceCDFRecord;
sourceCDF = unique(records.SourceCDF);
assert(numel(sourceCDF)==1);
obj = dataobj(char(sourceCDF));
v = getv(obj,'FHDU_SectoredRates');
u = getv(obj,'FHDU_SectoredRateUncertainties');
q = getv(obj,'FHDU_SectoredQuality');
e = getv(obj,'Epoch');
dt = getv(obj,'DeltaT');
labels = getv(obj,'Hydrogen_Channels_Label');
assert(strtrim(string(labels.data(10,:)))=="P1");
original = struct('SourceCDF',sourceCDF,'SourceCDFRecord',rows, ...
    'EpochNative',e.data(rows),'Rate',v.data(rows,10,:), ...
    'RateSigma',u.data(rows,10,:),'Quality',q.data(rows,10,:), ...
    'DeltaT',dt.data(rows));
raw = table(rows,records.EpochUTC,original.DeltaT, ...
    'VariableNames',{'SourceCDFRecord','EpochUTC','DeltaT_s_raw'});
for s=1:8
    raw.(sprintf('Rate_S%d_raw',s))=double(original.Rate(:,1,s));
    raw.(sprintf('RateSigma_S%d_raw',s))=double(original.RateSigma(:,1,s));
    raw.(sprintf('Quality_S%d_raw',s))=original.Quality(:,1,s);
end
[found,k]=ismember(dateshift(records.EpochUTC,'start','hour'),candidates.BinStartUTC);
assert(all(found));
raw.HourCandidateApplied=candidates.Applied(k);
raw.HourCandidateDecision=candidates.Decision(k);
raw.ContributesToDailyS1=records.Contributes_S1_to_S8(:,1);
raw.AppearsInVisibleHourlyPAD=ismember(dateshift(records.EpochUTC,'start','hour'), ...
    dateshift(hourly.EpochUTC(hourly.PADUsable),'start','hour'));
sector=(1:8)'; n=zeros(8,1); sumRate=n; meanRaw=n; meanHours=n; savedMean=n; savedFlux=n; pa=n;
for s=1:8
    valid=isfinite(rate(:,s)) & rate(:,s)>=0;
    n(s)=nnz(valid);sumRate(s)=sum(rate(valid,s));meanRaw(s)=sumRate(s)/n(s);
    use=candidates.SectorSampleCount(:,s)>0;
    meanHours(s)=sum(candidates.MeanRate(use,s).*candidates.SectorSampleCount(use,s))/sum(candidates.SectorSampleCount(use,s));
    savedMean(s)=daily.(sprintf('RawRate_S%d_1d',s));
    savedFlux(s)=daily.(sprintf('Flux_S%d_1d',s));
    pa(s)=daily.(sprintf('PA_S%d_deg',s));
end
assert(max(abs(meanRaw-savedMean))<1e-10);
assert(max(abs(meanHours-savedMean))<1e-10);
comparison=table(sector,n,sumRate,meanRaw,meanHours,savedMean,savedFlux,pa, ...
    'VariableNames',{'Sector','OriginalSampleCount','SumOriginalRate','MeanFromOriginalRecords','MeanFromAllHourlyBinsWeightedBySampleCount','SavedDailyMeanRate','SavedDailyFlux','SavedDailyPA_deg'});
spikes=ismember(rows,[3364;3368]);
assert(nnz(spikes)==2 && all(raw.Rate_S1_raw(spikes)==[350;438]));
report=struct('DateUTC','2019-05-28','DailyAuditFile',dailyFile,'HourlyAuditFile',hourlyFile, ...
    'SourceCDF',char(sourceCDF),'OriginalRecordCount',height(records), ...
    'S1DailySampleCount',n(1),'S1DailyMeanRate',savedMean(1),'S1DailyFlux',savedFlux(1), ...
    'S1DailyPA_deg',pa(1),'HourlyCandidateCount',height(candidates), ...
    'VisibleHourlyPADCount',nnz(hourly.PADUsable),'LargeS1RecordNumbers',rows(spikes)', ...
    'LargeS1Rates',raw.Rate_S1_raw(spikes)','LargeS1Quality',raw.Quality_S1_raw(spikes)', ...
    'TwoLargeRecordsFractionOfS1RateSum',sum(raw.Rate_S1_raw(spikes))/sumRate(1), ...
    'RawDailyIdentityMaxAbsError',max(abs(meanRaw-savedMean)), ...
    'AllHourWeightedIdentityMaxAbsError',max(abs(meanHours-savedMean)), ...
    'DailyContextDays',d.opts.ContextDays,'HourlyContextDays',h.opts.ContextDays, ...
    'BackgroundMode',d.opts.LECPBackgroundMode, ...
    'MissingDataPolicy','Preserve raw fill values in original table; use exact cached production contributor masks for verification. No new filtering.');
report.DailyColorLimits = existingColorLimits(d.pitchAngleTable,'1d',d.opts.ColorPercentiles);
report.HourlyColorLimits = existingColorLimits(h.pitchAngleTable,'1h',h.opts.ColorPercentiles);
writeCSV(raw,fullfile(out,'original_L1_records_20190528.csv'));
writeCSV(removevars(daily,'L1SourceRecords'),fullfile(out,'daily_actual_PAD_20190528.csv'));
writeCSV(removevars(hourly,'L1SourceRecords'),fullfile(out,'hourly_actual_PAD_20190528.csv'));
writeCSV(removevars(candidates,{'L1Rows','L2Rows','SourceRecords'}),fullfile(out,'hourly_all_candidates_20190528.csv'));
writeCSV(comparison,fullfile(out,'daily_mean_identity_by_sector.csv'));
save(fullfile(out,'daily_hourly_source_audit_20190528.mat'),'report','daily','hourly','candidates','records','original','comparison','dailyFile','hourlyFile');
fid=fopen(fullfile(out,'audit_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
disp(report);disp(comparison);
end
function writeCSV(t,file)
names=t.Properties.VariableNames;out=table;
for i=1:numel(names)
 v=t.(names{i});
 if isdatetime(v),v.Format='yyyy-MM-dd''T''HH:mm:ss.SSS''Z''';v=string(v);end
 if (isnumeric(v)||islogical(v)) && size(v,2)>1
  for s=1:size(v,2),out.(sprintf('%s_S%d',names{i},s))=v(:,s);end
 else,out.(names{i})=v;end
end
writetable(out,file,'Encoding','UTF-8');
end
function limits = existingColorLimits(t,suffix,pct)
% Reproduce the existing plot color limits for audit only.
names=arrayfun(@(s)sprintf('Flux_S%d_%s',s,suffix),1:7,'UniformOutput',false);
paNames=arrayfun(@(s)sprintf('PA_S%d_deg',s),1:7,'UniformOutput',false);
j=t{:,names};pa=t{:,paNames};
use=repmat(t.PADUsable,1,7) & isfinite(j) & j>0 & isfinite(pa);
x=sort(log10(j(use)));n=numel(x);
lo=x(max(1,round(pct(1)/100*n)));hi=x(min(n,max(1,round(pct(2)/100*n))));
lo=floor(lo*4)/4;hi=ceil(hi*4)/4;
if hi-lo<1,m=(lo+hi)/2;lo=m-.5;hi=m+.5;end
limits=[lo hi];
end
