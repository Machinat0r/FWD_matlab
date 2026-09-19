function source=Voyager_Supplement_V2_MAG(source)
% Real MAG 48-s CDF -> UTC hourly means -> daily means; only missing COHO days.
root='Z:/SPART-WORK/Data/Voyager/voyager2/mag/48s/reviewed_vim/';
e=dir(fullfile(root,'**','voyager2_48s_mag-vim_*_v*.cdf'));
files=sort(string(fullfile({e.folder},{e.name})).');
keys=regexp(cellstr(files),'_(\d{8})_v','tokens','once'); keys=string(cellfun(@(x)x{1},keys,'UniformOutput',false));
[~,keep]=unique(keys,'last'); files=files(keep);
d=source.Daily; begin=dateshift(d.EpochUTC(1),'start','day'); finish=dateshift(d.EpochUTC(end),'start','day')+days(1);
records=table; manifest=table;
for k=1:numel(files)
    token=regexp(char(files(k)),'_(\d{4})\d{4}_v','tokens','once'); yy=str2double(token{1});
    if yy<year(begin)||yy>year(finish), continue; end
    q=Voyager_Read_CDF_Product(files(k),'mag48s');
    t=q.Epoch(:); b=q.F1(:); use=t>=begin&t<finish;
    records=[records;table(t(use),b(use),repmat(k,nnz(use),1),find(use), ...
        'VariableNames',{'EpochUTC','F1_nT','SourceFileIndex','CDFRecord'})]; %#ok<AGROW>
    manifest=[manifest;table(k,files(k),string(Case1_File_SHA256(files(k))),{q.variable_meta.F1},{q.global_attributes}, ...
        'VariableNames',{'FileIndex','SourceFile','SHA256','F1Metadata','GlobalAttributes'})]; %#ok<AGROW>
end
records=sortrows(records,'EpochUTC');
% Only COHO-missing daily bins are candidates. Valid COHO days never use MAG.
[matched,dayIndex]=ismember(dateshift(records.EpochUTC,'start','day')+hours(12),d.EpochUTC); assert(all(matched));
[~,~,group]=unique(records.EpochUTC); multiplicity=accumarray(group,1);
candidate=~isfinite(d.BMean_nT(dayIndex));
unusedDuplicateRecords=records(multiplicity(group)>1 & ~candidate,:);
records=records(candidate,:);
assert(numel(unique(records.EpochUTC))==height(records),'Duplicate MAG timestamps on a missing COHO day require review.');
good=isfinite(records.F1_nT); t=records.EpochUTC(good); b=records.F1_nT(good);
[hour,~,bin]=unique(dateshift(t,'start','hour'));
hourMean=accumarray(bin,b,[],@mean); hourCount=accumarray(bin,1,[],@sum);
[matched,dayBin]=ismember(dateshift(hour,'start','day')+hours(12),d.EpochUTC); assert(all(matched));
n=height(d); meanMAG=accumarray(dayBin,hourMean,[n 1],@mean,NaN);
validHours=accumarray(dayBin,1,[n 1],@sum,0);
rawCount=accumarray(dayBin,hourCount,[n 1],@sum,0);
assert(sum(rawCount)==nnz(good));
originalB=d.BMean_nT; use=~isfinite(originalB)&isfinite(meanMAG);
d.BMean_nT(use)=meanMAG(use); d.MAGSampleCount(use)=validHours(use);
d.MAGDailySource=repmat("COHO",n,1); d.MAGDailySource(~isfinite(originalB))="missing";
d.MAGDailySource(use)="MAG48s_hourly_then_daily";
assert(isequaln(d.BMean_nT(isfinite(originalB)),originalB(isfinite(originalB))));
audit=struct('SourceFiles',manifest,'RawSelectedRecords',records,'UnusedDuplicateRecords',unusedDuplicateRecords, ...
    'Hourly',table(hour,hourMean,hourCount),'Daily',table(d.EpochUTC,originalB,meanMAG,validHours,rawCount,use), ...
    'Method','Original reviewed MAG F1 scalar; inherited CDF fill/valid screening. Mean of finite 48-s samples per UTC hour, then equal mean of finite hourly values per UTC day. Only wholly missing COHO daily values replaced. No interpolation, coverage threshold or modification of valid COHO daily values.', ...
    'CodeSHA256',Case1_File_SHA256([mfilename('fullpath'),'.m']));
source.Daily=d; source.MAGSupplement=audit;
source.Method.PreviousPanelA=source.Method.PanelA;
source.Method.PanelA='COHO hourly scalar magnitude daily mean; missing daily bins supplemented by reviewed MAG F1 48-s -> hourly mean -> daily mean.';
source.Method.Missing='No interpolation or synthetic filling; only measured MAG data supply missing COHO daily magnitude. Other variables retain inherited missing rules.';
source.Method.MAGSupplement=audit.Method;
valid=isfinite(d.BMean_nT); source.Coverage.ValidDays(1)=nnz(valid);
source.Coverage.FirstValidUTC(1)=dateshift(d.EpochUTC(find(valid,1)),'start','day');
source.Coverage.LastValidUTC(1)=dateshift(d.EpochUTC(find(valid,1,'last')),'start','day');
result=source;
save(fullfile(source.OutputFolder,'daily_with_MAG_supplement_audit.mat'),'result','-v7.3');
writetable(d,fullfile(source.OutputFolder,'daily_with_MAG_supplement.csv'));
fprintf('V2 MAG supplement: %d originally missing days supplied; %d still missing.\n',nnz(use),nnz(~isfinite(d.BMean_nT)));
end
