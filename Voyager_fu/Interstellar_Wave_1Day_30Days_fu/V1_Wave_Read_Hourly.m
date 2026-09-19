function [raw,manifest,duplicates] = V1_Wave_Read_Hourly(cfg)
% Read the latest local version of each original monthly COHO CDF.
base = fullfile(cfg.DataRoot,'voyager1','coho','1hr','l2','merged_mag_plasma');
entries = dir(fullfile(base,'**','voyager1_coho1hr_merged_mag_plasma_*_v*.cdf'));
files = string(fullfile({entries.folder},{entries.name})).';
keys = strings(size(files)); versions=zeros(size(files)); keep=false(size(files));
for k=1:numel(files)
    tok=regexp(char(files(k)),'_(\d{8})_v(\d+)\.cdf$','tokens','once');
    keys(k)=string(tok{1}); versions(k)=str2double(tok{2});
    date=datetime(tok{1},'InputFormat','yyyyMMdd','TimeZone','UTC');
    keep(k)=date<cfg.StopUTC && date+calmonths(1)>cfg.StartUTC;
end
files=files(keep); keys=keys(keep); versions=versions(keep);
selected=false(size(files)); uniqueKeys=unique(keys);
for k=1:numel(uniqueKeys)
    group=find(keys==uniqueKeys(k)); [~,j]=max(versions(group)); selected(group(j))=true;
end
files=sort(files(selected));
expected=(dateshift(cfg.StartUTC,'start','month'):calmonths(1):dateshift(cfg.StopUTC-days(1),'start','month')).';
expectedKeys=string(datestr(expected,'yyyymmdd'));
selectedKeys=keys(selected);
missingKeys=setdiff(expectedKeys,selectedKeys);
verifiedMissing=strings(0,1);
if isfield(cfg,'VerifiedMissingMonths'), verifiedMissing=string(cfg.VerifiedMissingMonths(:)); end
assert(all(ismember(missingKeys,verifiedMissing)), ...
    'Missing local monthly CDFs without official absence verification: %s',strjoin(missingKeys,', '));
assert(all(ismember(selectedKeys,expectedKeys)) && numel(unique(selectedKeys))==numel(files), ...
    'Unexpected monthly CDF selection.');
raw=table; manifest=table;
for k=1:numel(files)
    q=Voyager_Read_CDF_Product(files(k),'coho');
    assert(all(isfield(q,{'Epoch','BR','BT','BN','ABS_B'})),'Required magnetic variables absent.');
    names={'BR','BT','BN','ABS_B'};
    for j=1:numel(names)
        assert(strcmp(strtrim(q.variable_meta.(names{j}).attributes.UNITS),'nT'),'Unexpected field units.');
    end
    t=q.Epoch(:); b=[q.BR(:),q.BT(:),q.BN(:)]; n=numel(t);
    part=table(t,b,q.ABS_B(:),repmat(k,n,1),(1:n).', ...
        'VariableNames',{'EpochUTC','B_RTN_nT','Bmag_nT','SourceIndex','CDFRecordOneBased'});
    in=t>=cfg.StartUTC & t<cfg.StopUTC;
    raw=[raw;part(in,:)]; %#ok<AGROW>
    row=table(files(k),string(Case1_File_SHA256(files(k))),nnz(in), ...
        {q.variable_meta},{q.global_attributes}, ...
        'VariableNames',{'SourceFile','SHA256','RecordsInRange','VariableMetadata','GlobalAttributes'});
    manifest=[manifest;row]; %#ok<AGROW>
    if mod(k,12)==0, fprintf('Read %d / %d CDF files\n',k,numel(files)); end
end
raw=sortrows(raw,'EpochUTC');
[~,first,group]=unique(raw.EpochUTC,'stable'); ref=first(group);
assert(isequaln(raw.B_RTN_nT,raw.B_RTN_nT(ref,:)) && isequaln(raw.Bmag_nT,raw.Bmag_nT(ref)), ...
    'Conflicting duplicate source Epochs; no automatic choice is permitted.');
duplicates=raw(~ismember((1:height(raw)).',first),:); raw=raw(first,:);
dt=seconds(diff(raw.EpochUTC));
assert(all(dt>0),'Epoch must increase.');
assert(all(abs(dt/3600-round(dt/3600))<1e-8),'Unexpected nonhourly source cadence; inspect source.');
end

