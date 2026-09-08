function validation=Validate_V1_Daily_Overview
% Independent accumarray statistics and comparison with reference-event audit.
base='Z:/SPART-WORK/Data/Voyager/voyager1';
out = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview';
s=load(fullfile(out,'V1_daily_overview_audit.mat'));
d=s.result.Daily; raw=s.raw;
n=height(d); first=dateshift(d.EpochUTC(1),'start','day');
g=days(dateshift(raw.EpochUTC,'start','day')-first)+1;
fields={'B_nT','P1','P1'}; targets={'BMean_nT','P1Mean','P1Median'};
for k=1:3
    x=raw.(fields{k}); good=isfinite(x);
    fun=@mean; if k==3, fun=@median; end
    value=accumarray(g(good),x(good),[n 1],fun,NaN);
    assert(isequaln(value,d.(targets{k})),'Independent daily statistic mismatch.');
end
assert(all(diff(d.EpochUTC)==days(1)),'UTC grid is incomplete.');
assert(isequaln(sum(d.SectorDailyMean(:,[1 2 3 5 6 7]),2),d.P1SixSectorSum));
assert(all(isnan(d.P1SixSectorSum(~d.SixSectorsComplete))));
% Exact comparison to the user's reference event, when its audit is present.
e=dir(fullfile(base,'lecp','1d','derived','pitch_angle','2013-2021', ...
    'predicted_ck','V1_Case1-S01-L01_*_1d_nativeCDF_Epoch.mat'));
referenceCompared=0;
if numel(e)==1
    ref=load(fullfile(e.folder,e.name)); names=fieldnames(ref);
    for k=1:numel(names)
        t=ref.(names{k});
        if ~istable(t) || ~ismember('EpochUTC',t.Properties.VariableNames), continue, end
        cols=cell(1,6); sectors=[1 2 3 5 6 7];
        for j=1:6, cols{j}=sprintf('RawFlux_S%d_1d',sectors(j)); end
        if ~all(ismember(cols,t.Properties.VariableNames)), continue, end
        for j=1:height(t)
            row=find(dateshift(d.EpochUTC,'start','day')==dateshift(t.EpochUTC(j),'start','day'));
            assert(numel(row)==1);
            assert(isequaln(d.SectorDailyMean(row,sectors),t{j,cols}), ...
                'Reference-event sector flux mismatch.');
            referenceCompared=referenceCompared+1;
        end
    end
end
assert(referenceCompared>0,'Reference-event audit was not compared.');
archiveFile = 'Z:/SPART-WORK/Data/Voyager/source_verification/V1_daily_overview/official_archive_inventory.mat';
a=load(archiveFile,'audit');
assert(a.audit.VerifiedOnline);
assert(all(isfile(a.audit.Files.SourceFile)));
assert(s.result.Archive.VerifiedOnline); % 验证读取审计，不修改主结果
validation=struct('Passed',true,'DailyDays',n,'ReferenceRowsCompared',referenceCompared, ...
    'ValidatedUTC',datetime('now','TimeZone','UTC'));
save(fullfile(out,'validation.mat'),'validation'); disp(validation);
end



