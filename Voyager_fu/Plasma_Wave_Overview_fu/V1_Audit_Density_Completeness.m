function [density,audit] = V1_Audit_Density_Completeness(DataDir,PreviousDir,OutputDir,firstUTC,lastUTC)
% 核对现行官方密度 CSV、旧版源文件、旧图全部点与本次所用原始记录。
%% 明确源路径；CSV 例外已经用户批准
base=fullfile(DataDir,'pws','derived','electron_density','native');
currentFile=fullfile(base,'PDS_release_20260910','data','vg1-vlism-density-2012-2025.csv');
currentLabel=fullfile(base,'PDS_release_20260910','data','vg1-vlism-density-2012-2025.lblx');
liveFile=fullfile(base,'source_check_20260915','vg1-vlism-density-2012-2025.csv');
oldFile=fullfile(base,'PDS_annex_release_1','data','vg1-vlism-density-2012-2025.csv');
oldLabel=fullfile(base,'PDS_annex_release_1','data','vg1-vlism-density-2012-2025.lblx');
assert(strcmp(Case1_File_SHA256(currentFile),Case1_File_SHA256(liveFile)), ...
    'Official current CSV changed since preceding run; review new version.');
[density,currentRaw,currentText]=readDensity(currentFile,currentLabel);
[old,oldRaw,oldText]=readDensity(oldFile,oldLabel);
selected=density.EpochUTC>=firstUTC & density.EpochUTC<lastUTC;
assert(all(selected),'Some official density points lie outside current plot window.');
assert(all(isfinite(density.Density_cm3) & density.Density_cm3>0), ...
    'Source contains invalid densities; review official flags.');
density=sortrows(density,'EpochUTC');
% SCET + source 标识同类测量，保留同一 SCET 的不同 EPO/QTN 记录。
oldKey=string(oldRaw{:,1})+"|"+old.Source;
newKey=string(currentRaw{:,1})+"|"+lower(string(currentRaw{:,18}));
% 重复时间/来源键也按源行保留。版本比较使用 SCET + Source + N_e 的多重集合，
% 一条旧记录至多匹配一条新记录；此比较不改变任何绘图记录。
newNative=readDensity(currentFile,currentLabel);
oldScienceKey=oldKey+"|"+compose('%.17g',old.Density_cm3);
newScienceKey=newKey+"|"+compose('%.17g',newNative.Density_cm3);
usedOld=false(height(old),1); matchedOldRow=zeros(height(newNative),1);
for k=1:height(newNative)
    hit=find(~usedOld & oldScienceKey==newScienceKey(k),1,'first');
    if ~isempty(hit), usedOld(hit)=true; matchedOldRow(k)=hit; end
end
oldOnly=old(~usedOld,:);
newOnly=newNative(matchedOldRow==0,:);
common=intersect(oldKey,newKey,'stable');
densityChanges=table('Size',[0 5],'VariableTypes',{'string','string','string','double','double'}, ...
    'VariableNames',{'SCET_Source','OldDensityValues_cm3','CurrentDensityValues_cm3','OldRecordCount','CurrentRecordCount'});
for k=1:numel(common)
    a=old.Density_cm3(oldKey==common(k)); b=newNative.Density_cm3(newKey==common(k));
    if ~isequal(sort(unique(a)),sort(unique(b)))
        densityChanges(end+1,:)={common(k),join(compose('%.8g',unique(a)),','), ...
            join(compose('%.8g',unique(b)),','),numel(a),numel(b)};
    end
end
duplicateKeys=[table(old.EpochUTC,old.Source,old.Density_cm3,repmat("annex",height(old),1), ...
    'VariableNames',{'EpochUTC','Source','Density_cm3','Version'}); ...
    table(newNative.EpochUTC,newNative.Source,newNative.Density_cm3,repmat("current",height(newNative),1), ...
    'VariableNames',{'EpochUTC','Source','Density_cm3','Version'})];
dupKey=[oldKey+"|annex";newKey+"|current"];
[~,~,g]=unique(dupKey);counts=accumarray(g,1);
duplicateKeys=duplicateKeys(counts(g)>1,:);
writetable(duplicateKeys,fullfile(OutputDir,'density_duplicate_time_source_records.csv'));

writetable(oldOnly,fullfile(OutputDir,'density_old_version_only_records.csv'));
writetable(newOnly,fullfile(OutputDir,'density_current_version_only_records.csv'));
writetable(densityChanges,fullfile(OutputDir,'density_version_value_changes.csv'));
writetable(density,fullfile(OutputDir,'V1_PWS_density_native_records.csv'));

%% 检查此前 FIG 中 EPO 和 QTN 的全部坐标，确认是否漏画或轴限裁切
previousFile=fullfile(PreviousDir,'V1_1990_20250630_daily_PWS_5panels.fig');
f=openfig(previousFile,'invisible'); cleaner=onCleanup(@()close(f)); %#ok<NASGU>
names=["epo","qtn"]; count=zeros(2,1); expected=zeros(2,1); match=false(2,1); visible=false(2,1);
for k=1:2
    h=findobj(f,'Type','line','DisplayName',char(upper(names(k))));
    assert(numel(h)==1,'Previous density series not uniquely identified.');
    use=density.Source==names(k);
    count(k)=numel(h.XData);expected(k)=nnz(use);
    match(k)=isequaln(h.XData(:),datenum(density.EpochUTC(use))) && ...
        isequaln(h.YData(:),density.Density_cm3(use));
    ax=ancestor(h,'axes');
    visible(k)=all(h.XData>=ax.XLim(1)&h.XData<=ax.XLim(2)) && ...
        all(h.YData>=ax.YLim(1)&h.YData<=ax.YLim(2));
end
figureCheck=table(names(:),expected,count,match,visible,'VariableNames', ...
    {'Source','SourceRecords','PreviousFigurePoints','ExactCoordinateMatch','AllInsideAxisLimits'});
assert(all(match & visible),'Previous plot omitted or clipped official density records.');
writetable(figureCheck,fullfile(OutputDir,'density_previous_figure_check.csv'));

%% 完整时间轴逐年覆盖，以及所有连续无密度点日期段
years=(year(firstUTC):year(lastUTC)).';
yearCoverage=table(years,zeros(size(years)),zeros(size(years)),zeros(size(years)),zeros(size(years)), ...
    'VariableNames',{'Year','DensityRecords','DaysWithDensity','EPORecords','QTNRecords'});
for k=1:numel(years)
    use=year(density.EpochUTC)==years(k);
    yearCoverage.DensityRecords(k)=nnz(use);
    yearCoverage.DaysWithDensity(k)=numel(unique(dateshift(density.EpochUTC(use),'start','day')));
    yearCoverage.EPORecords(k)=nnz(use & density.Source=="epo");
    yearCoverage.QTNRecords(k)=nnz(use & density.Source=="qtn");
end
day=(firstUTC:days(1):lastUTC-days(1)).';
present=ismember(day,dateshift(density.EpochUTC,'start','day'));
d=diff([false;~present;false]);start=find(d==1);stop=find(d==-1)-1;
gaps=table(day(start),day(stop),stop-start+1, ...
    'VariableNames',{'FirstEmptyUTCDate','LastEmptyUTCDate','EmptyDays'});
writetable(yearCoverage,fullfile(OutputDir,'density_yearly_coverage.csv'));
writetable(gaps,fullfile(OutputDir,'density_all_empty_date_intervals.csv'));
audit=struct('CurrentCSV',currentFile,'CurrentLabel',currentLabel,'CurrentLabelText',currentText, ...
    'CurrentSHA256',Case1_File_SHA256(currentFile),'CurrentRawTable',currentRaw, ...
    'LiveVerificationCSV',liveFile,'OldCSV',oldFile,'OldLabel',oldLabel, ...
    'OldLabelText',oldText,'OldRawTable',oldRaw,'OldSHA256',Case1_File_SHA256(oldFile), ...
    'OfficialRecords',height(density),'OldRecords',height(old),'ValidDays',nnz(present), ...
    'FirstUTC',min(density.EpochUTC),'LastUTC',max(density.EpochUTC), ...
    'EPORecords',nnz(density.Source=="epo"),'QTNRecords',nnz(density.Source=="qtn"), ...
    'OldOnly',oldOnly,'CurrentOnly',newOnly,'ChangedDensity',densityChanges, ...
    'PreviousFigure',previousFile,'PreviousFigureCheck',figureCheck, ...
    'YearlyCoverage',yearCoverage,'EmptyDateIntervals',gaps, ...
    'VersionPolicy','Current official data release used in full. Annex version archived for comparison only. Official errata do not explain all density-record differences; no speculative merging.', ...
    'MissingPolicy','No interpolation, artificial daily density, density inversion or EPO/QTN source replacement.');
fprintf('Density: %d current records on %d UTC days; all previous points verified.\n',height(density),nnz(present));
fprintf('Version audit: old-only %d; current-only %d; changed density %d.\n',height(oldOnly),height(newOnly),height(densityChanges));
end

function [density,raw,label]=readDensity(file,labelFile)
label=fileread(labelFile);
token=regexp(label,'<md5_checksum>([^<]+)</md5_checksum>','tokens','once');
assert(strcmpi(fileMD5(file),token{1}),'Density source does not match official MD5.');
options=detectImportOptions(file,'VariableNamingRule','preserve','TextType','string');
options=setvartype(options,options.VariableNames{1},'string');
raw=readtable(file,options);
n=regexp(label,'<records>(\d+)</records>','tokens','once');
assert(height(raw)==str2double(n{1}) && width(raw)==22);
assert(strcmp(raw.Properties.VariableNames{12},'N_e (cm^-3)'));
t=datetime(raw{:,1},'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
t.Format="yyyy-MM-dd'T'HH:mm:ss.SSS'Z'";
source=lower(string(raw{:,18}));assert(all(ismember(source,["epo","qtn"])));
density=table(t,raw{:,12},raw{:,15},raw{:,16},source,(1:height(raw)).', ...
    'VariableNames',{'EpochUTC','Density_cm3','Minimum_cm3','Maximum_cm3','Source','OriginalCSVRow'});
end

function value=fileMD5(file)
fid=fopen(file,'rb');assert(fid>=0);cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');digest=java.security.MessageDigest.getInstance('MD5');
digest.update(typecast(bytes,'int8'));
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
