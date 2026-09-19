function result = Run_V2_Plasma_Overview(syncArchive,includePWSDensity)
% Voyager 2: 原 a/b/c 日均曲线 + PLS 质子密度、温度、总流速。
% 默认科学输入为原始网站 CDF；派生表仅用于结果交付和验证。
if nargin<1, syncArchive=true; end
if nargin<2, includePWSDensity=true; end
%% 路径与时间参数
CodeDir=fileparts(mfilename('fullpath')); ProjectCode=fileparts(CodeDir);
DataDir='Z:/SPART-WORK/Data/Voyager/voyager2/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/';
if includePWSDensity,OutputDir=fullfile(OutputDir,'With_Official_PWS_Density');end
OriginalDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V2/';
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
addpath(CodeDir,fullfile(ProjectCode,'Case1_PPT_VerticalLine_Events_7d_fu'), ...
    fullfile(ProjectCode,'Interstellar_Daily_Overview_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
assert(isfolder(DataDir),'Authorized Z: archive unavailable.');
if ~isfolder(OutputDir), mkdir(OutputDir); end
SourceOutput=fullfile(OutputDir,'source_recompute');
if ~isfolder(SourceOutput), mkdir(SourceOutput); end
firstUTC=datetime(1990,1,1,'TimeZone','UTC');
lastUTC=datetime(2025,1,1,'TimeZone','UTC'); % exclusive; previous V2 full overview
before=originalHashes(OriginalDir,includePWSDensity);
if syncArchive
    PowerShellExe='C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/powershell/pwsh.exe';
    command=sprintf('"%s" -NoProfile -NonInteractive -File "%s"',PowerShellExe, ...
        fullfile(CodeDir,'Sync_V2_Plasma_Overview_Archive.ps1'));
    status=system(command);assert(status==0,'Official source verification failed.');
    if includePWSDensity
        command=sprintf('"%s" -NoProfile -NonInteractive -File "%s"',PowerShellExe, ...
            fullfile(CodeDir,'Sync_V2_PWS_Density_Archive.ps1'));
        status=system(command);assert(status==0,'Official PWS density verification failed.');
    end
end
inventoryFile=fullfile(DataDir,'source_verification','plasma_overview_1990_2024','official_source_inventory.json');
inventory=jsondecode(fileread(inventoryFile));
assert(inventory.LatestOfficialCOHOYear==2024,'Official COHO time range changed; update plot end explicitly.');
items=inventory.Files;
files=string({items(strcmp({items.Product},'COHO')).File}).';
assert(~isempty(files) && all(isfile(files)));
pws=table; pwsAudit=struct;
if includePWSDensity
    [pws,pwsAudit]=V2_Read_Official_PWS_Density(DataDir);
    lastUTC=max(lastUTC,dateshift(max(pws.EpochUTC),'start','day')+days(1));
end
day=(firstUTC:days(1):lastUTC-days(1)).';
%% 原始 COHO CDF：磁场、LECP P1 及 PLS 参数
parts=cell(numel(files),1); cohoAudit=cell(numel(files),1);
for k=1:numel(files)
    q=Voyager_Read_CDF_Product(files(k),'coho');
    b=NaN(numel(q.Epoch),1); magnitudeVariable="missing";
    if isfield(q,'ABS_B') && any(isfinite(q.ABS_B))
        b=q.ABS_B(:); magnitudeVariable="ABS_B";
    elseif isfield(q,'F') && any(isfinite(q.F))
        b=q.F(:); magnitudeVariable="F";
    end
    variable=["protonFlux1_LECP","protonDensity","protonTemp","V"];
    values=NaN(numel(b),4);
    for c=1:4
        assert(isfield(q,variable(c)),'Expected COHO variable missing: %s',variable(c));
        values(:,c)=q.(variable(c))(:);
    end
    assert(contains(q.variable_meta.protonFlux1_LECP.attributes.FIELDNAM,'0.52-1.45'));
    assert(strcmp(strtrim(q.variable_meta.protonFlux1_LECP.attributes.UNITS),'1/(cm^2 sec ster MeV)'));
    assert(strcmpi(strtrim(q.variable_meta.protonDensity.attributes.UNITS),'n/cc'));
    assert(strcmpi(strtrim(q.variable_meta.protonTemp.attributes.UNITS),'deg k'));
    assert(strcmpi(strtrim(q.variable_meta.V.attributes.UNITS),'km/s'));
    part=table(q.Epoch(:),b,values(:,1),values(:,2),values(:,3),values(:,4), ...
        repmat(k,numel(b),1),(1:numel(b)).','VariableNames', ...
        {'EpochUTC','B_nT','P1','ProtonDensity_cm3','ProtonTemperature_K','ProtonSpeed_kms','FileIndex','CDFRecord'});
    use=part.EpochUTC>=firstUTC & part.EpochUTC<lastUTC; parts{k}=part(use,:);
    cohoAudit{k}=struct('File',files(k),'SHA256',Case1_File_SHA256(files(k)), ...
        'MagnitudeVariable',magnitudeVariable,'VariableMetadata',q.variable_meta, ...
        'GlobalAttributes',q.global_attributes,'SelectedRecords',nnz(use));
    if mod(k,50)==0 || k==numel(files), fprintf('COHO %d/%d\n',k,numel(files)); end
end
raw=sortrows(vertcat(parts{:}),'EpochUTC'); clear parts
[~,firstRecord,groups]=unique(raw.EpochUTC,'stable');reference=firstRecord(groups);
assert(isequaln(raw{:,2:6},raw{reference,2:6}),'Conflicting duplicate COHO records require review.');
duplicates=raw(~ismember((1:height(raw)).',firstRecord),:);raw=raw(firstRecord,:);
means=NaN(numel(day),5);counts=zeros(numel(day),5);
for c=1:5, [means(:,c),counts(:,c)]=dailyMean(raw.EpochUTC,raw{:,c+1},day); end
%% 同原图的已获批磁场补充：只补 COHO 整日缺测
source=struct;
source.OutputFolder=SourceOutput;
source.Daily=table(day+hours(12),means(:,1),counts(:,1), ...
    'VariableNames',{'EpochUTC','BMean_nT','MAGSampleCount'});
source.Method=struct('PanelA','Finite COHO hourly ABS_B/F daily mean.', ...
    'Missing','CDF fill/valid metadata; no interpolation.');
source.Coverage=table("MAG",firstUTC,lastUTC-days(1),nnz(isfinite(means(:,1))), ...
    'VariableNames',{'Panel','FirstValidUTC','LastValidUTC','ValidDays'});
source=Voyager_Supplement_V2_MAG(source);
means(:,1)=source.Daily.BMean_nT;counts(:,1)=source.Daily.MAGSampleCount;
%% 太阳黑子数及六面板日统计
solarRaw=readmatrix(SolarFile,'Delimiter',';','FileType','text');
solarEpoch=datetime(solarRaw(:,1),solarRaw(:,2),solarRaw(:,3),'TimeZone','UTC');
solarValue=solarRaw(:,5);solarValue(solarValue==-1)=NaN;
assert(numel(unique(solarEpoch))==numel(solarEpoch));
[found,solarRows]=ismember(day,solarEpoch);sunspot=NaN(numel(day),1);
sunspot(found)=solarValue(solarRows(found));
daily=table(day+hours(12),sunspot,means(:,1),means(:,2),means(:,3),means(:,4),means(:,5), ...
    counts(:,1),counts(:,2),counts(:,3),counts(:,4),counts(:,5),source.Daily.MAGDailySource,solarRows, ...
    'VariableNames',{'EpochUTC','SunspotNumber','B_nT','P1Mean','ProtonDensity_cm3', ...
    'ProtonTemperature_K','ProtonSpeed_kms','MAGSampleCount','P1SampleCount', ...
    'DensitySampleCount','TemperatureSampleCount','SpeedSampleCount','MAGDailySource','SolarSourceRow'});
% 旧 CSV 仅用于原 a/b/c 数值验证，不参与科学计算。
baseline=readtable(fullfile(OriginalDir,'all_daily_values.csv'));
assert(height(baseline)<=height(daily));
oldTime=baseline.EpochUTC;
if ~isdatetime(oldTime), oldTime=datetime(oldTime); end
oldTime.TimeZone='UTC';[foundOld,oldRows]=ismember(oldTime,daily.EpochUTC);assert(all(foundOld));
newValues=daily{oldRows,2:4};oldValues=[baseline.SunspotNumber baseline.B_nT baseline.P1Mean];
assert(isequal(isnan(newValues),isnan(oldValues)));
good=isfinite(newValues)&isfinite(oldValues);
assert(all(abs(newValues(good)-oldValues(good))<=1e-12*max(1,abs(oldValues(good)))));
preservation=table(["SunspotNumber";"B_nT";"P1Mean"], ...
    max(abs(newValues-oldValues),[],1,'omitnan').', ...
    max(abs(newValues-oldValues)./max(1,abs(oldValues)),[],1,'omitnan').', ...
    'VariableNames',{'Variable','MaximumAbsoluteDifference','MaximumRelativeDifference'});
fprintf('Previous a/b/c statistics match; MAG gap supplement preserved.\n');
%% 输出与审计
coverage=makeCoverage(daily);
if includePWSDensity
    assert(all(pws.EpochUTC>=firstUTC & pws.EpochUTC<lastUTC));
    coverage(end+1,:)={"PWSElectronDensity_cm3",min(pws.EpochUTC),max(pws.EpochUTC), ...
        numel(unique(dateshift(pws.EpochUTC,'start','day'))),0};
    writetable(pws,fullfile(OutputDir,'V2_PWS_density_native_records.csv'));
end
writetable(daily,fullfile(OutputDir,'V2_daily_six_panel_values.csv'));
writetable(coverage,fullfile(OutputDir,'coverage.csv'));
writetable(preservation,fullfile(OutputDir,'abc_preservation.csv'));
emptyIntervals=missingIntervals(daily);writetable(emptyIntervals,fullfile(OutputDir,'missing_date_intervals.csv'));
method=struct('StartUTC',firstUTC,'EndUTCExclusive',lastUTC, ...
    'PanelA','Original SILSO daily sunspot number; -1 missing, zero valid.', ...
    'PanelB',source.Method.PanelA, ...
    'PanelC','Finite arithmetic mean of original COHO hourly protonFlux1_LECP; 0.52-1.45 MeV.', ...
    'PanelD','Finite arithmetic mean of original COHO hourly protonDensity, supplied by PLS; proton number density in cm^-3.', ...
    'PanelE','Finite arithmetic mean of original COHO hourly protonTemp in K. Source derives T=60.5*Vth^2. Use source T without recomputing it from averaged Vth.', ...
    'PanelF','Finite arithmetic mean of original COHO hourly V in km/s; scalar bulk speed, no reconstruction from RTN components.', ...
    'Missing','Original CDF FILLVAL/VALIDMIN/VALIDMAX and existing reader handling. No added thresholds, interpolation, zero filling, smoothing or fit. Each parameter averaged independently.', ...
    'PLSProduct','COHO hourly plasma parameters retained as one consistent product. Independent high-resolution PLS files verified for archive completeness; not mixed into hourly daily means.', ...
    'EnergyLabel','V2 LECP P1: 0.52-1.45 MeV, following Decker et al. (1999), ICRC 26, vol. 6, p. 328, Fig. 1, and Rice et al. (2000), GRL 27, 509-512. Matches original COHO metadata; no flux rescaling.', ...
    'DensityType','Proton density from PLS; no electron/proton equivalence assumption and no PWS electron density substitution.', ...
    'Display','P1, proton density and temperature logarithmic; sunspot, B and speed linear. Thin lines; missing values break lines. Nonpositive values hidden only on log display.', ...
    'Time','Full V2 overview retains user-approved 1990 start and latest COHO year 2024. No data invented after source parameter coverage ends.');
if includePWSDensity
    method.PanelD='PLS proton density daily means plus all native official PWS electron density points; two source series retained separately in audit.';
    method.DensityType='PLS n_p and PWS n_e shown together as density, without converting species or numerically filling the PLS array.';
    method.PWSDensity=pwsAudit.Method;
    method.Time='Start 1990-01-01; end extended to the latest native PWS density point when later than COHO. No synthetic MAG, temperature, speed or flux values after their products end.';
    method.Display=[method.Display,' PWS native points without bridging gaps; no added source legend or processing footnotes, following the user request. Right plot margin accommodates the last point.'];
end
result=struct('Daily',daily,'RawCOHO',raw,'COHOSources',{cohoAudit},'COHOIdenticalDuplicates',duplicates, ...
    'MAGSupplement',source.MAGSupplement,'Method',method,'Coverage',coverage, ...
    'Preservation',preservation,'MissingIntervals',emptyIntervals,'OfficialInventory',inventory, ...
    'SolarFile',SolarFile,'SolarSHA256',Case1_File_SHA256(SolarFile),'ProtectedBefore',before);
result.PWSDensity=pws;result.PWSDensityAudit=pwsAudit;
result.IncludePWSDensity=includePWSDensity;
result.OutputFiles=drawFigure(daily,firstUTC,lastUTC,OutputDir,pws);
if includePWSDensity
    verificationFile=fullfile(DataDir,'source_verification','plasma_overview_density_extension','official_sources.json');
    result.PWSSourceInventory=jsondecode(fileread(verificationFile));
end
after=originalHashes(OriginalDir,includePWSDensity);assert(isequal(before,after),'Previous output changed unexpectedly.');
result.ProtectedAfter=after;result.CreatedUTC=datetime('now','TimeZone','UTC');result.MATLABVersion=version;
dependency=string({[mfilename('fullpath'),'.m'],which('Voyager_Read_CDF_Product'), ...
    which('Voyager_Supplement_V2_MAG'),which('dataobj'),fullfile(CodeDir,'Sync_V2_Plasma_Overview_Archive.ps1')}).';
if includePWSDensity
    dependency=[dependency;string(which('V2_Read_Official_PWS_Density')); ...
        string(fullfile(CodeDir,'Sync_V2_PWS_Density_Archive.ps1'))];
end
hash=strings(size(dependency));
for k=1:numel(dependency),hash(k)=Case1_File_SHA256(dependency(k));end
result.Code=table(dependency,hash,'VariableNames',{'File','SHA256'});
save(fullfile(OutputDir,'V2_plasma_overview_audit.mat'),'result','-v7');
disp(coverage);fprintf('Saved six-panel V2 overview: %s\n',result.OutputFiles(1));
end

function [avg,count]=dailyMean(time,value,day)
[matched,bin]=ismember(dateshift(time,'start','day'),day);good=matched&isfinite(value);n=numel(day);
avg=NaN(n,1);count=zeros(n,1);
if any(good)
    avg=accumarray(bin(good),value(good),[n 1],@mean,NaN);
    count=accumarray(bin(good),1,[n 1],@sum,0);
end
assert(sum(count)==nnz(good));
end

function h=originalHashes(folder,protectSix)
names=["V2_all_daily_5panels.png";"V2_all_daily_5panels.pdf";"V2_all_daily_5panels.fig";"all_daily_values.csv"];
paths=fullfile(folder,names);
if protectSix
    previous='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/';
    names=["V2_1990_20241231_daily_PLS_6panels.png";"V2_1990_20241231_daily_PLS_6panels.pdf"; ...
        "V2_1990_20241231_daily_PLS_6panels.fig";"V2_daily_six_panel_values.csv";"V2_plasma_overview_audit.mat"];
    paths=[paths;fullfile(previous,names)];
end
hash=strings(size(paths));
for k=1:numel(paths),hash(k)=Case1_File_SHA256(paths(k));end
h=table(paths,hash,'VariableNames',{'File','SHA256'});
end

function coverage=makeCoverage(daily)
names=string(daily.Properties.VariableNames(2:7)).';v=daily{:,2:7};
coverage=table(names,NaT(6,1,'TimeZone','UTC'),NaT(6,1,'TimeZone','UTC'),zeros(6,1),zeros(6,1), ...
    'VariableNames',{'Parameter','FirstValidUTC','LastValidUTC','ValidDays','NonpositiveDays'});
for k=1:6
    good=isfinite(v(:,k));t=daily.EpochUTC(good);
    if ~isempty(t)
        coverage.FirstValidUTC(k)=min(t);coverage.LastValidUTC(k)=max(t);
        coverage.ValidDays(k)=nnz(good);coverage.NonpositiveDays(k)=nnz(good & v(:,k)<=0);
    end
end
end

function gaps=missingIntervals(daily)
gaps=table('Size',[0 4],'VariableTypes',{'string','datetime','datetime','double'}, ...
    'VariableNames',{'Parameter','FirstEmptyUTCDate','LastEmptyUTCDate','EmptyDays'});
gaps.FirstEmptyUTCDate.TimeZone='UTC';gaps.LastEmptyUTCDate.TimeZone='UTC';
day=dateshift(daily.EpochUTC,'start','day');
for c=3:7
    missing=~isfinite(daily{:,c});change=diff([false;missing;false]);
    first=find(change==1);last=find(change==-1)-1;
    part=table(repmat(string(daily.Properties.VariableNames{c}),numel(first),1),day(first),day(last),last-first+1, ...
        'VariableNames',gaps.Properties.VariableNames);
    gaps=[gaps;part]; %#ok<AGROW>
end
end

function files=drawFigure(daily,firstUTC,lastUTC,out,pws)
%% 六个上下对齐 panel；不添加处理脚注
f=figure('Color','w','Position',[30 30 1600 1440],'Visible','off');
ax=gobjects(6,1);height=.1375;gap=.015;top=.95;
colors=[.70 .30 .06;.15 .15 .15;.48 .12 .62;.08 .36 .62;.75 .25 .12;.05 .50 .40];
labels={{'Daily sunspot number'},{'|B| daily mean','(nT)'}, ...
    {'LECP P1 daily mean','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'PLS proton density','Daily mean (cm^{-3})'}, ...
    {'PLS proton temperature','Daily mean (K)'},{'PLS bulk speed','Daily mean (km s^{-1})'}};
values=daily{:,2:7};t=datenum(daily.EpochUTC);isLog=[false false true true true false];
plotEnd=lastUTC;
if ~isempty(pws)
    labels{4}={'Density','(cm^{-3})'};
    plotEnd=lastUTC+days(45); % display margin only; no added data rows
end
for k=1:6
    ax(k)=axes(f,'Position',[.105 top-k*height-(k-1)*gap .87 height],'Tag',sprintf('panel_%d',k));
    y=values(:,k);if isLog(k),y(y<=0)=NaN;end
    h=plot(ax(k),t,y,'.-','Color',colors(k,:),'LineWidth',.4,'MarkerSize',3, ...
        'Tag',daily.Properties.VariableNames{k+1});
    assert(isequaln(h.YData(:),y));
    if k==4 && ~isempty(pws)
        hold(ax(k),'on');
        plot(ax(k),datenum(pws.EpochUTC),pws.ElectronDensity_cm3,'o', ...
            'Color',colors(k,:),'MarkerFaceColor',colors(k,:),'MarkerSize',4.5, ...
            'LineStyle','none','Tag','PWS_ElectronDensity_native');
        hold(ax(k),'off');
    end
    if isLog(k),set(ax(k),'YScale','log');end
    if k==4 && ~isempty(pws)
        limits=ylim(ax(k));ylim(ax(k),[limits(1) max(limits(2),1.5*max(pws.ElectronDensity_cm3))]);
    end
    ylabel(ax(k),labels{k});
    set(ax(k),'Tag',sprintf('panel_%d',k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',.12);
    xlim(ax(k),datenum([firstUTC plotEnd]));
    xticks(ax(k),datenum(datetime(1990:5:2025,1,1,'TimeZone','UTC')));
    if k<6,xticklabels(ax(k),{});else,xticklabels(ax(k),string(1990:5:2025));end
    text(ax(k),.008,.88,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
end
linkaxes(ax,'x');xlabel(ax(6),'UTC');
annotation(f,'textbox',[.105 .958 .87 .03],'String', ...
    sprintf('Voyager 2 | Daily | 1990-01-01 to %s | P1 0.52-1.45 MeV',datestr(lastUTC-days(1),'yyyy-mm-dd')), ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',15);
if isempty(pws)
    stem=fullfile(out,'V2_1990_20241231_daily_PLS_6panels');
else
    stem=fullfile(out,sprintf('V2_1990_%s_daily_6panels_PWS',datestr(lastUTC-days(1),'yyyymmdd')));
end
files=string({[stem,'.png'],[stem,'.pdf'],[stem,'.fig']});
exportgraphics(f,files(1),'Resolution',200);
exportgraphics(f,files(2),'ContentType','image','Resolution',200);
savefig(f,files(3));close(f);
end