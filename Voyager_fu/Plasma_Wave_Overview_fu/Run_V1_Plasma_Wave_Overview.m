function result = Run_V1_Plasma_Wave_Overview(syncArchive,useParallel)
% 保留 V1 原日均概览前 3 个 panel，增加官方 PWS 密度及波动电场频谱。
% 默认直接读取原始 CDF；密度使用用户 2026-09-14 批准的官方原始 CSV。
if nargin<1, syncArchive=true; end
if nargin<2, useParallel=true; end
%% 路径与时间参数
CodeDir=fileparts(mfilename('fullpath'));
ProjectCode=fileparts(CodeDir);
DataDir='Z:/SPART-WORK/Data/Voyager/voyager1/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Plasma_Wave_Overview/';
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
OriginalDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V1/';
addpath(CodeDir,fullfile(ProjectCode,'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
assert(isfolder(DataDir),'Authorized Z: data archive is unavailable.');
if ~isfolder(OutputDir), mkdir(OutputDir); end
firstUTC=datetime(1990,1,1,'TimeZone','UTC');
lastUTC=datetime(2025,7,1,'TimeZone','UTC'); % exclusive
day=(firstUTC:days(1):lastUTC-days(1)).'; nDay=numel(day);
before=originalHashes(OriginalDir);
if syncArchive
    PowerShellExe='C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/powershell/pwsh.exe';
    command=sprintf('"%s" -NoProfile -NonInteractive -File "%s" -StartDate 19900101 -StopDate 20250630 -Workers 8', ...
        PowerShellExe,fullfile(CodeDir,'Sync_V1_PWS_Archive.ps1'));
    status=system(command); assert(status==0,'Official source synchronization failed.');
end

%% 原始 COHO CDF：与原概览相同的日均总磁场和 P1
files=cohoFiles(fullfile(DataDir,'coho','1hr','l2','merged_mag_plasma'),firstUTC,lastUTC);
assert(~isempty(files),'Original COHO CDFs are unavailable.');
raw=table; cohoAudit=cell(numel(files),1);
for k=1:numel(files)
    q=Voyager_Read_CDF_Product(files(k),'coho');
    b=NaN(numel(q.Epoch),1); magnitudeVariable="missing";
    if isfield(q,'ABS_B') && any(isfinite(q.ABS_B))
        b=q.ABS_B(:); magnitudeVariable="ABS_B";
    elseif isfield(q,'F') && any(isfinite(q.F))
        b=q.F(:); magnitudeVariable="F";
    end
    j=NaN(size(b));
    if isfield(q,'protonFlux1_LECP'), j=q.protonFlux1_LECP(:); end
    part=table(q.Epoch(:),b,j,repmat(k,numel(b),1),(1:numel(b)).', ...
        'VariableNames',{'EpochUTC','B_nT','P1','FileIndex','CDFRecord'});
    use=part.EpochUTC>=firstUTC & part.EpochUTC<lastUTC;
    raw=[raw;part(use,:)]; %#ok<AGROW>
    cohoAudit{k}=struct('File',files(k),'SHA256',string(Case1_File_SHA256(files(k))), ...
        'MagnitudeVariable',magnitudeVariable,'VariableMetadata',q.variable_meta, ...
        'GlobalAttributes',q.global_attributes,'SelectedRecords',nnz(use));
    if mod(k,50)==0 || k==numel(files), fprintf('COHO read %d/%d\n',k,numel(files)); end
end
raw=sortrows(raw,'EpochUTC');
[~,firstRecord,groups]=unique(raw.EpochUTC,'stable');
reference=firstRecord(groups);
assert(isequaln(raw.B_nT,raw.B_nT(reference)) && isequaln(raw.P1,raw.P1(reference)), ...
    'Conflicting duplicate COHO records require review.');
duplicateAudit=raw(~ismember((1:height(raw)).',firstRecord),:);
raw=raw(firstRecord,:);
[meanB,countB]=dailyMean(raw.EpochUTC,raw.B_nT,day);
[meanP1,countP1]=dailyMean(raw.EpochUTC,raw.P1,day);
solarRaw=readmatrix(SolarFile,'Delimiter',';','FileType','text');
solarEpoch=datetime(solarRaw(:,1),solarRaw(:,2),solarRaw(:,3),'TimeZone','UTC');
solarValue=solarRaw(:,5); solarValue(solarValue==-1)=NaN;
assert(numel(unique(solarEpoch))==numel(solarEpoch));
[found,solarRows]=ismember(day,solarEpoch);
sunspot=NaN(nDay,1); sunspot(found)=solarValue(solarRows(found));
daily=table(day+hours(12),sunspot,meanB,meanP1,countB,countP1,solarRows, ...
    'VariableNames',{'EpochUTC','SunspotNumber','B_nT','P1Mean','MAGSampleCount','P1SampleCount','SolarSourceRow'});
% 旧 CSV 只用于核对已交付结果，绝不作为科学计算输入。
baselineFile=fullfile(OriginalDir,'all_daily_values.csv');
assert(isfile(baselineFile),'Previous overview values are needed for preservation check.');
baseline=readtable(baselineFile);
assert(height(baseline)==height(daily),'Baseline time window differs.');
columns={'SunspotNumber','B_nT','P1Mean'};
preservation=table(string(columns(:)),false(3,1),NaN(3,1), ...
    'VariableNames',{'Variable','ValuesMatch','MaximumAbsoluteDifference'});
for k=1:3
    x=daily.(columns{k}); y=baseline.(columns{k});
    assert(isequal(isnan(x),isnan(y)),'Baseline missing-value mask differs.');
    good=isfinite(x)&isfinite(y);
    delta=max(abs(x(good)-y(good)),[],'omitnan');
    if isempty(delta), delta=0; end
    preservation.MaximumAbsoluteDifference(k)=delta;
    preservation.ValuesMatch(k)=all(abs(x(good)-y(good))<=1e-12*max(1,abs(x(good))));
end
assert(all(preservation.ValuesMatch),'First three panels differ from previous output.');
fprintf('First three panels match previous overview values.\n');

%% 官方电子密度 CSV 与配套 PDS 标签：保留原始采样和来源
DensityDir=fullfile(DataDir,'pws','derived','electron_density','native','PDS_release_20260910','data');
densityFile=fullfile(DensityDir,'vg1-vlism-density-2012-2025.csv');
labelFile=fullfile(DensityDir,'vg1-vlism-density-2012-2025.lblx');
label=fileread(labelFile);
md5=regexp(label,'<md5_checksum>([^<]+)</md5_checksum>','tokens','once');
assert(~isempty(md5) && strcmpi(fileMD5(densityFile),md5{1}),'Density CSV MD5 does not match official label.');
options=detectImportOptions(densityFile,'VariableNamingRule','preserve','TextType','string');
options=setvartype(options,options.VariableNames{1},'string');
densityRaw=readtable(densityFile,options);
recordToken=regexp(label,'<records>(\d+)</records>','tokens','once');
assert(height(densityRaw)==str2double(recordToken{1}));
assert(width(densityRaw)==22 && strcmp(densityRaw.Properties.VariableNames{12},'N_e (cm^-3)'));
t=datetime(densityRaw{:,1},'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
ne=densityRaw{:,12}; source=lower(string(densityRaw{:,18}));
assert(all(ismember(source,["epo","qtn"])),'Unexpected density source category.');
density=table(t,ne,densityRaw{:,15},densityRaw{:,16},source,(1:height(densityRaw)).', ...
    'VariableNames',{'EpochUTC','Density_cm3','Minimum_cm3','Maximum_cm3','Source','OriginalCSVRow'});
density=density(density.EpochUTC>=firstUTC & density.EpochUTC<lastUTC,:);
density=sortrows(density,'EpochUTC'); % retain every original measurement
densityAudit=struct('CSVFile',densityFile,'LabelFile',labelFile, ...
    'CSV_SHA256',Case1_File_SHA256(densityFile),'Label_SHA256',Case1_File_SHA256(labelFile), ...
    'OfficialMD5',string(md5{1}),'RawTable',densityRaw,'LabelText',label, ...
    'UserApproval','2026-09-14: official source CSV allowed for electron density.', ...
    'Method','Plot original N_e at SCET separately for EPO/QTN; no averaging, interpolation, density inversion or source preference.');

%% 原始 PWS CDF：16 频道日平均电场
metadataDir=fullfile(DataDir,'pws','source_verification','overview_1990_20250630');
inventory=jsondecode(fileread(fullfile(metadataDir,'selected_official_files.json')));
downloadStatus=jsondecode(fileread(fullfile(metadataDir,'download_summary.json')));
assert(isempty(downloadStatus.Failures),'PWS download inventory contains failures.');
assert(strcmp(downloadStatus.Start,'19900101') && strcmp(downloadStatus.StopInclusive,'20250630'), ...
    'Download inventory window differs from the requested overview.');
mirrorComparison=jsondecode(fileread(fullfile(metadataDir,'SPDF_Iowa_comparison.json')));
assert(mirrorComparison.CompleteMatch,'Official mirror inventories differ; review before claiming completeness.');
pwsFiles=string({inventory.File}).';
assert(numel(pwsFiles)==downloadStatus.SelectedCDFs && all(isfile(pwsFiles)));
[wave,pwsAudit]=V1_Read_PWS_Daily_CDF(pwsFiles,firstUTC,lastUTC,useParallel);
assert(isequal(wave.EpochUTC,daily.EpochUTC));

%% 保存统计与可复现记录
coverage=makeCoverage(daily,density,wave);
method=struct('StartUTC',firstUTC,'EndUTCExclusive',lastUTC, ...
    'FirstThree','SILSO daily source and finite arithmetic means of original COHO hourly ABS_B/F and protonFlux1_LECP; unchanged definitions.', ...
    'PWS','Original calibrated electric_field (V/m), each frequency separately, finite native-record arithmetic mean per UTC day.', ...
    'Missing','Original CDF fill/valid ranges; missing remains NaN. No minimum sample threshold, interpolation, background subtraction or frequency integration.', ...
    'Density','Official original CSV N_e at original SCET, EPO and QTN retained as separate point series. No resampling.', ...
    'Display','P1 log axis; PWS logarithmic frequency and E color scale. Frequency-cell edges are geometric midpoints for display only.', ...
    'EColorLimits_Vm',[5e-7 1e-3], ...
    'ColorLimitsReason','Inherited SCALEMIN/SCALEMAX of representative calibrated PWS source; only display clipping, no removal of values.', ...
    'NoVelocityOrTemperature','User approved omission: V1 PLS stopped working in 1980.', ...
    'EnergyLabel','P1 0.57-1.78 MeV display only; preserve original source flux and metadata.');
writetable(daily,fullfile(OutputDir,'V1_first_three_daily_values.csv'));
writetable(density,fullfile(OutputDir,'V1_PWS_density_selected_records.csv'));
writetable(wave,fullfile(OutputDir,'V1_PWS_electric_field_daily_values.csv'));
writetable(coverage,fullfile(OutputDir,'coverage.csv'));
writetable(preservation,fullfile(OutputDir,'first_three_preservation.csv'));
frequency=wave.Properties.UserData.Frequency_Hz;
channelMap=table((1:16).',frequency(:),repmat("V/m",16,1), ...
    'VariableNames',{'Channel','Frequency_Hz','ElectricFieldUnit'});
writetable(channelMap,fullfile(OutputDir,'PWS_frequency_channels.csv'));
missingDays=wave(~wave.OfficialCDFPresent | ~any(isfinite(wave.EMean_Vm),2),{'EpochUTC','OfficialCDFPresent'});
writetable(missingDays,fullfile(OutputDir,'PWS_missing_days.csv'));
result=struct('Daily',daily,'Density',density,'Wave',wave,'Coverage',coverage, ...
    'Method',method,'COHOSources',{cohoAudit},'PWSSources',{pwsAudit}, ...
    'DensityAudit',densityAudit,'SolarFile',SolarFile,'SolarSHA256',Case1_File_SHA256(SolarFile), ...
    'COHOIdenticalDuplicates',duplicateAudit,'DownloadStatus',downloadStatus, ...
    'OfficialMirrorComparison',mirrorComparison,'PWSOfficialInventory',inventory,'ChannelMap',channelMap, ...
    'Preservation',preservation,'OriginalFilesBefore',before,'CreatedUTC',datetime('now','TimeZone','UTC'));
files=drawFigure(daily,density,wave,firstUTC,lastUTC,OutputDir,method.EColorLimits_Vm);
after=originalHashes(OriginalDir);
assert(isequal(before,after),'An original figure changed unexpectedly.');
result.OutputFiles=files; result.OriginalFilesAfter=after;
result.MainCodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
result.PWSReaderSHA256=Case1_File_SHA256(fullfile(CodeDir,'V1_Read_PWS_Daily_CDF.m'));
result.MATLABVersion=version;
dependencies=string({which('Voyager_Read_CDF_Product'),which('dataobj'),which('EpochTT'),fullfile(CodeDir,'Sync_V1_PWS_Archive.ps1'),which('V1_Format_Plasma_Wave_Figure')}).';
dependencyHash=strings(size(dependencies));
for k=1:numel(dependencies), dependencyHash(k)=Case1_File_SHA256(dependencies(k)); end
result.Dependencies=table(dependencies,dependencyHash,'VariableNames',{'File','SHA256'});
save(fullfile(OutputDir,'V1_plasma_wave_overview_audit.mat'),'result','-v7');
disp(coverage); fprintf('Saved new five-panel overview: %s\n',files(1));
end

function [avg,count]=dailyMean(time,value,day)
[matched,bin]=ismember(dateshift(time,'start','day'),day);
good=matched & isfinite(value); n=numel(day);
avg=NaN(n,1);count=zeros(n,1);
if any(good)
    avg=accumarray(bin(good),value(good),[n 1],@mean,NaN);
    count=accumarray(bin(good),1,[n 1],@sum,0);
end
assert(sum(count)==nnz(good));
end

function files=cohoFiles(folder,firstUTC,lastUTC)
listing=dir(fullfile(folder,'**','voyager1_coho1hr_merged_mag_plasma_*_v*.cdf'));
paths=sort(string(fullfile({listing.folder},{listing.name})).');
keys=strings(size(paths)); years=NaN(size(paths));
for k=1:numel(paths)
    token=regexp(char(paths(k)),'_(\d{8})_v','tokens','once');
    keys(k)=string(token{1}); years(k)=str2double(token{1}(1:4));
end
use=years>=year(firstUTC) & years<=year(lastUTC);
paths=paths(use);keys=keys(use);
[~,keep]=unique(keys,'last');files=paths(keep);
end

function h=originalHashes(folder)
names=["V1_all_daily_5panels.png";"V1_all_daily_5panels.pdf";"V1_all_daily_5panels.fig";"all_daily_values.csv"];
h=table(names,strings(4,1),'VariableNames',{'File','SHA256'});
for k=1:4, h.SHA256(k)=string(Case1_File_SHA256(fullfile(folder,names(k)))); end
end

function value=fileMD5(file)
fid=fopen(file,'rb');assert(fid>=0);cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('MD5');
digest.update(typecast(bytes,'int8'));
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function coverage=makeCoverage(daily,density,wave)
name=["Sunspot";"MAG";"P1";"PWS density";"PWS electric field (any channel)"];
coverage=table(name,NaT(5,1,'TimeZone','UTC'),NaT(5,1,'TimeZone','UTC'),zeros(5,1), ...
    'VariableNames',{'Parameter','FirstValidUTC','LastValidUTC','ValidDays'});
times={daily.EpochUTC,daily.EpochUTC,daily.EpochUTC,density.EpochUTC,wave.EpochUTC};
good={isfinite(daily.SunspotNumber),isfinite(daily.B_nT),isfinite(daily.P1Mean), ...
    isfinite(density.Density_cm3),any(isfinite(wave.EMean_Vm),2)};
for k=1:5
    t=times{k}(good{k});
    if ~isempty(t)
        coverage.FirstValidUTC(k)=min(t); coverage.LastValidUTC(k)=max(t);
        coverage.ValidDays(k)=numel(unique(dateshift(t,'start','day')));
    end
end
end

function files=drawFigure(daily,density,wave,firstUTC,lastUTC,out,colorLimits)
%% 五个对齐 panel，保留原图前三条日均曲线的细线、配色和单位
f=figure('Color','w','Position',[50 30 1600 1400],'Visible','off');
ax=gobjects(5,1);
for k=1:5
    ax(k)=axes(f,'Position',[0.085 0.065+(5-k)*0.18 0.83 0.16]);
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(k),datenum([firstUTC lastUTC]));
end
time=datenum(daily.EpochUTC);
values=[daily.SunspotNumber daily.B_nT daily.P1Mean];
colors=[.70 .30 .06;.15 .15 .15;.48 .12 .62];
labels={{'Daily sunspot number'},{'|B| daily mean','(nT)'}, ...
    {'P1 daily mean','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}};
for k=1:3
    y=values(:,k); if k==3, y(y<=0)=NaN; end
    h=plot(ax(k),time,y,'.-','Color',colors(k,:),'LineWidth',.4,'MarkerSize',3);
    assert(isequaln(h.YData(:),y));
    if k==3, set(ax(k),'YScale','log'); end
    ylabel(ax(k),labels{k});
end
hold(ax(4),'on');
sourceNames=["epo","qtn"]; sourceColors=[.08 .36 .62;.05 .55 .32];
for k=1:2
    use=density.Source==sourceNames(k)&isfinite(density.Density_cm3);
    plot(ax(4),datenum(density.EpochUTC(use)),density.Density_cm3(use),'.', ...
        'Color',sourceColors(k,:),'MarkerSize',5,'DisplayName',upper(sourceNames(k)));
end
ylabel(ax(4),{'PWS electron density','(cm^{-3})'});
legend(ax(4),'Location','northwest','Orientation','horizontal','Box','off');
freq=wave.Properties.UserData.Frequency_Hz;
logf=log10(freq);
frequencyEdges=10.^[logf(1)-(logf(2)-logf(1))/2, ...
    (logf(1:end-1)+logf(2:end))/2,logf(end)+(logf(end)-logf(end-1))/2];
timeEdges=datenum((firstUTC:days(1):lastUTC).');
z=log10(wave.EMean_Vm.');z(~isfinite(z))=NaN;
[X,Y]=meshgrid(timeEdges,frequencyEdges);
C=NaN(size(X));C(1:end-1,1:end-1)=z;
surface(ax(5),X,Y,zeros(size(X)),C,'FaceColor','flat','EdgeColor','none');
view(ax(5),2);set(ax(5),'YScale','log','Layer','top');
ylim(ax(5),[frequencyEdges(1) frequencyEdges(end)]);
yticks(ax(5),[10 100 1000 10000]);
ylabel(ax(5),{'PWS frequency','(Hz)'});
colormap(ax(5),turbo(256));clim(ax(5),log10(colorLimits));
bar=colorbar(ax(5));bar.Position=[.928 .065 .012 .16];
bar.Ticks=-6:-3;bar.TickLabels={'10^{-6}','10^{-5}','10^{-4}','10^{-3}'};
bar.Label.String='Daily mean E (V m^{-1})';bar.FontSize=10;
tickTime=datenum(datetime(1990:5:2025,1,1,'TimeZone','UTC'));
for k=1:5
    set(ax(k),'Position',[.085 .065+(5-k)*.18 .83 .16], ...
        'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',.12);
    xlim(ax(k),datenum([firstUTC lastUTC]));xticks(ax(k),tickTime);
    if k<5, xticklabels(ax(k),{}); else, xticklabels(ax(k),string(1990:5:2025)); end
    text(ax(k),.008,.89,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
end
linkaxes(ax,'x');xlabel(ax(5),'UTC');
sgtitle(f,'Voyager 1 | 1990-01-01 to 2025-06-30 | P1 0.57-1.78 MeV', ...
    'FontSize',15,'FontWeight','bold');
V1_Format_Plasma_Wave_Figure(f);
stem=fullfile(out,'V1_1990_20250630_daily_PWS_5panels');
files=string({[stem,'.png'],[stem,'.pdf'],[stem,'.fig']});
exportgraphics(f,files(1),'Resolution',220);
% PDF 采用同一真实图的光栅输出，避免日频谱 20 万色块造成超大文件。
exportgraphics(f,files(2),'ContentType','image','Resolution',220);
savefig(f,files(3));close(f);
end
