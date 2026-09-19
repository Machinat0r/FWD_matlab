function result = Run_V1_Density_Multienergy_Overview
% V1 日均多能道概览：LECP 3 能道、CRS 15 能道，保留原生 PWS 密度点。
% 所有磁场/粒子科学值直接读取原始 COHO CDF；密度使用已获批官方 CSV。
% 不计算波动电场。原五面板图和原始数据保持原样。
%% 路径与时间
CodeDir=fileparts(mfilename('fullpath'));
ProjectCode=fileparts(CodeDir);
DataDir='Z:/SPART-WORK/Data/Voyager/voyager1/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Density_Multienergy_Overview/';
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
OriginalDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V1/';
PreviousDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Plasma_Wave_Overview/';
addpath(CodeDir,fullfile(ProjectCode,'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
assert(isfolder(DataDir),'Authorized Z: archive unavailable.');
if ~isfolder(OutputDir), mkdir(OutputDir); end
firstUTC=datetime(1990,1,1,'TimeZone','UTC');
lastUTC=datetime(2025,7,1,'TimeZone','UTC'); % exclusive
day=(firstUTC:days(1):lastUTC-days(1)).';
before=protectedHashes(OriginalDir,PreviousDir);
channel=["protonFlux"+(1:3)+"_LECP","protonFlux"+(1:15)+"_CRS"];
energy=["0.57-1.78","3.4-17.6","22-31","3-4.6","4.6-6.2","6.2-7.7", ...
    "7.7-12.8","12.8-17.9","17.9-30","30-48","48-56","74.471-83.661", ...
    "132.834-154.911","154.911-174.866","174.866-187.713","187.713-220.475", ...
    "220.475-270.050","270.050-346.034"];
instrument=[repmat("LECP",1,3),repmat("CRS",1,15)];
channelMap=table(channel(:),instrument(:),energy(:), ...
    'VariableNames',{'CDFVariable','Instrument','Energy_MeV'});

%% 逐文件读取原始 CDF，保留能段元数据和样本出处
files=cohoFiles(fullfile(DataDir,'coho','1hr','l2','merged_mag_plasma'),firstUTC,lastUTC);
assert(~isempty(files));
parts=cell(numel(files),1); cohoAudit=cell(numel(files),1);
absent=cell(numel(files),1);
for k=1:numel(files)
    q=Voyager_Read_CDF_Product(files(k),'coho');
    b=NaN(numel(q.Epoch),1); magnitudeVariable="missing";
    if isfield(q,'ABS_B') && any(isfinite(q.ABS_B))
        b=q.ABS_B(:); magnitudeVariable="ABS_B";
    elseif isfield(q,'F') && any(isfinite(q.F))
        b=q.F(:); magnitudeVariable="F";
    end
    flux=NaN(numel(b),numel(channel)); absent{k}=strings(0,1);
    for c=1:numel(channel)
        if isfield(q,channel(c))
            flux(:,c)=q.(channel(c))(:);
            meta=q.variable_meta.(channel(c)).attributes;
            limits=regexp(char(meta.CATDESC),'([0-9]+\.?[0-9]*)\s*-\s*([0-9]+\.?[0-9]*)','tokens','once');
            wanted=split(energy(c),'-');
            assert(numel(limits)==2 && isequal(str2double(string(limits(:))),str2double(wanted(:))), ...
                'Energy metadata changed in %s: %s',files(k),channel(c));
            assert(strcmp(strtrim(meta.UNITS),'1/(cm^2 sec ster MeV)'), ...
                'Flux units require review: %s',files(k));
        else
            absent{k}(end+1,1)=channel(c);
        end
    end
    part=table(q.Epoch(:),b,flux,repmat(k,numel(b),1),(1:numel(b)).', ...
        'VariableNames',{'EpochUTC','B_nT','Flux','FileIndex','CDFRecord'});
    use=part.EpochUTC>=firstUTC & part.EpochUTC<lastUTC;
    parts{k}=part(use,:);
    cohoAudit{k}=struct('File',files(k),'SHA256',Case1_File_SHA256(files(k)), ...
        'MagnitudeVariable',magnitudeVariable,'VariableMetadata',q.variable_meta, ...
        'GlobalAttributes',q.global_attributes,'SelectedRecords',nnz(use),'AbsentChannels',absent{k});
    if mod(k,50)==0 || k==numel(files), fprintf('COHO %d/%d\n',k,numel(files)); end
end
raw=sortrows(vertcat(parts{:}),'EpochUTC'); clear parts
[~,firstRecord,groups]=unique(raw.EpochUTC,'stable');
reference=firstRecord(groups);
assert(isequaln(raw.B_nT,raw.B_nT(reference)) && isequaln(raw.Flux,raw.Flux(reference,:)), ...
    'Conflicting duplicate source Epochs require review.');
duplicateAudit=raw(~ismember((1:height(raw)).',firstRecord),:);
raw=raw(firstRecord,:);

%% 各能道分别按 UTC 日对有效值算术平均；零值保留，缺测保留 NaN
[meanB,countB]=dailyMean(raw.EpochUTC,raw.B_nT,day);
meanFlux=NaN(numel(day),numel(channel)); countFlux=zeros(size(meanFlux));
for c=1:numel(channel)
    [meanFlux(:,c),countFlux(:,c)]=dailyMean(raw.EpochUTC,raw.Flux(:,c),day);
end
solarRaw=readmatrix(SolarFile,'Delimiter',';','FileType','text');
solarEpoch=datetime(solarRaw(:,1),solarRaw(:,2),solarRaw(:,3),'TimeZone','UTC');
solarValue=solarRaw(:,5); solarValue(solarValue==-1)=NaN;
assert(numel(unique(solarEpoch))==numel(solarEpoch));
[found,solarRows]=ismember(day,solarEpoch);
sunspot=NaN(numel(day),1); sunspot(found)=solarValue(solarRows(found));
daily=table(day+hours(12),sunspot,meanB,countB,solarRows, ...
    'VariableNames',{'EpochUTC','SunspotNumber','B_nT','MAGSampleCount','SolarSourceRow'});
for c=1:numel(channel)
    daily.(channel(c)+"_Mean")=meanFlux(:,c);
    daily.(channel(c)+"_Count")=countFlux(:,c);
end
% 旧统计只用于验证，未参与上述科学计算。
baseline=readtable(fullfile(OriginalDir,'all_daily_values.csv'));
assert(height(baseline)==height(daily));
newValues=[sunspot meanB meanFlux(:,1)];
oldValues=[baseline.SunspotNumber baseline.B_nT baseline.P1Mean];
assert(isequal(isnan(newValues),isnan(oldValues)));
good=isfinite(newValues)&isfinite(oldValues);
assert(all(abs(newValues(good)-oldValues(good))<=1e-12*max(1,abs(oldValues(good)))));
preservation=table(["SunspotNumber";"B_nT";"P1Mean"], ...
    max(abs(newValues-oldValues),[],1,'omitnan').', ...
    'VariableNames',{'Variable','MaximumAbsoluteDifference'});

%% 密度源文件完整性、版本差异及原图点数审计
[density,densityAudit]=V1_Audit_Density_Completeness(DataDir,PreviousDir,OutputDir,firstUTC,lastUTC);

%% 保存独立结果，不覆盖旧图
coverage=makeCoverage(daily,meanFlux,channelMap,density);
writetable(daily,fullfile(OutputDir,'V1_all_proton_channels_daily.csv'));
writetable(channelMap,fullfile(OutputDir,'proton_channel_metadata.csv'));
writetable(coverage,fullfile(OutputDir,'coverage.csv'));
writetable(preservation,fullfile(OutputDir,'first_three_preservation.csv'));
method=struct('StartUTC',firstUTC,'EndUTCExclusive',lastUTC, ...
    'DailyMean','Finite arithmetic mean of original COHO records within each UTC day, separately for B and each proton channel.', ...
    'CRSCadence','Source metadata describe CRS as 6-hour data embedded in COHO hourly records. Counts describe valid COHO entries, not independent CRS measurements.', ...
    'Missing','Apply original CDF FILLVAL/VALIDMIN/VALIDMAX via existing IRFU reader. No interpolation, zero filling, background subtraction, new quality thresholds or energy combination.', ...
    'Zero','Source-valid zero included in arithmetic mean. Nonpositive daily flux hidden only on logarithmic display, retained in result tables.', ...
    'Density','All current official CSV records at native SCET, EPO/QTN separate points; no averaging, interpolation, old/new version merging or waveform re-inversion.', ...
    'Display','Thin connected daily lines for SSN, B and proton channels. Missing daily values break lines. No electric field panels.', ...
    'EnergyLabel','LECP P1 0.57-1.78 MeV retained; all 18 COHO energy labels and units verified against each input file.');
result=struct('Daily',daily,'RawCOHO',raw,'FluxMean',meanFlux,'FluxCount',countFlux, ...
    'Channels',channelMap,'COHOSources',{cohoAudit},'COHOIdenticalDuplicates',duplicateAudit, ...
    'Density',density,'DensityAudit',densityAudit,'Coverage',coverage, ...
    'Preservation',preservation,'Method',method,'SolarFile',SolarFile, ...
    'SolarSHA256',Case1_File_SHA256(SolarFile),'ProtectedBefore',before, ...
    'MATLABVersion',version,'CreatedUTC',datetime('now','TimeZone','UTC'));
groups={1:3,4:8,9:13,14:18};
names=["LECP_3channels","CRS_01_05","CRS_06_10","CRS_11_15"];
outFiles=strings(0,1);
for g=1:numel(groups)
    figFiles=drawGroup(daily,meanFlux,channelMap,density,groups{g},names(g),firstUTC,lastUTC,OutputDir);
    outFiles=[outFiles;figFiles(:)]; %#ok<AGROW>
end
after=protectedHashes(OriginalDir,PreviousDir);
assert(isequal(before,after),'A protected original file changed.');
result.ProtectedAfter=after; result.OutputFiles=outFiles;
dependency=string({[mfilename('fullpath'),'.m'],which('V1_Audit_Density_Completeness'), ...
    which('Voyager_Read_CDF_Product'),which('dataobj')}).';
hash=strings(size(dependency));
for k=1:numel(dependency), hash(k)=Case1_File_SHA256(dependency(k)); end
result.Code=table(dependency,hash,'VariableNames',{'File','SHA256'});
save(fullfile(OutputDir,'V1_density_multienergy_audit.mat'),'result','-v7');
disp(coverage); fprintf('Completed 4 figures: %s\n',OutputDir);
end

function [avg,count]=dailyMean(time,value,day)
[matched,bin]=ismember(dateshift(time,'start','day'),day);
good=matched & isfinite(value); n=numel(day);
avg=NaN(n,1); count=zeros(n,1);
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

function h=protectedHashes(original,previous)
paths=[fullfile(original,["V1_all_daily_5panels.png";"V1_all_daily_5panels.pdf"; ...
    "V1_all_daily_5panels.fig";"all_daily_values.csv"]); ...
    fullfile(previous,["V1_1990_20250630_daily_PWS_5panels.png"; ...
    "V1_1990_20250630_daily_PWS_5panels.pdf";"V1_1990_20250630_daily_PWS_5panels.fig"])];
hash=strings(size(paths));
for k=1:numel(paths), hash(k)=Case1_File_SHA256(paths(k)); end
h=table(paths,hash,'VariableNames',{'File','SHA256'});
end

function coverage=makeCoverage(daily,flux,channelMap,density)
names=["Sunspot";"MAG";channelMap.CDFVariable;"PWS electron density"];
n=numel(names);
coverage=table(names,NaT(n,1,'TimeZone','UTC'),NaT(n,1,'TimeZone','UTC'),zeros(n,1),zeros(n,1), ...
    'VariableNames',{'Parameter','FirstValidUTC','LastValidUTC','ValidDays','NonpositiveDays'});
v=[daily.SunspotNumber daily.B_nT flux];
for k=1:size(v,2)
    valid=isfinite(v(:,k)); t=daily.EpochUTC(valid);
    if ~isempty(t)
        coverage.FirstValidUTC(k)=min(t);coverage.LastValidUTC(k)=max(t);
        coverage.ValidDays(k)=nnz(valid);
        coverage.NonpositiveDays(k)=nnz(valid & v(:,k)<=0);
    end
end
t=density.EpochUTC(isfinite(density.Density_cm3));
coverage.FirstValidUTC(end)=min(t);coverage.LastValidUTC(end)=max(t);
coverage.ValidDays(end)=numel(unique(dateshift(t,'start','day')));
end

function files=drawGroup(daily,flux,channelMap,density,index,name,firstUTC,lastUTC,out)
%% 与原图相同的上下对齐布局，每个能道独立 logarithmic panel
nPanel=numel(index)+3;
f=figure('Color','w','Position',[30 30 1600 240*nPanel],'Visible','off');
ax=gobjects(nPanel,1); bottom=.05; top=.95; gap=.015;
panelHeight=(top-bottom-(nPanel-1)*gap)/nPanel;
for k=1:nPanel
    ax(k)=axes(f,'Position',[.11 top-k*panelHeight-(k-1)*gap .865 panelHeight],'Tag',sprintf('panel_%d',k));
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',.12);
    xlim(ax(k),datenum([firstUTC lastUTC]));
end
t=datenum(daily.EpochUTC);
plot(ax(1),t,daily.SunspotNumber,'.-','Color',[.70 .30 .06],'LineWidth',.4,'MarkerSize',3,'Tag','Sunspot');
ylabel(ax(1),'Daily sunspot number');
plot(ax(2),t,daily.B_nT,'.-','Color',[.15 .15 .15],'LineWidth',.4,'MarkerSize',3,'Tag','MAG');
ylabel(ax(2),{'|B| daily mean','(nT)'});
colors=[.48 .12 .62;.05 .55 .57;.78 .56 .05;.10 .38 .66;.68 .20 .29];
for k=1:numel(index)
    c=index(k); y=flux(:,c);y(y<=0)=NaN;
    h=plot(ax(k+2),t,y,'.-','Color',colors(k,:),'LineWidth',.4,'MarkerSize',3, ...
        'Tag',char(channelMap.CDFVariable(c)));
    assert(isequaln(h.YData(:),y));
    set(ax(k+2),'YScale','log');
    ylabel(ax(k+2),{sprintf('%s H %s MeV',channelMap.Instrument(c),channelMap.Energy_MeV(c)), ...
        'Daily mean J','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'});
end
hold(ax(end),'on');
sourceNames=["epo","qtn"]; sourceColors=[.08 .36 .62;.05 .55 .32];
for k=1:2
    use=density.Source==sourceNames(k)&isfinite(density.Density_cm3);
    plot(ax(end),datenum(density.EpochUTC(use)),density.Density_cm3(use),'.', ...
        'Color',sourceColors(k,:),'MarkerSize',5,'DisplayName',upper(sourceNames(k)), ...
        'Tag',char(sourceNames(k)));
end
ylabel(ax(end),{'PWS electron density','(cm^{-3})'});
lg=legend(ax(end),'Location','northwest','Orientation','horizontal','Box','off');
lg.Units='normalized'; lg.Position=[.16 bottom+panelHeight-.023 .12 .02];
tickTime=datenum(datetime(1990:5:2025,1,1,'TimeZone','UTC'));
for k=1:nPanel
    xlim(ax(k),datenum([firstUTC lastUTC]));
    xticks(ax(k),tickTime);
    if k<nPanel, xticklabels(ax(k),{}); else, xticklabels(ax(k),string(1990:5:2025)); end
    text(ax(k),.008,.88,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
end
linkaxes(ax,'x');xlabel(ax(end),'UTC');
titleText=sprintf('Voyager 1 | Daily | 1990-01-01 to 2025-06-30 | %s protons',channelMap.Instrument(index(1)));
annotation(f,'textbox',[.11 .958 .865 .03],'String',titleText,'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',15);
stem=fullfile(out,"V1_1990_20250630_daily_"+name);
files=stem+[".png",".pdf",".fig"];
exportgraphics(f,files(1),'Resolution',180);
exportgraphics(f,files(2),'ContentType','image','Resolution',180);
savefig(f,files(3));close(f);
end
