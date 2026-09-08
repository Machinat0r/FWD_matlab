function KH_audit_event_cdf_20260903(catalogPath,outputDir)
% MMS磁层顶KH事件的数据覆盖核查。只检查原始CDF，不进行FOTE计算。
% 输入：事件表(UTF-8 CSV)，输出目录(必须在MMS数据根目录内)。
% 输出：每事件/每卫星/每产品的CDF变量覆盖，以及四星共同覆盖时长。
% 覆盖按CDF原生Epoch及变量有效值计算；不插值，不填补缺口。
% B、Vi、Ve、E四类产品独立报告；共同覆盖不代表FOTE质量合格。
%% 1. 目录和输入
if nargin<1, catalogPath='C:\Users\Administrator\Documents\KH\MMS_KH_published_event_catalog.csv'; end
if nargin<2, outputDir='Z:\SPART-WORK\Data\MMS\derived\KH\catalog_audit_20260903'; end
dataRoot='Z:\SPART-WORK\Data\MMS';
addpath('C:\Users\Administrator\Documents\irfu-matlab-master\contrib\nasa_cdf_patch');
if ~isfolder(outputDir), mkdir(outputDir); end
cd(outputDir); % NASA CDF工具可能使用当前目录暂存
T=readtable(catalogPath,'TextType','string','VariableNamingRule','preserve');
cache=containers.Map('KeyType','char','ValueType','any');
products=["B","Vi","Ve","E"];
entries=struct([]);
summaries=struct([]);
%% 2. 按事件读取四星原始CDF
for ie=1:height(T)
    t1=datenum(char(T.StartUTC(ie)),'yyyy-mm-dd HH:MM:SS');
    t2=datenum(char(T.EndUTC(ie)),'yyyy-mm-dd HH:MM:SS');
    allSegments=cell(4,4);
    for ip=1:4
        for ic=1:4
            [folderParts,varName]=productSpec(products(ip),ic);
            files=eventFiles(dataRoot,folderParts,t1,t2);
            segments=[]; usedFiles={}; failures={}; validRecords=0;
            for jf=1:numel(files)
                fp=files{jf};
                if isKey(cache,fp)
                    rec=cache(fp);
                else
                    rec=readCoverage(fp,varName);
                    cache(fp)=rec;
                end
                if ~isempty(rec.error), failures{end+1}=rec.error; continue; end %#ok<AGROW>
                sg=rec.segments;
                if isempty(sg), continue; end
                sg=sg(sg(:,2)>=t1 & sg(:,1)<=t2,:);
                if isempty(sg), continue; end
                sg(:,1)=max(sg(:,1),t1); sg(:,2)=min(sg(:,2),t2);
                segments=[segments;sg]; %#ok<AGROW>
                usedFiles{end+1}=fp; %#ok<AGROW>
                validRecords=validRecords+rec.validRecords;
            end
            segments=mergeSegments(segments);
            allSegments{ip,ic}=segments;
            k=numel(entries)+1;
            entries(k).EventID=char(T.EventID(ie));
            entries(k).Spacecraft=ic;
            entries(k).Product=char(products(ip));
            entries(k).Variable=varName;
            entries(k).Files=numel(usedFiles);
            entries(k).CoverageSeconds=durationSeconds(segments);
            entries(k).FirstUTC=epochText(firstOrNaN(segments));
            entries(k).LastUTC=epochText(lastOrNaN(segments));
            entries(k).Segments=segments;
            entries(k).FilePaths=usedFiles;
            entries(k).ReadErrors=failures;
        end
    end
    summaries(ie).EventID=char(T.EventID(ie));
    for ip=1:4
        s=allSegments{ip,1};
        for ic=2:4, s=intersectSegments(s,allSegments{ip,ic}); end
        summaries(ie).([char(products(ip)) 'CommonSeconds'])=durationSeconds(s);
        summaries(ie).([char(products(ip)) 'SpacecraftCount'])=sum(cellfun(@(x)~isempty(x),allSegments(ip,:)));
        if ip==1, core=s; end
        if ip==2, core=intersectSegments(core,s); end
    end
    % 共同B+Vi原始有效数据时间段，仅报告数据基础，不增设拓扑判据。
    summaries(ie).BViCommonSeconds=durationSeconds(core);
    summaries(ie).BViCommonIntervals=segmentText(core);
    fprintf('%s: B=%d Vi=%d Ve=%d E=%d; common B+Vi %.1fs\n',...
        T.EventID(ie),summaries(ie).BSpacecraftCount,summaries(ie).ViSpacecraftCount,...
        summaries(ie).VeSpacecraftCount,summaries(ie).ESpacecraftCount,...
        summaries(ie).BViCommonSeconds);
    writeJson(fullfile(outputDir,'cdf_coverage.json'),struct('summaries',summaries,'entries',entries,...
        'complete',ie==height(T),'processed',ie,'total',height(T),...
        'method','原CDF Epoch与指定变量有限有效值；按原生采样周期的3倍分段；不插值，区段相交',...
        'checkedUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'))));
end
end

function [parts,v]=productSpec(p,ic)
sc=sprintf('mms%d',ic);
switch p
case "B", parts={sc,'fgm','brst','l2'}; v=[sc '_fgm_b_gse_brst_l2'];
case "Vi", parts={sc,'fpi','brst','l2','dis-moms'}; v=[sc '_dis_bulkv_gse_brst'];
case "Ve", parts={sc,'fpi','brst','l2','des-moms'}; v=[sc '_des_bulkv_gse_brst'];
case "E", parts={sc,'edp','brst','l2','dce'}; v=[sc '_edp_dce_gse_brst_l2'];
end
end

function files=eventFiles(root,parts,t1,t2)
% 读取各日目录；包含事件开始前最近文件，真实终止时间由CDF确定。
names={}; paths={}; times=[];
for day=floor(t1)-1:floor(t2)
    folder=fullfile(root,parts{:},datestr(day,'yyyy'),datestr(day,'mm'),datestr(day,'dd'));
    if ~isfolder(folder),continue;end
    d=dir(fullfile(folder,'*.cdf'));
    for j=1:numel(d)
        token=regexp(d(j).name,'_(\d{14})_v','tokens','once');
        if isempty(token),continue;end
        names{end+1}=d(j).name; paths{end+1}=fullfile(folder,d(j).name); %#ok<AGROW>
        times(end+1)=datenum(token{1},'yyyymmddHHMMSS'); %#ok<AGROW>
    end
end
if isempty(times),files={};return;end
[~,order]=sort(string(names)); names=names(order);paths=paths(order);times=times(order);
% 同一开始时刻只选最后版本；绝大多数目录已有唯一官方版本。
[~,last]=unique(times,'last');times=times(last);paths=paths(last);
[ts,order]=sort(times);paths=paths(order);
inside=find(ts>=t1 & ts<=t2);
prior=find(ts<t1,1,'last');
files=paths(unique([prior inside]));
end

function rec=readCoverage(path,varName)
rec=struct('segments',[],'validRecords',0,'error','');
try
    epochName='Epoch'; if contains(varName,'_edp_'),epochName=[varName(1:4) '_edp_epoch_brst_l2'];end
    info=spdfcdfinfo(path,'Variables',{epochName,varName});
    if size(info.Variables,1)<2 || isempty(info.Variables{2,1})
        rec.error=['missing variable: ' path ' : ' varName];return;
    end
    n=info.Variables{1,3};
    if n<1,return;end
    % 原生Epoch和向量直接读取，保存的仅是覆盖元数据。
    raw=spdfcdfread(path,'Variables',{epochName,varName},'CombineRecords',true);
    tt=raw{1}; vv=double(raw{2});tt=double(tt(:));
    if size(vv,1)~=numel(tt),vv=vv.';end
    if size(vv,1)~=numel(tt),error('Epoch/vector row mismatch');end
    valid=all(isfinite(vv),2) & all(abs(vv)<1e29,2) & isfinite(tt);
    delta=diff(tt);delta=delta(isfinite(delta)&delta>0);
    if isempty(delta),return;end
    dt=median(delta);
    % NaN/填充值记录及超过3个原生周期的采样间断切断区段。
    indices=find(valid);
    if isempty(indices),return;end
    breaks=[true;diff(indices)>1 | diff(tt(indices))>3*dt];
    starts=find(breaks);ends=[starts(2:end)-1;numel(indices)];
    rec.segments=[tt(indices(starts)),tt(indices(ends))];
    rec.validRecords=numel(indices);
catch ME
    rec.error=[path ' : ' ME.message];
end
end

function s=mergeSegments(s)
if isempty(s),s=zeros(0,2);return;end
s=sortrows(s,1);out=s(1,:);
for k=2:size(s,1)
    if s(k,1)<=out(end,2),out(end,2)=max(out(end,2),s(k,2));else,out(end+1,:)=s(k,:);end %#ok<AGROW>
end
s=out;
end
function out=intersectSegments(a,b)
out=zeros(0,2);i=1;j=1;
while i<=size(a,1)&&j<=size(b,1)
    x=max(a(i,1),b(j,1));y=min(a(i,2),b(j,2));
    if y>x,out(end+1,:)=[x y];end %#ok<AGROW>
    if a(i,2)<b(j,2),i=i+1;else,j=j+1;end
end
end
function d=durationSeconds(s)
d=sum(max(0,s(:,2)-s(:,1)))*86400;
end
function t=firstOrNaN(s)
if isempty(s),t=NaN;else,t=s(1,1);end
end
function t=lastOrNaN(s)
if isempty(s),t=NaN;else,t=s(end,2);end
end
function t=epochText(x)
if ~isfinite(x),t='';else,t=datestr(x,'yyyy-mm-dd HH:MM:SS.FFF');end
end
function text=segmentText(s)
text=cell(size(s,1),2);
for k=1:size(s,1),text{k,1}=epochText(s(k,1));text{k,2}=epochText(s(k,2));end
end
function writeJson(path,value)
fid=fopen(path,'w','n','UTF-8');c=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(value));
end

