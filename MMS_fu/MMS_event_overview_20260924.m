function MMS_event_overview_20260924(eventNumbers,spacecraft,previewOnly)
% 20个候选事件四星overview：原CDF -> IRFU -> PNG（供PPT插图）。
% 输入：事件编号1:20、卫星1:4；时间UTC，事件前后各10 min。
% 磁场/流速及MMS1-3电场用GSM；MMS4电场仅原生DSL XY。n:cm^-3, T:eV, V:km/s, E:mV/m, B:nT。
% burst按各变量真实记录优先；其余时段显示survey/fast，保留数据缺口。
% 无平滑、无跨缺口插值、无FOTE计算。密集曲线仅作保留极值的显示压缩。
%% 1. 路径与IRFU环境
if nargin<1, eventNumbers=1:15; end
if nargin<2, spacecraft=1:4; end
if nargin<3, previewOnly=false; end
codeRoot=fileparts(mfilename('fullpath'));
irfu='C:\Users\Administrator\Documents\irfu-matlab-master';
dataRoot='Z:\SPART-WORK\Data\MMS';
outRoot='C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924';
workRoot=fullfile(tempdir,'MMS_overview_20260924');
auditRoot=fullfile(dataRoot,'derived','event_overview_20260924');
if previewOnly,outRoot=fullfile(workRoot,'preview');end
for p={irfu,fullfile(irfu,'irf'),fullfile(irfu,'plots'),fullfile(irfu,'mission','mms'), ...
        fullfile(irfu,'mission','cluster'),fullfile(irfu,'contrib','nasa_cdf_patch')}
    addpath(p{1});
end
setenv('CDF_LEAPSECONDSTABLE',fullfile(irfu,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
if ~isfolder(outRoot), mkdir(outRoot); end
if ~isfolder(fullfile(workRoot,'pdf_pages')), mkdir(fullfile(workRoot,'pdf_pages')); end
events=jsondecode(fileread(fullfile(codeRoot,'MMS_event_overview_20260924_events.json')));
manifestPath=fullfile(codeRoot,'MMS_event_overview_20260924_manifest.json');
manifest=jsondecode(fileread(manifestPath));
if ~isfolder(auditRoot),mkdir(auditRoot);end
global MMS_DB
MMS_DB=mms_db;
MMS_DB.add_db(mms_local_file_db([dataRoot filesep]));
MMS_DB.cache.enabled=true;
MMS_DB.cache.cacheSizeMax=4096;
%% 2. 逐事件逐卫星读取；完成记录可用于断点重跑
for ie=eventNumbers
    ev=events(ie);
    eventStart=posixtime(datetime(ev.start,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'",'TimeZone','UTC'));
    eventEnd=posixtime(datetime(ev.end,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'",'TimeZone','UTC'));
    a=eventStart-600; b=eventEnd+600;
    tint=irf.tint([strrep(iso(a),' ','T') 'Z/' strrep(iso(b),' ','T') 'Z']);
    eventDir=fullfile(outRoot,[ev.id '_' ev.start(1:10)]);
    if ~isfolder(eventDir), mkdir(eventDir); end
    for ic=spacecraft
        base=sprintf('%s_%s_MMS%d_overview',ev.id,strrep(ev.start(1:10),'-',''),ic);
        doneFile=fullfile(auditRoot,[base '.json']);
        if previewOnly,doneFile=fullfile(workRoot,[base '_preview.json']);end
        pngFile=fullfile(eventDir,[base '.png']);
        pdfFile=fullfile(workRoot,'pdf_pages',[base '.pdf']);
        if previewOnly,pdfFile=fullfile(workRoot,[base '_preview.pdf']);end
        if ~previewOnly && isfile(doneFile) && isfile(pngFile)
            fprintf('SKIP %s already completed\n',base); continue
        end
        fprintf('START %s\n',base);
        try
            [D,sourceFiles,errors]=readEvent(ic,tint,a,b,ev.id,manifest,previewOnly);
            if ~isempty(errors),error('CDF read error: %s',strjoin(errors,' | '));end
            [fig,stats]=makeFigure(D,ev,ic,a,b,eventStart,eventEnd);
            exportgraphics(fig,pngFile,'Resolution',150,'BackgroundColor','white');
            layoutAudit=getappdata(fig,'layoutAudit');
            close(fig);
            result=struct('event',ev.id,'spacecraft',ic,'startUTC',iso(a),'endUTC',iso(b), ...
                'coordinateSystem',['B,V: GSM; E: ' D.Ecoord],'png',pngFile,'pdfPage',pdfFile,'sources',{sourceFiles}, ...
                'statistics',stats,'readErrors',{errors},'complete',true,'renderVersion',4,'layout',layoutAudit,'plotBase','Overview_download.m + Overview_download_mms4.m');
            fid=fopen(doneFile,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(result));fclose(fid);
            fprintf('DONE %s | panels with records %d/9 | read errors %d\n',base,stats.panelsAvailable,numel(errors));
        catch ME
            fprintf(2,'FAILED %s\n%s\n',base,getReport(ME,'extended','hyperlinks','off'));
            close all;
        end
        MMS_DB.cache.enabled=false; MMS_DB.cache.enabled=true;
    end
end
fprintf('OVERVIEW_BATCH_FINISHED\n');
end

function [D,sources,errors]=readEvent(ic,tint,a,b,eventId,manifest,previewOnly)
%% 3. 仪器各自原生时间轴；仅读取已确认下载的产品
names={'B','E','Vi','Ve','Ni','Ne','Ti','Te','Si','Se'};
for k=1:numel(names), D.(names{k})={[],[]}; end
D.Ecoord='GSM';if ic==4,D.Ecoord='DSL XY';end
D.Emeta={struct('records',0,'blocks',zeros(0,2)),struct('records',0,'blocks',zeros(0,2))};
sources={};errors={};
sc=sprintf('mms%d',ic);
for im=1:2
    mode='fast';bmode='srvy';if im==2,mode='brst';bmode='brst';end
    bp=[sc '_fgm_' bmode '_l2'];
    ep=[sc '_edp_' mode '_l2_dce'];
    ip=[sc '_fpi_' mode '_l2_dis-moms'];
    dp=[sc '_fpi_' mode '_l2_des-moms'];
    for product={bp,ep,ip,dp}
        prefix=product{1};
        ix=find(strcmp({manifest.candidate_coverage.event},eventId) & strcmp({manifest.candidate_coverage.product},prefix),1);
        if isempty(ix),continue;end
        candidates=manifest.candidate_coverage(ix).candidate_files;
        if isempty(candidates),continue;end
        if ischar(candidates),candidates={candidates};end
        usable={};
        for jf=1:numel(candidates)
            hit=find(strcmp({manifest.files.file_name},candidates{jf}),1);
            if ~isempty(hit) && isfile(manifest.files(hit).path)
                f=dir(manifest.files(hit).path);
                if f.bytes==manifest.files(hit).file_size,usable{end+1}=manifest.files(hit).path;end %#ok<AGROW>
            end
        end
        if isempty(usable)
            if previewOnly,continue;else,error('Required CDF not downloaded: %s',prefix);end
        end
        if numel(usable)~=numel(candidates)
            error('DownloadIncomplete:product','Incomplete required CDF download: %s',prefix);
        end
        sources=[sources usable]; %#ok<AGROW>
        try
            if strcmp(prefix,bp)
                D.B{im}=getSeries(bp,[sc '_fgm_b_gsm_' bmode '_l2'],tint);
                if ~isempty(D.B{im})
                    D.B{im}=D.B{im}(:,1:4);
                    D.B{im}(:,5)=sqrt(sum(D.B{im}(:,2:4).^2,2));
                end
            elseif strcmp(prefix,ep)
                evname=[sc '_edp_dce_gse_' mode '_l2'];if ic==4,evname=[sc '_edp_dce_dsl2d_' mode '_l2'];end
                [D.E{im},D.Emeta{im}]=readElectric(ep,evname,a,b,D.Ecoord);
            else
                particle='dis';key='i';if strcmp(prefix,dp),particle='des';key='e';end
                pre=[sc '_' particle '_'];
                velocity=getSeries(prefix,[pre 'bulkv_gse_' mode],tint);
                if ~isempty(velocity),velocity=irf_gse2gsm(velocity(:,1:4));end
                D.(['V' key]){im}=velocity;
                D.(['N' key]){im}=getSeries(prefix,[pre 'numberdensity_' mode],tint);
                para=getSeries(prefix,[pre 'temppara_' mode],tint);
                perp=getSeries(prefix,[pre 'tempperp_' mode],tint);
                if ~isempty(para) && isequal(size(para),size(perp)) && max(abs(para(:,1)-perp(:,1)))<1e-5
                    D.(['T' key]){im}=[para(:,1),(para(:,2)+2*perp(:,2))/3,para(:,2),perp(:,2)];
                end
                D.(['S' key]){im}=getSpectrum(prefix,[pre 'energyspectr_omni_' mode],tint,a,b);
            end
        catch ME
            errors{end+1}=[prefix ': ' ME.message]; %#ok<AGROW>
            fprintf(2,'READ ERROR %s\n',errors{end});
        end
    end
end
end


function [E,meta]=readElectric(prefix,var,a,b,coord)
% 电场按180秒分段读取，避免多小时burst一次占用大量内存。
% 每段GSM转换后只保留绘图所需的原始极值；原生覆盖与数量单独保存。
parts={};blocks=zeros(0,2);n=0;last=-Inf;
for t0=a:180:b-1e-6
    t1=min(t0+180,b);
    ti=irf.tint([strrep(iso(t0),' ','T') 'Z/' strrep(iso(t1),' ','T') 'Z']);
    A=getSeries(prefix,var,ti);
    if isempty(A),continue;end
    A=A(A(:,1)>last,:);if isempty(A),continue;end
    last=A(end,1);
    if strcmp(coord,'GSM'),A=irf_gse2gsm(A(:,1:4));else,A=A(:,1:3);end
    valid=all(isfinite(A(:,2:end)),2);
    n=n+sum(valid);blocks=[blocks;timeBlocks(A(:,1),valid)]; %#ok<AGROW>
    parts{end+1}=displayCurves(A,max(1200,ceil(30000*(t1-t0)/(b-a)))); %#ok<AGROW>
end
if isempty(parts),E=[];else,E=vertcat(parts{:});end
if ~isempty(blocks)
    blocks=sortrows(blocks,1);united=blocks(1,:);
    for j=2:size(blocks,1)
        if blocks(j,1)<=united(end,2)+1e-6
            united(end,2)=max(united(end,2),blocks(j,2));
        else,united(end+1,:)=blocks(j,:);end %#ok<AGROW>
    end
    blocks=united;
end
meta=struct('records',n,'blocks',blocks);
end

function A=getSeries(prefix,var,tint)
ts=mms.db_get_ts(prefix,var,tint);
if isempty(ts) || isempty(ts.data),A=[];return;end
A=irf.ts2mat(ts);A(:,2:end)=clean(A(:,2:end));
A=A(isfinite(A(:,1)),:);
[~,idx]=unique(A(:,1),'sorted');A=A(idx,:);
end

function S=getSpectrum(prefix,var,tint,a,b)
S=[];v=mms.db_get_variable(prefix,var,tint);
if isempty(v) || ~isfield(v,'data') || isempty(v.data),return;end
ts=mms.variable2ts(v);
t=ts.time.epochUnix; p=clean(double(ts.data));f=double(v.DEPEND_1.data);
keep=t>=a & t<=b;
if ~any(keep),return;end
if isvector(f),f=repmat(f(:)',numel(t),1);end
if size(f,1)~=numel(t) || ~isequal(size(f),size(p))
    error('Energy/time dimension mismatch for %s',var);
end
p(p<=0)=NaN;f(f<=0 | abs(f)>1e29)=NaN;
t=t(keep);p=p(keep,:);f=f(keep,:);
d=median(diff(t),'omitnan');if isempty(d)||~isfinite(d)||d<=0,d=4.5;end
S=struct('t',t,'p',p,'f',f,'dt',d/2,'p_label',{{'log_{10} DEF','keV/(cm^2 s sr keV)'}},'f_label','');
end

function x=clean(x)
x=double(x);x(~isfinite(x)|abs(x)>1e29)=NaN;
end

function [fig,stats]=makeFigure(D,ev,ic,a,b,eventStart,eventEnd)
% 将已核对的原生观测交给用户Overview_download.m复制版绘图。
fields={'B','Vi','Ve','E','Ni','Ne','Ti','Te'};P=struct;counts0=zeros(1,8);
for k=1:numel(fields)
    key=fields{k};parts={};
    for im=1:2
        A=D.(key){im};if isempty(A),continue;end
        if im==1
            if strcmp(key,'E'),A=maskLower(A,D.E{2},D.Emeta{2}.blocks);else,A=maskLower(A,D.(key){2});end
        end
        counts0(k)=counts0(k)+sum(any(isfinite(A(:,2:end)),2));
        if ~strcmp(key,'E'),A=displayCurves(A,18000);end
        parts{end+1}=[A;nan(1,size(A,2))]; %#ok<AGROW>
    end
    if isempty(parts)
        cols=2;if strcmp(key,'B'),cols=5;elseif ismember(key,{'Vi','Ve','E'}),cols=4;end
        if ismember(key,{'Ti','Te'}),cols=4;end
        if strcmp(key,'E')&&ic==4,cols=3;end
        P.(key)=[[a;b],nan(2,cols-1)];
    else,P.(key)=vertcat(parts{:});end
end
counts=[counts0(1:4),counts0(5)+counts0(6),counts0(7),counts0(8),0,0];
for k=1:2
    key='Si';if k==2,key='Se';end
    for im=1:2
        S=D.(key){im};if isempty(S),continue;end
        if im==1&&~isempty(D.(key){2})
            H=D.(key){2};blocks=timeBlocks(H.t,any(isfinite(H.p),2));
            for j=1:size(blocks,1),S.p(S.t>=blocks(j,1)&S.t<=blocks(j,2),:)=NaN;end
            D.(key){im}=S;
        end
        counts(7+k)=counts(7+k)+sum(any(isfinite(S.p),2));
    end
end
coverage=struct;fields={'B','Vi','Ve','E'};
for k=1:4
    meta={};if strcmp(fields{k},'E'),meta=D.Emeta;end
    coverage.(fields{k})=modeStats(D.(fields{k}),a,b,meta);
end
fig=Overview_download_events_20260924(P,D,ev,ic,a,b,eventStart,eventEnd,counts,coverage);
stats=struct('panelsAvailable',sum(counts>0),'panelRecords',counts,'coverage',coverage);
end

function A=maskLower(A,H,blocks)
if isempty(A)||isempty(H),return;end
if nargin<3,blocks=timeBlocks(H(:,1),all(isfinite(H(:,2:end)),2));end
for k=1:size(blocks,1),A(A(:,1)>=blocks(k,1)&A(:,1)<=blocks(k,2),2:end)=NaN;end
end

function blocks=timeBlocks(t,valid)
t=t(:);d=median(diff(t),'omitnan');if isempty(d)||~isfinite(d)||d<=0,d=1;end
ix=find(valid(:));
if isempty(ix),blocks=zeros(0,2);return;end
cuts=find(diff(t(ix))>3*d | diff(ix)>1);
start=[1;cuts+1];stop=[cuts;numel(ix)];
blocks=[t(ix(start))-d/2,t(ix(stop))+d/2];
end

function A=displayCurves(A,maxPoints)
% 先按真实时间缺口分段，再选每块每分量的原始极值点；不做平均。
if isempty(A),return;end
t=A(:,1);d=median(diff(t),'omitnan');if isempty(d)||~isfinite(d)||d<=0,d=1;end
cuts=find(diff(t)>3*d | diff(t)<=0);
starts=[1;cuts+1];ends=[cuts;size(A,1)];out=cell(numel(starts),1);
for seg=1:numel(starts)
    x=A(starts(seg):ends(seg),:);
    limit=max(200,ceil(maxPoints*size(x,1)/size(A,1)));
    if size(x,1)>limit
        block=max(2,ceil(size(x,1)/(limit/(2*(size(x,2)-1)+2))));
        gapEdges=find(any(diff(isfinite(x(:,2:end)),1,1)~=0,2));
        keep=[1;size(x,1);gapEdges;gapEdges+1];
        for q=1:block:size(x,1)
            rows=q:min(q+block-1,size(x,1));
            [~,lo]=min(x(rows,2:end),[],1,'omitnan');[~,hi]=max(x(rows,2:end),[],1,'omitnan');
            keep=[keep;rows(1);rows(end);reshape(rows(lo),[],1);reshape(rows(hi),[],1)]; %#ok<AGROW>
        end
        x=x(unique(keep),:);
    end
    out{seg}=[x;nan(1,size(x,2))];
end
A=vertcat(out{:});
end

function stats=modeStats(pair,a,b,meta)
if nargin<4,meta={};end
stats=struct;
for im=1:2
    key='surveyFast';if im==2,key='burst';end
    A=pair{im};
    if isempty(A),stats.(key)=struct('records',0,'seconds',0,'intervalsUTC',{{}},'blocks',zeros(0,2));continue;end
    blocks=timeBlocks(A(:,1),all(isfinite(A(:,2:end)),2));
    nNative=sum(all(isfinite(A(:,2:end)),2));
    if ~isempty(meta),blocks=meta{im}.blocks;nNative=meta{im}.records;end
    blocks(:,1)=max(blocks(:,1),a);blocks(:,2)=min(blocks(:,2),b);
    times=cell(size(blocks,1),2);
    for j=1:size(blocks,1),times{j,1}=iso(blocks(j,1));times{j,2}=iso(blocks(j,2));end
    stats.(key)=struct('records',nNative, ...
        'seconds',sum(max(0,blocks(:,2)-blocks(:,1))),'intervalsUTC',{times},'blocks',blocks);
end
end

function s=iso(t)
s=char(datetime(t,'ConvertFrom','posixtime','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'));
end
