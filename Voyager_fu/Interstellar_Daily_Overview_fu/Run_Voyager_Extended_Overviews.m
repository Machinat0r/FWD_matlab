function result=Run_Voyager_Extended_Overviews(spacecraft,syncArchive,rangeMode)
% 六图/航天器：2008起及1990起长时段，日/三日/月，原始CDF输入。
if nargin<1, spacecraft=1:2; end
if nargin<2, syncArchive=true; end
if nargin<3, rangeMode='extended'; end
rangeMode=validatestring(rangeMode,{'extended','all','heliopause','boundary_markers'});
%% 路径
CodeDir='C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu/';
DataDir='Z:/SPART-WORK/Data/Voyager/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/';
if strcmp(rangeMode,'heliopause'), OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Heliopause_Overviews/'; end
if strcmp(rangeMode,'boundary_markers')
    OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Boundary_Markers_1990/';
end
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
addpath(CodeDir); addpath(fullfile(fileparts(CodeDir(1:end-1)),'Case1_PPT_VerticalLine_Events_7d_fu'));
cfg=Case1_Config; Case1_Add_IRFU_Path(cfg.IRFURoot);
if ~isfolder(OutputDir), mkdir(OutputDir); end
solarRaw=readmatrix(SolarFile,'Delimiter',';','FileType','text');
solarTime=datetime(solarRaw(:,1),solarRaw(:,2),solarRaw(:,3),'TimeZone','UTC')+hours(12);
solarValue=solarRaw(:,5); solarValue(solarValue==-1)=NaN;
assert(numel(unique(solarTime))==numel(solarTime));
result=struct;
for sc=spacecraft
    out=fullfile(OutputDir,sprintf('V%d',sc)); if ~isfolder(out), mkdir(out); end
    if syncArchive, Voyager_Sync_Full_Overview(DataDir,sc); end
    hp=datetime(2012,8,25,'TimeZone','UTC'); if sc==2, hp=datetime(2018,11,5,'TimeZone','UTC'); end
    startLimit=[]; if strcmp(rangeMode,'heliopause'), startLimit=hp; end
    if strcmp(rangeMode,'boundary_markers'), startLimit=datetime(1990,1,1,'TimeZone','UTC'); end
    [source,hourly]=Voyager_Read_Full_Daily(sc,fullfile(out,'source_recompute'),startLimit);
    sourceAudit=fullfile(source.OutputFolder,'V1_daily_overview_audit.mat');
    if sc==2
        source=Voyager_Supplement_V2_MAG(source);
        sourceAudit=fullfile(source.OutputFolder,'daily_with_MAG_supplement_audit.mat');
    end
    daily=source.Daily;
    [found,index]=ismember(daily.EpochUTC,solarTime);
    solar=NaN(height(daily),1); solar(found)=solarValue(index(found));
    seven=sum(daily.SectorDailyMean(:,1:7),2);
    seven(~all(isfinite(daily.SectorDailyMean(:,1:7)),2))=NaN;
    audit=struct('Spacecraft',sc,'SourceAudit',sourceAudit, ...
        'VoyagerMethod',source.Method,'SolarSource',SolarFile,'SolarSHA256',Case1_File_SHA256(SolarFile), ...
        'SolarRows',index,'SolarOriginalColumns',solarRaw, ...
        'Method','SILSO raw daily values, -1 missing, zero retained. Existing arithmetic daily means; nonoverlapping 3-day windows anchored at plot start; calendar months clipped to interval. Three-day P1 median pools original hourly P1. Monthly P1 median panel averages daily medians. Daily sums include S4; other cadences exclude S4. S8 excluded throughout. No interpolation, threshold, lag or detrending.');
    ranges={'from2008','all'}; modes={'daily','three_day','monthly'};
    if strcmp(rangeMode,'all')
        ranges={'all'};
        % Preserve previously delivered 2008-window audits as output records only.
        % Every newly computed science value above comes directly from raw CDFs.
        priorFile=fullfile(out,'extended_overviews_audit.mat');
        if isfile(priorFile)
            previous=load(priorFile,'audit');
            for q=1:numel(modes)
                key=['from2008_',modes{q}];
                if isfield(previous.audit,key), audit.(key)=previous.audit.(key); end
            end
        end
    end
    if strcmp(rangeMode,'heliopause'), ranges={'heliopause'}; end
    boundaries=[];
    if strcmp(rangeMode,'boundary_markers')
        ranges={'all'}; modes={'daily'};
        boundaries=Voyager_Boundary_Dates(sc);
        audit.BoundaryMarkers=boundaries;
    end
    audit.RequestedRangeMode=rangeMode;
    audit.LongRangeStartUTC=datetime(1990,1,1,'TimeZone','UTC');
    for r=1:numel(ranges)
        first=source.Method.StartUTC;
        if strcmp(ranges{r},'from2008'), first=datetime(2008,1,1,'TimeZone','UTC'); end
        if strcmp(ranges{r},'all'), first=datetime(1990,1,1,'TimeZone','UTC'); end
        if strcmp(ranges{r},'heliopause'), first=hp; end
        last=source.Method.EndUTCExclusive;
        for m=1:numel(modes)
            tag=modes{m};
            if m==1
                nominal=(first:days(1):last-seconds(1)).'; nominalEnd=nominal+days(1);
                sectorSum=seven; sectors=1:7;
            elseif m==2
                nominal=(first:days(3):last-seconds(1)).'; nominalEnd=nominal+days(3);
                sectorSum=daily.P1SixSectorSum; sectors=[1 2 3 5 6 7];
            else
                nominal=(dateshift(first,'start','month'):calmonths(1):dateshift(last-seconds(1),'start','month')).';
                nominalEnd=nominal+calmonths(1); sectorSum=daily.P1SixSectorSum; sectors=[1 2 3 5 6 7];
            end
            begin=nominal; finish=nominalEnd; begin(begin<first)=first; finish(finish>last)=last;
            use=daily.EpochUTC>=first & daily.EpochUTC<last;
            inputs=[solar daily.BMean_nT daily.P1Mean daily.P1Median sectorSum];
            [values,count,bins]=groupValues(daily.EpochUTC(use),inputs(use,:),begin,finish,'mean');
            hourlyCount=[];
            if m==2
                h=hourly.EpochUTC>=first & hourly.EpochUTC<last;
                [values(:,4),hourlyCount]=groupValues(hourly.EpochUTC(h),hourly.P1(h),begin,finish,'median');
            end
            epoch=begin+(finish-begin)/2;
            w=table(epoch,begin,finish,days(finish-begin),begin~=nominal|finish~=nominalEnd, ...
                values(:,1),values(:,2),values(:,3),values(:,4),values(:,5),count, ...
                'VariableNames',{'EpochUTC','StartUTC','EndUTCExclusive','CalendarDays','PartialWindow', ...
                'SunspotNumber','B_nT','P1Mean','P1MedianPanel','SectorSum','ValidDays'});
            if m==2, w.P1MedianHourlySamples=hourlyCount; end
            files=drawOverview(w,sc,tag,ranges{r},sectors,out,boundaries);
            name=[ranges{r},'_',tag];
            writetable(w,fullfile(out,[name,'_values.csv']));
            audit.(name)=struct('Windows',w,'DailyBinIndex',bins,'OutputFiles',files,'Sectors',sectors);
            fprintf('V%d %s %s: %d points; %s to %s\n',sc,ranges{r},tag,height(w),string(first),string(last-days(1)));
        end
    end
    audit.CreatedUTC=datetime('now','TimeZone','UTC');
    audit.CodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
    save(fullfile(out,'extended_overviews_audit.mat'),'audit','-v7.3');
    writetable(source.Coverage,fullfile(out,'Voyager_daily_coverage.csv'));
    result.(sprintf('V%d',sc))=struct('OutputFolder',out,'Coverage',source.Coverage);
end
end

function [values,count,bins]=groupValues(t,x,begin,finish,operation)
% Exact shared half-open bins; no interpolation or minimum coverage.
edges=[begin;finish(end)]; bins=discretize(t,edges);
assert(all(isfinite(bins)) && all(t<finish(end)));
n=numel(begin); values=NaN(n,size(x,2)); count=zeros(n,size(x,2));
for p=1:size(x,2)
    good=isfinite(x(:,p));
    if any(good)
        values(:,p)=accumarray(bins(good),x(good,p),[n 1],str2func(operation),NaN);
        count(:,p)=accumarray(bins(good),1,[n 1],@sum,0);
    end
    assert(sum(count(:,p))==nnz(good));
end
% Independent direct checks include both partial endpoints and interior bins.
for k=unique(round(linspace(1,n,min(n,20))))
    rows=t>=begin(k)&t<finish(k);
    for p=1:size(x,2)
        v=x(rows & isfinite(x(:,p)),p); expected=NaN;
        if ~isempty(v), expected=feval(operation,v); end
        assert(isequaln(expected,values(k,p)));
    end
end
end

function files=drawOverview(w,sc,tag,range,sectors,out,boundaries)
%% 五面板，延续配色、细线及对数粒子轴
period='Daily'; cLabel='P1 daily median'; bLabel='P1 daily mean'; magLabel='|B| daily mean';
if strcmp(tag,'three_day')
    period='3-day'; cLabel='3-day P1 median'; bLabel='3-day mean P1'; magLabel='3-day mean |B|';
elseif strcmp(tag,'monthly')
    period='Monthly'; cLabel='Monthly mean of daily P1 medians'; bLabel='Monthly mean P1'; magLabel='Monthly mean |B|';
end
sectorLabel='P1 sum: S1-S7'; if numel(sectors)==6, sectorLabel='P1 sum: S1,S2,S3,S5,S6,S7'; end
units='(cm^{-2} s^{-1} sr^{-1} MeV^{-1})';
labels={{[period,' sunspot number']},{magLabel,'(nT)'},{bLabel,units},{cLabel,units},{sectorLabel,units}};
colors=[0.70 0.30 0.06;0.15 0.15 0.15;0.48 0.12 0.62;0.48 0.12 0.62;0.08 0.36 0.62];
x=[w.SunspotNumber w.B_nT w.P1Mean w.P1MedianPanel w.SectorSum];
f=figure('Color','w','Position',[50 30 1600 1400],'Visible','off');
layout=tiledlayout(f,5,1,'TileSpacing','compact','Padding','compact'); ax=gobjects(5,1);
for k=1:5
    ax(k)=nexttile(layout); y=x(:,k); if k>=3, y(y<=0)=NaN; end
    h=plot(ax(k),w.EpochUTC,y,'.-','Color',colors(k,:),'LineWidth',0.4,'MarkerSize',3);
    assert(isequaln(h.YData(:),y));
    if k>=3, set(ax(k),'YScale','log'); end
    ylabel(ax(k),labels{k});
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(k),[w.StartUTC(1) w.EndUTCExclusive(end)]);
    text(ax(k),0.008,0.89,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
    if k<5, ax(k).XTickLabel=[]; end
end
if ~isempty(boundaries)
    for k=1:5
        for j=1:2
            hmark=xline(ax(k),boundaries.TimeUTC(j),'--','Color',boundaries.Colors(j,:), ...
                'LineWidth',1.1,'Tag',char(boundaries.ShortName(j)));
            if k==1
                date=boundaries.TimeUTC(j); date.Format='yyyy-MM-dd';
                hmark.Label=sprintf('%s: %s',boundaries.Name(j),string(date));
                hmark.LabelOrientation='horizontal';
                hmark.LabelHorizontalAlignment='left';
                hmark.LabelVerticalAlignment='top';
                hmark.FontSize=10;
            end
        end
    end
end
linkaxes(ax,'x'); xlabel(ax(5),'UTC'); xtickformat(ax(5),'yyyy');
first=w.StartUTC(1); first.Format='yyyy-MM-dd'; last=w.EndUTCExclusive(end)-days(1); last.Format='yyyy-MM-dd';
title(layout,sprintf('Voyager %d | %s | %s to %s | P1 0.57-1.78 MeV',sc,period,string(first),string(last)), ...
    'FontSize',15,'FontWeight','bold');
stem=fullfile(out,sprintf('V%d_%s_%s_5panels',sc,range,tag));
if ~isempty(boundaries), stem=[stem,'_TS_HP']; end
files=string({[stem,'.png'],[stem,'.pdf'],[stem,'.fig']});
exportgraphics(f,files(1),'Resolution',220); exportgraphics(f,files(2),'ContentType','vector'); savefig(f,files(3)); close(f);
end
