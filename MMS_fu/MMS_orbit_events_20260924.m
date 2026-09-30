function MMS_orbit_events_20260924(eventNumbers,force)
% 20个事件的四星位置/构型图；基于用户 新建文件夹/MMS_orbit.m。
% 输入：事件编号1:20。时间段取中点，单一时刻直接使用；UTC，GSM。
% 图形算法使用mms.mms4_pl_conf原函数的兼容副本，位置单位km和RE=6372 km。
% MEC暂缺时读取NASA SSCWeb原始GSE星历，再由原图形函数转换为GSM。
%% 1. 路径与IRFU
if nargin<1,eventNumbers=1:20;end
if nargin<2,force=false;end
codeRoot=fileparts(mfilename('fullpath'));
irfu='C:\Users\Administrator\Documents\irfu-matlab-master';
dataRoot='Z:\SPART-WORK\Data\MMS';
outRoot='C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924\orbits';
auditRoot=fullfile(dataRoot,'derived','event_orbits_20260924');
for p={codeRoot,irfu,fullfile(irfu,'irf'),fullfile(irfu,'plots'),fullfile(irfu,'mission','mms'), ...
        fullfile(irfu,'mission','cluster'),fullfile(irfu,'contrib','nasa_cdf_patch')}
    addpath(p{1});
end
setenv('CDF_LEAPSECONDSTABLE',fullfile(irfu,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
if ~isfolder(outRoot),mkdir(outRoot);end
global MMS_DB
MMS_DB=mms_db;MMS_DB.add_db(mms_local_file_db([dataRoot filesep]));
MMS_DB.cache.enabled=true;
set(groot,'defaultFigureVisible','off');
events=jsondecode(fileread(fullfile(codeRoot,'MMS_event_overview_20260924_events.json')));
mec=jsondecode(fileread(fullfile(auditRoot,'mec_manifest.json')));
ssc=jsondecode(fileread(fullfile(auditRoot,'ssc_manifest.json')));
%% 2. 逐事件调用原轨道绘图
for ie=eventNumbers
    ev=events(ie);
    a=posixtime(datetime(ev.start,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'",'TimeZone','UTC'));
    b=posixtime(datetime(ev.end,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss'Z'",'TimeZone','UTC'));
    t=(a+b)/2;
    tint=irf.tint(EpochTT(EpochUnix(t)),1);
    base=[ev.id '_MMS1-4_orbit_GSM'];
    png=fullfile(outRoot,[base '.png']);
    recordFile=fullfile(auditRoot,[base '.json']);
    if ~force&&isfile(png)&&isfile(recordFile),fprintf('SKIP %s\n',ev.id);continue;end
    fprintf('START_ORBIT %s %s\n',ev.id,utc(t));
    try
        if ie<=15
            % 原MMS_orbit两行调用。副本只修正毫秒时间字符串兼容问题。
            h=MMS_orbit_pl_conf_20260924(tint);
            sourceType='MMS MEC L2 definitive ephemeris';
            files=mec.files(startsWith({mec.files.timetag},ev.start(1:10)));
            sources={files.path};
            sourceURL='https://lasp.colorado.edu/mms/sdc/public/about/how-to/';
        else
            item=ssc(strcmp({ssc.event},ev.id));
            R=readSSC(item.path,t);
            h=MMS_orbit_pl_conf_20260924(tint,R);
            sourceType='NASA SSCWeb ephemeris (GSE input)';
            sources={item.path};sourceURL=item.url;
        end
        h=MMS_orbit_pl_conf_20260924('gsm');
        fig=gcf;
        MMS_orbit_layout_20260924(fig,h);
        % 原程序允许修改坐标轴显示范围；扩展磁尾X范围以容纳所有事件。
        set(h([1 3 4]),'XLim',[-35 20]);
        set(h(1:7),'FontSize',18);
        % 将图上实际使用的位置保存为少量衍生审计数据，供PPT文字引用。
        D=get(fig,'UserData');positionKm=zeros(4,3);
        for ic=1:4
            C=D.r.(['C' num2str(ic)]);
            assert(all(isfinite(C(:)))&&C(1,1)<=t&&C(end,1)>=t,'Invalid position coverage');
            p=irf_resamp(C,t);
            positionKm(ic,:)=p(2:4);
        end
        modelText=findall(fig,'Type','text');
        strings=arrayfun(@(q)string(q.String),modelText,'UniformOutput',false);
        strings=string([strings{:}]);
        modelNotes=strings(contains(strings,'IMF'));
        exportgraphics(fig,png,'Resolution',150,'BackgroundColor','white');
        savefig(fig,fullfile(outRoot,[base '.fig']));
        result=struct('event',ev.id,'eventStartUTC',ev.start,'eventEndUTC',ev.end, ...
            'snapshotUTC',utc(t),'selection','midpoint of supplied interval','coordinateSystem','GSM', ...
            'positionKm',positionKm,'positionRE',positionKm/6372,'earthRadiusKm',6372, ...
            'sourceType',sourceType,'sources',{sources},'sourceURL',sourceURL, ...
            'plotter','MMS_orbit_pl_conf_20260924 (compatibility copy of mms.mms4_pl_conf)', ...
            'modelNotes',{cellstr(modelNotes)},'png',png,'complete',true);
        fid=fopen(recordFile,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(result));fclose(fid);
        close(fig);fprintf('DONE_ORBIT %s\n',ev.id);
    catch ME
        fprintf(2,'FAILED_ORBIT %s\n%s\n',ev.id,getReport(ME,'extended','hyperlinks','off'));
        close all;
    end
end
fprintf('ORBIT_BATCH_FINISHED\n');
end

function R=readSSC(file,t)
%% 3. 读取NASA官方JSON。只在原生样本包围目标时刻时线性对齐。
raw=unwrap(jsondecode(fileread(file)));
assert(strcmp(raw.Result.StatusCode,'SUCCESS'),'SSC status failed');
sats=raw.Result.Data;
common=(t-120:30:t+120)';
R.time=EpochTT(EpochUnix(common));
for ic=1:4
    sat=sats{find(cellfun(@(s)strcmp(s.Id,sprintf('mms%d',ic)),sats),1)};
    coords=sat.Coordinates{1};assert(strcmp(coords.CoordinateSystem,'GSE'),'Expected GSE');
    times=cellfun(@(s)s(1:23),sat.Time,'UniformOutput',false);
    native=posixtime(datetime(times,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS",'TimeZone','UTC'));
    xyz=[coords.X(:) coords.Y(:) coords.Z(:)];
    assert(all(isfinite(xyz(:)))&&native(1)<=common(1)&&native(end)>=common(end),'SSC does not bracket interval');
    assert(all(diff(native)>0)&&max(diff(native))<=90,'SSC native gap');
    mat=irf_resamp([native(:) xyz],common,'linear');
    R.(['gseR' num2str(ic)])=mat(:,2:4);
end
end

function x=unwrap(x)
% SSCWeb的Java类型包装去除，不改变任何数值或时间。
if iscell(x)
    if numel(x)==2&&ischar(x{1})&&(startsWith(x{1},'gov.')||startsWith(x{1},'java'))
        x=unwrap(x{2});
    else
        x=cellfun(@unwrap,x,'UniformOutput',false);
    end
elseif isstruct(x)
    f=fieldnames(x);
    for k=1:numel(f),x.(f{k})=unwrap(x.(f{k}));end
end
end

function s=utc(t)
s=char(datetime(t,'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
end
