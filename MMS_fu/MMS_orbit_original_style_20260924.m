%% 20个事件的四星轨道位置
% 参照原 MMS_orbit.m，直接调用mms.mms4_pl_conf。
% 单时刻事件用该时刻，时间段事件用中点；位置图使用GSM。
close all
clearvars -except EventList
clc
if ~exist('EventList','var'),EventList=1:20;end

%% 路径
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS\';
OutputDir='C:\Users\Administrator\Documents\KH\MMS_events_original_style_20260924\orbits';
RecordDir=[ParentDir 'derived\events_original_style_20260924'];
addpath(CodeDir,IRFDir);
irf('check_path');
mms.db_init('local_file_db',ParentDir);
if ~isfolder(OutputDir),mkdir(OutputDir);end
if ~isfolder(RecordDir),mkdir(RecordDir);end
set(groot,'defaultFigureVisible','off');
Events=jsondecode(fileread(fullfile(CodeDir,'MMS_event_overview_20260924_events.json')));

%% 原IRFU函数的日期显示兼容修正，仅去掉标题中的毫秒
% 原函数中的mm和mmm被时间格式解析器同时识别。备份原文件后只改两处格式字符串。
Source=which('mms.mms4_pl_conf');
s=fileread(Source);
old='utc_yyyy-mm-dd HH:MM:SS.mmm';
new='utc_yyyy-mm-dd HH:MM:SS';
if contains(s,old)
    Backup=fullfile(RecordDir,'mms4_pl_conf_before_date_format_fix.m');
    if ~isfile(Backup),copyfile(Source,Backup);end
    s=strrep(s,old,new);
    fid=fopen(Source,'w','n','UTF-8');fprintf(fid,'%s',s);fclose(fid);
    clear mms.mms4_pl_conf
end

%% 逐事件画图
for ie=EventList
    Event=Events(ie);
    t1=irf_time(Event.start,'utc>epoch');t2=irf_time(Event.end,'utc>epoch');
    t=(t1+t2)/2;
    tint=irf.tint(EpochTT(EpochUnix(t)),1);
    fprintf('START_ORBIT %s\n',Event.id);
    Result=struct('event',Event.id,'complete',false,'error','');
    try
        %% load position
        if ie<=15
            h=mms.mms4_pl_conf(tint);
            SourceType='MMS MEC L2 definitive ephemeris';
        else
            % MEC缺失的事件使用已归档的NASA SSCWeb原始GSE星历。
            % 直接读取JSON字段，不重采样、不另写解析或插值函数。
            File=[ParentDir 'ancillary\sscweb\2026\' Event.id '_gse.json'];
            Raw=jsondecode(fileread(File));
            assert(strcmp(Raw{2}.Result{2}.StatusCode,'SUCCESS'));
            Satellites=Raw{2}.Result{2}.Data{2};
            R=struct;
            for ic=1:4
                Satellite=Satellites{ic}{2};
                assert(strcmp(Satellite.Id,sprintf('mms%d',ic)));
                Coordinates=Satellite.Coordinates{2}{1}{2};
                assert(strcmp(Coordinates.CoordinateSystem,'GSE'));
                Times=Satellite.Time{2};
                NativeTime=zeros(numel(Times),1);
                for it=1:numel(Times)
                    UTC=Times{it}{2};
                    NativeTime(it)=irf_time([UTC(1:23) 'Z'],'utc>epoch');
                end
                if ic==1
                    CommonTime=NativeTime;
                    R.time=EpochTT(EpochUnix(CommonTime));
                else
                    assert(isequal(CommonTime,NativeTime),'SSC four spacecraft timestamps differ');
                end
                c_eval('R.gseR?=[Coordinates.X{2}(:) Coordinates.Y{2}(:) Coordinates.Z{2}(:)];',ic);
            end
            h=mms.mms4_pl_conf(tint,R);
            SourceType='NASA SSCWeb ephemeris (GSE input)';
        end
        h=mms.mms4_pl_conf('gsm');

        %% Init figure
        % 保留原函数全部8个panel，将排版改成两行，便于放入PPT。
        set(gcf,'Position',[0 0 1600 1050]);
        for i=1:8
            col=mod(i-1,4);row=floor((i-1)/4);
            set(h(i),'Units','normalized','Position',[0.055+0.245*col 0.56-0.47*row 0.185 0.36]);
        end
        set(h([1 3 4]),'XLim',[-35 20]);
        set(h(1:7),'FontSize',18);
        set(h(8),'Position',[0.78 0.09 0.21 0.36]);
        tx=findall(h(8),'Type','text');
        for it=1:numel(tx)
            str=string(tx(it).String);
            if any(str==["MMS1","MMS2","MMS3","MMS4"])
                ic=str2double(extractAfter(str,'MMS'));
                tx(it).Position=[0.04+0.55*mod(ic-1,2) 1.0-0.13*floor((ic-1)/2) 0];
                tx(it).FontSize=17;
            elseif contains(str,'MMS configuration')
                tx(it).Position=[0 0.71 0];tx(it).FontSize=17;
            elseif contains(str,'IMF')
                tx(it).Position=[0 0.43 0];tx(it).FontSize=17;
                tx(it).String=strrep(strrep(char(str),',By=',',\newline By='),',Bz=',',\newline Bz=');
            end
        end
        ln=findall(h(8),'Type','line');
        for il=1:numel(ln)
            if numel(ln(il).XData)==1
                ic=round(ln(il).XData/0.27)+1;
                if ic>=1&&ic<=4
                    ln(il).XData=0.55*mod(ic-1,2);ln(il).YData=1.0-0.13*floor((ic-1)/2);
                end
            end
        end

        %% 保存图和图上使用的位置
        D=get(gcf,'UserData');
        PositionKm=zeros(4,3);
        for ic=1:4
            c_eval('Position=irf_resamp(D.r.C?,t);',ic);
            PositionKm(ic,:)=Position(2:4);
        end
        set(gcf,'Renderer','painters');
        set(gcf,'paperpositionmode','auto');
        PNG=fullfile(OutputDir,[Event.id '_MMS1-4_orbit_GSM.png']);
        exportgraphics(gcf,PNG,'Resolution',150,'BackgroundColor','white');
        Result.complete=true;Result.png=PNG;
        Result.positionKm=PositionKm;Result.sourceType=SourceType;
        Result.snapshotUTC=irf_time(t,'epoch>utc');
        Result.plotter='mms.mms4_pl_conf';Result.coordinateSystem='GSM';
        close all;
    catch ME
        Result.error=getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'FAILED_ORBIT %s\n%s\n',Event.id,Result.error);
        close all;
    end
    fid=fopen(fullfile(RecordDir,[Event.id '_orbit.json']),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(Result));fclose(fid);
    fprintf('DONE_ORBIT %s complete=%d\n',Event.id,Result.complete);
end
fprintf('ORBIT_SCRIPT_FINISHED\n');
