%% MMS1磁尾overview：每4h两张图，2026-07-20至2026-08-15
% 仿照原Overview_download.m的直接脚本 / %%分区，复用此前overview绘图段。
% 图1：B、Vi、AE；图2：B、Vi、AE、Ni、离子全向能谱。
% 原始MMS CDF及CDAWeb OMNI 1min CDF直接由IRFU读取；无新增function。
% 时间UTC，矢量GSM；标题为时段中点位置，RE=6372km。
% burst优先，survey/fast补缺；AE保留原生1min采样，不平滑、不插值。
% 可在调用前设置WindowList=[1 2 ...]或SkipExisting=false。
close all
clearvars -except WindowList SkipExisting
clc
%% 路径及IRFU初始化
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS\';
RecordDir=[ParentDir 'derived\MMS1_tail_20260720_0815'];
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815';
DataDir=[ParentDir 'ancillary\omni\hro_1min\2026'];
addpath(CodeDir,IRFDir);irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
mms.db_init('local_file_db',ParentDir);
mms.db_init('db_cache_size_max',2048);
mms.db_init('db_cache_enabled',true);
DownloadState=jsondecode(fileread(fullfile(RecordDir,'science_complete.json')));
assert(DownloadState.complete,'请先补齐本批清单中的原始CDF，再运行绘图');
Windows=jsondecode(fileread(fullfile(RecordDir,'windows.json')));
if ~exist('WindowList','var'),WindowList=1:numel(Windows);end
if ~exist('SkipExisting','var'),SkipExisting=true;end
Folders={'B_Vi_AE','B_Vi_AE_Ni_Ei'};
for k=1:2
    if ~isfolder(fullfile(OutputDir,Folders{k})),mkdir(fullfile(OutputDir,Folders{k}));end
end
set(groot,'defaultFigureVisible','off');
%% AE：CDAWeb OMNI_HRO_1MIN原始CDF，AE_INDEX [nT]
Files={'omni_hro_1min_20260701_v01.cdf','omni_hro_1min_20260801_v01.cdf'};
SourceRoot='https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/2026/';
AE_all=[];AE_sources=cell(1,2);
for im=1:2
    FilePath=fullfile(DataDir,Files{im});URL=[SourceRoot Files{im}];
    if ~isfile(FilePath)
        websave([FilePath '.part'],URL,weboptions('Timeout',120));
        movefile([FilePath '.part'],FilePath);
    end
    dobj=dataobj(FilePath);
    AE_ts=get_ts(dobj,'AE_INDEX');A=irf.ts2mat(AE_ts);
    Raw=get_variable(dobj,'AE_INDEX');Fill=getfillval(dobj,'AE_INDEX');
    Expected=double(Raw.data);Expected(Expected==double(Fill))=NaN;
    % 整数型CDF变量：按原始FILLVAL在double数组中保留缺测。
    A(double(Raw.data)==double(Fill),2)=NaN;
    assert(isequaln(A(:,2),Expected(:)),'AE与原CDF不一致');
    assert(all(abs(diff(A(:,1))-60)<1e-3),'AE时间间隔非1min');
    AE_all=[AE_all;A]; %#ok<AGROW>
    AE_sources{im}=struct('path',FilePath,'url',URL,'variable','AE_INDEX','fillValue',double(Fill));
end
assert(all(abs(diff(AE_all(:,1))-60)<1e-3));
%% 分段、读取及绘图
ic=1;
for iw=WindowList
    Event=Windows(iw);RecordFile=fullfile(RecordDir,[Event.id '_overview.json']);
    if SkipExisting&&isfile(RecordFile)
        Previous=jsondecode(fileread(RecordFile));
        if Previous.complete&&all(startsWith(string(Previous.png),[OutputDir filesep]))&&all(isfile(string(Previous.png)))
            fprintf('SKIP %s already complete\n',Event.id);continue;
        end
    end
    tint=irf.tint([Event.startUTC '/' Event.endUTC]);
    tStart=tint.start.epochUnix;tStop=tint.stop.epochUnix;
    AE=irf_tlim(AE_all,tint.epochUnix);
    Result=struct('id',Event.id,'complete',false,'error','','png',{{}},...
        'startUTC',Event.startUTC,'endUTC',Event.endUTC,'spacecraft',1);
    fprintf('START %s %s / %s\n',Event.id,Event.startUTC,Event.endUTC);
    try
        %% load data
        % 保留B、gsmVi、Ni等原变量名称。
        c_eval('B?=cell(1,2);gsmVi?=cell(1,2);',ic);
        c_eval('Ni?=cell(1,2);',ic);
        energy_i=cell(1,2);
        specrec_i=cell(1,2);
        for im=1:2
            mode='brst';bmode='brst';
            if im==2,mode='fast';bmode='srvy';end
            % load B [GSM, nT]
            c_eval('B?_ts=mms.get_data([''B_gsm_'' bmode],tint,?);',ic);
            c_eval('if ~isempty(B?_ts), B?_ts=B?_ts.tlim(tint);Bt?_ts=B?_ts.abs;B?{im}=[irf.ts2mat(B?_ts) double(Bt?_ts.data)];end',ic);
            % load Vi [GSM, km/s]
            c_eval('Vi?_ts=mms.get_data([''Vi_gse_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Vi?_ts),gsmVi?_ts=irf_gse2gsm(Vi?_ts.tlim(tint));gsmVi?{im}=irf.ts2mat(gsmVi?_ts);end',ic);
            % load N [cm^-3]
            c_eval('Ni?_ts=mms.get_data([''Ni_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Ni?_ts),Ni?{im}=irf.ts2mat(Ni?_ts.tlim(tint));end',ic);
            % load FPI omnidirectional energy flux
            c_eval('energy_i{im}=mms.db_get_variable([''mms?_fpi_'' mode ''_l2_dis-moms''],[''mms?_dis_energyspectr_omni_'' mode],tint);',ic);
            for particle={'i'}
                eval(['energy=energy_' particle{1} '{im};']);
                if isempty(energy),continue;end
                % 用现有IRFU转换处理CDF填充值和FPI采样中心。
                energy_ts=mms.variable2ts(energy);
                specrec_p=struct('t',energy_ts.time.epochUnix);
                specrec_p.f=double(energy.DEPEND_1.data);
                specrec_p.p=double(energy_ts.data);
                rows=specrec_p.t>=tStart&specrec_p.t<=tStop;
                if size(specrec_p.f,1)==numel(rows),specrec_p.f=specrec_p.f(rows,:);end
                specrec_p.t=specrec_p.t(rows);specrec_p.p=specrec_p.p(rows,:);
                if isempty(specrec_p.t),continue;end
                specrec_p.p(specrec_p.p<=0)=NaN;
                specrec_p.dt=median(diff(specrec_p.t))/2;
                specrec_p.f_label='';
                specrec_p.p_label={' ','log10(keV/(cm^2 s sr keV))'};
                eval(['specrec_' particle{1} '{im}=specrec_p;']);
            end
        end

        %% burst优先，缺口保持空白
        % 只在有有效burst观测的连续区间隐藏survey/fast。
        % 以下只控制绘图显示；不平滑、不插值、不构造测量值。
        Names={'B','gsmVi','Ni'};
        NativeCounts=zeros(4,2);

        for iv=1:numel(Names)
            c_eval(['A=' Names{iv} '?;'],ic);
            for im=1:2
                if ~isempty(A{im}),NativeCounts(iv,im)=sum(any(isfinite(A{im}(:,2:end)),2));end
            end
            if ~isempty(A{1})
                dt=median(diff(A{1}(:,1)));valid=find(all(isfinite(A{1}(:,2:end)),2));
                cut=find(diff(valid)>1 | diff(A{1}(valid,1))>3*dt);
                first=[1;cut+1];last=[cut;numel(valid)];
                if isempty(valid),first=[];last=[];end
                for k=1:numel(first)
                    lo=A{1}(valid(first(k)),1)-dt/2;hi=A{1}(valid(last(k)),1)+dt/2;
                    if ~isempty(A{2}),A{2}(A{2}(:,1)>=lo&A{2}(:,1)<=hi,2:end)=NaN;end
                end
            end
            % 在实测时间缺口中插入NaN行，避免折线跨缺口连接。
            for im=1:2
                if size(A{im},1)<2,continue;end
                dt=median(diff(A{im}(:,1)));cut=find(diff(A{im}(:,1))>3*dt);
                A{im}=sortrows([A{im};[A{im}(cut,1)+dt nan(numel(cut),size(A{im},2)-1)]],1);
            end
            c_eval([Names{iv} '?=A;'],ic);
        end
        for particle={'i'}
            eval(['S=specrec_' particle{1} ';']);
            row=4;
            for im=1:2
                if ~isempty(S{im}),NativeCounts(row,im)=sum(any(isfinite(S{im}.p),2));end
            end
            if ~isempty(S{1})&&~isempty(S{2})
                dt=2*S{1}.dt;valid=find(any(isfinite(S{1}.p),2));
                cut=find(diff(valid)>1 | diff(S{1}.t(valid))>3*dt);
                first=[1;cut+1];last=[cut;numel(valid)];
                if isempty(valid),first=[];last=[];end
                for k=1:numel(first)
                    lo=S{1}.t(valid(first(k)))-dt/2;hi=S{1}.t(valid(last(k)))+dt/2;
                    S{2}.p(S{2}.t>=lo&S{2}.t<=hi,:)=NaN;
                end
            end
            eval(['specrec_' particle{1} '=S;']);
        end


        %% Init figure：同一段数据分别画3个和5个panel
        for FigureType=1:2
            n=3;if FigureType==2,n=5;end
            set(0,'DefaultAxesFontSize',12);set(0,'DefaultLineLineWidth',0.5);
            fn=figure('Visible','off','Color','w','Position',[10 10 1280 160*n+140]);
            set(fn,'UserData',struct('t_start_epoch',tStart));
            h=gobjects(1,n);
            for i=1:n,h(i)=irf_subplot(n,1,-i);end
            i=1;
        %% B plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=B?{im};',ic);
            if isempty(A),continue;end
            irf_plot([A(:,1) A(:,5)],'reduce','color','k','Linewidth',0.75);
            irf_plot([A(:,1) A(:,2)],'reduce','color','b','Linewidth',0.75);
            irf_plot([A(:,1) A(:,3)],'reduce','color','g','Linewidth',0.75);
            irf_plot([A(:,1) A(:,4)],'reduce','color','r','Linewidth',0.75);
        end
        irf_plot([tStart 0;tStop 0],'k--','Linewidth',0.75);
        grid off;ylabel('B [nT]','fontsize',12);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
        irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
        i=i+1;

        %% Vi plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=gsmVi?{im};',ic);
            if isempty(A),continue;end
            irf_plot([A(:,1) A(:,2)],'reduce','color','b','Linewidth',0.75);
            irf_plot([A(:,1) A(:,3)],'reduce','color','g','Linewidth',0.75);
            irf_plot([A(:,1) A(:,4)],'reduce','color','r','Linewidth',0.75);
        end
        irf_plot([tStart 0;tStop 0],'k--','Linewidth',0.75);
        grid off;ylabel('Vi [km/s]','fontsize',12);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
        irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
        i=i+1;


        %% AE plot：原生1min，缺测保留NaN
        axes(h(i));hold on;
        irf_plot(AE,'color','k','Linewidth',0.9);
        ylabel('AE [nT]','fontsize',12);grid off;
        i=i+1;
        if FigureType==2
        %% N plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=Ni?{im};',ic);
            if ~isempty(A),irf_plot(A,'color','r','Linewidth',0.75);end
        end
        grid off;ylabel('N [cm^{-3}]','fontsize',12);
        set(gca,'ColorOrder',[1 0 0]);
        irf_legend(gca,{'Ni'},[0.97 0.92]);
        i=i+1;

        %% plot ION energy spectrom
        axes(h(i));hold on;colormap(h(i),jet);
        for im=2:-1:1
            specrec_p_i=specrec_i{im};
            if isempty(specrec_p_i),continue;end
            irf_spectrogram(h(i),specrec_p_i,'log','donotshowcolorbar');
        end
        grid off;set(h(i),'yscale','log','ytick',[1e1 1e2 1e3 1e4],'fontsize',12);
        ylabel('Ei(ev)','fontsize',12);set(gca,'Ylim',[20 4e4]);
        if any(NativeCounts(4,:)>0)
            hcb=colorbar(h(i));ylabel(hcb,{'log10(keV/','(cm^2 s sr','keV))'},'fontsize',7);
        end
        i=i+1;


        end
        %% 时间轴和布局：紧凑panel、水平时间标签、无事件竖线
        irf_zoom(h,'x',tint);irf_plot_axis_align(h);
        Labels={'B [nT]','Vi [km/s]','AE [nT]','N [cm^{-3}]','Ei(ev)'};
        Available=[any(NativeCounts(1,:)),any(NativeCounts(2,:)),any(isfinite(AE(:,2))),any(NativeCounts(3,:)),any(NativeCounts(4,:))];
        Available=Available(1:n);
        for i=1:n
            PanelStep=0.81/n;
            set(h(i),'Position',[0.105 0.1+(n-i)*PanelStep 0.77 PanelStep-0.002]);
            set(h(i),'FontName','Arial','FontSize',12,'TickDir','out','Box','on');
            ylabel(h(i),Labels{i},'FontSize',12);
            if ~Available(i)
                cla(h(i));text(h(i),.5,.5,'No available L2 data','Units','normalized',...
                    'HorizontalAlignment','center','Color',[.5 .5 .5]);set(h(i),'YTick',[]);
            end
            if i<n,set(h(i),'XTickLabel',[]);xlabel(h(i),'');end
        end
        cb=findall(fn,'Type','colorbar');
        for k=1:numel(cb)
            pos=cb(k).Axes.Position;cb(k).Position=[0.89 pos(2) 0.009 pos(4)];cb(k).FontSize=10;
            cb(k).Label.FontSize=10;
        end
        irf_zoom(h,'x',tint);irf_zoom(h(1:min(n,4)),'y');
        if FigureType==2,set(h(5),'YLim',[20 4e4]);end
        for i=1:n,set(h(i).YLabel,'Units','normalized','Position',[-0.075 0.5 0]);end
        title(h(1),sprintf('X = %.2f  Y = %.2f  Z = %.2f',Event.midpointXYZ_RE),...
            'FontSize',15,'FontWeight','normal');
        set(h(1).Title,'Units','normalized','Position',[0.5 1.09 0]);
        axes(h(n));set(gca,"XTickLabelRotation",0);
        if ~strcmp(Event.startUTC(1:10),Event.endUTC(1:10))
            xlabel(h(n),[Event.startUTC(1:10) ' / ' Event.endUTC(1:10) ' UTC'],'FontSize',12);
        end
        set(h,'XGrid','off','XMinorGrid','off');drawnow;
        %% 出图保存部分
        set(fn,'Renderer','painters','paperpositionmode','auto');
        PNGFile=fullfile(OutputDir,Folders{FigureType},[Event.id '_MMS1_' Folders{FigureType} '.png']);
        exportgraphics(fn,PNGFile,'Resolution',150,'BackgroundColor','white');
        Result.png{FigureType}=PNGFile;
        close(fn);
        end
        %% 保存本段读取和绘图记录
        Result.coordinateSystem='GSM';Result.midpointXYZ_RE=Event.midpointXYZ_RE;
        Result.nativeCountNames={'B','Vi','Ni','Ei'};Result.nativeCounts=NativeCounts;
        Result.modeOrder={'burst','surveyFast'};Result.panelAvailable=Available;
        Result.aeRows=size(AE,1);Result.aeValid=sum(isfinite(AE(:,2)));
        Result.aeSources=AE_sources;Result.aeCadenceSeconds=60;
        Result.spectrumLimits=[20 4e4];Result.timeLabelRotation=0;
        Result.eventMarkers=false;Result.complete=true;
    catch ME
        Result.error=getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'FAILED %s\n%s\n',Event.id,Result.error);close all
    end
    fid=fopen(RecordFile,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(Result));fclose(fid);
    fprintf('DONE %s complete=%d\n',Event.id,Result.complete);
end
fprintf('TAIL_OVERVIEW_FINISHED\n');

