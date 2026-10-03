%% 四星磁场curlometer电流：2026-07-21 03:07-03:15 UTC，原始时间分辨率
% 参考用户Plot_TCS.m和原Overview_MMS1_tail_20260720_0815.m，直接脚本/%%分区。
% B/Vi/AE保留原MMS1绘图流程；在Vi与AE之间插入两栏。
% 电流直接截取已验证的四小时原生结果；其来源为四星FGM/MEC及IRFU c_4_j [nA/m^2]。
% 第三栏只画总电流模|J|；第四栏画有符号J_parallel和垂直电流幅值J_perp。
% J_parallel=J dot Bhat，正/负分别对应沿/反沿四星平均磁场。
% J_perp为三分量垂直矢量的模，非负；图例两分量不加绝对值符号。
% 同一时刻四星平均B逐点分解，保留原生磁场时间戳；不平滑、不做时间平均。
% 电流完整原生采样绘图，无reduce；真正数据缺口插入NaN以断开曲线。
% 原B/Vi/AE不做新增时间平均，磁场分解参考始终为逐时刻磁场。
% 既有计算按同一模式将四星磁场/位置以IRFU线性对齐到MMS1时刻；本图不重新计算。
% 对齐限于有限观测连续段，FGM连续性依据CDF bdeltahalf，MEC官方原生30s。
% 数值容差仅为epochUnix浮点精度；不跨观测缺口或在区间外补值。
% 同时有四星burst时优先用burst电流，其余时刻使用四星survey电流。
% 不加入divB/curlB、几何、距离等质量门槛；divB只留作诊断。
% 原B/Vi显示控制段含既有3*dt及reduce，保持与用户所给图一致；不用于新电流对齐。
% AE仍为原CDF AE_INDEX原生1min，FILLVAL保留NaN。
close all
clearvars -except WindowList
clc
%% 路径及IRFU初始化
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS\';
RecordDir=[ParentDir 'derived\MMS_current_20260721_0307_0315_native_20261002'];
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS_current_20260721_0307_0315_native_20261002';
PerpOpacity=0.5; % 垂直电流曲线不透明度：0全透明，1完全不透明。
DataDir=[ParentDir 'ancillary\omni\hro_1min\2026'];
addpath(CodeDir,IRFDir);irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
mms.db_init('local_file_db',ParentDir);
mms.db_init('db_cache_size_max',2048);
mms.db_init('db_cache_enabled',true);
SourceRecordDir=[ParentDir 'derived\MMS_current_20260721_0000_0400_native_20261001'];
DownloadState=jsondecode(fileread(fullfile(SourceRecordDir,'download_complete.json')));
assert(DownloadState.complete,'请先补齐本批四星原CDF');
Windows=jsondecode(fileread(fullfile(RecordDir,'windows.json')));
if ~exist('WindowList','var'),WindowList=1:numel(Windows);end
if ~isfolder(OutputDir),mkdir(OutputDir);end
Manifest=jsondecode(fileread(fullfile(SourceRecordDir,'manifest.json')));
% 四小时结果的原CDF已齐备，直接复用本地数据，不重新下载。
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
%% 逐事件读取和计算
for iw=WindowList
    Event=Windows(iw);
    CurrentLegendPosition=[0.97 0.92];
    tint=irf.tint([Event.startUTC '/' Event.endUTC]);
    tStart=tint.start.epochUnix;tStop=tint.stop.epochUnix;
    tintB=irf.tint(tint.start+(-1),tint.stop+1);
    tintR=irf.tint(tint.start+(-60),tint.stop+60);
    AE=irf_tlim(AE_all,tint.epochUnix);
    fprintf('START %s %s / %s\n',Event.id,Event.startUTC,Event.endUTC);
    Result=struct('id',Event.id,'complete',false,'error','','startUTC',Event.startUTC,...
        'endUTC',Event.endUTC,'coordinateSystem','GSM','currentUnit','nA/m^2');
    try
        %% load data：原有MMS1 B和Vi，独立保留原始时间轴
        ic=1;Bplot1=cell(1,2);gsmVi1=cell(1,2);
        for im=1:2
            mode='brst';bmode='brst';
            if im==2,mode='fast';bmode='srvy';end
            Bplot1_ts=mms.get_data(['B_gsm_' bmode],tint,1);
            if ~isempty(Bplot1_ts)
                Bplot1_ts=Bplot1_ts.tlim(tint);Bt1_ts=Bplot1_ts.abs;
                Bplot1{im}=[irf.ts2mat(Bplot1_ts) double(Bt1_ts.data)];
            end
            Vi1_ts=mms.get_data(['Vi_gse_fpi_' mode '_l2'],tint,1);
            if ~isempty(Vi1_ts)
                gsmVi1_ts=irf_gse2gsm(Vi1_ts.tlim(tint));
                gsmVi1{im}=irf.ts2mat(gsmVi1_ts);
            end
        end
        %% load R1-R4：原生MEC位置，只用于本图中点位置标题
        Rnative=cell(1,4);RSources=cell(1,4);
        for ic=1:4
            prefix=sprintf('mms%d_mec_srvy_l2_epht89d',ic);
            varname=sprintf('mms%d_mec_r_gsm',ic);
            R_ts=mms.db_get_ts(prefix,varname,tintR);
            assert(~isempty(R_ts),'缺少MMS%d MEC位置',ic);
            Rnative{ic}=irf.ts2mat(R_ts);
            RSources{ic}=mms.db_list_files(prefix,tintR);
        end
        %% 时段中点MMS1位置，仅用于原overview位置标题
        Rmid=irf_resamp(Rnative{1},(tStart+tStop)/2,'linear');
        Event.midpointXYZ_RE=Rmid(1,2:4)/6372;
        %% 原生电流直接截取：从已保存四小时结果选取03:07-03:15，不重新估计
        Previous=load(fullfile(SourceRecordDir,'EV01_20260721_0000_0400_native_current.mat'),...
            'J','Jabs','Jpar','Jperp','ReferenceB','CurrentHalf');
        SourceInfo=jsondecode(fileread(fullfile(SourceRecordDir,'EV01_20260721_0000_0400_current.json')));
        assert(SourceInfo.complete&&~SourceInfo.timeAveraging&&~SourceInfo.fixedBackgroundField);
        J=cell(1,2);Jabs=cell(1,2);Jpar=cell(1,2);Jperp=cell(1,2);
        ReferenceB=cell(1,2);CurrentHalf=cell(1,2);CurrentInfo=cell(1,2);
        ReferenceVerified=SourceInfo.instantaneousReferenceVerified;
        for im=1:2
            mode='brst';if im==2,mode='srvy';end
            CurrentInfo{im}=struct('mode',mode,'finiteCurrentRows',0);
            if isempty(Previous.J{im}),continue;end
            keep=Previous.J{im}(:,1)>=tStart&Previous.J{im}(:,1)<=tStop;
            J{im}=Previous.J{im}(keep,:);Jabs{im}=Previous.Jabs{im}(keep);
            Jpar{im}=Previous.Jpar{im}(keep);Jperp{im}=Previous.Jperp{im}(keep);
            ReferenceB{im}=Previous.ReferenceB{im}(keep,:);
            CurrentHalf{im}=Previous.CurrentHalf{im}(keep);
            assert(isequal(J{im}(:,1),ReferenceB{im}(:,1)),'电流和实时参考磁场时刻不一致');
            CurrentInfo{im}.finiteCurrentRows=sum(all(isfinite(J{im}(:,2:4)),2));
            CurrentInfo{im}.alignedRows=size(J{im},1);
            CurrentInfo{im}.sources=SourceInfo.current(im).sources;
            fprintf('EXTRACTED %s %s finite=%d rows=%d\n',Event.id,mode,...
                CurrentInfo{im}.finiteCurrentRows,size(J{im},1));
        end
        clear Previous
        %% 原始分辨率电流：保留时间戳和数值，缺口仅在显示数组加NaN断线
        CurrentPlot=cell(1,2);NativeSummary=cell(1,2);
        for im=1:2
            if isempty(J{im}),continue;end
            A=[J{im}(:,1) Jabs{im} Jpar{im} Jperp{im}];
            good=all(isfinite(A(:,2:4)),2);
            assert(all(diff(A(:,1))>0),'原始电流时间轴不递增');
            width=CurrentHalf{im};
            cut=find(diff(A(:,1))>width(1:end-1)+width(2:end)+8*eps(max(abs(A(:,1)))));
            GapRows=[(A(cut,1)+A(cut+1,1))/2 nan(numel(cut),3)];
            CurrentPlot{im}=sortrows([A;GapRows],1);
            mode='brst';if im==2,mode='srvy';end
            NativeSummary{im}=struct('mode',mode,'selectedRows',sum(good),...
                'medianCadenceSeconds',median(diff(A(:,1))),'displayGapRows',numel(cut),...
                'totalRange_nA_m2',[min(A(good,2)) max(A(good,2))],...
                'parallelRange_nA_m2',[min(A(good,3)) max(A(good,3))],...
                'perpendicularRange_nA_m2',[min(A(good,4)) max(A(good,4))]);
            assert(isequaln(CurrentPlot{im}(ismember(CurrentPlot{im}(:,1),A(:,1)),:),A),...
                '显示断线改变了原生电流数值');
            fprintf('NATIVE %s %s selected=%d dt=%.9f s\n',Event.id,mode,sum(good),median(diff(A(:,1))));
        end
        ic=1;
        %% burst优先，缺口保持空白
        % 只在有有效burst观测的连续区间隐藏survey/fast。
        % 以下只控制绘图显示；不平滑、不插值、不构造测量值。
        Names={'Bplot','gsmVi'};
        NativeCounts=zeros(2,2);

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
        %% Init figure：B/Vi/|J|/J_parallel和J_perp/AE，共5个panel
        n=5;set(0,'DefaultAxesFontSize',12);set(0,'DefaultLineLineWidth',0.5);
        fn=figure('Visible','off','Color','w','Position',[10 10 1280 160*n+140]);
        set(fn,'UserData',struct('t_start_epoch',tStart));
        h=gobjects(1,n);for i=1:n,h(i)=irf_subplot(n,1,-i);end
        i=1;
        %% B plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=Bplot?{im};',ic);
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


        %% 总电流强度：第三个panel，原生采样
        axes(h(i));hold on;
        for im=2:-1:1
            if ~isempty(CurrentPlot{im}),irf_plot(CurrentPlot{im}(:,[1 2]),'color','k','Linewidth',0.75);end
        end
        ylabel('|J| [nA/m^2]');grid off;
        set(gca,'ColorOrder',[0 0 0]);
        irf_legend(gca,{'|J|'},CurrentLegendPosition);i=i+1;
        %% 场向/垂直电流：第四个panel，实时磁场分解，J_parallel保留正负
        axes(h(i));hold on;
        for im=2:-1:1
            if isempty(CurrentPlot{im}),continue;end
            irf_plot(CurrentPlot{im}(:,[1 3]),'color','b','Linewidth',0.75);
            irf_plot(CurrentPlot{im}(:,[1 4]),'color','r','Linewidth',0.75);
        end
        irf_plot([tStart 0;tStop 0],'k--','Linewidth',0.75);
        ylabel('J [nA/m^2]');grid off;
        set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
        irf_legend(gca,{'J_{||}','J_{\perp}'},CurrentLegendPosition);i=i+1;
        %% AE plot：原生1min，缺测保留NaN
        axes(h(i));hold on;
        irf_plot(AE,'color','k','Linewidth',0.9);
        ylabel('AE [nT]','fontsize',12);grid off;
        i=i+1;
        %% 时间轴与版式：沿用原图中点位置标题和水平UTC标签
        irf_zoom(h,'x',tint);irf_plot_axis_align(h);
        Labels={'B [nT]','Vi [km/s]','|J| [nA/m^2]','J [nA/m^2]','AE [nT]'};
        currentAvailable=any(cellfun(@(A)~isempty(A)&&any(all(isfinite(A(:,2:4)),2)),CurrentPlot));
        Available=[any(NativeCounts(1,:)),any(NativeCounts(2,:)),currentAvailable,currentAvailable,any(isfinite(AE(:,2)))];
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
        irf_zoom(h,'x',tint);irf_zoom(h,'y');
        if currentAvailable,set(h(3),'YLim',[0 h(3).YLim(2)]);end
        for i=1:n,set(h(i).YLabel,'Units','normalized','Position',[-0.075 0.5 0]);end
        title(h(1),sprintf('X = %.2f  Y = %.2f  Z = %.2f',Event.midpointXYZ_RE),...
            'FontSize',15,'FontWeight','normal');
        set(h(1).Title,'Units','normalized','Position',[0.5 1.09 0]);
        axes(h(n));set(gca,"XTickLabelRotation",0);
        if ~strcmp(Event.startUTC(1:10),Event.endUTC(1:10))
            xlabel(h(n),[Event.startUTC(1:10) ' / ' Event.endUTC(1:10) ' UTC'],'FontSize',12);
        end
        set(h,'XGrid','off','XMinorGrid','off');drawnow;
        %% 垂直电流半透明：沿用IRFU绘出的顶点，只调整显示
        % 在irf_zoom完成后转换，保留原轴范围及原生采样顶点。
        % MATLAB官方透明曲线接口：patch的EdgeAlpha，末尾NaN防止闭合。
        % https://www.mathworks.com/help/matlab/creating_plots/add-transparency-to-graphics-objects.html
        set(h(4),'ClippingStyle','rectangle');
        PerpLines=findobj(h(4),'Type','line','Color',[1 0 0]);
        for ip=1:numel(PerpLines)
            Xperp=PerpLines(ip).XData;Yperp=PerpLines(ip).YData;
            patch(h(4),[Xperp(:).' NaN],[Yperp(:).' NaN],'r',...
                'FaceColor','none','EdgeColor','r','EdgeAlpha',PerpOpacity,'Clipping','on',...
                'LineWidth',PerpLines(ip).LineWidth);
        end
        delete(PerpLines);drawnow;
        %% 保存图件及可复现来源记录
        set(fn,'Renderer','painters','paperpositionmode','auto');
        Stem=[Event.id '_MMS1_B_Vi_J_AE'];
        PNGFile=fullfile(OutputDir,[Stem '.png']);PDFFile=fullfile(OutputDir,[Stem '.pdf']);
        % 先在临时目录完整导出，再替换正式图件。
        ExportTempDir=fullfile(tempdir,'MMS_current_20260721_0307_0315_native_20261002','export_alpha');
        if ~isfolder(ExportTempDir),mkdir(ExportTempDir);end
        TempPNG=fullfile(ExportTempDir,[Stem '.png']);TempPDF=fullfile(ExportTempDir,[Stem '.pdf']);
        exportgraphics(fn,TempPNG,'Resolution',150,'BackgroundColor','white');
        exportgraphics(fn,TempPDF,'ContentType','vector','BackgroundColor','white');
        movefile(TempPNG,PNGFile,'f');movefile(TempPDF,PDFFile,'f');
        Result.png=PNGFile;Result.pdf=PDFFile;Result.nativeCounts=NativeCounts;
        Result.current=CurrentInfo;Result.positionSources=RSources;
        Result.aeRows=size(AE,1);Result.aeValid=sum(isfinite(AE(:,2)));
        Result.aeSources=AE_sources;Result.midpointXYZ_RE=Event.midpointXYZ_RE;
        Result.panelAvailable=Available;Result.panelOrder={'B','Vi','|J|','Jparallel,Jperp','AE'};
        Result.parallelSigned=true;
        Result.parallelOpacity=1;Result.perpendicularOpacity=PerpOpacity;
        Result.perpendicularDefinition='norm of the three-component perpendicular vector; nonnegative';
        Result.parallelReference='instantaneous four-spacecraft mean magnetic field at each native current timestamp from c_4_j';
        Result.instantaneousReferenceVerified=ReferenceVerified;Result.fixedBackgroundField=false;
        CurrentDataFile=fullfile(RecordDir,[Event.id '_native_current.mat']);
        % 只保存本次请求的衍生电流及参考B，原始CDF仍直接读取，不转换全仪器数据。
        save(CurrentDataFile,'J','Jabs','Jpar','Jperp','ReferenceB','CurrentHalf',...
            'tStart','tStop','NativeSummary','-v7.3');
        Result.nativeCurrentData=CurrentDataFile;Result.nativeCurrentSummary=NativeSummary;
        Result.alignment='irf_resamp linear within finite continuous CDF-supported segments, no extrapolation';
        Result.timeAveraging=false;Result.currentAveragingSeconds=0;Result.smoothing=false;
        Result.originalPanelsTimeAveraging=false;
        Result.displayReduction='irf_plot reduce for B,Vi; current native timestamps without reduce; AE native 1min';
        Result.qualityFilter=false;
        Result.complete=currentAvailable;Result.previousNativeValuesMatch=true;
        Result.previousRecordDir=SourceRecordDir;
        close(fn);
    catch ME
        Result.error=getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'FAILED %s\n%s\n',Event.id,Result.error);close all
    end
    fid=fopen(fullfile(RecordDir,[Event.id '_current.json']),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(Result));fclose(fid);
    fprintf('DONE %s complete=%d\n',Event.id,Result.complete);
end
fprintf('NATIVE_CURRENT_OVERVIEW_FINISHED\n');
