%% 四星磁场curlometer电流：用户指定的五个事件（2026-10-01）
% 参考用户Plot_TCS.m和原Overview_MMS1_tail_20260720_0815.m，直接脚本/%%分区。
% B/Vi/AE保留原MMS1绘图流程；在Vi与AE之间插入两栏。
% 电流：四星FGM GSM磁场[nT]和MEC GSM位置[km]，IRFU c_4_j，输出[nA/m^2]。
% 第三栏只画总电流模|J|；第四栏画有符号J_parallel和垂直电流幅值J_perp。
% J_parallel=J dot Bhat，正/负分别对应沿/反沿四星平均磁场。
% J_perp为三分量垂直矢量的模，非负；图例两分量不加绝对值符号。
% 先用同一时刻四星平均B逐点分解，再分别对|J|、J_parallel、J_perp做独立60s算术平均。
% 每图起点开始分箱，末尾不足60s按实际观测平均；无观测的分钟保留NaN。
% 原B/Vi/AE不做新增时间平均，磁场分解参考始终为逐时刻磁场。
% 同一模式四星磁场/位置以irf_resamp(...,'linear')对齐到MMS1磁场时刻。
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
OriginalRecordDir=[ParentDir 'derived\MMS1_tail_20260720_0815'];
RecordDir=[ParentDir 'derived\MMS_current_5events_20261001'];
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS_current_5events_20261001';
PerpOpacity=0.5; % 垂直电流曲线不透明度：0全透明，1完全不透明。
DataDir=[ParentDir 'ancillary\omni\hro_1min\2026'];
addpath(CodeDir,IRFDir);irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
mms.db_init('local_file_db',ParentDir);
mms.db_init('db_cache_size_max',2048);
mms.db_init('db_cache_enabled',true);
DownloadState=jsondecode(fileread(fullfile(RecordDir,'download_complete.json')));
assert(DownloadState.complete,'请先补齐本批四星原CDF');
Windows=jsondecode(fileread(fullfile(RecordDir,'windows.json')));
if ~exist('WindowList','var'),WindowList=1:numel(Windows);end
if ~isfolder(OutputDir),mkdir(OutputDir);end
Manifest=jsondecode(fileread(fullfile(RecordDir,'manifest.json')));
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
    % 08-02右侧电流持续上升，图例放左侧避开曲线。
    if strcmp(Event.id,'W082_20260802115210'),CurrentLegendPosition=[0.03 0.92];end
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
        %% load R1-R4：直接读取各星MEC原生GSM位置[km]
        Rnative=cell(1,4);RSources=cell(1,4);
        for ic=1:4
            prefix=sprintf('mms%d_mec_srvy_l2_epht89d',ic);
            varname=sprintf('mms%d_mec_r_gsm',ic);
            R_ts=mms.db_get_ts(prefix,varname,tintR);
            assert(~isempty(R_ts),'缺少MMS%d MEC位置',ic);
            Rnative{ic}=irf.ts2mat(R_ts);
            RSources{ic}=mms.db_list_files(prefix,tintR);
        end
        %% 四星电流：分别计算burst和survey，统一到本模式MMS1磁场时刻
        J=cell(1,2);Jabs=cell(1,2);Jpar=cell(1,2);Jperp=cell(1,2);
        CurrentInfo=cell(1,2);CurrentHalf=cell(1,2);ReferenceVerified=false(1,2);
        for im=1:2
            bmode='brst';if im==2,bmode='srvy';end
            Bnative=cell(1,4);Bhalf=cell(1,4);BSources=cell(1,4);
            counts=zeros(1,4);cadence=zeros(1,4);
            for ic=1:4
                B_ts=mms.get_data(['B_gsm_' bmode],tintB,ic);
                if isempty(B_ts),continue;end
                B_ts=B_ts.tlim(tintB);Bnative{ic}=irf.ts2mat(B_ts);
                prefix=sprintf('mms%d_fgm_%s_l2',ic,bmode);
                varname=sprintf('mms%d_fgm_bdeltahalf_%s_l2',ic,bmode);
                D_ts=mms.db_get_ts(prefix,varname,tintB);
                assert(~isempty(D_ts),'缺少FGM原CDF采样宽度');
                assert(strcmp(D_ts.units,'s'),'FGM采样宽度单位非s');
                [found,index]=ismember(B_ts.time.ttns,D_ts.time.ttns);
                assert(all(found),'FGM时间与bdeltahalf不对应');
                Bhalf{ic}=double(D_ts.data(index));
                counts(ic)=sum(all(isfinite(Bnative{ic}(:,2:4)),2));
                cadence(ic)=median(diff(Bnative{ic}(:,1)));
                BSources{ic}=mms.db_list_files(prefix,tintB);
            end
            CurrentInfo{im}=struct('mode',bmode,'nativeBCounts',counts,'medianCadenceSeconds',cadence,...
                'sources',{BSources},'finiteCurrentRows',0);
            if any(cellfun(@isempty,Bnative)),continue;end
            T=Bnative{1}(:,1);inside=T>=tStart&T<=tStop;
            CurrentHalf{im}=Bhalf{1}(inside);T=T(inside);
            if isempty(T),continue;end
            %% 时间对齐：在CDF支持的连续段内直接调用IRFU，不外推
            Aligned=cell(2,4);
            for iv=1:2
                for ic=1:4
                    if iv==1
                        A=Bnative{ic};HalfWidth=Bhalf{ic};
                    else
                        A=Rnative{ic};HalfWidth=15*ones(size(A,1),1);
                        % 本批官方MEC原生30s；15s是相邻轨道时刻间距的一半。
                    end
                    assert(all(diff(A(:,1))>0),'原数据时间必须严格递增');
                    valid=find(all(isfinite(A(:,2:4)),2)&isfinite(HalfWidth)&HalfWidth>0);
                    out=[T nan(numel(T),3)];
                    if ~isempty(valid)
                        tol=8*eps(max(abs(A(:,1))));
                        cut=find(diff(valid)>1 | diff(A(valid,1))>...
                            HalfWidth(valid(1:end-1))+HalfWidth(valid(2:end))+tol);
                        first=[1;cut+1];last=[cut;numel(valid)];
                        for k=1:numel(first)
                            rows=valid(first(k):last(k));
                            target=T>=A(rows(1),1)&T<=A(rows(end),1);
                            if numel(rows)>=2&&any(target)
                                out(target,:)=irf_resamp(A(rows,:),T(target),'linear');
                            elseif numel(rows)==1
                                exact=T==A(rows,1);out(exact,2:4)=repmat(A(rows,2:4),sum(exact),1);
                            end
                        end
                    end
                    Aligned{iv,ic}=out;
                end
            end
            c_eval('B?=Aligned{1,?};R?=Aligned{2,?};',1:4);
            validInput=true(numel(T),1);
            for ic=1:4
                validInput=validInput&all(isfinite(Aligned{1,ic}(:,2:4)),2)&all(isfinite(Aligned{2,ic}(:,2:4)),2);
            end
            %% c_4_j：完全复用用户Plot_TCS.m的8参数调用
            [J_B,divB,Bmean]=c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);
            J_B(~validInput,2:4)=NaN;Bmean(~validInput,2:4)=NaN;
            J_B(:,2:4)=J_B(:,2:4)*1e9;
            % 核对参考磁场：同一时刻的四星空间平均，没有时间背景平均。
            Binstant=0.25*B1(:,2:4)+0.25*B2(:,2:4)+0.25*B3(:,2:4)+0.25*B4(:,2:4);
            assert(isequal(Bmean(:,1),J_B(:,1)),'磁场与电流时刻不一致');
            assert(isequaln(Bmean(validInput,2:4),Binstant(validInput,:)),...
                'c_4_j参考磁场与逐时刻四星平均不一致');
            ReferenceVerified(im)=true;
            [J_parallel,J_perp]=irf_dec_parperp(Bmean,J_B);
            J{im}=J_B;Jabs{im}=irf_abs(J_B,1);
            Jpar{im}=J_parallel(:,2);Jperp{im}=irf_abs(J_perp,1);
            good=all(isfinite(J_B(:,2:4)),2);
            CurrentInfo{im}.parallelRange_nA_m2=[min(Jpar{im}(good)) max(Jpar{im}(good))];
            CurrentInfo{im}.parallelNegativeRows=sum(good&Jpar{im}<0);
            CurrentInfo{im}.parallelPositiveRows=sum(good&Jpar{im}>0);
            residual=abs(Jabs{im}.^2-Jpar{im}.^2-Jperp{im}.^2);
            assert(all(residual(good)<=1e-10*max(1,Jabs{im}(good).^2)),'场向/垂直分解恒等式失败');
            CurrentInfo{im}.finiteCurrentRows=sum(good);
            CurrentInfo{im}.alignedRows=numel(T);
            CurrentInfo{im}.minMeanB_nT=min(irf_abs(Bmean(good,:),1));
            CurrentInfo{im}.decompositionMaxResidual=max(residual(good));
            CurrentInfo{im}.divBoverCurlB=abs(divB(good,2)*1e9)./Jabs{im}(good);
            % 保存诊断分位数；不按诊断值删除电流。
            diag=CurrentInfo{im}.divBoverCurlB;
            CurrentInfo{im}.divBoverCurlB=prctile(diag(isfinite(diag)),[50 95 100]);
            fprintf('CURRENT %s %s valid=%d / %d\n',Event.id,bmode,sum(good),numel(T));
        end
        %% 电流模式显示：四星有效burst覆盖处隐藏survey
        if ~isempty(J{1})&&~isempty(J{2})
            valid=find(all(isfinite(J{1}(:,2:4)),2));
            if ~isempty(valid)
                width=CurrentHalf{1};
                cut=find(diff(valid)>1|diff(J{1}(valid,1))>...
                    width(valid(1:end-1))+width(valid(2:end))+8*eps(max(abs(J{1}(:,1)))));
                first=[1;cut+1];last=[cut;numel(valid)];
                for k=1:numel(first)
                    lo=J{1}(valid(first(k)),1);hi=J{1}(valid(last(k)),1);
                    hide=J{2}(:,1)>=lo&J{2}(:,1)<=hi;
                    J{2}(hide,2:4)=NaN;Jabs{2}(hide)=NaN;Jpar{2}(hide)=NaN;Jperp{2}(hide)=NaN;
                end
            end
        end
        %% 电流1min平均：先逐点分解，再按每图起点分箱；每60s输出一个均值
        % burst有效时沿用burst优先，其余用survey；每个有效采样点等权算术平均。
        % IRFU平均窗口为(中心-30s,中心+30s]，相邻窗口不重叠，thresh=0不剔除极值。
        % 仅使用原图时窗内的有效观测；无观测窗口保留NaN，不插值补分钟。
        BinLeft=tStart+60*(0:ceil((tStop-tStart)/60)-1)';
        BinRight=min(BinLeft+60,tStop);BinCenter=BinLeft+30;
        PlotTime=(BinLeft+BinRight)/2;
        MinuteInput=[];MinuteCountByMode=zeros(numel(BinLeft),2);
        for im=1:2
            if isempty(J{im}),continue;end
            A=[J{im}(:,1) Jabs{im} Jpar{im} Jperp{im}];
            A=A(A(:,1)>tStart&A(:,1)<=tStop&all(isfinite(A(:,2:4)),2),:);
            if isempty(A),continue;end
            ib=discretize(A(:,1),[BinLeft;BinLeft(end)+60],'IncludedEdge','right');
            MinuteCountByMode(:,im)=accumarray(ib,1,[numel(BinLeft) 1],@sum,0);
            MinuteInput=[MinuteInput;A]; %#ok<AGROW>
        end
        MinuteInput=sortrows(MinuteInput,1);
        MinuteCurrent=[PlotTime nan(numel(PlotTime),3)];
        if ~isempty(MinuteInput)
            assert(all(diff(MinuteInput(:,1))>0),'优先模式合并后时间不递增');
            MeanCurrent=irf_resamp(MinuteInput,BinCenter,'mean','window',60,'thresh',0);
            assert(isequal(MeanCurrent(:,1),BinCenter),'分钟平均时间轴改变');
            MinuteCurrent(:,2:4)=MeanCurrent(:,2:4);
        end
        MinuteCount=sum(MinuteCountByMode,2);
        assert(sum(MinuteCount)==size(MinuteInput,1),'分钟分箱点数不守恒');
        assert(isequal(all(isfinite(MinuteCurrent(:,2:4)),2),MinuteCount>0),...
            '分钟缺测与有效输入点数不一致');
        MinuteInfo=struct('windowSeconds',60,'binOriginUTC',Event.startUTC,...
            'binConvention','(left,right]','method','irf_resamp mean window 60 thresh 0',...
            'sampleWeighting','equal weight per finite preferred-mode sample',...
            'decompositionBeforeAveraging',true,'rows',size(MinuteCurrent,1),...
            'finiteRows',sum(MinuteCount>0),'partialLastBinSeconds',BinRight(end)-BinLeft(end),...
            'columns',{{'epochUnix','mean_absJ_nA_m2','mean_Jparallel_nA_m2','mean_Jperp_nA_m2'}},...
            'binLeftEpoch',BinLeft,'binRightEpoch',BinRight,...
            'validInputCountByModeBurstSurvey',MinuteCountByMode,'values',MinuteCurrent);
        fprintf('MINUTE %s bins=%d valid=%d lastDuration=%.6f s\n',...
            Event.id,size(MinuteCurrent,1),sum(MinuteCount>0),BinRight(end)-BinLeft(end));
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


        %% 总电流强度：第三个panel仅画|J|的1min均值
        axes(h(i));hold on;
        irf_plot(MinuteCurrent(:,[1 2]),'color','k','Linewidth',0.75);
        ylabel('|J| [nA/m^2]');grid off;
        set(gca,'ColorOrder',[0 0 0]);
        irf_legend(gca,{'|J|'},CurrentLegendPosition);i=i+1;
        %% 场向/垂直电流：第四个panel，分别平均逐时刻分量，J_parallel保留正负
        axes(h(i));hold on;
        irf_plot(MinuteCurrent(:,[1 3]),'color','b','Linewidth',0.75);
        irf_plot(MinuteCurrent(:,[1 4]),'color','r','Linewidth',0.75);
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
        currentAvailable=any(all(isfinite(MinuteCurrent(:,2:4)),2));
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
        % 在irf_zoom完成后转换，保留原轴范围及1min均值顶点。
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
        ExportTempDir=fullfile(tempdir,'MMS_current_5events_20261001','export_alpha');
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
        Result.minuteAverage=MinuteInfo;
        Result.alignment='irf_resamp linear within finite continuous CDF-supported segments, no extrapolation';
        Result.timeAveraging=true;Result.currentAveragingSeconds=60;Result.smoothing=false;
        Result.originalPanelsTimeAveraging=false;
        Result.displayReduction='irf_plot reduce for B,Vi; current 1min means without reduce; AE native 1min';
        Result.qualityFilter=false;
        Baseline=jsondecode(fileread(fullfile(OriginalRecordDir,[Event.id '_overview.json'])));
        assert(isequal(NativeCounts,Baseline.nativeCounts(1:2,:)),'B/Vi原始点数与附件基线不一致');
        Result.originalBViCountsMatch=true;Result.complete=currentAvailable;
        close(fn);
    catch ME
        Result.error=getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'FAILED %s\n%s\n',Event.id,Result.error);close all
    end
    fid=fopen(fullfile(RecordDir,[Event.id '_current.json']),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(Result));fclose(fid);
    fprintf('DONE %s complete=%d\n',Event.id,Result.complete);
end
fprintf('FIVE_CURRENT_OVERVIEWS_FINISHED\n');
