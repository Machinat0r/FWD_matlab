%% EV01 MMS1 overview：保留前9个panel，去掉FEEPS和HPCA
% 参照原 Overview_download.m / Overview_download_mms4.m 的直接脚本写法。
% EventList和Spacecraft可在运行前指定，例如 EventList=4; Spacecraft=1;
% 原始CDF -> mms.get_data / mms.db_get_ts -> IRFU坐标变换 -> 原overview绘图段。
% 无自定义function；各cell第1项为burst，第2项为survey/fast。
close all
clearvars -except EventList Spacecraft
clc
if ~exist('EventList','var'), EventList=1; end
if ~exist('Spacecraft','var'), Spacecraft=1; end

%% 路径
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS\';
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS_EV01_MMS1_20260929';
RecordDir=[ParentDir 'derived\EV01_MMS1_9panels_20260929'];
addpath(CodeDir,IRFDir);
irf('check_path'); % 使用IRFU自带的路径初始化
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
if ~isfolder(OutputDir),mkdir(OutputDir);end
if ~isfolder(RecordDir),mkdir(RecordDir);end
mms.db_init('local_file_db',ParentDir);
mms.db_init('db_cache_size_max',2048); % 现成IRFU缓存限制为2GB
mms.db_init('db_cache_enabled',true);
Events=jsondecode(fileread(fullfile(CodeDir,'MMS_event_overview_20260924_events.json')));
set(groot,'defaultFigureVisible','off');

%% 事件和卫星
for ie=EventList
    Event=Events(ie);
    tEvent=[irf_time(Event.start,'utc>epoch');irf_time(Event.end,'utc>epoch')];
    tint=irf.tint(EpochTT(EpochUnix(tEvent(1)-600)),diff(tEvent)+1200);
    tStart=tint.start.epochUnix;tStop=tint.stop.epochUnix;
    for ic=Spacecraft
        fprintf('START %s MMS%d\n',Event.id,ic);
        Result=struct('event',Event.id,'spacecraft',ic,'complete',false,'error','');
        try
        %% load data
        % 保留B、gsmVi、gsmVe、Ni、Ne、Ti、Te等原变量名称。
        c_eval('B?=cell(1,2);gsmVi?=cell(1,2);gsmVe?=cell(1,2);',ic);
        c_eval('Ni?=cell(1,2);Ne?=cell(1,2);Ti?=cell(1,2);Te?=cell(1,2);',ic);
        energy_i=cell(1,2);energy_e=cell(1,2);
        specrec_i=cell(1,2);specrec_e=cell(1,2);
        for im=1:2
            mode='brst';bmode='brst';
            if im==2,mode='fast';bmode='srvy';end
            % load B [GSM, nT]
            c_eval('B?_ts=mms.get_data([''B_gsm_'' bmode],tint,?);',ic);
            c_eval('if ~isempty(B?_ts), B?_ts=B?_ts.tlim(tint);Bt?_ts=B?_ts.abs;B?{im}=[irf.ts2mat(B?_ts) double(Bt?_ts.data)];end',ic);
            % load Vi / Ve [GSM, km/s]
            c_eval('Vi?_ts=mms.get_data([''Vi_gse_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Vi?_ts),gsmVi?_ts=irf_gse2gsm(Vi?_ts.tlim(tint));gsmVi?{im}=irf.ts2mat(gsmVi?_ts);end',ic);
            c_eval('Ve?_ts=mms.get_data([''Ve_gse_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Ve?_ts),gsmVe?_ts=irf_gse2gsm(Ve?_ts.tlim(tint));gsmVe?{im}=irf.ts2mat(gsmVe?_ts);end',ic);
            % load N [cm^-3]
            c_eval('Ni?_ts=mms.get_data([''Ni_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Ni?_ts),Ni?{im}=irf.ts2mat(Ni?_ts.tlim(tint));end',ic);
            c_eval('Ne?_ts=mms.get_data([''Ne_fpi_'' mode ''_l2''],tint,?);',ic);
            c_eval('if ~isempty(Ne?_ts),Ne?{im}=irf.ts2mat(Ne?_ts.tlim(tint));end',ic);
            % load T [eV]，(Tparallel+2*Tperpendicular)/3沿用原overview
            c_eval('Ti_para?_ts=mms.db_get_ts([''mms?_fpi_'' mode ''_l2_dis-moms''],[''mms?_dis_temppara_'' mode],tint);',ic);
            c_eval('Ti_perp?_ts=mms.db_get_ts([''mms?_fpi_'' mode ''_l2_dis-moms''],[''mms?_dis_tempperp_'' mode],tint);',ic);
            c_eval('if ~isempty(Ti_para?_ts)&&~isempty(Ti_perp?_ts),Ti_para?=irf.ts2mat(Ti_para?_ts.tlim(tint));Ti_perp?=irf.ts2mat(Ti_perp?_ts.tlim(tint));assert(isequal(Ti_para?(:,1),Ti_perp?(:,1)));Ti?{im}=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3,Ti_para?(:,2),Ti_perp?(:,2)];end',ic);
            c_eval('Te_para?_ts=mms.db_get_ts([''mms?_fpi_'' mode ''_l2_des-moms''],[''mms?_des_temppara_'' mode],tint);',ic);
            c_eval('Te_perp?_ts=mms.db_get_ts([''mms?_fpi_'' mode ''_l2_des-moms''],[''mms?_des_tempperp_'' mode],tint);',ic);
            c_eval('if ~isempty(Te_para?_ts)&&~isempty(Te_perp?_ts),Te_para?=irf.ts2mat(Te_para?_ts.tlim(tint));Te_perp?=irf.ts2mat(Te_perp?_ts.tlim(tint));assert(isequal(Te_para?(:,1),Te_perp?(:,1)));Te?{im}=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3,Te_para?(:,2),Te_perp?(:,2)];end',ic);
            % load FPI omnidirectional energy flux
            c_eval('energy_i{im}=mms.db_get_variable([''mms?_fpi_'' mode ''_l2_dis-moms''],[''mms?_dis_energyspectr_omni_'' mode],tint);',ic);
            c_eval('energy_e{im}=mms.db_get_variable([''mms?_fpi_'' mode ''_l2_des-moms''],[''mms?_des_energyspectr_omni_'' mode],tint);',ic);
            for particle={'i','e'}
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
        Names={'B','gsmVi','gsmVe','Ni','Ne','Ti','Te'};
        NativeCounts=zeros(9,2);

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
        for particle={'i','e'}
            eval(['S=specrec_' particle{1} ';']);
            row=8;if strcmp(particle{1},'e'),row=9;end
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

        %% Init figure
        n=9;i=1;
        set(0,'DefaultAxesFontSize',10);
        set(0,'DefaultLineLineWidth',0.5);
        fn=figure('Visible','off','Color','w','Position',[10 10 1280 max(960,90*n+120)]);
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
        grid off;ylabel('B [nT]','fontsize',10);
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
        grid off;ylabel('Vi [km/s]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
        irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
        i=i+1;

        %% Ve plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=gsmVe?{im};',ic);
            if isempty(A),continue;end
            irf_plot([A(:,1) A(:,2)],'reduce','color','b','Linewidth',0.75);
            irf_plot([A(:,1) A(:,3)],'reduce','color','g','Linewidth',0.75);
            irf_plot([A(:,1) A(:,4)],'reduce','color','r','Linewidth',0.75);
        end
        irf_plot([tStart 0;tStop 0],'k--','Linewidth',0.75);
        grid off;ylabel('Ve [km/s]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
        irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.97 0.92]);
        i=i+1;

        %% Electric field
        % 只读取GSE电场并调用irf_gse2gsm；仅有DSL时不画电场。
        % 直接调用现有mms.db_get_ts。burst电场按600s读取。
        % 显示压缩直接调用IRFU包自带reduce_to_width，不另写压缩函数。
        axes(h(i));hold on;Ecounts=zeros(1,2);Eblocks=[];
        mode='fast';
        c_eval('E?_ts=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_gse_fast_l2'',tint);',ic);
        c_eval('E?=[];',ic);
        c_eval('if ~isempty(E?_ts),E?=irf.ts2mat(E?_ts.tlim(tint));end',ic);
        c_eval('if ~isempty(E?),E?=irf_gse2gsm(E?);end',ic);
        c_eval('Efast=E?;',ic);
        if ~isempty(Efast),Ecounts(2)=sum(any(isfinite(Efast(:,2:end)),2));end
        c_eval('Efiles=mms.db_list_files(''mms?_edp_brst_l2_dce'',tint);',ic);
        for t0=tStart:600:tStop-1e-6
            if isempty(Efiles),break;end
            fprintf('  E burst %s MMS%d +%.0fs\n',Event.id,ic,t0-tStart);
            tintE=irf.tint(EpochTT(EpochUnix(t0)),min(600,tStop-t0));
            c_eval('E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tintE,?);',ic);
            c_eval('E?=[];if ~isempty(E?_ts),E?=irf.ts2mat(E?_ts.tlim(tintE));end',ic);
            c_eval('A=E?;',ic);
            if isempty(A),continue;end
            A=A(A(:,1)>=t0&A(:,1)<min(t0+600,tStop),:);
            if isempty(A),continue;end
            A=irf_gse2gsm(A);
            Ecounts(1)=Ecounts(1)+sum(any(isfinite(A(:,2:end)),2));
            dt=median(diff(A(:,1)));valid=find(all(isfinite(A(:,2:end)),2));
            cut=find(diff(valid)>1|diff(A(valid,1))>3*dt);
            first=[1;cut+1];last=[cut;numel(valid)];
            if isempty(valid),first=[];last=[];end
            for k=1:numel(first)
                lo=A(valid(first(k)),1)-dt/2;hi=A(valid(last(k)),1)+dt/2;
                Eblocks=[Eblocks;lo hi]; %#ok<AGROW>
                if ~isempty(Efast),Efast(Efast(:,1)>=lo&Efast(:,1)<=hi,2:end)=NaN;end
                rows=valid(first(k):last(k));
                for j=2:size(A,2)
                    [tPlot,EPlot]=reduce_to_width(A(rows,1),A(rows,j),1600,[A(rows(1),1) A(rows(end),1)]);
                    colors='bgr';
                    irf_plot([tPlot EPlot],'color',colors(j-1),'Linewidth',0.75);
                end
            end
        end
        % fast只填补没有burst的时段；分段绘制避免跨数据缺口连线。
        if ~isempty(Efast)
            dt=median(diff(Efast(:,1)));valid=find(all(isfinite(Efast(:,2:end)),2));
            cut=find(diff(valid)>1|diff(Efast(valid,1))>3*dt);
            first=[1;cut+1];last=[cut;numel(valid)];
            if isempty(valid),first=[];last=[];end
            for k=1:numel(first)
                rows=valid(first(k):last(k));
                for j=2:size(Efast,2)
                    [tPlot,EPlot]=reduce_to_width(Efast(rows,1),Efast(rows,j),1600,[Efast(rows(1),1) Efast(rows(end),1)]);
                    colors='bgr';
                    irf_plot([tPlot EPlot],'color',colors(j-1),'Linewidth',0.75);
                end
            end
        end
        irf_plot([tStart 0;tStop 0],'k--','Linewidth',0.75);
        grid off;ylabel('E [mV/m]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
        irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
        i=i+1;

        %% N plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=Ne?{im};',ic);
            if ~isempty(A),irf_plot(A,'color','b','Linewidth',0.75);end
            c_eval('A=Ni?{im};',ic);
            if ~isempty(A),irf_plot(A,'color','g','Linewidth',0.75);end
        end
        grid off;ylabel('N [cm^{-3}]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
        irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
        i=i+1;

        %% Ti plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=Ti?{im};',ic);
            if isempty(A),continue;end
            irf_plot([A(:,1) A(:,2)],'color','k','Linewidth',0.75);
            irf_plot([A(:,1) A(:,3)],'color','b','Linewidth',0.75);
            irf_plot([A(:,1) A(:,4)],'color','r','Linewidth',0.75);
        end
        grid off;ylabel('Ti [eV]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
        irf_legend(gca,{'Ti','T_/_/','T_⊥'},[0.97 0.92]);
        i=i+1;

        %% Te plot
        axes(h(i));hold on;
        for im=2:-1:1
            c_eval('A=Te?{im};',ic);
            if isempty(A),continue;end
            irf_plot([A(:,1) A(:,2)],'color','k','Linewidth',0.75);
            irf_plot([A(:,1) A(:,3)],'color','b','Linewidth',0.75);
            irf_plot([A(:,1) A(:,4)],'color','r','Linewidth',0.75);
        end
        grid off;ylabel('Te [eV]','fontsize',10);
        set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
        irf_legend(gca,{'Te','T_/_/','T_⊥'},[0.97 0.92]);
        i=i+1;

        %% plot ION energy spectrom
        axes(h(i));hold on;colormap(h(i),jet);
        for im=2:-1:1
            specrec_p_i=specrec_i{im};
            if isempty(specrec_p_i),continue;end
            irf_spectrogram(h(i),specrec_p_i,'log','donotshowcolorbar');
        end
        grid off;set(h(i),'yscale','log','ytick',[1e1 1e2 1e3 1e4],'fontsize',10);
        ylabel('Ei(ev)','fontsize',10);set(gca,'Ylim',[1 4e4]);
        if any(NativeCounts(8,:)>0)
            hcb=colorbar(h(i));ylabel(hcb,{'log10(keV/','(cm^2 s sr','keV))'},'fontsize',7);
        end
        i=i+1;

        %% plot e energy spectrom
        axes(h(i));hold on;colormap(h(i),jet);
        for im=2:-1:1
            specrec_p_e=specrec_e{im};
            if isempty(specrec_p_e),continue;end
            irf_spectrogram(h(i),specrec_p_e,'log','donotshowcolorbar');
        end
        grid off;set(h(i),'yscale','log','ytick',[1e1 1e2 1e3 1e4],'fontsize',10);
        ylabel('Ee(ev)','fontsize',10);set(gca,'Ylim',[1 4e4]);
        if any(NativeCounts(9,:)>0)
            hcb=colorbar(h(i));ylabel(hcb,{'log10(keV/','(cm^2 s sr','keV))'},'fontsize',7);
        end
        i=i+1;

        %% 统一时间轴：沿用原overview现成函数
        irf_zoom(h,'x',tint);
        irf_plot_axis_align(h);
        % irf_subplot(n,1,-i)已采用紧凑panel；仅统一左右边界与页边距。
        for i=1:n
            PanelStep=0.84/n;
            set(h(i),'Position',[0.105 0.07+(n-i)*PanelStep 0.77 PanelStep-0.0015]);
            set(h(i),'FontName','Arial','FontSize',10,'TickDir','out','Box','on');
            if i<n,set(h(i),'XTickLabel',[]);xlabel(h(i),'');end
        end
        cb=findall(fn,'Type','colorbar');
        for k=1:numel(cb)
            ax=cb(k).Axes;pos=ax.Position;
            cb(k).Position=[0.89 pos(2) 0.009 pos(4)];cb(k).FontSize=8;
        end

        Coord='GSM';
        title(h(1),sprintf('%s  MMS%d  %s  UTC',Event.id,ic,Event.start(1:10)),'FontSize',13,'FontWeight','normal');
        Available=[sum(NativeCounts(1,:)),sum(NativeCounts(2,:)),sum(NativeCounts(3,:)),sum(Ecounts),sum(NativeCounts(4:5,:),'all'),sum(NativeCounts(6,:)),sum(NativeCounts(7,:)),sum(NativeCounts(8,:)),sum(NativeCounts(9,:))]>0;
        for i=find(~Available)
            cla(h(i));text(h(i),.5,.5,'No available L2 data','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
            set(h(i),'YTick',[]);
        end
        % cla保留坐标轴标签，缺失panel仍保留原物理量标注。
        Labels={'B [nT]','Vi [km/s]','Ve [km/s]','E [mV/m]','N [cm^{-3}]','Ti [eV]','Te [eV]','Ei(ev)','Ee(ev)'};
        for i=1:n,ylabel(h(i),Labels{i},'FontSize',10);end
        irf_zoom(h,'x',tint);
        irf_zoom(h(1:7),'y'); % IRFU现成函数避免panel边界刻度文字重叠
        for i=1:n,set(h(i).YLabel,'Units','normalized','Position',[-0.075 0.5 0]);end
        set(h(1).Title,'Units','normalized','Position',[0.5 1.1 0]);
        % 沿用原 Overview_download.m 的时间标签设置。
        axes(h(n));set(gca,"XTickLabelRotation",0);
        set(h,'XGrid','off','XMinorGrid','off');
        if ~Available(4)
            cla(h(4));text(h(4),.5,.5,'No GSE/GSM E data','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
            set(h(4),'YTick',[]);ylabel(h(4),'E [mV/m]','FontSize',10);
            set(h(4).YLabel,'Units','normalized','Position',[-0.075 0.5 0]);
        end
        drawnow;

        %% 能谱纵轴：原下限低于10 eV时改为20 eV，上限保持原值
        SpectrumPanels=[];SpectrumLimitsBefore=[];SpectrumLimitsAfter=[];
        for i=1:n
            if ~any(contains(string(Labels{i}),{'Ei(ev)','Ee(ev)'})),continue;end
            SpectrumPanels(end+1)=i;
            EnergyLimits=get(h(i),'YLim');SpectrumLimitsBefore(end+1,:)=EnergyLimits;
            if EnergyLimits(1)<10
                set(h(i),'YLim',[20 EnergyLimits(2)]);
            end
            SpectrumLimitsAfter(end+1,:)=get(h(i),'YLim');
        end

        %% 出图保存部分
        set(gcf,'Renderer','painters');
        set(gcf,'paperpositionmode','auto');
        Base=sprintf('%s_MMS%d_overview',Event.id,ic);
        Result.startUTC=char(tint.start.utc);Result.endUTC=char(tint.stop.utc);
        Result.coordinateSystem=Coord;Result.nativeCounts=NativeCounts;
        Result.timeLabelRotation=h(n).XTickLabelRotation;Result.eventMarkers=false;
        Result.electricSource='GSE L2 converted with irf_gse2gsm; DSL-only data excluded';
        Result.nativeCountNames=[Names {'energy_i','energy_e'}];
        Result.modeOrder={'burst','surveyFast'};Result.electricCounts=Ecounts;
        Result.electricBurstIntervals=Eblocks;Result.panelAvailable=Available;
        Result.panelCount=n;
        Result.feepsHpcaPanels=false;
        Result.spectrumPanels=SpectrumPanels;
        Result.spectrumLimitsBefore=SpectrumLimitsBefore;Result.spectrumLimitsAfter=SpectrumLimitsAfter;
        Result.spectrumLowerLimitRule="below 10 eV -> 20 eV; upper unchanged";
        Result.status='no L2 data';Result.png='';
        if any(Available)
            Result.png=fullfile(OutputDir,[Base '.png']);
            exportgraphics(fn,Result.png,'Resolution',150,'BackgroundColor','white');
            Result.status='plotted';
        end
        Result.complete=true;
        close(fn);
        catch ME
            Result.error=getReport(ME,'extended','hyperlinks','off');
            fprintf(2,'FAILED %s MMS%d\n%s\n',Event.id,ic,Result.error);
            close all
        end
        fid=fopen(fullfile(RecordDir,sprintf('%s_MMS%d_overview.json',Event.id,ic)),'w','n','UTF-8');
        fprintf(fid,'%s',jsonencode(Result));fclose(fid);
        fprintf('DONE %s MMS%d complete=%d\n',Event.id,ic,Result.complete);
    end
end
fprintf('OVERVIEW_SCRIPT_FINISHED\n');
