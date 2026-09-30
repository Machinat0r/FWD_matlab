function fn=Overview_download_events_20260923(P,D,ev,ic,a,b,eventStart,eventEnd,counts,coverage)
% 由 Overview_download.m 复制并保留所需绘图段；Electric field来自Overview_download_mms4.m。
% 2026-09-23：供20事件/四星批处理调用。原脚本保留不动。
% 改动：时间/卫星/路径参数化；删除无关分析段和事件特定轴限；缺测留空。
% B、V为GSM；MMS1-3 E为GSM，MMS4 E仅DSL XY。n:cm^-3，T/energy:eV。
% 输入P为原生观测的绘图数组（burst优先，缺口保留），D保留能谱及原生模式覆盖。
%% 与原脚本相同的变量命名
c_eval('B?=P.B(:,1:4); Bt?=P.B(:,[1 5]);',ic);
c_eval('gsmVi?=P.Vi; gsmVe?=P.Ve; E?=P.E;',ic);
c_eval('Ni?=P.Ni; Ne?=P.Ne; Ti?=P.Ti; Te?=P.Te;',ic);
%% Init figure：沿用原脚本分区、irf_subplot与c_eval绘图
n=9;i=1;
fn=figure('Visible','off','Color','w','Units','pixels','Position',[40 40 1500 1800]);
set(fn,'UserData',struct('t_start_epoch',a));
annotation(fn,'textbox',[.105 .961 .83 .028],'String',sprintf('%s  |  MMS%d  |  %s',ev.id,ic,ev.start(1:10)), ...
    'EdgeColor','none','FontSize',18,'FontWeight','bold','Interpreter','none');
annotation(fn,'textbox',[.105 .938 .83 .023],'String',sprintf('%s to %s UTC | B,V: GSM; E: %s', ...
    datestr(datetime(a,'ConvertFrom','posixtime'),'yyyy-mm-dd HH:MM:SS'),datestr(datetime(b,'ConvertFrom','posixtime'),'yyyy-mm-dd HH:MM:SS'),D.Ecoord), ...
    'EdgeColor','none','FontSize',10,'Interpreter','none');
annotation(fn,'textbox',[.105 .917 .83 .021],'String',['User event note: ' ev.label], ...
    'EdgeColor','none','FontSize',10,'Interpreter','none');
h=gobjects(9,1);
for j=1:9,h(j)=irf_subplot(n,1,-j);end
for j=1:8,set(h(j),'Position',[.105 .82-(j-1)*.100 .74 .083]);end
%% B plot：复制原B plot段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',10);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% Vi plot：沿用原矢量绘图段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
ylabel('Vi [km/s]','fontsize',8);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.05 0.92]);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% Ve plot：复制原Ve plot段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
ylabel('Ve [km/s]','fontsize',8);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.05 0.92]);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% Electric field：复制原MMS4 overview电场段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
if ic~=4, c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on; end
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
if ic==4,irf_legend(gca,{'E_x','E_y'},[0.97 0.92]);else,irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);end
ylabel(['E ' D.Ecoord ' [mV/m]'],'fontsize',10);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% N plot：复制原N plot段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
if ic==4,set(gca,'ColorOrder',[0 1 0]);irf_legend(gca,{'Ni'},[0.97 0.92]);else, irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);end
ylabel('N [cm^{-3}]','fontsize',8);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% T plot：按原N plot双曲线写法显示原程序标量温度
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
c_eval("irf_plot([Te?(:,1) Te?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Ti?(:,1) Ti?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
if ic==4,set(gca,'ColorOrder',[0 1 0]);irf_legend(gca,{'Ti'},[0.97 0.92]);else, irf_legend(gca,{'Te','Ti'},[0.97 0.92]);end
ylabel('T [eV]','fontsize',8);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% plot ION energy spectrom：沿用原能谱段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
hold(h(i),'on');colormap(h(i),jet)
for im=1:2
    specrec_p_ION=D.Si{im};
    if isempty(specrec_p_ION),continue;end
    irf_spectrogram(h(i),specrec_p_ION,'log','donotshowcolorbar');
end
grid off;set(h(i),'yscale','log');
set(h(i),'ytick',[1 1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ion energy [eV]','fontsize',10);
set(h(i),'Ylim',[1 4e4]);
hcb=colorbar(h(i));hcb.Position=[.862 .82-(i-1)*.100 .013 .083];
hcb.Label.String={'log_{10} DEF','keV/(cm^2 s sr keV)'};hcb.Label.FontSize=8;hcb.FontSize=8;
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% plot e energy spectrom：沿用原能谱段
axes(h(i));
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
if counts(i)>0
hold(h(i),'on');colormap(h(i),jet)
for im=1:2
    specrec_p_e=D.Se{im};
    if isempty(specrec_p_e),continue;end
    irf_spectrogram(h(i),specrec_p_e,'log','donotshowcolorbar');
end
grid off;set(h(i),'yscale','log');
set(h(i),'ytick',[1 1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Electron energy [eV]','fontsize',10);
set(h(i),'Ylim',[1 4e4]);
hcb=colorbar(h(i));hcb.Position=[.862 .82-(i-1)*.100 .013 .083];
hcb.Label.String={'log_{10} DEF','keV/(cm^2 s sr keV)'};hcb.Label.FontSize=8;hcb.FontSize=8;
set(h(i),'Position',[.105 .82-(i-1)*.100 .74 .083]);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;
%% 模式覆盖：新增，用于标明实际burst覆盖和survey/fast背景
axes(h(9));set(h(9),'Position',[.105 .055 .74 .04]);hold(h(9),'on');
fields={'B','Vi','Ve','E'};
for k=1:4
    for im=1:2
        A=D.(fields{k}){im};if isempty(A),continue;end
        mode='surveyFast';col=[.65 .69 .74];if im==2,mode='burst';col=[.96 .46 .08];end
        % 直接使用读取阶段已核对的连续覆盖秒数，避免从显示点重新估计覆盖。
        blocks=coverage.(fields{k}).(mode).blocks;
        for j=1:size(blocks,1),plot(h(9),blocks(j,:)-a,[5-k 5-k],'Color',col,'LineWidth',5);end
    end
end
ylim(h(9),[.5 4.5]);yticks(h(9),1:4);yticklabels(h(9),{'E','V_e','V_i','B'});ylabel(h(9),'Mode');
text(h(9),1.02,.72,'gray: survey/fast','Units','normalized','FontSize',8);
text(h(9),1.02,.24,'orange: burst','Units','normalized','Color',[.85 .36 .05],'FontSize',8);
%% 统一时间轴与保存前样式，保留原图逐panel组织
span=b-a;if span<=1800,step=300;elseif span<=5400,step=600;else,step=1800;end
tick=(ceil(a/step)*step:step:floor(b/step)*step)-a;
tickLabels=cellstr(datetime(tick+a,'ConvertFrom','posixtime','TimeZone','UTC','Format','HH:mm'));
ylabels={'B [nT]','Vi [km/s]','Ve [km/s]',['E ' D.Ecoord ' [mV/m]'],'N [cm^{-3}]','T [eV]','Ion energy [eV]','Electron energy [eV]'};
for k=1:8,ylabel(h(k),ylabels{k},'FontSize',10);end
for k=1:9
    if k<=6,ylim(h(k),'auto');end
    xlim(h(k),[0 span]);xticks(h(k),tick);grid(h(k),'off');
    set(h(k),'FontName','Arial','FontSize',10,'TickDir','out','Box','on','Layer','top');hold(h(k),'on');
    for bound=unique([eventStart,eventEnd]),xline(h(k),bound-a,'--','Color',[.35 .35 .35],'LineWidth',.8);end
    if k<9,xticklabels(h(k),[]);xlabel(h(k),'');else,xticklabels(h(k),tickLabels);xlabel(h(k),'UTC');end
end
for k=1:4
    name=fields{k};panel=[1 2 3 4];s='No released data';
    lo=coverage.(name).surveyFast.records>0;hi=coverage.(name).burst.records>0;
    if lo&&hi,s='burst + survey/fast';elseif hi,s='burst';elseif lo,s='survey/fast';end
    text(h(panel(k)),1.015,.38,s,'Units','normalized','FontSize',8,'Interpreter','none');
end
annotation(fn,'textbox',[.105 .002 .85 .021], ...
    'String','L2 CDF | Adapted from Overview_download.m | Dashed lines: supplied event bounds | Gaps remain blank', ...
    'EdgeColor','none','FontSize',8,'Interpreter','none');
end
