function fn=Overview_download_events_20260924(P,D,ev,ic,a,b,eventStart,eventEnd,counts,coverage)
% 由 Overview_download.m 复制并保留所需绘图段；Electric field来自Overview_download_mms4.m。
% 2026-09-24：供20事件/四星批处理调用。原脚本保留不动。
% 按用户要求：删除模式时间条，紧凑排列panel，恢复原程序ylabel与能谱色标。
% Ti/Te分开，直接复用原温度三曲线（标量、平行、垂直）；缺测保留。
% B、V为GSM；MMS1-3 E为GSM，MMS4 E仅DSL XY。n:cm^-3，T/energy:eV。
% 输入P为原生观测的绘图数组（burst优先，缺口保留），D保留能谱及原生模式覆盖。
%% 与原脚本相同的变量命名
c_eval('B?=P.B(:,1:4); Bt?=P.B(:,[1 5]);',ic);
c_eval('gsmVi?=P.Vi; gsmVe?=P.Ve; E?=P.E;',ic);
c_eval('Ni?=P.Ni; Ne?=P.Ne; Ti?=P.Ti; Te?=P.Te;',ic);
%% Init figure：沿用原脚本分区、irf_subplot与c_eval绘图
n=9;i=1;
fn=figure('Visible','off','Color','w','Units','pixels','Position',[40 40 1600 1250]);
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
% 负编号irf_subplot沿用原程序紧凑模式，统一边距后仅留0.002图高的间隙。
panelPos=zeros(n,4);step=.83/n;gap=.002;
for j=1:n,panelPos(j,:)=[.09 .90-j*step .77 step-gap];set(h(j),'Position',panelPos(j,:));end
%% B plot：复制原B plot段
axes(h(i));
set(h(i),'Position',panelPos(i,:));
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
set(h(i),'Position',panelPos(i,:));
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
set(h(i),'Position',panelPos(i,:));
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
set(h(i),'Position',panelPos(i,:));
if counts(i)>0
c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
if ic~=4, c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on; end
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
if ic==4,irf_legend(gca,{'E_x','E_y'},[0.97 0.92]);else,irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);end
ylabel('E [mV/m]','fontsize',10);
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% N plot：复制原N plot段
axes(h(i));
set(h(i),'Position',panelPos(i,:));
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

%% Ti plot：复用原Overview_download_mms4.m温度段
axes(h(i));set(h(i),'Position',panelPos(i,:));
if counts(i)>0
c_eval("irf_plot([Ti?(:,1) Ti?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Ti?(:,1) Ti?(:,3)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Ti?(:,1) Ti?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(gca,{'Ti','T_/_/','T_⊥'},[0.97 0.92]);
ylabel('Ti [eV]','fontsize',8);
else
text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% Te plot：复用原Overview_download.m温度段
axes(h(i));set(h(i),'Position',panelPos(i,:));
if counts(i)>0
c_eval("irf_plot([Te?(:,1) Te?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te?(:,1) Te?(:,3)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te?(:,1) Te?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(gca,{'Te','T_/_/','T_⊥'},[0.97 0.92]);
ylabel('Te [eV]','fontsize',8);
else
text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% plot ION energy spectrom：沿用原能谱段
axes(h(i));
set(h(i),'Position',panelPos(i,:));
if counts(i)>0
hold(h(i),'on');colormap(h(i),jet)
for im=1:2
    specrec_p_ION=D.Si{im};
    if isempty(specrec_p_ION),continue;end
    irf_spectrogram(h(i),specrec_p_ION,'log','donotshowcolorbar');
end
grid off;set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ei(ev)','fontsize',8);
set(h(i),'Ylim',[1 4e4]);
hcb=colorbar(h(i));hcb.Position=[.875 panelPos(i,2) .010 panelPos(i,4)];
hcb.Label.String={' ','log10(keV/(cm^2 s sr keV))'};hcb.Label.FontSize=7;hcb.FontSize=8;
set(h(i),'Position',panelPos(i,:));
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;

%% plot e energy spectrom：沿用原能谱段
axes(h(i));
set(h(i),'Position',panelPos(i,:));
if counts(i)>0
hold(h(i),'on');colormap(h(i),jet)
for im=1:2
    specrec_p_e=D.Se{im};
    if isempty(specrec_p_e),continue;end
    irf_spectrogram(h(i),specrec_p_e,'log','donotshowcolorbar');
end
grid off;set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ee(ev)','fontsize',8);
set(h(i),'Ylim',[1 4e4]);
hcb=colorbar(h(i));hcb.Position=[.875 panelPos(i,2) .010 panelPos(i,4)];
hcb.Label.String={' ','log10(keV/(cm^2 s sr keV))'};hcb.Label.FontSize=7;hcb.FontSize=8;
set(h(i),'Position',panelPos(i,:));
else
    text(.5,.5,'No available L2 measurements','Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5]);
end
i=i+1;
%% 统一时间轴与保存前样式，保留原图逐panel组织
span=b-a;if span<=1800,step=300;elseif span<=5400,step=600;else,step=1800;end
tick=(ceil(a/step)*step:step:floor(b/step)*step)-a;
tickLabels=cellstr(datetime(tick+a,'ConvertFrom','posixtime','TimeZone','UTC','Format','HH:mm'));
ylabels={'B [nT]','Vi [km/s]','Ve [km/s]','E [mV/m]','N [cm^{-3}]','Ti [eV]','Te [eV]','Ei(ev)','Ee(ev)'};
for k=1:n,ylabel(h(k),ylabels{k},'FontSize',10);end
for k=1:9
    if k<=7,ylim(h(k),'auto');end
    xlim(h(k),[0 span]);xticks(h(k),tick);grid(h(k),'off');
    set(h(k),'FontName','Arial','FontSize',10,'TickDir','out','TickLength',[.004 .004],'Box','on','Layer','top');hold(h(k),'on');
    for bound=unique([eventStart,eventEnd]),xline(h(k),bound-a,'--','Color',[.35 .35 .35],'LineWidth',.8);end
    if k<9,xticklabels(h(k),[]);xlabel(h(k),'');else,xticklabels(h(k),tickLabels);xlabel(h(k),'UTC');end
end
% 紧凑panel间距下，去掉恰好落在上下边缘的y刻度标签，避免相邻数字重叠。
for k=1:7
    yl=ylim(h(k));yt=yticks(h(k));
    yt=yt(yt>yl(1)+.025*diff(yl)&yt<yl(2)-.025*diff(yl));
    if ~isempty(yt),yticks(h(k),yt);end
end
for k=8:9
    set(h(k),'YScale','log','YLim',[1 4e4],'YTick',[1e1 1e2 1e3 1e4],'YMinorTick','off');
end
for k=1:n
    if counts(k)==0,yticks(h(k),[]);end
end
drawnow;
positions=zeros(n,4);
for k=1:n,positions(k,:)=h(k).Position;end
setappdata(fn,'layoutAudit',struct('axisCount',n,'ylabels',{ylabels}, ...
    'positions',positions,'modeTimeline',false,'temperatureColumns',{{'scalar','parallel','perpendicular'}}));
end
