function THEMISPlot(ic,tint,Dir)
%% load data
c_eval("B_? = th_read_l2_change_by_fwd('th?_fgl_gsm',tint);",ic);
c_eval("Vi_? = th_read_l2_change_by_fwd('th?_peir_velocity_gsm',tint);",ic);
% c_eval("E_? = th_read_l2_change_by_fwd('th?_eff_dot0_gsm',tint);",ic);
c_eval("E_? = th_read_l2_change_by_fwd('th?_eff_dot0_gsm',tint);",ic);
c_eval("Vires_? = irf_resamp(Vi_?,B_?);",ic);
c_eval("Eres_? = irf_resamp(E_?,B_?);",ic);
c_eval("Evxb_? = Eres_? + irf_cross(Vires_?,B_?)/1000;",ic);
c_eval("Evxb_?(:,1) = Eres_?(:,1);",ic);
%% Init figure
n=6;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Bx plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B_?(:,2))])-1 max([max(B_?(:,2))])+1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1]]);
% % % irf_legend(gca,{'B_x'},[0.97 0.92]);
% % % ylabel('B [nT]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Bx & Bz plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(B_?(:,2)) min(B_?(:,4))])-3 max([max(B_?(:,2)) max(B_?(:,4))])+3]);",ic);
% set(gca,'Ylim',[-35 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% By plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B_?(:,1) B_?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min(B_?(:,3))-1 max(B_?(:,3))+1]);",ic);
%     set(gca,'Ylim',[-12 10], 'ytick',[-12 -8 -4 0 4 8],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
irf_legend(gca,{'B_y'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%%  Bz plot
% % % h(i)=irf_subplot(n,1,-i);
% % % % c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[0 max([max(B_?(:,4))])+1]);",ic);
% % % % c_eval("set(gca,'Ylim',[min([min(B_?(:,4))])-1 max([max(B_?(:,4))])+1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[1 0 0]]);
% % % irf_legend(gca,{'B_z'},[0.97 0.92]);
% % % ylabel('B [nT]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Vix plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Vi_?(:,1) Vi_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,2))])/10)*10-15 fix(max([max(Vi_?(:,2))])/10)*10+15],'fontsize',9);",ic);
% set(gca,'Ylim',[-220 220], 'ytick',[-200 -100 0 100 200],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1]]);
irf_legend(gca,{'Vi_x'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% Viy plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,3))])/10)*10-15 fix(min([max(Vi_?(:,3))])/10)*10+15],'fontsize',9);",ic);
% set(gca,'Ylim',[-150 220], 'ytick',[-100 0 100 200],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 1 0]]);
irf_legend(gca,{'Vi_y'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% Viz plot
h(i)=irf_subplot(n,1,-i);
% c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,4))])/10)*10-10 fix(max([max(Vi_?(:,4))])/10)*10+10],'fontsize',9);",ic);
% set(gca,'Ylim',[-70 120], 'ytick',[-50 0 50 100],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[1 0 0]]);
irf_legend(gca,{'Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% E+vxb plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,3)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,4)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-100 100], 'ytick',[-100 -50 0 50 100],'fontsize',9);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evxb_?(:,2)) min(Evxb_?(:,3)) min(Evxb_?(:,4))])-2 max([max(Evxb_?(:,2)) max(Evxb_?(:,3)) max(Evxb_?(:,4))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.12]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%%  出图保存部分
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
% Dir =   'C:\Matlab\bin\新建文件夹\fwd\Figures\';
name = datestr(datenum(1970,1,1,0,0,0)+mean(tint)/86400,'yyyymmdd HH:MM:SS');
name = strrep(name,':','');
figname = [Dir strrep(name,' ','') '-' ic];
print(gcf, '-dpdf', [figname '.pdf']);    
clf
end