clear;
clc;
mms.db_init('local_file_db','D:\MMS\');
Units=irf_units;
me=Units.me;
% tint=irf.tint('2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z');
tint=irf.tint('2018-09-08T14:50:43.00Z/2018-09-08T14:54:00.00Z');



for ic=1:2
% load B


c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);

c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['Vi?_gsm=irf_gse2gsm(Vi?);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);


end


 %% Init figure
n=8;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


%% B plot
h(i)=irf_subplot(n,1,-i);
% irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('B [nT]','fontsize',10);
% set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
set(gca,'Ylim',[-5 25], 'ytick',[ 0 5 10 15 20]);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z'},[0.05 0.92]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

  %% Vi plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Vi1(:,1) Vi1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Vi_x Vi_y [km/s]','fontsize',10);
set(gca,'Ylim',[fix(min([min(Vi1_gsm(:,2)) min(Vi1_gsm(:,3)) min(Vi1_gsm(:,4))])/10)*10-10 fix(max(Vit1(:,2))/10)*10+10]);
%set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y'},[0.05 0.92]);
i=i+1;


%% Vi plot
%mms1
 h(i)=irf_subplot(n,1,-i);
ylabel('Vi [km/s]','fontsize',10);
a(:,1)=0:0.1498:60;
Vi1_gsm(:,5)=0;
b(:,1)=Vi1_gsm(:,5);
kk = length(a(:,1))-1;
quiver(a(1:6:kk,1),b(1:6:kk,1),Vi1_gsm(1:6:kk,2),Vi1_gsm(1:6:kk,3),0.3,'color','k');
% ylabel(h(2),'Vix Viy','fontsize',10);
set(gca,'ColorOrder',[0 0 0]);
set(gca,'xticklabel',[0 0.1 0.2])
% irf_legend(gca,{'mms1'},[0.05 0.05]);
set(gca,'ytick',[-2 2])
i=i+1;

% mms2
% h(i)=irf_subplot(n,1,-i);
% ylabel('Vi [km/s]','fontsize',10);
% c(:,1)=Vi2(:,2)';
% d(:,1)=Vi2(:,3)';
% % a(:,1)=0:0.1498:89.8;
% Vi2(:,5)=0
% b(:,1)=Vi2(:,5)
% kk = length(a(:,1));
% hi=quiver(a(1:3:kk,1),b(1:3:kk,1),c(1:3:kk,1),d(1:3:kk,1),'b');
% set(hi,'maxheadsize',0);
% ylabel(h(3),'Vix Viy','fontsize',10);
% set(gca,'ColorOrder',[0 0 0]);
% % irf_legend(gca,{'mms2'},[0.05 0.92]);
% set(gca,'ytick',[-1:1:2])
% i=i+1;

% mms3
% h(i)=irf_subplot(n,1,-i);
% ylabel('Vi [km/s]','fontsize',10);
% e(:,1)=Vi3(:,2)';
% f(:,1)=Vi3(:,3)';
% % a(:,1)=0:0.1498:110;
% Vi3(:,5)=0
% g(:,1)=Vi3(:,5)
% kk = length(a(:,1));
% quiver(a(1:3:kk,1),g(1:3:kk,1),f(1:3:kk,1),e(1:3:kk,1),'g');
% ylabel(h(4),'Vix Viy','fontsize',10);
% set(gca,'ColorOrder',[0 0 0]);
% irf_legend(gca,{'mms3'},[0.05 0.92]);
% set(gca,'ytick',[-2:2:4])
% i=i+1;

% mms4
% h(i)=irf_subplot(n,1,-i);
% ylabel('Vi [km/s]','fontsize',10);
% c(:,1)=Vi4(:,2)';
% d(:,1)=Vi4(:,3)';
% % a(:,1)=0:0.1498:110;
% Vi4(:,5)=0
% b(:,1)=Vi4(:,5)
% kk = length(a(:,1));
% quiver(a(1:3:kk,1),b(1:3:kk,1),d(1:3:kk,1),c(1:3:kk,1),'b');
% ylabel(h(5),'Vix Viy','fontsize',10);
% set(gca,'ColorOrder',[0 0 0]);
% irf_legend(gca,{'mms4'},[0.05 0.92]);
% set(gca,'ytick',[-1:1:3])
% i=i+1;




%%
% n=1;
% i=1;
% set(0,'DefaultAxesFontSize',8);
% set(0,'DefaultLineLineWidth', 0.5);
% fn=figure(1);clf;
% set(gcf,'PaperUnits','centimeters')
% xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% h(i)=irf_subplot(n,1,-i);
% a(:,1)=0:0.1498:110;
% b(:,1)=Vi1(:,1)*0+2;
% c(:,1)=Vi1(:,1)*0-14;
% quiver3(a(1:3:kk,1),b(1:3:kk,1),c(1:3:kk,1),Vi1(1:3:kk,4),Vi1(1:3:kk,2),Vi1(1:3:kk,3),'k');hold on;
% xlabel('z');
% ylabel('x');
% zlabel('y');
% b(:,1)=Vi2(:,1)*0-12;
% c(:,1)=Vi2(:,1)*0+6;
% quiver3(a(1:3:kk,1),b(1:3:kk,1),c(1:3:kk,1),Vi2(1:3:kk,4),Vi2(1:3:kk,2),Vi2(1:3:kk,3),'r');hold on;
% b(:,1)=Vi3(:,1)*0+3;
% c(:,1)=Vi3(:,1)*0+3;
% quiver3(a(1:3:kk,1),b(1:3:kk,1),c(1:3:kk,1),Vi3(1:3:kk,4),Vi3(1:3:kk,2),Vi3(1:3:kk,3),'g');hold on;
% b(:,1)=Vi4(:,1)*0+6;
% c(:,1)=Vi4(:,1)*0+1;
% quiver3(a(1:3:kk,1),b(1:3:kk,1),c(1:3:kk,1),Vi4(1:3:kk,4),Vi4(1:3:kk,2),Vi4(1:3:kk,3),'b');hold off;




%   set(h(1:n),'fontsize',8);
  irf_zoom(tint,'x',h(1:3));
%irf_adjust_panel_position;
  irf_plot_axis_align(h)
  

 %% 
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['vxyz0724s'];
% print(gcf, '-dpdf', [figname '.pdf']);