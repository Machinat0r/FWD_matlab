function Plot_BV(tint,IC,Name,flagTime)
% see also SDCFilenames,SDCFilesDownload,SDCDataMove
%% 底下就是原来的overview程序
global OutputDir

ic=IC;
% load B
c_eval(['B?_ts=mms.get_data(''B_gsm_srvy'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
      
% load FPI
% c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_fast_l2'', tint, ?);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['Vi?_gsm=irf_gse2gsm(Vi?);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);

%% Init figure
n=2;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 40; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% B plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)], 'k--', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
%irf_plot([Bbf1(:,1) Bbf1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))]))-1 fix(max(Bt1(:,2)))+1]);
%     set(gca,'Ylim',[-12 22], 'ytick',[-5 0 5 10 15 20 25],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)

i=i+1; 
%% Vi plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[fix(min([min(Vi1_gsm(:,2)) min(Vi1_gsm(:,3)) min(Vi1_gsm(:,4))])/10)*10-10 fix(max(Vit1(:,2))/10)*10+10],'fontsize',9);
%set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%%
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)

irf_pl_mark(h(1:n),flagTime,'k')

%%  出图保存部分
    set(gcf,'render','painters');
    set(gcf,'paperpositionmode','auto')
    % figname = [OutputDir,'OverviewFig\',Name(2:end-2)];    
    mkdir([OutputDir,'Figs\']);
    figname = [OutputDir,'Figs\',Name];  
    colormap(jet)
    print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
    close all
end