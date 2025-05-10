clc,clear;
tint=irf.tint('2008-02-22T06:25:30Z/2008-02-22T06:26:30Z');
%% load
tem=dataobj('E:\data_MMS\T20080222\thd_l2_esa_20080222_v01.cdf')
tem2=dataobj('E:\data_MMS\T20080222\thd_l2_fgm_20080222_v01.cdf')
tem3=dataobj('E:\data_MMS\T20080222\thd_l2_efi_20080222_v01.cdf')

thd_fgs_gsm=getmat(tem2,'thd_fgs_gsm');
thd_peer_density=getmat(tem,'thd_peer_density');
thd_peir_velocity_gsm=getmat(tem,'thd_peir_velocity_gsm');
thd_peer_avgtemp=getmat(tem,'thd_peer_avgtemp');
thd_peir_en_eflux_yaxis=getmat(tem,'thd_peir_en_eflux_yaxis');
thd_peir_en_eflux=getmat(tem,'thd_peir_en_eflux');
thd_eff_dot0_dsl=getmat(tem3,'thd_eff_dot0_dsl');
 %% Init figure
n_subplots=6;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize =18; ySize = 10; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% plot Bxyz
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
 irf_plot([thd_fgs_gsm(:,1) thd_fgs_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([thd_fgs_gsm(:,1) thd_fgs_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([thd_fgs_gsm(:,1) thd_fgs_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([thd_fgs_gsm(:,1) thd_fgs_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Bz [nT]','fontsize',10);
set(gca,'Ylim',[-20 30], 'ytick',[-15 0 15]);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
%% plot N
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([thd_peer_density(:,1) thd_peer_density(:,2)], 'color','k', 'Linewidth',0.75);hold on;
grid off;
ylabel('N [cm^{-3}]','fontsize',10);
 set(gca,'Ylim',[0.1 0.7], 'ytick',[ 0.2 0.4 0.6]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)

%% plot T
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([thd_peer_avgtemp(:,1) thd_peer_avgtemp(:,2)], 'color','k', 'Linewidth',0.75);hold on;
grid off;
ylabel('Te [ev]','fontsize',10);
set(gca,'Ylim',[0 5000], 'ytick',[ 1000 3000 5000]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)




%% plot Vi
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([thd_peir_velocity_gsm(:,1) thd_peir_velocity_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([thd_peir_velocity_gsm(:,1) thd_peir_velocity_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([thd_peir_velocity_gsm(:,1) thd_peir_velocity_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([thd_peir_velocity_gsm(:,1) thd_peir_velocity_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Vi [km/s]','fontsize',10);
 set(gca,'Ylim',[-400 400], 'ytick',[-600 -400 -200 0 200 400 600]);
irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.1 0.82]);

%% plot i energy spectrom
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(5),jet)
specrec_p_i=struct('t',thd_peir_en_eflux_yaxis(:,1));
specrec_p_i.f=transpose(thd_peir_en_eflux_yaxis(1,:));%energy levels
specrec_p_i.p=thd_peir_en_eflux;%data matrix
specrec_p_i.f_label='';
specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
[h(5), hcb5]=irf_spectrogram(h(5),specrec_p_i);
grid off;
ylabel('Ei(ev)','fontsize',10)
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
set(gca,'Ylim',[5e1 1e4]);
caxis(gca,[3 6])
poscbar5=get(hcb5,'pos');
poscbar5(3)=poscbar5(3)*0.5;
set(hcb5,'pos',poscbar5);
irf_legend(gca,'b',[0.99 0.98],'color','w','fontsize',12)

%% plot Exyz
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
 irf_plot([thd_eff_dot0_dsl(:,1) thd_eff_dot0_dsl(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([thd_eff_dot0_dsl(:,1) thd_eff_dot0_dsl(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([thd_eff_dot0_dsl(:,1) thd_eff_dot0_dsl(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([thd_eff_dot0_dsl(:,1) thd_eff_dot0_dsl(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('E [nT]','fontsize',10);
set(gca,'Ylim',[-20 30], 'ytick',[-15 0 15]);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)

set(h(1:5),'fontsize',7);
  irf_adjust_panel_position;
  irf_zoom(tint,'x',h(1:6));
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['0222d'];
print(gcf, '-dpdf', [figname '.pdf']);
