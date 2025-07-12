clear;clc

TT = '2024-05-10T12:00:00Z\2024-05-11T00:00:00Z';
Tsta = strsplit(TT,'\');
Tend = Tsta{2};Tsta = Tsta{1};
%% load data
tint = [iso2epoch(Tsta),iso2epoch(Tend)];
thb = dataobj('/Users/fwd/Documents/Ti~mor~/M/3-Superstorm_GRL/thb_l2_esa_20240510_v01.cdf');
N_thb=getmat(thb,'thb_peef_density');
Ni_thb=getmat(thb,'thb_peif_density');
V_thb=getmat(thb,'thb_peif_velocity_gse');
V_thb=irf_abs(V_thb);
N_thb = irf_resamp(N_thb,V_thb);
units = irf_units;
Pv_thb=[N_thb(:,1)   1e9*1e6*0.5*units.mp.*N_thb(:,2).*(1e3.*V_thb(:,5)).^2];
N_omni= irf_get_data_omni_modified(tint,'n','omni_min');

%% Plot
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

n = 3;i = 1;
h(i)=irf_subplot(n,1,-i);
irf_plot([N_thb(:,1) N_thb(:,2)],'b');hold on;
irf_plot([Ni_thb(:,1) Ni_thb(:,2)],'g');hold on;
irf_plot([N_omni(:,1) N_omni(:,2)],'r');
ylim([0 80])
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'THB-Ne','THB-Ni','OMNI-Np'},[0.05 0.92]);
ylabel('N [cm^{-3}]')
i = i+1;

h(i)=irf_subplot(n,1,-i);
irf_plot([V_thb(:,1) V_thb(:,5)],'k');
ylim([300 900])
ylabel('Vi [km/s]')
i = i+1;

h(i)=irf_subplot(n,1,-i);
irf_plot(Pv_thb,'k');
ylim([0 40])
ylabel('Pdyn [nPa]')
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)

set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')