clear;close all
clc

TT = '2018-12-27T21:00:00.000Z/2018-12-27T23:00:00.000Z';
Tsta = strsplit(TT,'/');
Tend = Tsta{2};Tsta = Tsta{1};
%% load data
tint = [iso2epoch(Tsta),iso2epoch(Tend)];
Bx_omni= irf_get_data_omni_modified(tint,'Bx','omni_min');
By_omni= irf_get_data_omni_modified(tint,'By','omni_min');
Bz_omni= irf_get_data_omni_modified(tint,'Bz','omni_min');
N_omni= irf_get_data_omni_modified(tint,'n','omni_min');
Vx_omni= irf_get_data_omni_modified(tint,'Vx','omni_min');
Vy_omni= irf_get_data_omni_modified(tint,'Vy','omni_min');
Vz_omni= irf_get_data_omni_modified(tint,'Vz','omni_min');
P_omni= irf_get_data_omni_modified(tint,'P','omni_min');
%%
B = [Bx_omni By_omni(:,2) Bz_omni(:,2)];
Vi = [Vx_omni Vy_omni(:,2) Vz_omni(:,2)];
%% plot
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

n = 4;i = 1;

h(i)=irf_subplot(n,1,-i);
irf_plot([B(:,1) B(:,2)],'b');hold on;
irf_plot([B(:,1) B(:,3)],'g');hold on;
irf_plot([B(:,1) B(:,4)],'r');
grid off
ylim([-10 10])
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'Bx','By','Bz'},[0.05 0.92]);
ylabel('B [nT]')
i = i+1;

h(i)=irf_subplot(n,1,-i);
irf_plot([Vi(:,1) Vi(:,2)],'b');hold on;
irf_plot([Vi(:,1) Vi(:,3)],'g');hold on;
irf_plot([Vi(:,1) Vi(:,4)],'r');
grid off
ylim([-400 100])
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'Vx','Vy','Vz'},[0.05 0.92]);
ylabel('Vi [km/s]')
i = i+1;

h(i)=irf_subplot(n,1,-i);
irf_plot(N_omni,'k');
% ylim([0 40])
grid off
ylabel('N [cm^{-3}]')
i=i+1;

h(i)=irf_subplot(n,1,-i);
irf_plot(P_omni,'k');
grid off
% ylim([0 40])
ylabel('Pdyn [nPa]')


irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)

set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')