clear;clc;close all
filepath = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/omni_hro_1min_20190801_v01.cdf';
TT = '2015-10-18T00:00:00.00Z/2015-10-19T00:00:00.00Z';
tint = irf.tint(TT);

dst= irf_get_data_omni_modified(tint,'bsnx','omni_min');
% dst(dst==999999)=0;
h = irf_plot(dst);
ylabel('ts [Re]','fontsize',8);
grid off
irf_zoom(tint,'x',h);

set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')