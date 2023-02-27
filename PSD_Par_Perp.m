clear;
clc;
close all

mms.db_init('local_file_db','D:\MMS\');

ic=1;
Tintr = irf.tint('2021-07-21T13:19:55.00Z/2021-07-21T13:20:15.00Z');
% Tintr = irf.tint('2017-08-23T15:38:20Z/2017-08-23T15:39:15Z');
% Tintr=irf.tint('2017-07-06T03:08:50.00Z/2017-07-06T03:09:23.00Z');
% Tintr = irf.tint('2020-07-09T16:49:35.00Z/2020-07-09T16:49:50.00Z');

%% Load data
%
tic;
c_eval('diste = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',Tintr);',ic);
%theta=diste.DEPEND_2.data;
c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',Tintr);',ic);
c_eval('energy=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy_brst'',Tintr);',ic);

% c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
c_eval('phi=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_phi_brst'',Tintr);',ic);
% c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);

c_eval('stepTable=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_steptable_parity_brst'',Tintr);',ic);
toc;

Bxyz=mms.get_data('B_gse_brst',Tintr,ic);
% Bxyz = Bxyz.resample(diste);
%
% c_eval('diste = mms.db_get_variable(''mms?_fpi_fast_l2_des-dist'',''mms?_des_dist_fast'',Tintr);',ic);
% theta=diste.DEPEND_2.data;
% c_eval('diste = mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_dist_fast'',Tintr);',ic);
% c_eval('energy=mms.db_get_variable(''mms?_fpi_fast_l2_des-dist'',''mms?_des_energy_fast'',Tintr);',ic);
% 
% % c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
% % c_eval('phi=mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_phi_fast'',Tintr);',ic);
% % c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);
% 
% c_eval('stepTable=mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_steptable_parity_fast'',Tintr);',ic);
% toc;
% 
% Bxyz=mms.get_data('B_gse_brst',Tintr,ic);

%% Produce a single PAD at a selected time
% tint1 = irf_time('2021-07-21T13:19:58.00Z','utc>epochTT');
tint1=irf.tint('2021-07-21T13:19:59.200Z/2021-07-21T13:19:59.400Z');
[paddist01,thetapad,energypad,tintpad1] = mms_get_pitchangledist_my_change(diste,Bxyz,tint1); 
paddist01 = paddist01*1e30; 
% [paddist1,thetapad,energypad,tintpad1] = mms_get_pitchangledist_my_change(diste,Bxyz,tint1); 
% paddist1 = paddist1*1e30; 

paddist1=irf.nanmean(paddist01.data);
% paddist1=irf.nanmean(paddist01);
paddist1=squeeze(paddist1);
% % % %
% % % paddist1=paddist01;
for ii=1:32
    sum_paddist1(ii,:)=nansum(paddist1(ii,:));
    paddist1_0(ii,:) = nansum(paddist1(ii,1:6)); %0-33
    paddist1_180(ii,:) = nansum(paddist1(ii,25:30)); %147-180
    paddist1_90(ii,:) = nansum(paddist1(ii,13:18)); %75-105
end




%% Plot PAD

h(1)=irf_subplot(1,1,1);

loglog(energypad(1,:),paddist1_0,'b');hold on
loglog(energypad(1,:),paddist1_90,'g');hold on
loglog(energypad(1,:),paddist1_180,'r');hold on


% set(h(1),'Xlim',[1e1 1e5]);
% 
% set(h(1),'Ylim',[1e-20 1e6]);
set(h(1),'yscale','log');set(h(1),'xscale','log');
% end

ylabel('PSD (s^3 km^{-6})');
xlabel('E (ev)');

set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'0','90','180'},[0.05 0.92]);

% % % ylabel(h(1),' PSD (s^3 km^{-6})','FontSize',15);
% % % xlabel(h(1),'\theta (deg.)','FontSize',15)
% % % set(h(1),'yscale','log');
% % % set(h(1),'Ylim',[1e-4 4*1e0],'FontSize',12);
% % % % set(h(1),'Ylim',[3e-6 3e1],'FontSize',12);
irf_zoom(h(1),'x',[0 180]);
box on;
% irf_zoom(h(2),'y',yrange);
tintutc1 = tintpad1.utc;
 
% title(h(1),strcat('MMS',num2str(ic),'-',tintutc1(12:19),'UT'));
title(h(1),strcat('MMS',num2str(ic),'-',tintutc1(23:2:38),'UT'));

%%

set(gcf,'color','w');
set(gcf,'render','painters');
set(h,'box','on')

figname=['W3'];

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf, '-dpdf', [figname '.pdf']);


