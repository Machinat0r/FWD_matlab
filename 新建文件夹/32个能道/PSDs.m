% Script to plot electron PSD around pitch angles 0, 90, and 180 deg 
% and PSD versus pitch angle L1b brst data 
%
% Written by D. B. Graham

ic = 1; % Spacecraft number

% cd D:\MATLAB\mms_db
mms.db_init('local_file_db','D:\Matlab\xy-matlab\MMS\mms_db\data')
Tintr = irf.tint('2017-06-11T01:59:25Z/2017-06-11T01:59:55Z');

%% Load data
c_eval('fpiFile? = dataobj(''D:\Matlab\xy-matlab\MMS\mms_db\data\mms1\fpi\brst\l2\des-dist\2017\06\11\mms1_fpi_brst_l2_des-dist_20170611015853_v3.2.0.cdf'');',ic);
c_eval('diste_struct = get_variable(fpiFile?,''mms?_des_dist_brst'');',ic);
theta=diste_struct.DEPEND_2.data;
energy0=diste_struct.DEPEND_3.data(1,1:end);
energy1=diste_struct.DEPEND_3.data(2,1:end);
tic;
c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',Tintr);',ic);
% c_eval('energy0=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy0_brst'',Tintr);',ic);
% c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
c_eval('phi=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_phi_brst'',Tintr);',ic);
% c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);
c_eval('stepTable=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_steptable_parity_brst'',Tintr);',ic);
toc;
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tintr);',ic);

%% Produce a single PAD at a selected time
% c_eval('tint? = irf_time(''2015-12-14T00:59:03.10Z'',''utc>epochTT'');',ic);
% tint1 = irf_time('2015-12-14T00:59:03.000000Z','utc>epochTT');
% [paddist1,thetapad1,energypad1,tintpad1] = mms.get_pitchangledist(diste,phi,theta,stepTable,energy0,energy1,Bxyz,tint1); 
% 
% paddist1 = paddist1*1e30; %convert to commonly used s^3 km^-6

%% Plot PAD

fn=figure;
set(fn,'Position',[10 10 600 250])
    h(1)=axes('position',[0.08 0.08 0.2 0.8]); 
%     h(2)=axes('position',[0.58 0.08 0.2 0.8]); 
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

% ymin = 10^-4;
% ymax = ceil(max(max(log10(paddist1))));
% yrange = [ymin 10^ymax];
% plot(h(1),energypad,paddist(:,1),'k',energypad,mean(paddist(:,[6 7]),2),'r',energypad,paddist(:,12),'b');
% ylabel(h(1),'f_e (s^3 km^{-6})');
% xlabel(h(1),'E (eV)')
% set(h(1),'yscale','log');
% set(h(1),'xscale','log');
% irf_zoom(h(1),'y',yrange);
% irf_zoom(h(1),'x',[10 3e4]);
% irf_legend(h(1),{'0 deg'},[0.91 0.92],'color','k')
% irf_legend(h(1),{'90 deg'},[0.91 0.84],'color','r')
% irf_legend(h(1),{'180 deg'},[0.91 0.76],'color','b')   
%% 
tint1 = irf_time('2017-06-11T01:59:40Z','utc>epochTT');
[paddist1,thetapad1,energypad1,tintpad1] = mms.get_pitchangledist(diste,phi,theta,stepTable,energy0,energy1,Bxyz,tint1); 

paddist1 = paddist1*1e30; %convert to commonly used s^3 km^-6

ymin = 10^-4;
ymax = ceil(max(max(log10(paddist1))));
yrange = [ymin 10^ymax];

jetcolor = colormap('jet');
lll = length(jetcolor(:,1));
vcolors = floor(lll/length(energypad1));
vcolors = [1:length(energypad1)]*vcolors;
    
c_eval('plot(h(1),thetapad1,squeeze(paddist1(?,:)),''color'',jetcolor(vcolors(?),:));',1);
hold(h(1),'on');
c_eval('plot(h(1),thetapad1,squeeze(paddist1(?,:)),''color'',jetcolor(vcolors(?),:));',[2:32]);
hold(h(1),'off')
ylabel(h(1),'f_e (s^3 km^{-6})');
% xlabel(h(1),'\theta (deg.)')
set(h(1),'yscale','log');
irf_zoom(h(1),'x',[0 180]);
irf_zoom(h(1),'y',yrange);
set(h(1),'xtick',[]);
tintutc = tintpad1.utc;
title(h(1),strcat(tintutc(12:23),'UT'));
for i=1:length(energypad1)
    a=roundn(energypad1(1,:)/1000,-3);
c_eval('irf_legend(h(1),''? keV'',[1.01 1-i*0.03],''color'',jetcolor(vcolors(i),:));',a(i));
end

%% 
% set(h(2:6),'ytick',[]);
% set(h(2),'ytick',[]);
set(h(1),'xtick',[0 90 180]);
set(gcf,'color','w');
% print(gcf, '-dpng', ['case20151214_figure71' '.png']);