clear;
clc;
close all

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
mms.db_init('local_file_db', ParentDir);

ic=1;

Tintr = irf.tint('2021-08-22T06:40:00.00Z/2021-08-22T06:42:00.00');

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
% tint1 = irf_time('2017-06-11T01:59:25.00Z','utc>epochTT');
tint1=irf.tint('2021-08-22T06:40:33.00Z/2021-08-22T06:40:36.00');
[paddist01,thetapad,energypad,tintpad1] = mms_get_pitchangledist_my_change(diste,Bxyz,tint1); 
paddist01 = paddist01*1e30; 
%
paddist1=irf.nanmean(paddist01.data);
paddist1=squeeze(paddist1);
%

for ii=1:32,
    sum_paddist1(ii,:)=nansum(paddist1(ii,:));
end


% tint2 = irf_time('2017-06-19T03:52:10.00Z','utc>epochTT');
tint2=irf.tint('2021-08-22T06:41:13.00Z/2021-08-22T06:41:16.00');
[paddist02,thetapad,energypad,tintpad2] = mms_get_pitchangledist_my_change(diste,Bxyz,tint2); 
paddist02 = paddist02*1e30;
%
paddist2=irf.nanmean(paddist02.data);
paddist2=squeeze(paddist2);
%
for ii=1:32,
    sum_paddist2(ii,:)=nansum(paddist2(ii,:));
end
%
% tint3 = irf_time('2017-06-19T03:52:15.00Z','utc>epochTT');
tint3=irf.tint('2021-08-22T06:41:43.00Z/2021-08-22T06:41:46.00');
[paddist03,thetapad,energypad,tintpad3] = mms_get_pitchangledist_my_change(diste,Bxyz,tint3); 
paddist03 = paddist03*1e30; 
%
paddist3=irf.nanmean(paddist03.data);
paddist3=squeeze(paddist3);
%

for ii=1:32,
    sum_paddist3(ii,:)=nansum(paddist3(ii,:));
end

%
% tint4 = irf_time('2017-06-19T03:52:19.00Z','utc>epochTT');
% % % tint4=irf.tint('2017-06-11T01:59:50Z/2017-06-11T01:59:53Z');
% % % [paddist04,thetapad,energypad,tintpad4] = mms_get_pitchangledist_my_change(diste,Bxyz,tint4); 
% % % paddist04 = paddist04*1e30; 
% % % %
% % % paddist4=irf.nanmean(paddist04.data);
% % % paddist4=squeeze(paddist4);
% % % %
% % % 
% % % for ii=1:32,
% % %     sum_paddist4(ii,:)=nansum(paddist4(ii,:));
% % % end

% tint5 = irf_time('2017-06-19T03:52:19.00Z','utc>epochTT');
% tint5=irf.tint('2017-06-11T01:59:50Z/2017-06-11T01:59:55Z');
% [paddist05,thetapad,energypad,tintpad5] = mms_get_pitchangledist_my_change(diste,Bxyz,tint5); 
% paddist05 = paddist05*1e30; 
% %
% paddist5=irf.nanmean(paddist05.data);
% paddist5=squeeze(paddist5);
% %
% 
% for ii=1:32,
%     sum_paddist5(ii,:)=nansum(paddist5(ii,:));
% end

paddist_1=paddist1(17:32,:);
paddist_2=paddist2(17:32,:);
paddist_3=paddist3(17:32,:);
% paddist_4=paddist4(17:32,:);

% paddist_5=paddist5(17:32,:);
%% bi-Maxwellian function
ne1t=4.0;
ne2t=15;
ne3t=1.5;
% ne4t=0.05    ;
vth1=2.5995e8;
vth2=4.3831e8;
vth3=1.3660e8;  %1.366明显错了
% vth4=9.2e8;


Tperpoverpar1t=0.4595;
Tperpoverpar2t=1.2582;
Tperpoverpar3t=1.1212;
% Tperpoverpar4t=1.1280;
Me=9.3e-31;  
E1=energypad(1:21);
vpar1=sqrt(2*E1*1.6e-19/Me)*100;    
vd1=-8e7;
vd2=3e7;   %to Scott 3e7
vd3=3e7;   %to Scott 3e7
% vd4=-7e7;
% c_eval('vdnorm?=vd?/vth?;',1:4);
% c_eval('temp?=(vth?/4.19e7)^2;',1:4);
% c_eval('temppar?=3*temp?/(1+2*Tperpoverpar?t);',1:4);
% %% 180°
% f180=ne1t/(sqrt(pi)*vth1).^3*exp(-((-vpar1-vd1)/vth1).^2)*(Tperpoverpar1t)^-1+ne2t/(sqrt(pi)*vth2).^3*exp(-((-vpar1-vd2)/vth2).^2)*(Tperpoverpar2t)^-1+ne3t/(sqrt(pi)*vth3).^3*exp(-((-vpar1-vd3)/vth3).^2)*(Tperpoverpar3t)^-1;
% 
% % f180=ne1t/(sqrt(pi)*vth1).^3*exp(-((-vpar1-vd1)/vth1).^2)*(Tperpoverpar1t)^-1+ne2t/(sqrt(pi)*vth2).^3*exp(-((-vpar1-vd2)/vth2).^2)*(Tperpoverpar2t)^-1+ne3t/(sqrt(pi)*vth3).^3*exp(-((-vpar1-vd3)/vth3).^2)*(Tperpoverpar3t)^-1++ne4t/(sqrt(pi)*vth4).^3*exp(-((-vpar1-vd4)/vth4).^2)*(Tperpoverpar4t)^-1;
% f180=f180*1e30;
% %% 0°
% f0=ne1t/(sqrt(pi)*vth1).^3*exp(-((vpar1-vd1)/vth1).^2)*(Tperpoverpar1t)^-1+ne2t/(sqrt(pi)*vth2).^3*exp(-((vpar1-vd2)/vth2).^2)*(Tperpoverpar2t)^-1+ne3t/(sqrt(pi)*vth3).^3*exp(-((vpar1-vd3)/vth3).^2)*(Tperpoverpar3t)^-1;
% 
% % f0=ne1t/(sqrt(pi)*vth1).^3*exp(-((vpar1-vd1)/vth1).^2)*(Tperpoverpar1t)^-1+ne2t/(sqrt(pi)*vth2).^3*exp(-((vpar1-vd2)/vth2).^2)*(Tperpoverpar2t)^-1+ne3t/(sqrt(pi)*vth3).^3*exp(-((vpar1-vd3)/vth3).^2)*(Tperpoverpar3t)^-1+ne4t/(sqrt(pi)*vth4).^3*exp(-((vpar1-vd4)/vth4).^2)*(Tperpoverpar4t)^-1;
% f0=f0*1e30;
% %% 90°
% vperp1=vpar1;
% 
% f90=ne1t/(sqrt(pi)*vth1).^3*exp(-(-vd1/vth1).^2)*(Tperpoverpar1t)^-1*exp(-(vperp1.^2/(vth1).^2*(Tperpoverpar1t)^-1))+ne2t/(sqrt(pi)*vth2).^3*exp(-(-vd2/vth2).^2)*(Tperpoverpar2t)^-1*exp(-(vperp1.^2/(vth2).^2*(Tperpoverpar2t)^-1))+ne3t/(sqrt(pi)*vth3).^3*exp(-(-vd3/vth3).^2)*(Tperpoverpar3t)^-1*exp(-(vperp1.^2/(vth3).^2*(Tperpoverpar3t)^-1));
% 
% % f90=ne1t/(sqrt(pi)*vth1).^3*exp(-(-vd1/vth1).^2)*(Tperpoverpar1t)^-1*exp(-(vperp1.^2/(vth1).^2*(Tperpoverpar1t)^-1))+ne2t/(sqrt(pi)*vth2).^3*exp(-(-vd2/vth2).^2)*(Tperpoverpar2t)^-1*exp(-(vperp1.^2/(vth2).^2*(Tperpoverpar2t)^-1))+ne3t/(sqrt(pi)*vth3).^3*exp(-(-vd3/vth3).^2)*(Tperpoverpar3t)^-1*exp(-(vperp1.^2/(vth3).^2*(Tperpoverpar3t)^-1))+ne4t/(sqrt(pi)*vth4).^3*exp(-(-vd4/vth4).^2)*(Tperpoverpar4t)^-1*exp(-(vperp1.^2/(vth4).^2*(Tperpoverpar4t)^-1));
% f90=f90*1e30;
%%
for k=1:32
E=energypad(k);
j=thetapad;
vpar=sqrt(2*E*1.6e-19/Me)*cosd(j)*100;
vperp=sqrt(2*E*1.6e-19/Me)*sind(j)*100;
c_eval('f?=ne1t/(sqrt(pi)*vth1).^3.*exp(-((vpar-vd1)/vth1).^2).*(Tperpoverpar1t)^-1.*exp(-(vperp.^2/(vth1).^2.*(Tperpoverpar1t)^-1))+ne2t/(sqrt(pi)*vth2).^3.*exp(-((vpar-vd2)/vth2).^2).*(Tperpoverpar2t)^-1.*exp(-(vperp.^2/(vth2).^2.*(Tperpoverpar2t)^-1))+ne3t/(sqrt(pi)*vth3).^3.*exp(-((vpar-vd3)/vth3).^2).*(Tperpoverpar3t)^-1.*exp(-(vperp.^2/(vth3).^2*(Tperpoverpar3t)^-1));',k);

% c_eval('f?=ne1t/(sqrt(pi)*vth1).^3.*exp(-((vpar-vd1)/vth1).^2).*(Tperpoverpar1t)^-1.*exp(-(vperp.^2/(vth1).^2.*(Tperpoverpar1t)^-1))+ne2t/(sqrt(pi)*vth2).^3.*exp(-((vpar-vd2)/vth2).^2).*(Tperpoverpar2t)^-1.*exp(-(vperp.^2/(vth2).^2.*(Tperpoverpar2t)^-1))+ne3t/(sqrt(pi)*vth3).^3.*exp(-((vpar-vd3)/vth3).^2).*(Tperpoverpar3t)^-1.*exp(-(vperp.^2/(vth3).^2*(Tperpoverpar3t)^-1))+ne4t/(sqrt(pi)*vth4).^3.*exp(-((vpar-vd4)/vth4).^2).*(Tperpoverpar4t)^-1.*exp(-(vperp.^2/(vth4).^2*(Tperpoverpar4t)^-1));',k);
c_eval('f?=f?*1e30;',k);
end


%% Plot PAD

fn=figure;
set(fn,'Position',[100 20 1000 600])
    h(1)=axes('position',[0.06 0.19 0.18 0.65]); 
    h(2)=axes('position',[0.40 0.19 0.18 0.65]); 
    h(3)=axes('position',[0.74 0.19 0.18 0.65]); 
%     h(4)=axes('position',[0.72 0.19 0.18 0.65]);
    
%     h(5)=axes('position',[0.84 0.19 0.18 0.65]); 
    
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',1.5); 
%%
% set(h(1),'defaultLineLineWidth',0.5); 
% ymin = 10^-1;
% ymax = ceil(max(max(log10(sum_paddist1))));
% yrange = [ymin 10^ymax];
% 
% loglog(h(1),energypad,sum_paddist1,'k.-');
% hold(h(1),'on');
% loglog(h(1),energypad,sum_paddist2,'r.-');
% 
% % hold(h(1),'on');
% % plot(h(1),E1,f0,'k');
% % plot(h(1),E1,f180,'b');
% % plot(h(1),E1,f90,'r');
% hold(h(1),'off');
%  
% ylabel(h(1),'f_e (s^3 km^{-6})');
% xlabel(h(1),'E (eV)');
% % set(h(1),'yscale','log');
% % set(h(1),'xscale','log');
% % irf_zoom(h(1),'y',yrange);
% irf_zoom(h(1),'x',[0 4e4]);
% tintutc1 = tintpad1.utc;
% tintutc2 = tintpad2.utc;
% irf_legend(h(1),strcat(tintutc1(23:2:38),'UT'),[0.96, 0.98],'color','k');
% irf_legend(h(1),strcat(tintutc2(23:2:38),'UT'),[0.96, 0.93],'color','r');
%%
jetcolor = colormap('jet');

% c_eval('plot(h(2),thetapad,paddist(?,:),''color'',jetcolor(?+9,:));',1);
hold(h(1),'on');
% c_eval('plot(h(2),j,f?,''color'',jetcolor(?+2,:));',[2:18]);
% paddist_1=paddist1(13:32,:);
c_eval('plot(h(1),thetapad,paddist_1(?,:),''color'',jetcolor(15*?+12,:));',1:16);
% hold(h(1),'on');
% c_eval('plot(h(1),thetapad,paddist_2(?,:),''--'',''color'',jetcolor(3*?+12,:));',1:16);
% for i=1;32,
% plot(h(2),thetapad,paddist(i,:),'color',jetcolor(2*i,:) );
% end
% c_eval('plot(h(2),thetapad,paddist(?,:));',[1:32]);
hold(h(1),'off');

% 
% for i=1:length(energypad(1:32))
%     a=roundn(energypad(1,:)/1000,-3);
% c_eval('irf_legend(h(2),''? keV'',[1.02 1.03-i*0.033],''color'',jetcolor(2*i,:));',a(i));
% end

ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'\theta (deg.)')
set(h(1),'yscale','log');
set(h(1),'Ylim',[1e-4 5*1e1]);
irf_zoom(h(1),'x',[0 180]);
box on;
% irf_zoom(h(2),'y',yrange);
tintutc1 = tintpad1.utc;
 
% title(h(1),strcat('MMS',num2str(ic),'-',tintutc1(12:23),'UT'));
title(h(1),strcat('MMS',num2str(ic),'-',tintutc1(23:2:38),'UT'));

%%
jetcolor = colormap('jet');
% set(h(2),'defaultLineLineWidth',1.5); 
% c_eval('plot(h(2),thetapad,paddist(?,:),''color'',jetcolor(?+9,:));',1);
hold(h(2),'on');
% c_eval('plot(h(2),j,f?,''color'',jetcolor(?+2,:));',[2:18]);
% paddist_1=paddist1(13:32,:);
c_eval('plot(h(2),thetapad,paddist_2(?,:),''color'',jetcolor(15*?+12,:));',1:16);
% hold(h(1),'on');
% c_eval('plot(h(2),thetapad,paddist_4(?,:),''--'',''color'',jetcolor(3*?+12,:));',1:16);
% for i=1;32,
% plot(h(2),thetapad,paddist(i,:),'color',jetcolor(2*i,:) );
% end
% c_eval('plot(h(2),thetapad,paddist(?,:));',[1:32]);
hold(h(2),'off');

% 
% for i=1:length(energypad(1:32))
%     a=roundn(energypad(1,:)/1000,-3);
% c_eval('irf_legend(h(2),''? keV'',[1.02 1.03-i*0.033],''color'',jetcolor(2*i,:));',a(i));
% end

% ylabel(h(2),'f_e (s^3 km^{-6})');
xlabel(h(2),'\theta (deg.)')
set(h(2),'yscale','log');
set(h(2),'Ylim',[1e-4 5*1e1]);
irf_zoom(h(2),'x',[0 180]);
box on;
% irf_zoom(h(2),'y',yrange);
tintutc2 = tintpad2.utc;
 
% title(h(2),strcat('MMS',num2str(ic),'-',tintutc2(12:23),'UT'));
title(h(2),strcat('MMS',num2str(ic),'-',tintutc2(23:2:38),'UT'));
%%
jetcolor = colormap('jet');
hold(h(3),'on');
%  set(h(3),'defaultLineLineWidth',1.5); 
c_eval('plot(h(3),thetapad,paddist_3(?,:),''color'',jetcolor(15*?+12,:));',1:16);

hold(h(3),'off');


for i=1:length(energypad(17:32))
    a=roundn(energypad(1,17:32)/1000,-3);
c_eval('irf_legend(h(3),''? keV'',[1.02 0.86-i*0.045],''color'',jetcolor(15*i+12,:));',a(i));
end

% ylabel(h(3),'f_e (s^3 km^{-6})');
xlabel(h(3),'\theta (deg.)')
set(h(3),'yscale','log');
set(h(3),'Ylim',[1e-4 5*1e1]);
irf_zoom(h(3),'x',[0 180]);
% irf_zoom(h(2),'y',yrange);
tintutc3 = tintpad3.utc;
% title(h(3),strcat('MMS',num2str(ic),'-',tintutc3(12:23),'UT'));
title(h(3),strcat('MMS',num2str(ic),'-',tintutc3(23:2:38),'UT'));

%%
% % % jetcolor = colormap('jet');
% % % 
% % % % c_eval('plot(h(2),thetapad,paddist(?,:),''color'',jetcolor(?+9,:));',1);
% % % hold(h(4),'on');
% % % % c_eval('plot(h(2),j,f?,''color'',jetcolor(?+2,:));',[2:18]);
% % % % paddist_1=paddist1(13:32,:);
% % % c_eval('plot(h(4),thetapad,paddist_4(?,:),''color'',jetcolor(3*?+12,:));',1:16);
% % % % for i=1;32,
% % % % plot(h(2),thetapad,paddist(i,:),'color',jetcolor(2*i,:) );
% % % % end
% % % % c_eval('plot(h(2),thetapad,paddist(?,:));',[1:32]);
% % % hold(h(4),'off');
% % % 
% % % for i=1:length(energypad(17:32))
% % %     a=roundn(energypad(1,17:32)/1000,-3);
% % % c_eval('irf_legend(h(4),''? keV'',[1.02 0.86-i*0.045],''color'',jetcolor(3*i+12,:));',a(i));
% % % end
% % % 
% % % 
% % % % ylabel(h(4),'f_e (s^3 km^{-6})');
% % % xlabel(h(4),'\theta (deg.)')
% % % set(h(4),'yscale','log');
% % % set(h(4),'Ylim',[1e-4 5*1e1]);
% % % irf_zoom(h(4),'x',[0 180]);
% % % box on;
% % % % irf_zoom(h(2),'y',yrange);
% % % tintutc4 = tintpad4.utc;
% % %  
% % % % title(h(4),strcat('MMS',num2str(ic),'-',tintutc4(12:23),'UT'));
% % % title(h(4),strcat('MMS',num2str(ic),'-',tintutc4(23:2:38),'UT'));
%%

set(gcf,'color','w');
set(gcf,'render','painters');
set(h,'box','on')

figname=['W3-P2'];

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf, '-dpdf', [figname '.pdf']);


