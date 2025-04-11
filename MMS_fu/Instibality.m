close all
clear;clc

global ParentDir 
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TT = '2015-11-19T12:40:00.000Z/2015-11-19T14:20:00.000Z';
TT = '2015-11-19T14:07:30.000Z/2015-11-19T14:09:30.000Z';
% TT = '2019-03-29T13:23:00.000Z/2019-03-29T13:23:30.000Z';
% TT = '2017-07-18T13:04:50.000Z/2017-07-18T13:05:00.000Z';
% TT = '2015-11-19T13:20:00.000Z/2015-11-19T13:25:00.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); 
filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms, dis-moms');
filenames = [filenames1, filenames2, filenames3, filenames4];

[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
% SDCFilesDownload(filenames,TempDir)

% SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 32, 'CheckSize', 0)
% SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);
units = irf_units;

% load B
c_eval(['B?_ts=mms.get_data(''B_gse_srvy'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);

% load N
c_eval('Ne?_ts = mms.get_data(''Ne_fpi_fast_l2'',tint,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_fast_l2'',tint,?);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);

%load T
c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_fast_l2_des-moms'',''mms?_des_temppara_fast'',tint);',ic);
c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_fast_l2_des-moms'',''mms?_des_tempperp_fast'',tint);',ic);
c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);

c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_temppara_fast'',tint);',ic);
c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_tempperp_fast'',tint);',ic);
c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);

%load V
c_eval('Ve?_ts = mms.get_data(''Ve_gse_fpi_fast_l2'',tint,?);',ic)
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
% c_eval(['gsmVe?_ts=irf_gse2gsm(Ve?_ts);'],ic);
% c_eval(['Ve?=irf.ts2mat(gsmVe?_ts);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);

c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_fast_l2'',tint,?);',ic); 
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
% c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
% c_eval(['Vi?=irf.ts2mat(gsmVi?_ts);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);

% load energy spectra
c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_fast_l2_des-moms'',''mms?_des_pitchangdist_lowen_fast'',tint);',ic)
c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_fast_l2_des-moms'',''mms?_des_pitchangdist_miden_fast'',tint);',ic)
c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_fast_l2_des-moms'',''mms?_des_pitchangdist_highen_fast'',tint);',ic)
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_fast_l2_des-moms'',''mms?_des_energyspectr_omni_fast'',tint);',ic)
c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_energyspectr_omni_fast'',tint);',ic)

 energy_e=energy_e1; energy_e.data=energy_e1.data;  energy_e.nrec=energy_e1.nrec;

% load R
Pos = mms.get_data('R_gse',tint);
c_eval('R? = Pos.gseR?;',ic)
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];',ic)
c_eval('R? = irf_resamp(R?,Bt?);',ic)

%% move mean data
% % % dfB = mean(1/(B1(2,1)-B1(1)));
% % % c_eval('Bmean? = B?;',ic);
% % % dt = dfB*5*60; % 5 min
% % % c_eval('Bmean?(:,2) = movmean(B?(:,2), dt);');
% % % c_eval('Bmean?(:,3) = movmean(B?(:,3), dt);');
% % % c_eval('Bmean?(:,4) = movmean(B?(:,4), dt);');
% % % c_eval('Bmeant? = irf_abs(Bmean?);');
% % % c_eval('Bmeant? = Bmeant?(:,[1,5]);');
% % % % 
% % % % dfe = mean(1/(Ne1(2,1)-Ne1(1)));
% % % % dt = dfe*5*60;
% % % % c_eval('Ne?(:,2) = movmean(Ne?(:,2), dt);');
% % % % c_eval('Te?(:,2) = movmean(Te?(:,2), dt);');
% % % % 
% % % % dfi = mean(1/(Ni1(2,1)-Ni1(1)));
% % % % dt = dfi*5*60;
% % % % c_eval('Ni?(:,2) = movmean(Ni?(:,2), dt);');
% % % % c_eval('Ti?(:,2) = movmean(Ti?(:,2), dt);');
% % % 
% % % % resample
% % % c_eval('Ne? = irf_resamp(Ne?, B?);',ic);
% % % c_eval('Ni? = irf_resamp(Ni?, B?);',ic);
% % % c_eval('Te? = irf_resamp(Te?, B?);',ic);
% % % c_eval('Ti? = irf_resamp(Ti?, B?);',ic);
%% ICI
% % c_eval('Pm? = [Bt?(:,1) 10^(-9)*Bt?(:,2).^2 / (2*units.mu0)];',ic); %nPa
% % c_eval("Pte? = [Ne?(:,1) units.e*1e15*Ne?(:,2).*Te?(:,2)];",ic);
% % c_eval("Pti? = [Ni?(:,1) units.e*1e15*Ni?(:,2).*Ti?(:,2)];",ic);
% % c_eval('beta? = [Bt?(:,1) Pti?(:,2) ./ Pm?(:,2)];',ic);
% % 
% % % 2*deltaB/B
% % c_eval('deltaB? = abs(Bt?(:,2)-Bmeant?(:,2));',ic);
% % c_eval('ICI? = [Bt?(:,1) 2 * deltaB? ./ Bmeant?(:,2)];',ic);
%% cone angel
Bx_omni= irf_get_data_omni_modified(tint,'Bx','omni_min');
By_omni= irf_get_data_omni_modified(tint,'ByGSM','omni_min');
Bz_omni= irf_get_data_omni_modified(tint,'BzGSM','omni_min');
B_omni = [Bx_omni, By_omni(:,2), Bz_omni(:,2)];

B1_res = irf_resamp(B1, B_omni);
B1_res(:,2:4) = 0*B1_res(:,2:4); 
B1_res(:,4) = B1_res(:,4)+1;

B_dot = irf_dot(B_omni, B1_res);
% theta = acosd(B_dot(:,2)./vecnorm(B_omni(:,2:4),2,2)./vecnorm(B1_res(:,2:4),2,2)); %cone angle
theta = atan2d(B_omni(:,3), B_omni(:,4));
theta = [B_omni(:,1), theta];
%% Sound Speed
% Cs = sqrt(gamma*p/rho)
%% Init figure
n=8;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% B plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
set(gca,'Ylim',[-20 60], 'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% clock angle plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([theta(:,1) theta(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;

%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
set(gca,'Ylim',[-180 180], 'ytick',[-180 -90 0 90 180]);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('\theta [deg]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% ICI
% % h(i)=irf_subplot(n,1,-i);
% % c_eval("irf_plot([ICI?(:,1) ICI?(:,2)], 'color','k', 'Linewidth',0.75);",1); hold on;
% % c_eval("irf_plot([ICI?(:,1) ICI?(:,2)], 'color','b', 'Linewidth',0.75);",2); hold on;
% % c_eval("irf_plot([ICI?(:,1) ICI?(:,2)], 'color','g', 'Linewidth',0.75);",3); hold on;
% % c_eval("irf_plot([ICI?(:,1) ICI?(:,2)], 'color','r', 'Linewidth',0.75);",4); hold on;
% % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
% % % set(gca,'Ylim',[-30 80], 'fontsize',9);
% % pos1=get(gca,'pos');
% % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % ylabel('ICI','fontsize',10);
% % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % i=i+1;
%% Bx plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,2))])-0.1 max(B?(:,2))+0.1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Bx [nT]','fontsize',8);
% % % i=i+1;
%% By plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,3))])-0.1 max(B?(:,3))+0.1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('By [nT]','fontsize',8);
% % % i=i+1;
%% Bz plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,4))])-0.1 max(B?(:,4))+0.1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Bz [nT]','fontsize',8);
% % % i=i+1;
%% Btotal plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt2(:,1) Bt2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt3(:,1) Bt3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt4(:,1) Bt4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(Bt?(:,2))])-0.1 max(Bt?(:,2))+0.1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Btotal [nT]','fontsize',8);
% % % i=i+1;
%% N plot
h(i)=irf_subplot(n,1,-i);

%滤波
%     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
%     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;

%非滤波
c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2)) min(Ni?(:,2))])-0.02]) max([max(Ne?(:,2)) max(Ni?(:,2))])+0.02]);",ic)
    set(gca,'Ylim',[20 42],'fontsize',9);
% pos1=get(h(1),'pos');
%  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
%  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
 irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
ylabel('N [cm^{-3}]','fontsize',8);
i=i+1;
%% Ve plot
% % % h(i)=irf_subplot(n,1,-i);
% % % dspan = 4;
% % % c_eval('gsmVe? = [smooth(gsmVe?(:,1), dspan), smooth(gsmVe?(:,2), dspan), smooth(gsmVe?(:,3), dspan), smooth(gsmVe?(:,4), dspan)];',ic);
% % % c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vebf?(:,1) Vebf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vebf?(:,1) Vebf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([Vebf?(:,1) Vebf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("quiver(gsmVe?(:,1),0*gsmVe?(:,1),gsmVe?(:,2),gsmVe?(:,3));",ic);hold on;
% % % % irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % 
% % % grid off;
% % % ylabel('Ve [km/s]','fontsize',8);
% % % % c_eval("set(gca,'Ylim',[fix(min([min(gsmVe?(:,2)) min(gsmVe?(:,3)) min(gsmVe?(:,4))])/10)*10-10 fix(max(Vet?(:,2))/10)*10+10]);",ic);
% % % 
% % % set(gca,'Ylim',[-1500 2000]);
% % % % set(gca,'Ylim',[-1000 1000], 'ytick',[-600 -400 -200 0 200 400 600]);
% % % %irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Ve_N','Ve_M','Ve_L','|Ve|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.05 0.92]);
% % % i=i+1;
%% Electric field
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% % % % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('E [mV/m]','fontsize',8)
% % % i=i+1;
%% Vi plot
h(i)=irf_subplot(n,1,-i);
% c_eval("irf_plot([Vi?(:,1) Vi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Vi?(:,1) Vi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Vi?(:,1) Vi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;

irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi?(:,1) Vi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",ic);
set(gca,'Ylim',[0 250]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
ylabel('|Vi| [km/s]','fontsize',8);
i=i+1;
%% Te plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Te_para?(:,1) (Te_para?(:,2)+2*Te_perp?(:,2))/3], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te_para?(:,1) Te_para?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te_perp?(:,1) Te_perp?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[fix(min([min(Te_para?(:,2)) min(Te_para?(:,2)) min(Te_perp?(:,2))])/10)*10-10 fix(max([max(Te_para?(:,2)) max(Te_para?(:,2)) max(Te_perp?(:,2))])/10)*10+10],'fontsize',9);",ic);
set(gca,'Ylim',[20 40]);
% irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(gca,{'Te','T_/_/','T_⊥'},[0.97 0.92]);
ylabel('Te [eV]','fontsize',8);
i=i+1;
%% Ti plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Ti_para?(:,1) (Ti_para?(:,2)+2*Ti_perp?(:,2))/3], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Ti_para?(:,1) Ti_para?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Ti_perp?(:,1) Ti_perp?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold off;
grid off;
ylabel('Ti [eV]','fontsize',8);
% c_eval("set(gca,'Ylim',[fix(min([min(Ti_para?(:,2)) min(Ti_para?(:,2)) min(Ti_perp?(:,2))])/10)*10-10 fix(max([max(Ti_para?(:,2)) max(Ti_para?(:,2)) max(Ti_perp?(:,2))])/10)*10+10]);",ic);
set(gca,'Ylim',[0 5000]);
% irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(gca,{'Ti','Tipara','Tiperp'},[0.97 0.92]);
i=i+1;
%% Pressure
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Pb?(:,1) Pb?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Pt?(:,1) Pt?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Pt?(:,1) Pt?(:,2)+Pb?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vit?(:,1) Pd?], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Pthe_para(:,1) Pthe_para(:,2)], 'color','g', 'Linewidth',0.75);hold on;
% % % % irf_plot([Pthe_perp(:,1) Pthe_perp(:,2)], 'color','y', 'Linewidth',0.75);hold on;
% % % grid off;
% % % % set(h(i),'yscale','log');
% % % % set(h(i),'ytick',[0 0.25 0.5],'fontsize',9);
% % % % set(gca,'Ylim',[0 0.5]);
% % % c_eval("set(gca,'Ylim',[0 max(Pt?(:,2)+Pb?(:,2))+0.01]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0];[0 1 0]]);
% % % irf_legend(gca,{'Pm','Pthe','Ptotal'},[0.97 0.92]);
% % % % pos3=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('P [nPa]','fontsize',8)
% % % i=i+1; 
%% plot low e pad
% % % %     %0-200eV
% % % h(i)=irf_subplot(n,1,-i);
% % % % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % colormap(h(i),jet)
% % % specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
% % % specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
% % % specrec_p_elow.p=energy_low.data;%data matrix
% % % specrec_p_elow.f_label='';
% % % specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
% % % [h(i), hcb6]=irf_spectrogram(h(i),specrec_p_elow);
% % % ylabel('PA low','fontsize',8)
% % % % set(gca,'yscale','log');
% % % set(h(i),'ytick',[0 90 180]);
% % % % caxis(gca,[7 7.7]);
% % % %irf_legend(h(i),'g',[0.99 0.98],'color','w','fontsize',12);
% % % poscbar6=get(hcb6,'pos');
% % % poscbar6(3)=poscbar6(3)*0.5;
% % % set(hcb6,'pos',poscbar6);
% % % i=i+1;
%% plot mid e pad
% % % %     %200-2000eV
% % % h(i)=irf_subplot(n,1,-i);
% % % %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % colormap(h(i),jet)
% % % 
% % % specrec_p_emid=struct('t',irf_time(energy_mid.DEPEND_0.data,'ttns>epoch'));
% % % specrec_p_emid.f=transpose(energy_mid.DEPEND_1.data(1,1:30));%energy levels
% % % specrec_p_emid.p=energy_mid.data;%data matrix
% % % specrec_p_emid.f_label='';
% % % specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
% % % [h(i), hcb7]=irf_spectrogram(h(i),specrec_p_emid);
% % % ylabel('PA mid','fontsize',8)
% % % %set(gca,'yscale','log');
% % % set(h(i),'ytick',[0 90 180]);
% % % % clim(gca,[6 8]);
% % % %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
% % % poscbar7=get(hcb7,'pos');
% % % poscbar7(3)=poscbar7(3)*0.5;
% % % set(hcb7,'pos',poscbar7);
% % % i=i+1;
%% plot high e pad
%2k-30keV
% % % h(i)=irf_subplot(n,1,-i);
% % % % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % colormap(h(i),jet)
% % % 
% % % specrec_p_ehigh=struct('t',irf_time(energy_high.DEPEND_0.data,'ttns>epoch'));
% % % specrec_p_ehigh.f=transpose(energy_high.DEPEND_1.data(1,1:30));%energy levels
% % % specrec_p_ehigh.p=energy_high.data;%data matrix
% % % specrec_p_ehigh.f_label='';
% % % specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
% % % [h(i), hcb6]=irf_spectrogram(h(i),specrec_p_ehigh);
% % % ylabel('PA high','fontsize',8)
% % % 
% % % set(h(i),'ytick',[0 90 180]);
% % % % clim(gca,[6 7.6]);
% % % irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
% % % poscbar6=get(hcb6,'pos');
% % % poscbar6(3)=poscbar6(3)*0.5;
% % % set(hcb6,'pos',poscbar6);
% % % i=i+1;
%% plot e energy spectrom
h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)

specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
specrec_p_e.p=energy_e.data;%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ee(ev)','fontsize',8)
set(gca,'Ylim',[1e1 3e4]);
caxis(gca,[5 9])

% irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
poscbar8=get(hcb8,'pos');
poscbar8(3)=poscbar8(3)*0.5;
%poscbar6(1)=poscbar6(1)*0.5;
set(hcb8,'pos',poscbar8);
i=i+1;
%% plot ION energy spectrom
h(i)=irf_subplot(n,1,-i);

c_eval("specrec_p_i=struct('t',irf_time(energy_i?.DEPEND_0.data,'ttns>epoch'));",ic);
c_eval("specrec_p_i.f=transpose(energy_i?.DEPEND_1.data(1,1:32));",ic);%energy levels
c_eval("specrec_p_i.p=energy_i?.data;",ic);%data matrix
specrec_p_i.f_label='';
specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb7]=irf_spectrogram(h(i),specrec_p_i);
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ei(ev)','fontsize',8)
set(gca,'Ylim',[1e1 3e4]);
caxis(gca,[4.5 8])


% irf_legend(gca,'f',[0.99 0.98],'colo6.4r','k','fontsize',12);
poscbar7=get(hcb7,'pos');
poscbar7(3)=poscbar7(3)*0.5;
set(hcb7,'pos',poscbar7);
% % % load('/Users/fwd/Documents/MATLAB/Code/zwz/color_rwb.mat');
% % % cm = othercolor('PuBu8');
% % % % cm = othercolor('Greys9');
% % % cm = flip(cm);
% % % % cm(:,3) = linspace(0.15,0.851, 256);
% % % colormap(color_rwb)
i=i+1;
%% pressure & entropy & interchange instability 
% % % ic = 1:2;
% % % 
% % % c_eval('Bt? = irf_resamp(Bt?, Ni?);',ic)
% % % c_eval('R? = irf_resamp(R?, Ni?);',ic)
% % % c_eval('B? = irf_resamp(B?, Ni?);',ic)
% % % c_eval('Ni?=irf_resamp(Ni?,Bt?);',ic)
% % % c_eval('Ne?=irf_resamp(Ne?,Bt?);',ic)
% % % c_eval('Te?=irf_resamp(Te?,Bt?);',ic)
% % % c_eval('Ti?=irf_resamp(Ti?,Bt?);',ic)
% % % c_eval('Pb?=[Bt?(:,1) ((Bt?(:,2).^2))/(2*units.mu0)*1e-9];',ic);%nPa
% % % c_eval('Pti? = irf_multiply(11604.505*units.kB*1e6*1e9,[Ni?(:,1) Ni?(:,2)],1,[Ti?(:,1) Ti?(:,2)],1);',ic);
% % % c_eval('Pte? = irf_multiply(11604.505*units.kB*1e6*1e9,[Ne?(:,1) Ne?(:,2)],1,[Te?(:,1) Te?(:,2)],1);',ic);
% % % c_eval('Pt? = [Pti?(:,1) Pti?(:,2)+Pte?(:,2)];',ic);
% % % 
% % % c_eval('v?=10^0.7368.*sqrt(R?(:,2).^2+R?(:,3).^2).^0.7634.*B?(:,4).^(-0.3059)./sqrt(B?(:,4).^2+2.*400.*pi.*Pt?(:,2));',ic);
% % % c_eval('etr?=Pt?(:,2).*v?.^(5/3);',ic);
% % % 
% % % c_eval('beta? = [Bt?(:,1) Pt?(:,2)./Pb?(:,2)];',ic);
% % % c_eval('pe_pb? = [Bt?(:,1) Pte?(:,2)./Pb?(:,2)];',ic);
% % % 
% % % deltaR = irf_abs(R2(:,2:4) - R1(:,2:4));
% % % deltaR = 1./deltaR(:,4);
% % % c=(deltaR.*(v1-v2)./v2-beta2(:,2).*(deltaR.*(Pt1(:,2)-Pt2(:,2)))./(2.*Pt2(:,2)));
% % % d=deltaR.*(etr1-etr2);
% % % f=c.*d;
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % % c_eval("irf_plot([beta?(:,1) beta?(:,2)], 'color','k', 'Linewidth',0.75);",1); hold on;
% % % c_eval("irf_plot([Bt?(:,1) f], 'color','k', 'Linewidth',0.75); hold on;",1);
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(h(i),'yscale','log');
% % % % set(h(i),'ytick',[1 2 3 4],'fontsize',9);
% % % % set(gca,'Ylim',[-5 5]);
% % % % set(gca,'Ylim',[round(min(S1))-1 round(max(S1))]);
% % % set(gca,'ColorOrder',[0 0 1]);
% % % % irf_legend(gca,{'/beta'},[0.97 0.92]);
% % % % pos3=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % % ylabel('\beta','fontsize',12)
% % % ylabel('Criterion','fontsize',12)
% % % i=i+1; 
%% plot waves
% % % c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
% % % c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
% % % c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% % % % Bscm=Bscm{1};            %Bscm??cell
% % % c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
% % % magB = Bxyz.abs;
% % % 
% % % %gse2gsm
% % % c_eval(['Egse=irf.ts2mat(Exyz);'],ic);
% % % c_eval(['Egsm=irf_gse2gsm(Egse);'],ic);
% % % Exyz.data=Egsm(:,2:4);
% % % try
% % % c_eval(['Bscmgse=irf.ts2mat(Bscm);'],ic);
% % % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
% % % Bscm.data=Bscmgsm(:,2:4);
% % % 
% % % % Rotate E and B into field-aligned coordinates
% % % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
% % % Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);
% % % % Bandpass filter E and B waveforms
% % % dfE = 1/median(diff(Exyz.time.epochUnix));
% % % dfB = 1/median(diff(Bscm.time.epochUnix));
% % % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
% % % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
% % % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
% % % catch
% % % % % 当Bscm发生bug时其会变为{1,2}的cell，点进去发现两部分是一样的，有时候重启matlab会好使有时候不好使就用下面这部分（到wave transforms之前）
% % % c_eval(['Bscmgse=irf.ts2mat(Bscm{1,1});'],ic);
% % % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
% % % Bscm{1,1}.data=Bscmgsm(:,2:4);
% % % 
% % % % Rotate E and B into field-aligned coordinates
% % % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
% % % Bscmfac = irf_convert_fac(Bscm{1,1},Bxyz,[1 0 0]);
% % % % Bandpass filter E and B waveforms
% % % dfE = 1/median(diff(Exyz.time.epochUnix));
% % % dfB = 1/median(diff(Bscm{1,1}.time.epochUnix));
% % % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
% % % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
% % % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
% % % end
% % % 
% % % % Wavelet transforms
% % % nf = 100;
% % % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[1 4000]);
% % % % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 50000]);
% % % Bwavelet = irf_wavelet(Bscmfac,'nf',nf,'f',[1 4000]);
% % % 
% % % %compress wavelet transform data 10 point average
% % % nc = 20;
% % % idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
% % % Ewavelettimes = Ewavelet.t(idx);
% % % Ewaveletx = zeros(length(idx),nf);
% % % Ewavelety = zeros(length(idx),nf);
% % % Ewaveletz = zeros(length(idx),nf);
% % % for ii = [1:length(idx)];
% % %         Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % % end
% % % specperpE=struct('t',Ewavelettimes);
% % % specperpE.f=Ewavelet.f;
% % % specperpE.p=Ewaveletx+Ewavelety;
% % % specperpE.f_label='';
% % % specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % specparE=struct('t',Ewavelettimes);
% % % specparE.f=Ewavelet.f;
% % % specparE.p=Ewaveletz;
% % % specparE.f_label='';
% % % specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % specE=struct('t',Ewavelettimes);
% % % specE.f=Ewavelet.f;
% % % specE.p=Ewaveletx+Ewavelety+Ewaveletz;
% % % specE.f_label='';
% % % specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % 
% % % idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
% % % Bwavelettimes = Bwavelet.t(idx);
% % % Bwaveletx = zeros(length(idx),nf);
% % % Bwavelety = zeros(length(idx),nf);
% % % Bwaveletz = zeros(length(idx),nf);
% % % for ii = [1:length(idx)];
% % %         Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % % end
% % % specB=struct('t',Bwavelettimes);
% % % specB.f=Bwavelet.f;
% % % specB.p=Bwaveletx+Bwavelety+Bwaveletz;
% % % specB.f_label='';
% % % specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
% % % 
% % % 
% % % % Compute characteristic frequencies
% % % Units=irf_units; % read in standard units
% % % Me=Units.me;
% % % Mp=Units.mp;
% % % e=Units.e;
% % % epso=Units.eps0;
% % % mu0=Units.mu0;
% % % Mp_Me = Mp/Me;
% % % B_SI=magB.data*1e-9;
% % % Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
% % % Wce = e*B_SI/Me;
% % % Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
% % % Fce = Wce/2/pi;
% % % Fce01=Fce*0.1;
% % % Fce05=Fce*0.5;
% % % Fpe = Wpe/2/pi;
% % % Fcp = Fce/Mp_Me;
% % % Fpp = Wpp/2/pi;
% % % Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
% % % Fpe = irf.ts_scalar(magB.time,Fpe);
% % % Fce = irf.ts_scalar(magB.time,Fce);
% % % Flh = irf.ts_scalar(magB.time,Flh);
% % % Fpp = irf.ts_scalar(magB.time,Fpp);
% % % Fce01=irf.ts_scalar(magB.time,Fce01);
% % % Fce05=irf.ts_scalar(magB.time,Fce05);
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % colormap(h(i),jet)
% % % [a8,b8]=irf_spectrogram(h(i),specE,'log');
% % % 
% % % hold(h(i),'on');
% % % % irf_plot(h(i),Fpe,'color','k','LineWidth',1)
% % % irf_plot(h(i),Flh,'color','k','LineWidth',1)
% % % % % % irf_plot(h(i),Fce,'color','r','LineWidth',1)
% % % % % % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
% % % % % % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
% % % hold(h(i),'off');
% % % 
% % % % irf_legend(h(i),'(h)',[0.99 0.97],'color','w','fontsize',12)
% % % caxis(h(i),[-7 -1]);
% % % set(h(i),'yscale','log');
% % % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
% % % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
% % % ylabel(b8,{'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',10);
% % % grid(h(i),'off');
% % % poscbar8=get(b8,'pos');
% % % poscbar8(3)=poscbar8(3)*0.5;
% % % set(b8,'pos',poscbar8);
% % % i=i+1;
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % colormap(h(i),jet)
% % % [a9,b9]=irf_spectrogram(h(i),specB,'log');
% % % %[h(i),b9]=irf_spectrogram(h(i),specB,'log');
% % % 
% % % hold(h(i),'on');
% % % irf_plot(h(i),Flh,'color','k','LineWidth',1)
% % % % % % irf_plot(h(i),Fce,'color','r','LineWidth',1)
% % % % % % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
% % % % % % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
% % % hold(h(i),'off');
% % % 
% % % caxis(h(i),[-8 1]);
% % % set(h(i),'yscale','log');
% % % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
% % % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
% % % ylabel(b9,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',10);
% % % grid(h(i),'off');
% % % poscbar9=get(b9,'pos');
% % % poscbar9(3)=poscbar9(3)*0.5;
% % % set(b9,'pos',poscbar9);
% % % i=i+1;
% 
%   set(h(1:n),'fontsize',8);
% %   irf_zoom(tint,'x',h(1:n));



%   irf_pl_mark(h(1:n),[iso2epoch('2017-07-24T12:56:52.00Z')],'g');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.159Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.399Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.589Z')],'b');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.789Z')],'k'); 
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30Z')],'g');
%   irf_pl_marak(h(1:8),[iso2epoch('2015-10-16T13:04:30.209Z')],'k');
%  add_position(gca,gseR1), xlabel(gca,'')
%  irf_zoom(tintlmn,'x',h(4:7))

%%  出图保存部分
% c_eval("irf_pl_mark(h(?),B1(tempidx_B1,1),'k');",1:n)
irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
irf_plot_axis_align(h)
colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

% % % cd  /Users/fwd/Documents/Ti~mor~/M/DF&MP/1-Comparison&Implication/Figures/
% % % % rmdir(TempDir,'s'); 
% % % figname = 'mp';
% % %     print(gcf, '-dpdf', [figname '.pdf']);    