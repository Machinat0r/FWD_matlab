close all
clear;clc

global ParentDir 
ParentDir = 'D:/MMS/'; 
DownloadDir = 'D:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


% TT = '2020-08-03T01:45:23.000Z/2020-08-03T01:45:38.000Z'; % case 16-short
TT = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z'; 
% TT = '2021-07-17T17:25:10.000Z/2021-07-17T17:25:30.000Z'; % case 15

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:3;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms');
filenames = [filenames1, filenames2, filenames3, filenames4];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% % % [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
[filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 64, 'CheckSize', 0)
% SDCFilesDownload(filenames,TempDir)
% % % 
% % % SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 32, 'CheckSize', 0)
SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 64, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
c_eval(['Bt?_ts=B?_ts.abs;'],iic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],iic);
 % c_eval(['B?=irf_gse2gsm(B?);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],iic);
% lvbo
c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',iic);
c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,3);',iic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],iic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);
c_eval('E?_err_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);
c_eval('E?_err=irf.ts2mat(E?_err_ts);',ic);

% load FPI
c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',iic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],iic);

% Te
c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);

% c_eval('dfTe_para? = 1/median(diff(Te_para?_ts.time.epochUnix));',ic);
% c_eval('Te_parabf? = Te_para?_ts.filt(0,1.5,dfE?,5);',ic);
% c_eval(['Te_parabf?=irf.ts2mat(Te_parabf?);'],ic);

c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',iic);
c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],iic);
c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',iic);
c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],iic);
c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],iic);

c_eval('Ve?_ts = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
c_eval(['gsmVe?_ts=irf_gse2gsm(Ve?_ts);'],ic);
c_eval(['gsmVe?=irf.ts2mat(gsmVe?_ts);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);

c_eval('dfVe? = 1/median(diff(gsmVe?_ts.time.epochUnix));',ic);
c_eval('Vebf? = gsmVe?_ts.filt(0,1,dfVe?,5);',ic);
% c_eval(['Vebf?=irf.ts2mat(Vebf?);'],ic);



c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',iic); 
% c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],iic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],iic);
c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],iic);
c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],iic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],iic);

%% R
units = irf_units;
Pos = mms.get_data('R_gsm',tint);
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
c_eval('R? = irf_resamp(R?,Bt?);',iic)
%% Maxwell Tensor
% f = 1/miu0 * (delta x B) x B
% [div_T,~] =c_4_grad('R?','B?','div_T');
% div_T(:,2:end) = div_T(:,2:end) * 1e-7; % nPa/m

[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
div_T = [divTshear(:,1) divTshear(:,2:end) - divPb(:,2:end)];
%% pressure
%Pressure
%Pt
c_eval('Ni? = irf_resamp(Ni?,B?);');
c_eval('Ti? = irf_resamp(Ti?,B?);');
c_eval("Pti? = [Ni?(:,1) 1.6*10^(-19)*Ni?(:,2).*Ti?(:,2)];");
c_eval('Rho? = [Ni?(:,1) 1e6 * Ni?(:,2) * units.mp];');
gradn = c_4_grad('R?','Rho?','grad');
gradp = c_4_grad('R?','Pti?','grad');

crosstemp = irf_cross(gradn, gradp);
Baroclinic = [crosstemp(:,1) crosstemp(:,2:4)./Rho1(:,2).^2];
Baroclinic(:,2:end) = Baroclinic(:,2:end); %s-2 (?)
Baroclinic = irf_abs(Baroclinic);

Rho1 = irf_resamp(Rho1, B1);
div_T(:,2:end) = div_T(:,2:end)./Rho1(:,2)*units.me/units.mp;
div_T = irf_abs(div_T);
%% E+VxB
c_eval('Eres? = irf_resamp(E?, B?);',2);
c_eval('Veres? = irf_resamp(gsmVe?, B?);',2);
c_eval('Vires? = irf_resamp(gsmVi?, B?);',2);

c_eval('Eres?(:,1) = zeros(size(Eres?(:,1)));',2);
c_eval('EVixB = Eres? + 1e-3*irf_cross(Vires?, B?);',2);
c_eval('EVexB = Eres? + 1e-3*irf_cross(Veres?, B?);',2);
%% J & JdotE
%J
 c_eval('Je?_ts = -units.e*Ne?_ts*gsmVe?_ts*1e3*1e6*1e9;',ic);
c_eval('Ji?_ts = units.e*Ne?_ts*gsmVi?_ts.resample(Ne?_ts.time)*1e3*1e6*1e9;',ic);
c_eval('J?_ts = (Je?_ts+Ji?_ts);',ic);
c_eval('J? = irf.ts2mat(J?_ts);',ic);

%JdotE
c_eval('E_resJ = irf_resamp(E?,J?);',ic);
c_eval('JdotEtemp = [E_resJ(:,1) irf_dot(J?(:,2:4),E_resJ(:,2:4))];',ic);
c_eval('JdotE? = irf_resamp(JdotEtemp,Ne?);',ic);

% J dot E'
[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
J_B(:,2:4) = 1e9*J_B(:,2:4);
c_eval('E_resJ? = irf_resamp(E?,J_B);',ic);
c_eval('JdotE_B? = [E_resJ?(:,1) irf_dot(J_B(:,2:4),E_resJ?(:,2:4))];',ic);

% c_eval('gsmVe? = irf_resamp(gsmVe?,B?);',ic)
% c_eval('Eplus? = [E_resJ?(:,1) E_resJ?(:,2:4) + 1e3*irf_cross(1e3*gsmVe?(:,2:4),1e-9*B?(:,2:4))];',ic)
% c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J_B(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('B?_resE = irf_resamp(B?,E?);',ic)
c_eval('gsmVe?_resE = irf_resamp(gsmVe?,E?);',ic)
c_eval('J_B = irf_resamp(J_B,E?);',ic)
c_eval('Eplus? = [E?(:,1) E?(:,2:4) + 1e3*irf_cross(1e3*gsmVe?_resE(:,2:4),1e-9*B?_resE(:,2:4))];',ic)
% c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J_B(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('J? = irf_resamp(J?,E?);',ic);
c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J?(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('JdotEint = [JdotEplus?(:,1) cumsum(JdotEplus?(:,2))];',ic);

%% div and PI

[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
J_B(:,2:4) = 1e9*J_B(:,2:4);
c_eval('E_resJ = irf_resamp(E?,J_B);',ic);
c_eval('JdotE_B = [E_resJ(:,1) irf_dot(J_B(:,2:4),E_resJ(:,2:4))];',ic);

%% Init figure
ic = 2;
n=7;
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
c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-3 max(Bt?(:,2))+3]);",ic);
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
set(gca,'Ylim',[-10 20])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% E+VixB
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B?(:,1) EVixB(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([B?(:,1) EVixB(:,3)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([B?(:,1) EVixB(:,4)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?_err(1:2:end,1) E?_err(1:2:end,2)], 'm--', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?_err(1:2:end,1) E?_err(2:2:end,2)], 'm--', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[-30 30]);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% c_eval("set(gca,'Ylim',[min([min(Evixb?(:,1)) min(Evixb?(:,2)) min(Evixb?(:,3))])-2 max([max(Evixb?(:,1)) max(Evixb?(:,2)) max(Evixb?(:,3))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'xtick',[])
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E+VixB [mV/m]','fontsize',12)
i=i+1;
%% E+VexB
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B?(:,1) EVexB(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([B?(:,1) EVexB(:,3)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([B?(:,1) EVexB(:,4)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?_err(1:2:end,1) E?_err(1:2:end,2)], 'm--', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?_err(1:2:end,1) E?_err(2:2:end,2)], 'm--', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% c_eval("set(gca,'Ylim',[min([min(Evexb?(:,1)) min(Evexb?(:,2)) min(Evexb?(:,3))])-2 max([max(Evexb?(:,1)) max(Evexb?(:,2)) max(Evexb?(:,3))])+2]);",ic);
set(gca,'Ylim',[-30 30]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VexB_x','E+VexB_y','E+VexB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'xtick',[])
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%% x-
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Baroclinic(:,1) Baroclinic(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([div_T(:,1) div_T(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-3 max(Bt?(:,2))+3]);",ic);
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
set(gca,'Ylim',[-1e3 1e3])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'a_{baroclinic}','a_{Lorentz}'},[0.97 0.92]);
ylabel('a [s^{-2}]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% y-
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Baroclinic(:,1) Baroclinic(:,3)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([div_T(:,1) div_T(:,3)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-3 max(Bt?(:,2))+3]);",ic);
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
set(gca,'Ylim',[-1e3 1e3])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'a_{baroclinic}','a_{Lorentz}'},[0.97 0.92]);
ylabel('a [s^{-2}]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% z-
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Baroclinic(:,1) Baroclinic(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([div_T(:,1) div_T(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-3 max(Bt?(:,2))+3]);",ic);
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
set(gca,'Ylim',[-1e3 1e3])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'a_{baroclinic}','a_{Lorentz}'},[0.97 0.92]);
ylabel('a [s^{-2}]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% total-
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Baroclinic(:,1) Baroclinic(:,5)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([div_T(:,1) div_T(:,5)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-3 max(Bt?(:,2))+3]);",ic);
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
set(gca,'Ylim',[1e-1 1e3])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
set(gca,'yscale','log')
irf_legend(gca,{'a_{baroclinic}','a_{Lorentz}'},[0.97 0.92]);
ylabel('a [s^{-2}]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%%
% c_eval("irf_pl_mark(h(?),B1(tempidx_B1,1),'k');",1:n)
irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)

%%  出图保存部分
colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

% % % cd  /Users/fwd/Documents/Ti~mor~/M/DF&MP/1-Comparison&Implication/Figures/
% % % % rmdir(TempDir,'s'); 
% % % figname = 'mp';
% % %     print(gcf, '-dpdf', [figname '.pdf']);    