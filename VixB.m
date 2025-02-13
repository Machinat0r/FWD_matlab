close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TT = '2021-08-15T03:35:15.00Z/2021-08-15T03:35:30.00Z';
% TT = '2021-08-22T06:39:30.00Z/2021-08-22T06:43:00.00Z';
% TT = '2018-02-06T13:29:00.00Z/2018-02-06T13:30:30.00Z';
% % % TT = '2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z';
% TT = '2015-11-04T04:34:00.00Z/2015-11-04T04:37:00.00Z';
% TT = '2018-07-03T15:50:10.00Z/2018-07-03T15:50:25.00Z';
% TT = '2018-08-19T18:24:30.00Z/2018-08-19T18:26:20.00Z';
% TT = '2019-08-16T01:03:33.00Z/2019-08-16T01:05:13.00Z';
% % TT='2017-08-23T15:38:30.00Z/2017-08-23T15:39:15.00Z';
% TT = '2021-07-10T12:41:23.00Z/2021-07-10T12:42:23.00Z';
% TT = '2018-07-03T15:50:00.00Z/2018-07-03T15:51:00.00Z';
% TT = '2017-08-20T02:01:30.00Z/2017-08-20T02:03:00.00Z';
% TT = '2017-08-07T16:01:00.00Z/2017-08-07T16:02:00.00Z';
% TT = '2021-07-21T12:46:20.00Z/2021-07-21T12:46:40.00Z';
% TT = '2017-05-05T20:06:30.00Z/2017-05-05T20:07:10.00Z';
% TT = '2022-08-18T23:53:00.00Z/2022-08-18T23:54:00.00Z';
% TT = '2022-08-19T01:13:40.00Z/2022-08-19T01:14:40.00Z';
% TT = '2020-08-02T16:56:10.00Z/2020-08-02T16:56:25.00Z';
% TT = '2017-06-25T05:06:58.00Z/2017-06-25T05:07:02.00Z';
TT = '2017-06-11T17:54:00.000Z/2017-06-11T17:57:00.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames = [filenames1,filenames2,filenames3,filenames4];

[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [filenames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
[filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
% SDCFilesDownload_NAS(filenames_fast,TempDir)
SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
c_eval(['Bt?_ts=B?_ts.abs;'],iic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],iic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],iic);
% lvbo
c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',iic);
c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,3);',iic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],iic);

% c_eval('Bbff? = B?_ts.filt(0,0.8,dfB?,3);',ic);
% c_eval(['Bbff?=irf.ts2mat(Bbff?);'],ic);

%         c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);

c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

c_eval(['Vexb?=irf_cross(E?,B?);'],ic);
c_eval(['Vexb?(:,2:4)=1e3*Vexb?(:,2:4)./[Bt?_res(:,2).^2 Bt?_res(:,2).^2 Bt?_res(:,2).^2];'],ic);%km/s


c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
c_eval('Ebf? = E?_gsm.filt(0.8,1.1,dfE?,3);',ic);
c_eval(['Ebf?=irf.ts2mat(Ebf?);'],ic);

% load Ve
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

% load Vi
c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
% c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);

%Vifac
c_eval(['Bt?_resVi=irf_resamp(Bt?,Vi?);'],ic);
c_eval(['Vifac?=irf_convert_fac(gsmVi?,B?,[1,0,0]);'],ic);

%VixB
c_eval('B?_res = irf_resamp(B?,gsmVi?);',ic);
c_eval('E?_res = irf_resamp(E?,gsmVi?);',ic);
c_eval('Evxb? = 1e3*cross(1000*gsmVi?(:,2:4),1e-9*B?_res(:,2:4));',ic); %mV/m
c_eval('Evxb_andE? = Evxb? + E?_res(:,2:4);',ic);

%VexB
c_eval('B?_resVe = irf_resamp(B?,gsmVe?);',ic);
c_eval('E?_resVe = irf_resamp(E?,gsmVe?);',ic);
c_eval('Evxb?_Ve = 1e3*cross(1000*gsmVe?(:,2:4),1e-9*B?_resVe(:,2:4));',ic); %mV/m
c_eval('Evxb_andE?_Ve = Evxb?_Ve + E?_resVe(:,2:4);',ic);
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
c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

%% Electric field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on;

c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',8)
i=i+1;
%% E res
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([E?_res(:,1) E?_res(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?_res(:,1) E?_res(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?_res(:,1) E?_res(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',8)
i=i+1;
%% Ve plot
% % % h(i)=irf_subplot(n,1,-i);
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
% % % c_eval("set(gca,'Ylim',[fix(min([min(gsmVe?(:,2)) min(gsmVe?(:,3)) min(gsmVe?(:,4))])/10)*10-10 fix(max(Vet?(:,2))/10)*10+10]);",ic);
% % % 
% % % % c_eval("set(gca,'Ylim',[-700 700]);",ic);
% % % % set(gca,'Ylim',[-1000 1000], 'ytick',[-600 -400 -200 0 200 400 600]);
% % % %irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Ve_N','Ve_M','Ve_L','|Ve|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.05 0.92]);
% % % i=i+1;


%% V_ExB plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([Vexb1(:,1) Vexb1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% % % grid off;
% % % ylabel('Vexb [km/s]','fontsize',10);
% % % set(gca,'Ylim',[fix(min([min(Vexb1(:,2)) min(Vexb1(:,3)) min(Vexb1(:,4))]))-1 fix(max(max(Vexb1(:,2:4))))+1]);
% % % % set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Vexb_x','Vexb_y','Vexb_z'},[0.05 0.92]);
% % % i=i+1;
%% Vi plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Vit?(:,1) Vit?(:,2)], 'color','w', 'Linewidth',0.5);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",ic);
%set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',8);
i=i+1;
%% VixB field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVi?(:,1) Evxb?(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVi?(:,1) Evxb?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVi?(:,1) Evxb?(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evxb?(:,1)) min(Evxb?(:,2)) min(Evxb?(:,3))])-2 max([max(Evxb?(:,1)) max(Evxb?(:,2)) max(Evxb?(:,3))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'VixB_x','VixB_y','VixB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%% E+VixB
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVi?(:,1) Evxb_andE?(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVi?(:,1) Evxb_andE?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVi?(:,1) Evxb_andE?(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evxb_andE?_Ve(:,1)) min(Evxb_andE?(:,2)) min(Evxb_andE?(:,3))])-2 max([max(Evxb_andE?(:,1)) max(Evxb_andE?(:,2)) max(Evxb_andE?(:,3))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%% VexB field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVe?(:,1) Evxb?_Ve(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVe?(:,1) Evxb?_Ve(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVe?(:,1) Evxb?_Ve(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evxb?_Ve(:,1)) min(Evxb?_Ve(:,2)) min(Evxb?_Ve(:,3))])-2 max([max(Evxb?_Ve(:,1)) max(Evxb?_Ve(:,2)) max(Evxb?_Ve(:,3))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'VexB_x','VexB_y','VexB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%% E+VexB
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVe?(:,1) Evxb_andE?_Ve(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVe?(:,1) Evxb_andE?_Ve(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([gsmVe?(:,1) Evxb_andE?_Ve(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evxb_andE?_Ve(:,1)) min(Evxb_andE?_Ve(:,2)) min(Evxb_andE?_Ve(:,3))])-2 max([max(Evxb_andE?_Ve(:,1)) max(Evxb_andE?_Ve(:,2)) max(Evxb_andE?_Ve(:,3))])+2]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VexB_x','E+VexB_y','E+VexB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%%
irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)

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
colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

% cd  C:\Matlab\bin\新建文件夹\fwd\
% rmdir(TempDir,'s'); 
% figname = 'overview';
%     print(gcf, '-dpdf', [figname '.pdf']);    