close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


% TT = '2017-06-19T03:54:16.000Z/2017-06-19T03:54:20.000Z';
% TT = '2017-07-06T01:44:55.000Z/2017-07-06T01:45:04.000Z';
% % % TT = '2017-07-06T17:31:58.000Z/2017-07-06T17:31:59.500Z';
% TT = '2019-07-29T16:06:15.000Z/2019-07-29T16:06:19.000Z';
% TT = '2017-07-12T11:54:33.600Z/2017-07-12T11:54:35.400Z';
% TT = '2017-05-25T04:40:29.500Z/2017-05-25T04:40:33.000Z';
% % % TT = '2017-05-28T03:55:56.300Z/2017-05-28T03:55:57.300Z';
% TT = '2017-06-25T04:41:51.500Z/2017-06-25T04:41:52.500Z'; %vey too slow
% TT = '2017-06-28T02:45:22.000Z/2017-06-28T02:45:25.000Z';
% TT = '2017-07-03T23:11:10.000Z/2017-07-03T23:11:25.000Z'; %peak not coincide
% TT = '2017-07-06T01:32:42.660Z/2017-07-06T01:32:43.300Z'; % case 8
% TT = '2018-07-06T12:37:04.500Z/2018-07-06T12:37:05.200Z'; % case 12
% TT = '2017-07-06T08:40:50.750Z/2017-07-06T08:40:51.800Z';
% TT = '2017-07-09T10:46:32.500Z/2017-07-09T10:46:34.500Z'; % case 10
% TT = '2017-07-17T09:43:57.700Z/2017-07-17T09:43:58.200Z'; % case 11
% TT = '2019-08-16T08:25:06.500Z/2019-08-16T08:25:09.000Z'; %case 13
% TT = '2021-09-05T00:43:44.700Z/2021-09-05T00:43:49.000Z'; % case 14
% TT = '2021-07-17T17:25:19.200Z/2021-07-17T17:25:21.500Z'; % case 15 , 1
% % % TT = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z'; % case 16, 1/3
 % % % TT = '2017-06-27T22:44:35.000Z/2017-06-27T22:44:40.000Z';
 % TT = '2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z';
 % TT = '2017-06-28T17:04:20.000Z/2017-06-28T17:04:45.000Z';
 % % % TT = '2022-08-19T03:52:22.000Z/2022-08-19T03:52:28.000Z';
 % TT = '2020-07-24T09:42:07.000Z/2020-07-24T09:42:15.000Z';
 TT = '2019-09-24T02:17:06.500Z/2019-09-24T02:17:10.000Z';
 


tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1;
% % % filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% % % filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
% % % filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% % % filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% % % % filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% % % % filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms');
% % % filenames = [filenames1, filenames2];
% % % % % 
% % % [filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% % % % % [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% % % % % [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);
% % % 
% % % SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 64, 'CheckSize', 1)
% % % % SDCFilesDownload(filenames,TempDir)
% % % % % 
% % % % % SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % % % SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval(['B1_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
c_eval(['Bt1_ts=B1_ts.abs;'],iic); 
c_eval(['B1=irf.ts2mat(B1_ts);'],iic);
 % c_eval(['B?=irf_gse2gsm(B?);'],ic);
c_eval(['Bt1=irf.ts2mat(Bt1_ts);'],iic);

% c_eval('B? = irf_resamp(B?, B1);',iic);
% c_eval('Bt? = irf_resamp(Bt?, Bt1);',iic);
% lvbo
% c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',iic);
% c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,3);',iic);
% c_eval(['Bbf?=irf.ts2mat(Bbf?);'],iic);

% c_eval('Bbff? = B?_ts.filt(0,0.8,dfB?,3);',ic);
% c_eval(['Bbff?=irf.ts2mat(Bbff?);'],ic);

% L=[-0.25,-0.50,0.83];%最大变化方向
% N=[0.34,-0.85,-0.40];%外法向
% M=cross(N,L);
% c_eval('B1=irf_newxyz(B1,L,M,N);',ic);

% load E
c_eval(['E1_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et1_ts=E1_ts.abs;'],ic); 
c_eval(['E1=irf_gse2gsm(E1_ts);'],ic);
c_eval(['E1=irf.ts2mat(E1_ts);'],ic);
c_eval(['Et1=irf.ts2mat(Et1_ts);'],ic);
c_eval(['E1_resamp=irf_resamp(E1,B1);'],ic);

c_eval(['Bt1_res=irf_resamp(Bt1,Et1);'],ic);

c_eval(['Efac1=irf_convert_fac(E1,B1,[1,0,0]);'],ic);

c_eval(['Vexb1=irf_cross(E1_resamp,B1);'],ic);
c_eval(['Vexb1(:,2:4)=1e3*Vexb1(:,2:4)./[Bt1(:,2).^2 Bt1(:,2).^2 Bt1(:,2).^2];'],ic);%km/s

% % % c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
% % % c_eval('Ebf? = E?_gsm.filt(0.8,1.1,dfE?,3);',ic);
% % % c_eval(['Ebf?=irf.ts2mat(Ebf?);'],ic);

% c_eval('Ebff? = E?_gsm.filt(0,0.8,dfE?,3);',ic);
% c_eval(['Ebff?=irf.ts2mat(Ebf?);'],ic);

%         c_eval('Elmn?=irf_newxyz(Ebf1,L,M,N);',ic);


% load FPI
c_eval('Ne1_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne1=irf.ts2mat(Ne1_ts);'],ic);
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
% c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
% % % c_eval('dfNe? = 1/median(diff(Ne?_ts.time.epochUnix));',ic);
% % % c_eval('Nebf? = Ne?_ts.filt(0,1,dfNe?,5);',ic);
% % % c_eval(['Nebf?=irf.ts2mat(Nebf?);'],ic);

% % % c_eval('dfNi? = 1/median(diff(Ni?_ts.time.epochUnix));',ic);
% % % c_eval('Nibf? = Ni?_ts.filt(0,1,dfNi?,5);',ic);
% % % c_eval(['Nibf?=irf.ts2mat(Nibf?);'],ic);

c_eval('Ve1_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet1_ts=Ve1_ts.abs;'],ic); 
c_eval(['Ve1=irf.ts2mat(Ve1_ts);'],ic);
c_eval(['gsmVe1_ts=irf_gse2gsm(Ve1_ts);'],ic);
c_eval(['gsmVe1=irf.ts2mat(gsmVe1_ts);'],ic);
c_eval(['Vet1=irf.ts2mat(Vet1_ts);'],ic);

% c_eval('dfVe? = 1/median(diff(gsmVe?_ts.time.epochUnix));',ic);
% c_eval('Vebf? = gsmVe?_ts.filt(0,1,dfVe?,5);',ic);
% c_eval(['Vebf?=irf.ts2mat(Vebf?);'],ic);

%% R
% units = irf_units;
% Pos = mms.t_data('R_gsm',tint);
% c_eval('R? = Pos.gsmR?;')
% c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
% c_eval('R? = irf_resamp(R?,Bt?);',1:iic)

%% Smooth
dspan = 5;
Ve1 = [smooth(gsmVe1(:,1),dspan),smooth(gsmVe1(:,2),dspan),...
    smooth(gsmVe1(:,3),dspan),smooth(gsmVe1(:,4),dspan)];

Vexb1 = irf_resamp(Vexb1, Ve1);
Vexb1 = [smooth(Vexb1(:,1),dspan),smooth(Vexb1(:,2),dspan),...
    smooth(Vexb1(:,3),dspan),smooth(Vexb1(:,4),dspan)];
% Ve1(:,2:4) = Ve1(:,2:4) - Vexb1(:,2:4);
% Ve1(:,2) = Ve1(:,2) + 75;

c_eval('gsmVe? = [smooth(gsmVe?(:,1),dspan),smooth(gsmVe?(:,2),dspan),smooth(gsmVe?(:,3),dspan),smooth(gsmVe?(:,4),dspan)];',1)
dd=1;
xx = Ve1(1:dd:end,1);
B1_smooth = [smooth(B1(:,1),dspan),smooth(B1(:,2),dspan),...
    smooth(B1(:,3),dspan),smooth(B1(:,4),dspan)];
dd_B = 10;
BB = B1(1:dd_B:end,1);
% xx=xx';
%%
% % h(1) = irf_subplot(3,1,-1);
% % quiver(xx,0*xx,Ve_smooth(1:dd:end,2),Ve_smooth(1:dd:end,3))
% % 
% % h(2) = irf_subplot(3,1,-2);
% % quiver(xx,0*xx,Ve_smooth(1:dd:end,2),Ve_smooth(1:dd:end,4))

h(1) = irf_subplot(1,1,-1);
c_eval("irf_plot([Ve1(:,1) Ve1(:,2)], 'color','b', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([Ve1(:,1) Ve1(:,3)], 'color','g', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([Ve1(:,1) Ve1(:,4)], 'color','r', 'Linewidth',0.75);",1); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Ve1(:,1) Ve1(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[fix(min([min(Ve?(:,2)) min(Ve?(:,3)) min(Ve?(:,4))])/10)*11 fix(max(Vet?(:,2))/10)*10.5],'fontsize',9);",ic);
% set(gca,'Ylim',[-250 250]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Ve_x','Ve_y','Ve_z','|Ve|'},[0.97 0.92]);
ylabel('Ve [km/s]','fontsize',8);
i=i+1;
%%
figure(2)
h(1)=irf_subplot(4,1,-1);
plot(B1(:,1), B1(:,2), 'color','b', 'Linewidth',0.75); hold on;
plot(B1(:,1), B1(:,3), 'color','g', 'Linewidth',0.75); hold on;
plot(B1(:,1), B1(:,4), 'color','r', 'Linewidth',0.75); hold on;
plot(B1(:,1), B1(:,2)*0, 'color','k', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[-10 20]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z'},[0.97 0.92]);
ylabel('B [nT]','fontsize',8);

h(2) = irf_subplot(4,1,-2);

q=quiver(BB,0*BB,B1_smooth(1:dd_B:end,2),B1_smooth(1:dd_B:end,3),0.5);
% ylim([-0.4 0.4])
q.ShowArrowHead = 'off';

h(3) = irf_subplot(4,1,-3); % Vex Vey
q=quiver(xx,0*xx,Ve1(1:dd:end,2),Ve1(1:dd:end,3),0.2);
% ylim([-15 5])
q.ShowArrowHead = 'off';

h(4) = irf_subplot(4,1,-4); % Vey Vez
quiver(xx,0*xx,Ve1(1:dd:end,3),Ve1(1:dd:end,4),0.5)

irf_zoom(tint,'x',h(1:4));
irf_plot_axis_align(h)
irf_timeaxis(h(1:4));
set(h(1:3),'XTickLabel',[]);  
%%
% figure(3)
% h(1)=irf_subplot(4,1,-1);
% plot(B1(:,1), B1(:,2), 'color','b', 'Linewidth',0.75); hold on;
% plot(B1(:,1), B1(:,3), 'color','g', 'Linewidth',0.75); hold on;
% plot(B1(:,1), B1(:,2)*0, 'color','k', 'Linewidth',0.75); hold on;
% set(gca,'Ylim',[-5 10]);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y'},[0.97 0.92]);
% ylabel('B [nT]','fontsize',8);
% 
% h(2) = irf_subplot(4,1,-2);
% irf_plot([gsmVe1(:,1), gsmVe1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe2(:,1), gsmVe2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe3(:,1), gsmVe3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), gsmVe4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), 0*gsmVe4(:,2)], 'k--', 'Linewidth',0.75); hold on;
% set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% ylabel('Vex [km/s]','fontsize',8);
% grid off
% set(h(2),'Ylim',[-300 2200]);
% 
% 
% h(3) = irf_subplot(4,1,-3);
% irf_plot([gsmVe1(:,1), gsmVe1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe2(:,1), gsmVe2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe3(:,1), gsmVe3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), gsmVe4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), 0*gsmVe4(:,2)], 'k--', 'Linewidth',0.75); hold on;
% set(h(3),'Ylim',[-1000 600]);
% grid off
% 
% h(4) = irf_subplot(4,1,-4);
% irf_plot([gsmVe1(:,1), gsmVe1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe2(:,1), gsmVe2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe3(:,1), gsmVe3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), gsmVe4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([gsmVe4(:,1), 0*gsmVe4(:,2)], 'k--', 'Linewidth',0.75); hold on;
% set(h(4),'Ylim',[-500 1500]);
% grid off
% 
% irf_zoom(tint,'x',h(1:4));
% irf_plot_axis_align(h)
%%
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')