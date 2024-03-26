clear;
clc;

mms.db_init('local_file_db','D:\MMS\');

ic=1;

Tsta='2016-09-27T01:18:21.000Z';   
Tend='2016-09-27T01:18:21.700Z';
Tint=irf.tint(Tsta,Tend);
% Tint=irf.tint('2019-07-22T17:09:45.00Z/2019-07-22T17:11:00.00Z');
%% Load Data 
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);
magB = Bxyz.abs;
B=irf.ts2mat(Bxyz);
Bt=irf.ts2mat(magB);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
% E_temp=irf.ts2mat(Exyz_gse);
% Exyz=irf_gse2gsm(Exyz_gse);
E=irf.ts2mat(Exyz);
% Exyz = TSeries(Exyz_gse.time,[Exyz_gsm(:,2:4)]);
% Exyz = irf.ts_vec_xyz(Exyz_gse.time,Exyz_gsm(:,2:4));
c_eval('Bscm_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);
% Bscm_ts = Bscm_ts{1,1};
% Bscm1 = irf_gse2gsm(Bscm_ts);
% Bscm1 = irf.ts2mat(Bscm_ts);

% for ii=1:length(Bscm_cell)
%     c_eval('Bscm?=Bscm_cell{?};',ii); 
%     
%     c_eval('time?=Bscm?.time;',ii); 
% end
% t_Bscm=time1;
% data_Bscm=Bscm1.data;
% for ii=2:length(Bscm_cell)
%     
%     c_eval('t_Bscm.epoch((length(t_Bscm.epoch(:,1))+1):(length(t_Bscm.epoch(:,1))+length(time?.epoch(:,1))),1)=time?.epoch;',ii); 
%     c_eval('data_Bscm((length(data_Bscm(:,1))+1):(length(data_Bscm(:,1))+length(Bscm?.data(:,1))),1:end)=Bscm?.data;',ii); 
% end
% Bscm_gse=irf.ts_vec_xyz(t_Bscm.epoch,data_Bscm);
% Bscm=irf_gse2gsm(Bscm_gse);
% % % Bscm1=irf_gse2gsm(Bscm_ts{1,1});
% % % Bscm2=irf_gse2gsm(Bscm_ts{1,2});
% Bscm{1,2}=irf_gse2gsm(Bscm_cell{1,2});

% 
% Bscm_mat1=irf.ts2mat(Bscm_gse1);
% Bscm_mat2=irf.ts2mat(Bscm_gse2);
% Bscm_mat=[Bscm_mat1;Bscm_mat2];
%Bscm=irf_gse2gsm(Bscm_gse);
% Bscm = irf.ts_vec_xyz(irf_time(Bscm_mat(:,1),'epoch>epochtt'),Bscm_mat_gsm(:,2:4));
% Bscm=Bscm{1};            %Bscm是cell
% c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);


% % L=[-0.676137 0 -0.736776];
% % M=[-0.367706 -0.874235 0.328266];
% % N=[-0.644115 0.485502 0.591103];
% % 
% % % B_LMN=irf_newxyz(Bxyz,L,M,N);
% % % Bxyzmag = TSeries(Bxyz.time,[B_LMN.data magB.data]);
% % Bscmfac1=irf_newxyz(Bscm_ts,L,M,N);
% % % Bxyzmag = TSeries(Bscm_ts.time,B_LMN.data);
% % Exyzfac=irf_newxyz(Exyz,L,M,N);
% % Exyzfac=TSeries(Exyz.time,Exyzfac.data);
%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac1 = irf_convert_fac(Bscm_ts,Bxyz,[0 1 0]);

% Bscmfac2 = irf_convert_fac(Bscm2,Bxyz,[1 0 0]);
%% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB1 = 1/median(diff(Bscm_ts.time.epochUnix));
% dfB2 = 1/median(diff(Bscm2.time.epochUnix));

lf=100;
Exyzfachf = Exyzfac.filt(lf,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,lf,dfE,5);
Bscmfachf1 = Bscmfac1.filt(lf,0,dfB1,5);
% Bscmfachf2 = Bscmfac2.filt(lf,0,dfB2,5);

Efachf=irf.ts2mat(Exyzfachf);
Bfachf1=irf.ts2mat(Bscmfachf1);
% Bfachf2=irf.ts2mat(Bscmfachf2);
% % % %% units
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
%% Init figure
n_subplots=3;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

Bfachf1 = Bfachf1(1721:1804,:);
%% Bfachf plot 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bfachf1(:,1) Bfachf1(:,3)], 'color','b', 'Linewidth',0.75); hold on;

% % % irf_plot([Bfachf1(:,1) Bfachf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf1(:,1) Bfachf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;

% % % irf_plot([Bfachf2(:,1) Bfachf2(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaB [nT]');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{\perp 1}','\deltaB_{\perp 2}','\deltaB_{||}'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 0];[0 0 0];[0 0 0]]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);
%% Bfachf plot 2
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([Bfachf1(:,1) Bfachf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Bfachf1(:,1) Bfachf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;

% % % irf_plot([Bfachf2(:,1) Bfachf2(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaB [nT]');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{\perp 1}','\deltaB_{\perp 2}','\deltaB_{||}'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 0];[0 0 0];[0 0 0]]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);
%% Bfachf plot 3
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % irf_plot([Bfachf1(:,1) Bfachf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf1(:,1) Bfachf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;

% % % irf_plot([Bfachf2(:,1) Bfachf2(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bfachf2(:,1) Bfachf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaB [nT]');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{\perp 1}','\deltaB_{\perp 2}','\deltaB_{||}'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 0];[0 0 0];[0 0 0]]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);