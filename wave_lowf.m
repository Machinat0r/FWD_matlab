clear;
clc;

mms.db_init('local_file_db','D:\MMS\');

ic=3;
% Tint=irf.tint('2018-08-22T15:34:32.00Z/2018-08-22T15:34:36.00Z');
% Tsta='2018-06-02T18:28:50Z';   
% Tend='2018-06-02T18:31:50Z';
% Tsta='2017-08-23T15:38:30.00Z';   
% Tend='2017-08-23T15:39:13.00Z';

% Tsta='2020-08-03T02:36:04.00Z';
% Tend='2020-08-03T02:36:30.00Z';   
% Tsta = '2017-08-20T02:02:00.00Z';
% Tend = '2017-08-20T02:03:00.00Z';
% Tsta = '2020-07-05T00:31:00.00Z';
% Tend = '2020-07-05T00:32:00.00Z';
Tsta='2019-08-05T16:24:00.00Z';   
Tend='2019-08-05T16:25:00.00Z';
% Tsta='2021-08-22T06:40:45.000Z';   
% Tend='2021-08-22T06:41:45.000Z';
% Tsta='2021-07-22T12:44:30.00Z';   
% Tend='2021-07-22T12:45:30.00Z';
Tint=irf.tint(Tsta,Tend);
% Tint=irf.tint('2019-07-22T17:09:45.00Z/2019-07-22T17:11:00.00Z');
%% Load Data 
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tint);',ic);
magB = Bxyz.abs;
B=irf.ts2mat(Bxyz);
Bt=irf.ts2mat(magB);
c_eval('Exyz_gse=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
% E_temp=irf.ts2mat(Exyz_gse);
Exyz=irf_gse2gsm(Exyz_gse);
E=irf.ts2mat(Exyz);
% Exyz = TSeries(Exyz_gse.time,[Exyz_gsm(:,2:4)]);
% Exyz = irf.ts_vec_xyz(Exyz_gse.time,Exyz_gsm(:,2:4));
c_eval('Bscm_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);

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
flag = 1;
try
    Bscm1=irf_gse2gsm(Bxyz);
%     Bscm1=irf_gse2gsm(Bscm_ts);
catch
    Bscm1=irf_gse2gsm(Bscm_ts{1,1});
    Bscm2=irf_gse2gsm(Bscm_ts{1,2});
    if length(Bscm1)~=length(Bscm2)
    flag = 2;
    end
end
% Bscm2=irf_gse2gsm(Bscm_ts{1,2});
% Bscm{1,2}=irf_gse2gsm(Bscm_cell{1,2});

% 
% Bscm_mat1=irf.ts2mat(Bscm_gse1);
% Bscm_mat2=irf.ts2mat(Bscm_gse2);
% Bscm_mat=[Bscm_mat1;Bscm_mat2];
%Bscm=irf_gse2gsm(Bscm_gse);
% Bscm = irf.ts_vec_xyz(irf_time(Bscm_mat(:,1),'epoch>epochtt'),Bscm_mat_gsm(:,2:4));
% Bscm=Bscm{1};            %Bscm是cell
% c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);

% L=-[-0.562 0.719 -0.410];
% M=-[-0.555 -0.695 -0.458];
% N=[0.614 0.030 -0.789];

% B_LMN=irf_newxyz(Bxyz,L,M,N);
% Bxyzmag = TSeries(Bxyz.time,[B_LMN.data magB.data]);
% E_LMN=irf_newxyz(Exyz,L,M,N);
% E_LMN=TSeries(Exyz_gse.time,E_LMN(:,2:4));
%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac1 = irf_convert_fac(Bscm1,Bxyz,[1 0 0]);
lf = 50;
if flag == 2
    Bscmfac2 = irf_convert_fac(Bscm2,Bxyz,[1 0 0]);
    dfB2 = 1/median(diff(Bscm2.time.epochUnix));
    Bscmfachf2 = Bscmfac2.filt(0,lf,dfB2,5);
%     Bscmfachf2 = Bscmfac2.filt(40,300,dfB2,5);
    Bfachf2=irf.ts2mat(Bscmfachf2);
end
%% Bandpass filter E and B waveforms
% dfE = 1/median(diff(Exyz.time.epochUnix));
dfB1 = 1/median(diff(Bscm1.time.epochUnix));
dfE = 128;


% 
Exyzfachf = Exyzfac.filt(0.1,lf,dfE,5);
% Exyzfaclf = Exyzfac.filt(0.1,lf,dfE,5);
Bscmfachf1 = Bscmfac1.filt(0.1,lf,dfB1,5);
% Bscmfachf1 = Bscmfac1.filt(0.02,0,128,5);
% Bscmfachf1 = Bscmfac1.filt(40,300,dfB1,5);
% % % 
% % % hf = 40;
% % % Exyzfachf = Exyzfac.filt(lf,hf,dfE,5);
% % % Bscmfachf1 = Bscmfac1.filt(lf,hf,dfB1,5);


Efachf=irf.ts2mat(Exyzfachf);
Bfachf1=irf.ts2mat(Bscmfachf1);

%% Wavelet transforms
nf = 100;
% Ewavelet = irf_wavelet(Exyzfachf,'nf',nf,'f',[lf 1]);
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[0.1 lf]);
Bwavelet1 = irf_wavelet(Bscmfachf1,'nf',nf,'f',[0.1 lf]);
% Bwavelet1 = irf_wavelet(Bscmfachf1,'nf',nf,'f',[lf 1]);
if flag == 2
Bwavelet2 = irf_wavelet(Bscmfachf2,'nf',nf,'f',[0.1 lf]);
end

%compress wavelet transform data 10 point average
nc = 20;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specperpE=struct('t',Ewavelettimes);
specperpE.f=Ewavelet.f;
specperpE.p=Ewaveletx+Ewavelety;
specperpE.f_label='';
specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};

specparE=struct('t',Ewavelettimes);
specparE.f=Ewavelet.f;
specparE.p=Ewaveletz;
specparE.f_label='';
specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};


idx1 = [nc/2:nc:length(Bwavelet1.t)-nc/2];
Bwavelettimes1 = Bwavelet1.t(idx1);
Bwaveletx1 = zeros(length(idx1),nf);
Bwavelety1 = zeros(length(idx1),nf);
Bwaveletz1 = zeros(length(idx1),nf);
for ii = [1:length(idx1)];
        Bwaveletx1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,1}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
        Bwavelety1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,2}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
        Bwaveletz1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,3}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
end
specB1=struct('t',Bwavelettimes1);
specB1.f=Bwavelet1.f;
specB1.p=Bwaveletx1+Bwavelety1+Bwaveletz1;
specB1.f_label='';
specB1.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};

if flag ==2
idx2 = [nc/2:nc:length(Bwavelet2.t)-nc/2];
Bwavelettimes2 = Bwavelet2.t(idx2);
Bwaveletx2 = zeros(length(idx2),nf);
Bwavelety2 = zeros(length(idx2),nf);
Bwaveletz2 = zeros(length(idx2),nf);
for ii = [1:length(idx2)];
        Bwaveletx2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,1}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
        Bwavelety2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,2}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
        Bwaveletz2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,3}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
end
specB2=struct('t',Bwavelettimes2);
specB2.f=Bwavelet2.f;
specB2.p=Bwaveletx2+Bwavelety2+Bwaveletz2;
specB2.f_label='';
specB2.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
end
%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fce01=Fce*0.1;
Fce05=Fce*0.5;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fcp = irf.ts_scalar(magB.time,Fcp);
Fpe = irf.ts_scalar(magB.time,Fpe);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fce01=irf.ts_scalar(magB.time,Fce01);
Fce05=irf.ts_scalar(magB.time,Fce05);
%%
%% Init figure 2
n_subplots=7;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B(:,1) B(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt(:,1) Bt(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('B [nT]','fontsize',12);
set(h(1),'Ylim',[fix(min([min(B(:,2)) min(B(:,3)) min(B(:,4))])/10)*10-10 fix(max(Bt(:,2))/10)*10+10]);
% set(gca,'Ylim',[-10 25]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.1 0.12]);
%% E plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([E(:,1) E(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Et(:,1) Et(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('E [mV/m]','fontsize',12);
% set(h(2),'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
set(gca,'Ylim',[-20 20]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.1 0.12]);

%% Efachf plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Efachf(:,1) Efachf(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaE [mV m^{-1}]');
% set(gca,'Ylim',[-10 10], 'ytick',[-9:3:9]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaE_{\perp 1}','\deltaE_{\perp 2}','\deltaE_{||}'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 0];[0 0 0];[0 0 0]]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);
%% Especperp
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(4), hcb]=irf_spectrogram(h(4),specperpE); hold on;
irf_plot(h(4),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Fpe,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(4),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(4),Fce05,'color','c','LineWidth',1.5);hold off;
grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e1 1e2 1e3 1e4]);
set(gca,'ylim',[0.1 10]);
caxis(h(4),[-7 -2]);
ylabel(h(4),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);
%% Especpar
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(5), hcb]=irf_spectrogram(h(5),specparE); hold on;
irf_plot(h(5),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Fpe,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(5),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(5),Fce05,'color','c','LineWidth',1.5);hold off;

grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e1 1e2 1e3 1e4]);
set(gca,'ylim',[0.1 10]);
caxis(h(5),[-7 -2]);
ylabel(h(5),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);
%% Bfachf plot 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)*0],'k--', 'Linewidth',0.75);hold on;

if flag ==2
irf_plot([Bfachf2(:,1) Bfachf2(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,2)*0],'k--', 'Linewidth',0.75);
end
hold off; grid off;
ylabel('\deltaB [nT]');
set(gca,'ylim',[min(min(Bfachf1(:,2:4))) max(max(Bfachf1(:,2:4)))]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{\perp 1}','\deltaB_{\perp 2}','\deltaB_{||}'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 0];[0 0 0];[0 0 0]]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);
%% Bspec 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(7), hcb]=irf_spectrogram(h(7),specB1); hold on;
if flag ==2
[h(7), hcb]=irf_spectrogram(h(7),specB2); hold on;
end
irf_plot(h(7),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(7),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(7),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(7),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(7),Fce05,'color','c','LineWidth',1.5);hold off;
grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e-1 1e0 1e1 1e2 1e3 1e4]);
set(gca,'ylim',[0.1 10]);
% caxis(h(7),[-3.5 1.5]);
ylabel(h(7),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);

%%
colormap(jet);

irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',Tint);

% irf_pl_mark(h(1:7),[iso2epoch('2015-10-16T13:04:26Z')],'k');
% irf_pl_mark(h(1:7),[iso2epoch('2015-12-14T00:59:04Z')],'k');

set(h(1:7),'fontsize',12);
set(gcf,'paperpositionmode','auto');
irf_adjust_panel_position;
figname = 'wave';
% figname=['C:\Users\fwd\Desktop\Ti~mor~\M\Formation of the rolling-pin distribution of suprathermal electrons behind dipolarization fronts\obs\wave'];
% print(gcf, '-dpdf', [figname '.pdf']);