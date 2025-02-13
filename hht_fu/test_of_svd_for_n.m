%------written by Wending Fu, Mar.2024 in Beijing------------
clear;clc;close;warning('off','MATLAB:singularMatrix')
%% Data Load
global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
TT = '2016-01-07T09:34:27.00Z/2016-01-07T09:34:37.00Z';
% % % TT = '2015-09-19T07:43:29.00Z/2015-09-19T07:43:32.00Z';
Tint = irf.tint(TT);
Tintlong = Tint+[-60 60];

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 3;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames = [filenames1,filenames2,filenames3,filenames4];

[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [filenames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

c_eval('B?=mms.get_data(''B_gse_fgm_brst_l2'',Tintlong,?);',ic);

c_eval('Bxyz_brst?=mms.get_data(''B_gsm_brst'',Tint,?);',ic);
% % % c_eval('Bxyz_srvy?=mms.get_data(''B_gsm_srvy'',Tint,?);',ic);
c_eval('Bscm?=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);
c_eval('Bscm? = irf_gse2gsm(Bscm?);',ic);

c_eval('Exyz_brst?=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',ic);
c_eval('Exyz_brst? = irf_gse2gsm(Exyz_brst?);',ic);
% c_eval('E?_gse_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
% c_eval('Exyz_brst?=Exyz_brst?.resample(Bxyz_brst?);',ic);

R  = mms.get_data('R_gsm',Tintlong);

c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?(:,1:3));',ic);
c_eval('magB=Bxyz_brst?.abs;',ic);
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);

%% Magnetic Field Background
c_eval('B0 = Bxyz_brst?.filt(0,1,128,5);',ic);
c_eval('E0 = Exyz_brst?.filt(0,1,128,5);',ic);
%% Test EM wave
% % % B0 = [1, 1.5, 5]; E0 = [0.5, 0.2, 0.3];
% % % B_amp = 2; E_amp = 3;
% % % 
% % % t = transpose(0:1/8192:1);
% % % 
% % % f = transpose(linspace(500,2000,size(t,1))); 
% % % % f = transpose(linspace(1000,1000,size(t,1))); 
% % % omega = 2 * pi * f; angle_deg = 5;
% % % angle_rad = deg2rad(angle_deg);
% % % 
% % % B_direction = null(B0);
% % % B_direction = transpose(B_direction(:,1));
% % % B_direction = B_direction / norm(B_direction);
% % % k = [sin(angle_rad), 0, cos(angle_rad)];
% % % E_direction = cross(B_direction, k);
% % % 
% % % E = bsxfun(@times, E_amp * cos(cumtrapz(t, omega)), E_direction) + E0;
% % % B = bsxfun(@times, B_amp * sin(cumtrapz(t, omega)), B_direction) + B0;
% % % 
% % % noise_std = 0.1;
% % % E_noise = noise_std * randn(size(E));
% % % B_noise = noise_std * randn(size(B));
% % % 
% % % E = E + E_noise;
% % % B = B + B_noise;
% % % 
% % % Bxyz_brst1 = [t,B0.*ones(size(t))];
% % % E0 = [t,E0.*ones(size(t))];
% % % Bscm1 = [t, B];
% % % Exyz_brst1 = [t,E];
% % % B0 = Bxyz_brst1;
% Rxyz1 = zeros(size(Bxyz_brst1));
%% polarization from svd-EM
fmin=2e2;
fmax=4e3;
% c_eval('fmax=0.5/(Bxyz_brst?.time(2).epochUnix-Bxyz_brst?.time(1).epochUnix);',ic);%奈奎斯特频率
c_eval('[k,theta,t_wavelet,freq_wavelet, B_denoise, E_denoise] = svd_EM(Bxyz_brst?, E0, Bscm?, Exyz_brst?, [fmin,fmax]);',ic);
k_norm = vecnorm(permute(k,[3,1,2]));
k_norm = reshape(k_norm, [size(k_norm,2), size(k_norm, 3)]);
freq_time_mat = repmat(freq_wavelet',[size(k,1),1]);
wave_velocity = 2*pi*freq_time_mat./k_norm;
theta(theta>90) = 180 - theta(theta>90);
%% wavelet
% nT^2/Hz
% c_eval('Bsum = irf_wavelet(Bscm?,''returnpower'',1,''cutedge'',0,''nf'',100,''f'',[fmin fmax],''fs'',8192);',ic);
c_eval('Bsum = irf_wavelet(B_denoise,''returnpower'',1,''cutedge'',0,''nf'',100,''f'',[fmin fmax],''fs'',8192);',ic);
Bwavelet = Bsum.p{1}+Bsum.p{2}+Bsum.p{3};
% (mV/m)^2/Hz
% c_eval('Esum = irf_wavelet(Exyz_brst?,''returnpower'',1,''cutedge'',0,''nf'',100,''f'',[fmin fmax],''fs'',8192);',ic);
c_eval('Esum = irf_wavelet(E_denoise,''returnpower'',1,''cutedge'',0,''nf'',100,''f'',[fmin fmax],''fs'',8192);',ic);
Ewavelet = Esum.p{1}+Esum.p{2}+Esum.p{3};
%% polarization from svd-old
% % % [~,planarity,thetak,elliptict]=irf_wavepolarize_magneticSVD(Bscm1,B0);
c_eval('polarization=irf_ebsp(Exyz_brst?,Bscm?,Bxyz_brst?,B0,Rxyz?,[fmin fmax],''polarization'',''fac'');',ic);
frequency = polarization.f;
time = polarization.t;
thetak = polarization.k_tp(:,:,1);
% % % Bsum = polarization.bb_xxyyzzss(:,:,4);
% % % Esum2D = polarization.ee_ss;
% % % Esum = polarization.ee_xxyyzzss(:,:,4);
% % % ellipticity = polarization.ellipticity;
% % % dop = polarization.dop;
% % % thetak = polarization.k_tp(:,:,1);
% % % planarity = polarization.planarity;
% % % t=polarization.pf_rtp(:,:,2);
% % % cost=cosd(t);
% % % v_ph=sqrt(Esum./Bsum)*1e6;
% % % Bsumthres = 1e-7;
% % % removepts = find(Bsum < Bsumthres);
% % % ellipticity(removepts) = NaN;
% % % thetak(removepts) = NaN;
% % % dop(removepts) = NaN;
% % % planarity(removepts) = NaN;
% % % v_ph(removepts) = NaN;
% % % cost(removepts) = NaN;
%% Compute characteristic frequencies
% % % Units=irf_units; % read in standard units
% % % Me=Units.me;
% % % Mp=Units.mp;
% % % e=Units.e;
% % % epso=Units.eps0;
% % % mu0=Units.mu0;
% % % Mp_Me = Mp/Me;
% % % B_SI=magB.data*1e-9;
% % % Wpe = sqrt(ne.resample(Bxyz_brst1).data*1e6*e^2/Me/epso);
% % % Wce = e*B_SI/Me;
% % % Wpp = sqrt(ne.resample(Bxyz_brst1).data*1e6*e^2/Mp/epso);
% % % Fce = Wce/2/pi;
% % % Fce01=Fce*0.1;
% % % Fce05=Fce*0.5;
% % % Fpe = Wpe/2/pi;
% % % Fcp = Fce/Mp_Me;
% % % Fpp = Wpp/2/pi;
% % % Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
% % % Fce = irf.ts_scalar(magB.time,Fce);
% % % Flh = irf.ts_scalar(magB.time,Flh);
% % % Fpp = irf.ts_scalar(magB.time,Fpp);
% % % Fce01=irf.ts_scalar(magB.time,Fce01);
% % % Fce05=irf.ts_scalar(magB.time,Fce05);
%% Init figure
n_subplots=4;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);
%% E wave
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
specrec=struct('t',Esum.t,'p_label',['mV^{2}m^{-2}Hz^{-1}']);
  specrec.f=Esum.f;
  specrec.p=Ewavelet;
  specrec.f_label='';
  specrec.p_label={'log_{10}E^{2}','mV^{2}m^{-2}Hz^{-1}'};
irf_spectrogram(h(1),specrec,'log','donotfitcolorbarlabel');
% % % irf_spectrogram(h(1),Esum,'log','donotfitcolorbarlabel');
hold(h(1),'on');
% % % irf_plot(h(1),Flh,'color','k','LineWidth',1.5)
% % % irf_plot(h(1),Fce,'color','r','LineWidth',1.5)
% % % irf_plot(h(1),Fce01,'color','w','LineWidth',1.5)
% % % irf_plot(h(1),Fce05,'color','c','LineWidth',1.5)
hold(h(1),'off');
  irf_legend(h(1),'(a)',[0.99 0.98],'color','w','fontsize',12)
set(h(1),'yscale','log');
% set(h(1),'ytick',[1e1 1e2 1e3 1e4]);
set(h(1),'ytick',[1e1 1e2 1e3]);
% caxis(h(5),[-6 -2.5]);
% clim(h(1),[-6 -1.3]);
ylabel(h(1),'f (Hz)','fontsize',12);
set(h(1),'XTickLabel',[])
%% B wave
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
specrec=struct('t',Bsum.t);
specrec.f=Bsum.f;
specrec.p=Bwavelet;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
irf_spectrogram(h(2),specrec,'log','donotfitcolorbarlabel');
% % % irf_spectrogram(h(2),Bsum,'log','donotfitcolorbarlabel');
hold(h(2),'on');
% % % irf_plot(h(2),Flh,'color','k','LineWidth',1.5)
% % % irf_plot(h(2),Fce,'color','r','LineWidth',1.5)
% % % irf_plot(h(2),Fce01,'color','w','LineWidth',1.5)
% % % irf_plot(h(2),Fce05,'color','c','LineWidth',1.5)
hold(h(2),'off');
irf_legend(h(2),'(b)',[0.99 0.98],'color','w','fontsize',12)
set(h(2),'yscale','log');
% set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
set(h(2),'ytick',[1e1 1e2 1e3]);
% clim(h(2),[-7.5 -3]);
ylabel(h(2),'f (Hz)','fontsize',12);
set(h(2),'XTickLabel',[])
%% theta from svd-old
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=thetak;
specrec.f_label='';
specrec.p_label={'\theta_{k}'};
irf_spectrogram(h(3),specrec,'lin','donotfitcolorbarlabel');
% % % irf_spectrogram(h(3),thetak,'lin','donotfitcolorbarlabel');
hold(h(3),'on');
% % % irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
% % % irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
% % % irf_plot(h(3),Fce01,'color','w','LineWidth',1.5)
% % % irf_plot(h(3),Fce05,'color','c','LineWidth',1.5)
hold(h(3),'off');
  irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
set(h(3),'yscale','log');
% set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
set(h(3),'ytick',[1e1 1e2 1e3]);
clim(h(3),[0, 90])
ylabel(h(3),'f (Hz)','fontsize',12);
set(h(3),'XTickLabel',[])
%% wave velocity
% % % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % specrec=struct('t',Bsum.t);
% % % specrec.f=freq_wavelet;
% % % specrec.p=wave_velocity/1e3;
% % % specrec.f_label='';
% % % specrec.p_label={'log V [km/s]'};
% % % irf_spectrogram(h(3),specrec,'log','donotfitcolorbarlabel');
% % % % % % irf_spectrogram(h(3),thetak,'lin','donotfitcolorbarlabel');
% % % hold(h(3),'on');
% % % % % % irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
% % % % % % irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
% % % % % % irf_plot(h(3),Fce01,'color','w','LineWidth',1.5)
% % % % % % irf_plot(h(3),Fce05,'color','c','LineWidth',1.5)
% % % hold(h(3),'off');
% % %   irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
% % % set(h(3),'yscale','log');
% % % % set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
% % % set(h(3),'ytick',[1e1 1e2 1e3]);
% % % % clim(h(3),[0, 90])
% % % ylabel(h(3),'f (Hz)','fontsize',12);
% % % set(h(3),'XTickLabel',[])
%% theta from svd-EM
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);
% % % n_mod = sum(n,3);
  specrec=struct('t',Bsum.t);
    specrec.f=freq_wavelet;
    specrec.p=theta;
    % specrec.p = n_mod;
    specrec.f_label='';
    specrec.p_label={'\theta_{k}'};
    irf_spectrogram(h(4),specrec,'lin','donotfitcolorbarlabel');
    hold(h(4),'on');
% % % irf_plot(h(4),Flh,'color','k','LineWidth',1.5)
% % % irf_plot(h(4),Fce,'color','r','LineWidth',1.5)
% % % irf_plot(h(4),Fce01,'color','w','LineWidth',1.5)
% % % irf_plot(h(4),Fce05,'color','c','LineWidth',1.5)
hold(h(4),'off');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)
set(h(4),'yscale','log');
% set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
set(h(4),'ytick',[1e1 1e2 1e3]);
clim(h(4),[0, 90])
ylabel(h(4),'f (Hz)','fontsize',12);
set(h(4),"XTickLabelRotation",0)
%%
% % % irf_zoom(h(1:n_subplots),'x',Tint);
colormap('jet')
% % % set(gcf,'render','painters');
% % % set(gcf,'paperpositionmode','auto')