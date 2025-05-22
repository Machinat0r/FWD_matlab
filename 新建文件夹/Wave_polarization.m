clc;
clear;
close all;
global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
TT = '2016-01-07T09:34:27.00Z/2016-01-07T09:34:37.00Z';
Tint = irf.tint(TT);
Tintlong = Tint+[-60 60];
mms.db_init('local_file_db',ParentDir);

ic=1:4;
c_eval('B?=mms.get_data(''B_gse_fgm_brst_l2'',Tintlong,?);',ic);

ic=1;
c_eval('Bxyz_brst?=mms.get_data(''B_gsm_brst'',Tint,?);',ic);
c_eval('Bxyz_srvy?=mms.get_data(''B_gsm_srvy'',Tint,?);',ic);
c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);
c_eval('Bscm = irf_gse2gsm(Bscm);',ic);

c_eval('Exyz_brst?=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',ic);
% c_eval('E?_gse_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
% c_eval('Exyz_brst?=Exyz_brst?.resample(Bxyz_brst?);',ic);

R  = mms.get_data('R_gsm',Tintlong);

c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?(:,1:3));',ic);
% c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?);',ic);
c_eval('magB=Bxyz_brst?.abs;',ic);
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);
%% calculate B0

B0 = Bxyz_brst1.filt(0,1,128,5);
%% polarization
fmin=2e2;
fmax=4e3;
% c_eval('fmax=0.5/(Bxyz_brst?.time(2).epochUnix-Bxyz_brst?.time(1).epochUnix);',ic);%奈奎斯特频率
c_eval('polarization=irf_ebsp(Exyz_brst?,Bscm,Bxyz_brst?,B0,Rxyz1,[fmin fmax],''polarization'',''fac'');',ic);
% irf_pl_ebsp(polarization);

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(Bxyz_brst1).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(Bxyz_brst1).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fce01=Fce*0.1;
Fce05=Fce*0.5;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fce01=irf.ts_scalar(magB.time,Fce01);
Fce05=irf.ts_scalar(magB.time,Fce05);

% save polarization;
% load polarization;
%% plot
frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Esum2D = polarization.ee_ss;
Esum = polarization.ee_xxyyzzss(:,:,4);
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
t=polarization.pf_rtp(:,:,2);
cost=cosd(t);
v_ph=sqrt(Esum./Bsum)*1e6;


Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
v_ph(removepts) = NaN;
cost(removepts) = NaN;

%% Init figure 2
n_subplots=10;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% 磁场Bx
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        irf_plot(h(1),B1.x,'color','k');hold(h(1),'on');
        irf_plot(h(1),B2.x,'color','r');hold(h(1),'on');
        irf_plot(h(1),B3.x,'color','g');hold(h(1),'on');
        irf_plot(h(1),B4.x,'color','b');hold(h(1),'on');
        hold(h(1),'on');
        hold(h(1),'off');
        grid(h(1),'off');
        set(h(1), 'Ylim',[-3.2 -0.3]); 
        ylabel(h(1),{'Bx (nT)'},'fontsize',9,'Interpreter','tex');
        set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
        irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.02 9.9]);
%         irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
        set(h(1),'FontSize',10);  

%% 磁场By
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        irf_plot(h(2),B1.y,'color','k');hold(h(2),'on');
        irf_plot(h(2),B2.y,'color','r');hold(h(2),'on');
        irf_plot(h(2),B3.y,'color','g');hold(h(2),'on');
        irf_plot(h(2),B4.y,'color','b');hold(h(2),'on');
        hold(h(2),'on');
        hold(h(2),'off');
        grid(h(2),'off');
        set(h(2), 'Ylim',[-2.3 1.2]); 
        ylabel(h(2),{'By (nT)'},'fontsize',9,'Interpreter','tex');
%         irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.02 7]);
%         irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
        set(h(2),'FontSize',10);  
%% 磁场Bz
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        irf_plot(h(3),B1.z,'color','k');hold(h(3),'on');
        irf_plot(h(3),B2.z,'color','r');hold(h(3),'on');
        irf_plot(h(3),B3.z,'color','g');hold(h(3),'on');
        irf_plot(h(3),B4.z,'color','b');hold(h(3),'on');
        hold(h(3),'on');
        hold(h(3),'off');
        grid(h(3),'off');
        set(h(3), 'Ylim',[3 8.2]); 
        ylabel(h(3),{'Bz (nT)'},'fontsize',9,'Interpreter','tex');
%         irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.02 7]);
%         irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
        set(h(3),'FontSize',10); 
%% 磁场|B|
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        irf_plot(h(4),B1.abs,'color','k');hold(h(4),'on');
        irf_plot(h(4),B2.abs,'color','r');hold(h(4),'on');
        irf_plot(h(4),B3.abs,'color','g');hold(h(4),'on');
        irf_plot(h(4),B4.abs,'color','b');hold(h(4),'on');
        hold(h(4),'on');
        hold(h(4),'off');
        grid(h(4),'off');
        set(h(4), 'Ylim',[3.3 8]); 
        ylabel(h(4),{'|B| (nT)'},'fontsize',9,'Interpreter','tex');
        set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%         irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.02 7]);
%         irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
        set(h(4),'FontSize',10); 
        
%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time,'p_label',['mV^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Esum;
    specrec.f_label='';
    specrec.p_label={'log_{10}E^{2}','mV^{2}m^{-2}Hz^{-1}'};
    irf_spectrogram(h(5),specrec,'log','donotfitcolorbarlabel');
hold(h(5),'on');
irf_plot(h(5),Flh,'color','k','LineWidth',1.5)
irf_plot(h(5),Fce,'color','r','LineWidth',1.5)
irf_plot(h(5),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(5),Fce05,'color','c','LineWidth',1.5)
hold(h(5),'off');
  irf_legend(h(5),'(a)',[0.99 0.98],'color','w','fontsize',12)
set(h(5),'yscale','log');
% set(h(1),'ytick',[1e1 1e2 1e3 1e4]);
set(h(5),'ytick',[1e1 1e2 1e3]);
% caxis(h(5),[-6 -2.5]);
caxis(h(5),[-6 -1.3]);
ylabel(h(5),'f (Hz)','fontsize',12);

%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=Bsum;
    specrec.f_label='';
    specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
    irf_spectrogram(h(6),specrec,'log','donotfitcolorbarlabel');
    hold(h(6),'on');
irf_plot(h(6),Flh,'color','k','LineWidth',1.5)
irf_plot(h(6),Fce,'color','r','LineWidth',1.5)
irf_plot(h(6),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(6),Fce05,'color','c','LineWidth',1.5)
hold(h(6),'off');
  irf_legend(h(6),'(b)',[0.99 0.98],'color','w','fontsize',12)
set(h(6),'yscale','log');
% set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
set(h(6),'ytick',[1e1 1e2 1e3]);
caxis(h(6),[-7.5 -3]);
ylabel(h(6),'f (Hz)','fontsize',12);

% h(3)=irf_panel('vph');
%   specrec=struct('t',time);
%     specrec.f=frequency;
%     specrec.p=v_ph;
%     specrec.f_label='';
%     specrec.p_label={'log_{10}|E|/|B|','m s^{-1}'};
%     irf_spectrogram(h(3),specrec,'log','donotfitcolorbarlabel');
%     hold(h(3),'on');
% irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
% irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
% irf_plot(h(3),Fce01,'color','w','LineWidth',1.5)
% irf_plot(h(3),Fce05,'color','c','LineWidth',1.5)
% hold(h(3),'off');
%   irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
% set(h(3),'yscale','log');
% set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
% caxis(h(3),[5, 9])
% ylabel(h(3),'f (Hz)','fontsize',12);

%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=ellipticity;
    specrec.f_label='';
    specrec.p_label={'Ellipticity'};
    irf_spectrogram(h(7),specrec,'lin','donotfitcolorbarlabel');
    hold(h(7),'on');
irf_plot(h(7),Flh,'color','k','LineWidth',1.5)
irf_plot(h(7),Fce,'color','r','LineWidth',1.5)
irf_plot(h(7),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(7),Fce05,'color','c','LineWidth',1.5)
hold(h(7),'off');
  irf_legend(h(7),'(c)',[0.99 0.98],'color','k','fontsize',12)
set(h(7),'yscale','log');
% set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
set(h(7),'ytick',[1e1 1e2 1e3]);
% caxis(h(7),[-1, 1])
ylabel(h(7),'f (Hz)','fontsize',12);

%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=thetak;
    specrec.f_label='';
    specrec.p_label={'\theta_{k}'};
    irf_spectrogram(h(8),specrec,'lin','donotfitcolorbarlabel');
    hold(h(8),'on');
irf_plot(h(8),Flh,'color','k','LineWidth',1.5)
irf_plot(h(8),Fce,'color','r','LineWidth',1.5)
irf_plot(h(8),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(8),Fce05,'color','c','LineWidth',1.5)
hold(h(8),'off');
  irf_legend(h(8),'(d)',[0.99 0.98],'color','k','fontsize',12)
set(h(8),'yscale','log');
% set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
set(h(8),'ytick',[1e1 1e2 1e3]);
caxis(h(8),[0, 90])
ylabel(h(8),'f (Hz)','fontsize',12);

% h(5)=irf_panel('dop');
%   specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
%     specrec.f=frequency;
%     specrec.p=dop;
%     specrec.f_label='';
%     specrec.p_label={'DOP'};
%     irf_spectrogram(h(5),specrec,'lin','donotfitcolorbarlabel');
%     hold(h(5),'on');
% irf_plot(h(5),Flh,'color','k','LineWidth',1.5)
% irf_plot(h(5),Fce,'color','r','LineWidth',1.5)
% irf_plot(h(5),Fce01,'color','w','LineWidth',1.5)
% irf_plot(h(5),Fce05,'color','c','LineWidth',1.5)
% hold(h(5),'off');
%   irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)
% set(h(5),'yscale','log');
% set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
% caxis(h(5),[0, 1])
% ylabel(h(5),'f (Hz)','fontsize',12);

%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=planarity;
    specrec.f_label='';
    specrec.p_label={'planarity'};
    irf_spectrogram(h(9),specrec,'lin','donotfitcolorbarlabel');
    hold(h(9),'on');
irf_plot(h(9),Flh,'color','k','LineWidth',1.5)
irf_plot(h(9),Fce,'color','r','LineWidth',1.5)
irf_plot(h(9),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(9),Fce05,'color','c','LineWidth',1.5)
hold(h(9),'off');
  irf_legend(h(9),'(e)',[0.99 0.98],'color','k','fontsize',12)
set(h(9),'yscale','log');
% set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
set(h(9),'ytick',[1e1 1e2 1e3]);
caxis(h(9),[0, 1])
ylabel(h(9),'f (Hz)','fontsize',12);


%% 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  specrec=struct('t',time);
    specrec.f=frequency;
    specrec.p=cost;
    specrec.f_label='';
    specrec.p_label={'Par PFlux'};
    irf_spectrogram(h(10),specrec,'lin','donotfitcolorbarlabel');
    hold(h(10),'on');
irf_plot(h(10),Flh,'color','k','LineWidth',1.5)
irf_plot(h(10),Fce,'color','r','LineWidth',1.5)
irf_plot(h(10),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(10),Fce05,'color','c','LineWidth',1.5)
hold(h(10),'off');
  irf_legend(h(10),'(f)',[0.99 0.98],'color','w','fontsize',12)
set(h(10),'yscale','log');
% set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
set(h(10),'ytick',[1e1 1e2 1e3]);
 caxis(h(10),[-1,1])
ylabel(h(10),'f (Hz)','fontsize',12);

% set(h(5:10),'ytick', [0.5 1 2 5 10]);
set(h(5:10),'ytick', [5 10 20 50 100 200]);
% irf_adjust_panel_position
irf_plot_axis_align(h);
% irf_pl_number_subplots(h,[0.02, 0.98]);
irf_zoom(h,'x',Tint);
set(h(1:10),'xgrid','off','ygrid','off')

% n=256;
%  rgb = flipud(seismicColorMap(n));
%  colormap(h(i_subplot-1),rgb);
%%
colormap(h(5),jet);
colormap(h(6),jet);
colormap(h(7),jet);
colormap(h(8),jet);
colormap(h(9),jet);
colormap(h(10),jet);
% colormap(h(7),jet);
% colormap(h(8),jet);

Tsta='2021-07-22T05:20:55Z';   
Tend='2021-07-22T05:20:58Z';
%%
tint2=irf.tint('2021-07-22T05:20:56Z/2021-07-22T05:20:57');
for ii = 1:6
    c_eval('irf_pl_mark(h(?),tint2,[0.74 0.74 0.74],''edgecolor'',''k'',''linestyle'',''--'')',ii);
end
%%
% color=[192 192 192]./255;
% set(h(1:6),'color',color);
% set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
% irf_adjust_panel_position;
figname=['wave_polarization' '_' Tsta(1:4) Tsta(6:7) Tsta(9:10) '-' Tsta(12:13) ...
    Tsta(15:16) Tsta(18:19) '_' Tend(12:13) Tend(15:16) Tend(18:19)];
% print(gcf, '-dpng', [figname '.png']);