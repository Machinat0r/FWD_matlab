clear
clc

mms.db_init('local_file_db','D:\Works\mms_db\data');
%tint = irf.tint('2018-09-11T08:05:38.00Z/2018-09-11T08:05:42.00Z');
tint=irf.tint('2019-03-29T13:22:50.00Z/2019-03-29T13:23:30.00Z');
%tint = irf.tint('2018-09-08T13:23:30.00Z/2018-09-08T13:24:30.00Z');
ic=1;

%% Load Data 
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% Bscm=Bscm{1};            %Bscm??cell
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);


% Bscm=mms.db_get_variable('mms1_scm_brst_l2_scb','mms1_scm_acb_gse_scb_brst_l2',tint);


Bscm_t=Bscm.time;
Bscm_data=Bscm.data;
% Bscm_data=[Bscm{1, 1}.data;Bscm{1, 2}.data];
Bscm=irf.ts_vec_xyz(Bscm_t, Bscm_data);




% % merge data/time 
% c_eval('Bscm_1 = dataobj(''D:\MATLAB\xy-matlab\MMS\mms_db\data\mms1\scm\brst\l2\scb\2017\07\18\mms1_scm_brst_l2_scb_20170718130353_v2.2.0.cdf'');',ic);
% c_eval('Bscm1 = get_variable(Bscm_1,''mms?_scm_acb_gse_scb_brst_l2'');',ic);
% 
% c_eval('Bscm_2 = dataobj(''D:\MATLAB\xy-matlab\MMS\mms_db\data\mms1\scm\brst\l2\scb\2017\07\18\mms1_scm_brst_l2_scb_20170718130503_v2.2.0.cdf'');',ic);
% c_eval('Bscm2 = get_variable(Bscm_2,''mms?_scm_acb_gse_scb_brst_l2'');',ic);
% % data merge
% data1=Bscm1.data; data2=Bscm2.data; data=[data1;data2];Bscm=Bscm1; Bscm.data=data; Bscm.nrec=Bscm1.nrec+Bscm2.nrec;
% 
% % time merge
% data1=Bscm1.DEPEND_0.data;data2=Bscm2.DEPEND_0.data; data=[data1;data2]; Bscm.DEPEND_0.data=data;






magB = Bxyz.abs;

%gse2gsm
c_eval(['Egse=irf.ts2mat(Exyz);'],ic);
c_eval(['Egsm=irf_gse2gsm(Egse);'],ic);
Exyz.data=Egsm(:,2:4);

c_eval(['Bscmgse=irf.ts2mat(Bscm);'],ic);
c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
Bscm.data=Bscmgsm(:,2:4);


%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Efacmat=irf.ts2mat(Exyzfac)
Eperp=sqrt(Efacmat(:,2).^2+Efacmat(:,3).^2);

Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);
%% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm.time.epochUnix));
% Exyzfachf = Exyzfac.filt(4000,1,dfE,5);
% Exyzfaclf = Exyzfac.filt(1,4000,dfE,5);
% Bscmfachf = Bscmfac.filt(4000,1,dfB,5);

%% Wavelet transforms
nf = 100;
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[1 4000]);
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[1 4000]);
Bwavelet = irf_wavelet(Bscmfac,'nf',nf,'f',[1 4000]);

%compress wavelet transform data 10 point average
nc = 20;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
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

specE=struct('t',Ewavelettimes);
specE.f=Ewavelet.f;
specE.p=Ewaveletx+Ewavelety+Ewaveletz;
specE.f_label='';
specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

nc = 20;
idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
Bwavelettimes = Bwavelet.t(idx);
Bwaveletx = zeros(length(idx),nf);
Bwavelety = zeros(length(idx),nf);
Bwaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
        Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specB=struct('t',Bwavelettimes);
specB.f=Bwavelet.f;
specB.p=Bwaveletx+Bwavelety+Bwaveletz;
specB.f_label='';
specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};


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
Fpe = irf.ts_scalar(magB.time,Fpe);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fce01=irf.ts_scalar(magB.time,Fce01);
Fce05=irf.ts_scalar(magB.time,Fce05);

%% Init figure
n_subplots=3;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 5; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%%
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([Efacmat(:,1) Efacmat(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Efacmat(:,1) Eperp], 'color','b', 'Linewidth',0.75); hold off;
% set(gca,'ColorOrder',[[0 1 0];[0 0 1]]);
% irf_legend(gca,{'E||     ','Eกอ'},[0.05 0.95],'fontsize',5);
% ylabel({'E';'[mV/m]'},'fontsize',5);



% h=irf_plot(5,'newfigure');
% set(gcf,'position',[100 100 650 800]);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
mmsColors=[0 0 1; 0 1 0; 1 0 0; 0 0 0];
set(h(1),'ColorOrder',mmsColors);
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
hold(h(1),'on');
irf_plot(h(1),B1, 'Linewidth',0.75);
irf_plot(h(1),magB,'color','k', 'Linewidth',0.75);
irf_plot(h(1),[B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.5);
hold(h(1),'off');
ylabel(h(1),{'B','(nT)'},'Interpreter','tex','fontsize',12);
set(h(1),'ytick',[5 10 15]);
% irf_zoom(h(1),'y',[-60 60]);
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}','|B|'},[0.1 0.92],'fontsize',15)
irf_legend(h(1),'a)',[0.01 0.95],'color','k','fontsize',12)
grid(h(1),'off');

% h(2)=irf_panel('ELMN');
% mmsColors=[1 0 0; 0 0 1 ; 0 1 0];
% set(h(2),'ColorOrder',mmsColors)
% hold(h(2),'on');
% irf_plot(h(2),E_LMN);
% hold(h(2),'off');
% ylabel(h(2),{'E (mV m^{-1})'},'Interpreter','tex','fontsize',12);
% irf_legend(h(2),{'E_{L}','E_{M}','E_{N}'},[0.98 0.12],'fontsize',15)
% irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',12)
% set(h(2),'ylim',[-20 20],'ytick',[-10:10:10])
% grid(h(2),'off');


% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% mmsColors=[0 0 1; 0 1 0; 1 0 0];
% set(h(1),'ColorOrder',mmsColors)
% hold(h(1),'on');
% irf_plot(h(1),Exyzfac);
% hold(h(1),'off');
% ylabel(h(1),{'E','(mV/m)'},'Interpreter','tex','fontsize',8);
% irf_legend(gca,{'E_{\perp1}        ','E_{\perp2}        ','E_{||}'},[0.1 0.92],'fontsize',8);

% h(3)=irf_panel('Especperp');
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
[a2,b2]=irf_spectrogram(h(2),specE,'log');
hold(h(2),'on');
% irf_plot(h(2),Fpe,'color','k','LineWidth',1.5)
irf_plot(h(2),Flh,'color','k','LineWidth',1.5)
irf_plot(h(2),Fce,'color','r','LineWidth',1.5)
irf_plot(h(2),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(2),Fce05,'color','c','LineWidth',1.5)
hold(h(2),'off');
% irf_legend(h(2),'(h)',[0.99 0.97],'color','w','fontsize',12)
% caxis(h(2),[-6 0]);
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(2),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
ylabel(b2,{'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',10);
grid(h(2),'off');
colormap(h(2),jet); 
% % h(4)=irf_panel('Especpar');
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% [a2,b2]=irf_spectrogram(h(2),specparE,'log');
% hold(h(2),'on');
% irf_plot(h(2),Flh,'color','k','LineWidth',1.5)
% irf_plot(h(2),Fce,'color','r','LineWidth',1.5)
% irf_plot(h(2),Fce01,'color','w','LineWidth',1.5)
% irf_plot(h(2),Fce05,'color','c','LineWidth',1.5)
% hold(h(2),'off');
% % irf_legend(h(2),'c)',[0.01 0.95],'color','w','fontsize',12)
% %caxis(h(2),[-6 0]);
% set(h(2),'yscale','log');
% set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
% ylabel(h(2),{'f','(Hz)'},'fontsize',8,'Interpreter','tex');
% ylabel(b2,{'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',8);
% grid(h(2),'off');
% colormap(h(2),jet); 


% % h(4)=irf_panel('Especpar');
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% [a3,b3]=irf_spectrogram(h(3),specperpE,'log');
% hold(h(3),'on');
% irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
% irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
% irf_plot(h(3),Fce01,'color','w','LineWidth',1.5)
% irf_plot(h(3),Fce05,'color','c','LineWidth',1.5)
% hold(h(3),'off');
% % irf_legend(h(3),'d)',[0.01 0.95],'color','w','fontsize',12)
% %caxis(h(3),[-6 0]);
% set(h(3),'yscale','log');
% set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
% ylabel(h(3),{'f','(Hz)'},'fontsize',8,'Interpreter','tex');
% ylabel(b3,{'log_{10} E_{กอ}^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',8);
% grid(h(3),'off');
% colormap(h(3),jet); 


% h(6)=irf_panel('Bscmhf');
% mmsColors=[0 0 1; 0 1 0; 1 0 0];
% set(h(6),'ColorOrder',mmsColors)
% hold(h(6),'on');
% irf_plot(h(6),Bscmfachf);
% hold(h(6),'off');
% ylabel(h(6),{'\delta B (nT)'},'Interpreter','tex','fontsize',12);
% irf_legend(h(6),{'B_{\perp 1}','B_{\perp 2}','B_{||}'},[0.98 0.12],'fontsize',15)
% irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',12)
% irf_legend(h(6),'f > 10 Hz',[0.1 0.1],'color','k','fontsize',12)
% set(h(6),'ylim',[-2 2],'ytick',[-1:1:1])
% grid(h(6),'off');


%logB
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
[a3,b3]=irf_spectrogram(h(3),specB,'log');
hold(h(3),'on');
irf_plot(h(3),Flh,'color','k','LineWidth',1.5)
irf_plot(h(3),Fce,'color','r','LineWidth',1.5)
irf_plot(h(3),Fce01,'color','w','LineWidth',1.5)
irf_plot(h(3),Fce05,'color','c','LineWidth',1.5)
hold(h(3),'off');
% irf_legend(h(4),'e)',[0.01 0.95],'color','w','fontsize',12)
 %caxis(h(4),[-9 -4]);
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(3),{'f','(Hz)'},'fontsize',8,'Interpreter','tex');
ylabel(b3,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',8);
grid(h(3),'off');
colormap(h(3),jet); 



c_eval('title(h(1),''MMS?'',''fontsize'',12)',ic);
irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',tint);
set(h(1:3),'fontsize',8);
set(gcf,'render','painters');
% figname=['wave1'];
% print(gcf, '-dpdf', [figname '.pdf']);
% % set(gcf,'paperpositionmode','auto')
% % irf_adjust_panel_position;
%  set(gcf,'render','painters');
% set(h,'box','on');
% 
 figname=['waveMP'];
% 
 print(gcf, '-dpng', [figname '.png']);
