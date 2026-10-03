


clear
%--------------------------------------------

% Tsta='2003-10-09T02:05:00Z';
% Tend='2003-10-09T02:30:00Z';

Tsta='2003-10-09T02:21:20Z';
Tend='2003-10-09T02:28:20Z';

Tevent='2003-10-09T02:24:32Z'; 
%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
Tevt=iso2epoch(Tevent);
%--------------------------------------------

% caa_download(tint,'C*_CP_FGM_FULL')
% caa_download(tint,'C*_CP_EFW_L2_E3D_GSE')
% caa_download(tint,'C*_CP_STA_CWF_GSE'); 

caa_load FGM
caa_load EFW
caa_load STA


ic=4
%background magnetic field
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic); varname=irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
dobjname=irf_ssub('C?_CP_EFW_L2_E3D_INERT',ic);  varname=irf_ssub('E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT',ic);
c_eval(['E?=getmat(' dobjname ',''' varname ''');'],ic);
c_eval(['B?=irf_tlim(B?,tint);'],ic);
c_eval(['E?=irf_tlim(E?,tint);'],ic);

%low hybrid frequency
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic);  Bmag=irf_ssub('B_mag__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['BmagC?=getmat(' dobjname ',''' Bmag ''');'],ic);
c_eval(['BmagC?=irf_tlim(BmagC?,tint);'],ic);
c_eval(['fce_C?=irf_multiply(28,BmagC?,1,BmagC?,0);'],ic);
c_eval(['fci_C?=irf_multiply(0.015,BmagC?,1,BmagC?,0);'],ic);
c_eval(['flh_C?=irf_multiply(1,fce_C?,0.5,fci_C?,0.5);'],ic);

c_eval(['Btst=B?;'],ic);
c_eval(['Etst=E?;'],ic);
c_eval(['fci_tst=fci_C?;'],ic);
c_eval(['flh_tst=flh_C?;'],ic);


%% fft 
%nopfft=4096*16;  %number of point in fft
nopfft=1024*2; 
steplength = nopfft/2;   %step length

%B fft
fs=1/(Btst(2,1)-Btst(1,1));   %sampling rate in Hz
[Bzcomplex, freqB, timeB]=irf_wavefft(Btst(:,4), 'hamming', steplength, nopfft, fs);
timeB=linspace(Btst(1,1),Btst(end,1),length(timeB));
Bz_psd=Bzcomplex.' .* Bzcomplex';

indB_evt=find(timeB>=Tevt); indB_evt=indB_evt(1);
Bz_psd_evt=Bz_psd(indB_evt,:);


%E fft
fs=1/(Etst(2,1)-Etst(1,1));   %sampling rate in Hz
[Eycomplex, freqE, timeE]=irf_wavefft(Etst(:,3), 'hamming', steplength, nopfft, fs);
timeE=linspace(Etst(1,1),Etst(end,1),length(timeE));
Ey_psd=Eycomplex.' .* Eycomplex';

indE_evt=find(timeE>=Tevt); indE_evt=indE_evt(1);
Ey_psd_evt=Ey_psd(indE_evt,:);


%% save as structure format
Bz_spec=struct('t',timeB,'f',freqB','p',Bz_psd,'f_unit','Hz');
Ey_spec=struct('t',timeE,'f',freqE','p',Ey_psd,'f_unit','Hz');



%% Init figure
set(0,'defaultLineLineWidth', 0.5);
set(0,'defaultAxesFontSize', 10);
set(0,'defaultTextFontSize', 10);
set(0,'defaultAxesFontUnits', 'pixels');
figNumber=figure( ...
          'Name','Dataset coverage', ...
          'Tag','XYXgsm');clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = 1; yTop = -5;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])



%% Bz_psd plot
h(1)=axes('position',[0.15 0.83 0.6 0.15]); % [x y dx dy]
irf_spectrogram(gca, Bz_spec); hold on;
irf_plot(fci_tst,'k', 'LineWidth',0.75); hold on;
irf_plot(flh_tst,'k', 'LineWidth',0.75); hold off;
grid off;

ylabel(gca,'f Bz [Hz]');
set(gca,'Yscale','lin'); 


%% Ey_psd plot
h(2)=axes('position',[0.15 0.675 0.6 0.15]); % [x y dx dy]
irf_spectrogram(gca, Ey_spec); hold on;
irf_plot(fci_tst,'k', 'LineWidth',0.75); hold on;
irf_plot(flh_tst,'k', 'LineWidth',0.75); hold off;
grid off;

ylabel(gca,'f Ey [Hz]');
set(gca,'Yscale','lin'); 


%% PSD power law 1
h(3)=axes('position',[0.15 0.4 0.3 0.18]); % [x y dx dy]
plot(gca, freqB,Bz_psd_evt, 'k', 'LineWidth',1.3); hold on;
plot(gca, freqE,Ey_psd_evt, 'r', 'LineWidth',1.3); hold off;
grid on;

Tlab=epoch2iso(timeB(indB_evt));
irf_legend(gca,{[Tlab(12:22) 'UT']},[0.35, 1.02], 'color','k', 'FontSize',11);
ylabel(gca,'PSD, nT^2/Hz, (mV/m)^2/Hz');

  
%% PSD power law 2
h(4)=axes('position',[0.5 0.4 0.3 0.18]); % [x y dx dy]
plot(gca, freqB,Bz_psd_evt, 'k', 'LineWidth',1.3); hold on;
plot(gca, freqE,Ey_psd_evt, 'r', 'LineWidth',1.3); hold off;
grid on;

Tlab=epoch2iso(timeB(indB_evt));
irf_legend(gca,{[Tlab(12:22) 'UT']},[0.35, 1.02], 'color','k', 'FontSize',11);


%% refine plot
set(h(1:2),'ylim',[0.1 12]);  
set(h(3:4),'xscale','log', 'yscale','log');
set(h(3:4),'xlim',[1e-3 1e2], 'xtick',[1e-2 1e-1 1e0 1e1 1e2]);
set(h(3:4),'ylim',[1e-8 1e6], 'ytick',[1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4 1e6]);  

    
%% adjust colorbar position
colormap(hs_jet);
hcol1=colorbar('peer',h(1),'EastOutside');  
posFig=get(h(1),'Position'); 
left=posFig(1)+posFig(3)+0.015; low=posFig(2); width=0.015; height=posFig(4);
set(hcol1,'Position',[left low width height]);
ylabel(hcol1,'nT^2/Hz');
hcol2=colorbar('peer',h(2),'EastOutside');  
posFig=get(h(2),'Position'); 
left=posFig(1)+posFig(3)+0.015; low=posFig(2); width=0.015; height=posFig(4);
set(hcol2,'Position',[left low width height]);
ylabel(hcol2,'(mV/m)^2/Hz');



%% Annotation
Tsta='2003-10-09T02:24:24Z'; Tend='2003-10-09T02:24:40Z'; tint=[iso2epoch(Tsta) iso2epoch(Tend)];
irf_zoom(tint,'x',h(1:2)) 

%% save
set(gcf,'render','painters');
date=[Tlab(1:4) Tlab(6:7) Tlab(9:10) '_' Tlab(12:13) Tlab(15:16) Tlab(18:19)]; 
figname=['turbulence_power_law_fft_' date];
print(gcf, '-dpdf', [figname '.pdf']);

% set(gcf,'paperpositionmode','auto'); % to get the same on paper as on screen
% set(gcf,'render','painters');
% date=[Tlab(1:4) Tlab(6:7) Tlab(9:10) '_' Tlab(12:13) Tlab(15:16) Tlab(18:19)]; 
% figname=['turbulence_power_law_fft_' date];
% print(gcf, '-depsc2', '-loose', [figname '.eps']);  % -loose is necessary for MATLAB2014 and later

%In terminal window, run the following commands
%epstool --copy --bbox delme.eps delme_crop.eps
%ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None delme_crop.eps
  
  
  


