


clear
%--------------------------------------------
% Tsta='2003-10-09T02:05:00Z';
% Tend='2003-10-09T02:30:00Z';

Tsta='2003-10-09T02:05:00Z';
Tend='2003-10-09T02:15:00Z';

Tsta='2003-10-09T02:20:00Z';
Tend='2003-10-09T02:30:00Z';

Tsta='2003-10-09T02:24:31Z';
Tend='2003-10-09T02:24:33Z';

Tevent='2003-10-09T02:24:10Z'; 
%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
Tevt=iso2epoch(Tevent);
%--------------------------------------------

% caa_download(tint,'C4_CP_FGM_FULL_ISR2')
% caa_download(tint,'C4_CP_EFW_L2_E3D_INERT')
% caa_download(tint,'C4_CP_STA_CWF_ISR2'); 

caa_load FGM
caa_load EFW
caa_load STA


ic=4
%background magnetic field
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic); varname=irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
dobjname=irf_ssub('C?_CP_EFW_L2_E3D_INERT',ic);  varname=irf_ssub('E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT',ic);
c_eval(['E?=getmat(' dobjname ',''' varname ''');'],ic);
B=irf_tlim(B4,tint);
E=irf_tlim(E4,tint);

%low hybrid frequency
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic);  Bmag=irf_ssub('B_mag__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['BmagC?=getmat(' dobjname ',''' Bmag ''');'],ic);
Bmag=irf_tlim(BmagC4,tint);
fce=irf_multiply(28,Bmag,1,Bmag,0);
fci=irf_multiply(0.015,Bmag,1,Bmag,0);
flh=irf_multiply(1,fce,0.5,fci,0.5);



%% fft 
fs=1/(B(2,1)-B(1,1));   %sampling rate in Hz
[Bz_psd,f_B] = pmtm(B(:,4),3.5,[],fs);

%E fft
fs=1/(E(2,1)-E(1,1));   %sampling rate in Hz
[Ey_psd,f_E] = pmtm(E(:,3),3.5,[],fs);


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
xLeft = 3; yTop = -1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


%% PSD power law 1
h(1)=axes('position',[0.15 0.4 0.35 0.18]); % [x y dx dy]
plot(gca, f_B,Bz_psd, 'k', 'LineWidth',0.5); hold on;
plot(gca, f_E,Ey_psd, 'r', 'LineWidth',0.5); hold off;
grid on;

irf_legend(gca,{[Tsta(12:19) '--' Tend(12:19) 'UT']},[0.35, 1.02], 'color','k', 'FontSize',11);
ylabel(gca,'PSD, nT^2/Hz, (mV/m)^2/Hz');


%% refine plot
set(h(1),'xscale','log', 'yscale','log');
set(h(1),'xlim',[1e-3 1e2], 'xtick',[1e-2 1e-1 1e0 1e1 1e2]);
set(h(1),'ylim',[1e-8 1e4], 'ytick',[1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4 1e6]);  



%% save
set(gcf,'render','painters');
date=[Tsta(1:4) Tsta(6:7) Tsta(9:10) '_' Tsta(12:13) Tsta(15:16) Tsta(18:19) '-' Tend(12:13) Tend(15:16) Tend(18:19)]; 
figname=['turbulence_power_law_pmtm_' date];
print(gcf, '-dpdf', [figname '.pdf']);


% set(gcf,'paperpositionmode','auto'); % to get the same on paper as on screen
% set(gcf,'render','painters');
% date=[Tsta(1:4) Tsta(6:7) Tsta(9:10) '_' Tsta(12:13) Tsta(15:16) Tsta(18:19) '-' Tend(12:13) Tend(15:16) Tend(18:19)]; 
% figname=['turbulence_power_law_pmtm_' date];
% print(gcf, '-depsc2', '-loose', [figname '.eps']);  % -loose is necessary for MATLAB2014 and later

%In terminal window, run the following commands
%epstool --copy --bbox delme.eps delme_crop.eps
%ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None delme_crop.eps
  




