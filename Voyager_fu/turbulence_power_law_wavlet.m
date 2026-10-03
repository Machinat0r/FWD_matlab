


clear
%--------------------------------------------
Tsta='2003-10-09T02:05:00Z';
Tend='2003-10-09T02:30:00Z';

Tsta='2003-10-09T02:05:00Z';
Tend='2003-10-09T02:17:00Z';

Tevtst='2003-10-09T02:05:00Z';
Tevted='2003-10-09T02:15:00Z';

%--------------------------------------------
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
tint_evt=[iso2epoch(Tevtst) iso2epoch(Tevted)];
%--------------------------------------------

% caa_download(tint,'C3_CP_FGM_FULL_ISR2')
% caa_download(tint,'C3_CP_EFW_L2_E3D_INERT')
% caa_download(tint,'C3_CP_STA_CWF_ISR2'); 

caa_load FGM
caa_load EFW
caa_load STA

ic=3
%background magnetic field
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic); varname=irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
dobjname=irf_ssub('C?_CP_EFW_L2_E3D_INERT',ic);  varname=irf_ssub('E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT',ic);
c_eval(['E?=getmat(' dobjname ',''' varname ''');'],ic);
B=irf_tlim(B3,tint);
E=irf_tlim(E3,tint);

%low hybrid frequency
dobjname=irf_ssub('C?_CP_FGM_FULL_ISR2',ic);  Bmag=irf_ssub('B_mag__C?_CP_FGM_FULL_ISR2',ic);
c_eval(['BmagC?=getmat(' dobjname ',''' Bmag ''');'],ic);
Bmag=irf_tlim(BmagC3,tint);
fce=irf_multiply(28,Bmag,1,Bmag,0);
fci=irf_multiply(0.015,Bmag,1,Bmag,0);
flh=irf_multiply(1,fce,0.5,fci,0.5);

B=irf_resamp(B,E);
fs=1/(E(2,1)-E(1,1));   %sampling rate in Hz

Bz_spec=irf_wavelet(B(:,[1 4]));
Ey_spec=irf_wavelet(E(:,[1 3]));

pBtmp=Bz_spec.p{:}; tBtmp=Bz_spec.t; fB=Bz_spec.f; pBtmp2=[tBtmp pBtmp];
pEtmp=Ey_spec.p{:}; tEtmp=Ey_spec.t; fE=Ey_spec.f; pEtmp2=[tEtmp pEtmp];

pBtmp3=irf_tlim(pBtmp2,tint_evt); pB=pBtmp3(:,2:end);
pEtmp3=irf_tlim(pEtmp2,tint_evt); pE=pEtmp3(:,2:end);

BzPSD=nanmean(pB,1);
EyPSD=nanmean(pE,1);



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
irf_plot(fci,'k', 'LineWidth',0.75); hold on;
irf_plot(flh,'k', 'LineWidth',0.75); hold off;
grid off;

ylabel(gca,'f Bz [Hz]');
hcb=colorbar('peer',gca);
ylabel(hcb,'nT^2/Hz');


%% Ey_psd plot
h(2)=axes('position',[0.15 0.675 0.6 0.15]); % [x y dx dy]
irf_spectrogram(gca, Ey_spec); hold on;
irf_plot(fci,'k', 'LineWidth',0.75); hold on;
irf_plot(flh,'k', 'LineWidth',0.75); hold off;
grid off;

ylabel(gca,'f Ey [Hz]');
hcb=colorbar('peer',gca);
ylabel(hcb,'(mV/m)^2/Hz');


%% PSD power law 1
h(3)=axes('position',[0.3 0.4 0.3 0.18]); % [x y dx dy]
plot(gca, fB,BzPSD, 'k', 'LineWidth',0.75); hold on;
plot(gca, fE,EyPSD, 'r', 'LineWidth',0.75); hold off;
grid on;

irf_legend(gca,{[Tevtst(12:19) '-' Tevted(12:19) 'UT']},[0.15, 1.02], 'color','k', 'FontSize',11);
ylabel(gca,'PSD, nT^2/Hz, (mV/m)^2/Hz');


%% refine plot
set(h(1:2),'ylim',[0 12]);
set(h(3),'xscale','log', 'yscale','log');
set(h(3),'xlim',[1e-3 1e2], 'xtick',[1e-2 1e-1 1e0 1e1 1e2]);
set(h(3),'ylim',[1e-8 1e4], 'ytick',[1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4 1e6]);


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
irf_zoom(h(1:2),'x',tint);


%% save
set(gcf,'renderer','opengl');
date=[Tevtst(1:4) Tevtst(6:7) Tevtst(9:10) '_' Tevtst(12:13) Tevtst(15:16) Tevtst(18:19) '-' Tevted(12:13) Tevted(15:16) Tevted(18:19)]; 
figname=['turbulence_power_law_wavlet_' date];
print(gcf, '-dpng','-r300',[figname '.png']);



