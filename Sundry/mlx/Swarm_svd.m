clear;clc
SwarmA = importdata('/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/Sundry/mlx/SwarmA_NED.dat');
Bfull = SwarmA.data(:,1:3);
Bdelta = SwarmA.data(:,4:6);

f = 50;
t = transpose((0:size(Bfull,1)-1)/f);
E = [t, nan(size(Bfull))];
Bfull = [t,Bfull];
Bdelta = [t,Bdelta];
%%
lf = 0.1; hf = 16;
res = irf_ebsp(E,Bdelta,[],Bfull,[],[lf, hf],'polarization','fac');
spec_t=res.t;
spec_f=res.f;

B_power_spectrum=res.bb_xxyyzzss(:,:,4);
specBt=struct('t',spec_t);
specBt.f=spec_f;
specBt.p=B_power_spectrum;
specBt.f_label='';
specBt.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};


theta_k=res.k_tp(:,:,1);


specBthe=struct('t',spec_t);
specBthe.f=spec_f;
specBthe.p=theta_k;
specBthe.f_label='';
specBthe.p_label={'theta','degree'};

idx=isnan(theta_k);

planarity_k=res.planarity;
specBpla=struct('t',spec_t);
specBpla.f=spec_f;
% planarity_k(idx)=NaN;
specBpla.p=planarity_k;
specBpla.f_label='';
specBpla.p_label={'planarity',''};

ellipticity_k=res.ellipticity;
specBell=struct('t',spec_t);
specBell.f=spec_f;
% ellipticity_k(idx)=NaN;
specBell.p=(ellipticity_k);
specBell.f_label='';
specBell.p_label={'ellipticity',''};
%% Init figure
n_subplots=5;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;           
set(gcf,'PaperUnits','centimeters');
xSize = 19; ySize = 23; coef=floor(min(600/xSize,600/ySize));
xLeft = (20-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[20 20 700 800]);
colormap(jet);


%% B_total
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bdelta(:,1) Bdelta(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bdelta(:,1) Bdelta(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bdelta(:,1) Bdelta(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bdelta(:,1) Bdelta(:,2).*0], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([B_bg(:,1) B_bg(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([B_bg(:,1) B_bg(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([B_bg(:,1) B_bg(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([B_bg(:,1) B_bg(:,2).*0], 'color','k', 'Linewidth',0.75); hold on;
grid off;
ylabel({'B [nT]'},'fontsize',8);
% set(gca,'ColorOrder',[0 0 0]);
set(gca,'Ylim',[-150 150]);
% irf_legend(gca,{'|B|'},[0.1 0.12],'fontsize',8);



%% B_psd_total
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(2), hcb2]=irf_spectrogram(h(2),specBt,'log');
hold on
ylabel({'f','(Hz)'},'fontsize',12);
poscbarx=get(hcb2,'pos');
set(gca,'yscale','lin');
% set(gca,'yscale','log');
set(gca,'Ylim',[lf hf]);
% set(gca,'Ylim',[1 16]);
caxis(h(2),[-4 0]);
% caxis(h(1),[-3.5 -2]);
poscbarx(3)=poscbarx(3)*0.5;
set(hcb2,'pos',poscbarx);
set(hcb2,'fontsize',8);

%% theta
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(3), hcb3]=irf_spectrogram(h(3),specBthe,'lin');
hold on
% irf_plot(gca,Flh,'color','w','LineWidth',0.75);
ylabel({'f','(Hz)'},'fontsize',12);
poscbarx=get(hcb3,'pos');
set(gca,'yscale','lin');
% set(gca,'yscale','log');
set(gca,'Ylim',[lf hf]);
% caxis(h(3),[-2 0]);
% caxis(h(1),[-3.5 -2]);
poscbarx(3)=poscbarx(3)*0.5;
set(hcb3,'pos',poscbarx);
set(hcb3,'fontsize',8);

%% pla
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(4), hcb4]=irf_spectrogram(h(4),specBpla,'lin');
hold on
% irf_plot(gca,Flh,'color','w','LineWidth',0.75);
ylabel({'f','(Hz)'},'fontsize',12);
poscbarx=get(hcb4,'pos');
set(gca,'yscale','lin');
% set(gca,'yscale','log');
set(gca,'Ylim',[lf hf]);
% set(gca,'Ylim',[1 16]);
caxis(h(4),[0 1]);
% caxis(h(1),[-3.5 -2]);
poscbarx(3)=poscbarx(3)*0.5;
set(hcb4,'pos',poscbarx);
set(hcb4,'fontsize',8);

%% ell
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(5), hcb5]=irf_spectrogram(h(5),specBell,'lin');
hold on
% irf_plot(gca,Flh,'color','w','LineWidth',0.75);
ylabel({'f','(Hz)'},'fontsize',12);
poscbarx=get(hcb5,'pos');
set(gca,'yscale','lin');
% set(gca,'yscale','log');
set(gca,'Ylim',[lf hf]);
% set(gca,'Ylim',[1 16]);
% caxis(h(3),[-2 0]);
% caxis(h(1),[-3.5 -2]);
poscbarx(3)=poscbarx(3)*0.5;
set(hcb5,'pos',poscbarx);
set(hcb5,'fontsize',8);


for ii=1:i_subplot-1
    pospanel(ii,:)=get(h(ii),'pos');
end
for ii=1:i_subplot-2
    pospanel(ii,2)=pospanel(ii,2)+(i_subplot-1-ii)*0.01*0.8;
    set(h(ii),'pos',pospanel(ii,:));
end

  set(h(1:end),'fontsize',8);
  irf_adjust_panel_position;
  %% save figure
  set(gcf,'render','painters');%矢量图
%   set(gcf,'visible','off');
%   figname=['MS_WAVE――1026'];
%   print(gcf, '-dpdf', [figname '.pdf'])
%   print(gcf, '-dpng', [figname '.png'])

%%
% irf_zoom(Tint,'x',h(1:n_subplots));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)
set(gca,"XTickLabelRotation",0)
%%
[Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_means(Bdelta,Bfull);
%%
figure(2)
n_subplots=2;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;           
set(gcf,'PaperUnits','centimeters');
xSize = 19; ySize = 29; coef=floor(min(600/xSize,600/ySize));
xLeft = (20-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[20 20 400 800]);
colormap(jet);
%%
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_spectrogram(h(1),Bpsd,'log');
caxis(h(1),[-7 -2]);
set(h(1),'yscale','log');
set(h(1),'Ylim',[lf,hf]);
% set(h(1),'ytick', [ 0.02 0.1 1 4 ]);
grid(h(1),'off');
ylabel(h(1),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
set(h(1),'FontSize',10); 
colormap(h(1), jet);
%%
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
specWaveangle=struct('t',waveangle.t);
specWaveangle.f=waveangle.f;
specWaveangle.p=waveangle.p;
specWaveangle.f_label='';
specWaveangle.p_label={' ','k-B angle'};
[h(2) hcb2]=irf_spectrogram(h(2),specWaveangle,'lin');
ylabel({'f','(Hz)'},'fontsize',12);
set(gca,'ylim',frange);
set(gca,'yscale','log');

%%
% irf_zoom(Tint,'x',h(1:n_subplots));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)
set(gca,"XTickLabelRotation",0)

