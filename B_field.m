
%Notice: load the CAA data first (CAA_LOAD) 
cd E:\Cluster
Tsta='2002-08-21T08:18:30Z';
Tend='2002-08-21T08:20:30Z';
tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH

ic=4;
% %    Magnetic fields
c_eval("caa_download(tint,'C?_CP_FGM_FULL')",ic);
% caa_download(tint,'C2_CP_FGM_FULL');
% caa_download(tint,'C4_CP_FGM_FULL');

caa_load C  %load data from datebase form C

%background magnetic field

dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
c_eval('B?_gsm = irf_gse2gsm(B?_gse);',ic);
c_eval('B?=irf_abs(B?_gsm);',ic);

c_eval('B = B?;',ic);

%% Init figure
n_subplots=2;
i_subplot=1;
set(0,'DefaultAxesFontSize',9);%set property for font
set(0,'DefaultLineLineWidth', 1); %set property for line width
fn=figure(66);clf;
set(gcf,'PaperUnits','centimeters') %gcf=get handle to current figure
xSize = 13.7; ySize = 20; coef=floor(min(800/xSize,800/ySize)); % pevious 60
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Bz
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B(:,1) B(:,4)], 'color','k', 'Linewidth',0.75);hold on;
irf_plot([B(:,1) B(:,4)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('B_z [nT]')
set(gca,'Ylim',[-10 7], 'ytick',[-10 -5 0 5]);
grid off;

%% By plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B(:,1) B(:,3)], 'color','k', 'Linewidth',0.75);hold on;
irf_plot([B(:,1) B(:,3)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('B_y  [nT]')
set(gca,'Ylim',[-4.5 15], 'ytick',[-5 0 5 10 15]);
grid off;

% %% Bx plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([B(:,1) B(:,2)], 'color','k', 'Linewidth',0.75); 
% ylabel('B_x  [nT]')
% set(gca,'Ylim',[3 36], 'ytick',[5 15 25 35]);
% grid off

% %% B tot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([B(:,1) B(:,5)], 'color','k', 'Linewidth',0.75); 
% set(gca,'Ylim',[0 25], 'ytick',[0 5 10 15 20]);
% grid off;

% set(gca,'Ylim',[0.0 0.9])
%% Adjust the position
irf_adjust_panel_position;


%% Annotation
% irf_zoom(tint,'x',h)
% set(h(1:end-1),'XTickLabe','');
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;

% leg = 'abcdefghi';
% for ii=1:i_subplot-1
%     set(h(ii),'ColorOrder',[0 0 0]);
%     irf_legend(h(ii),{['(' leg(ii) ')']},[0.02, 0.05]);
% end

% set(gcf,'render','painters');
% figname=['B_field' num2str(ic)];
% print(gcf, '-dpdf', [figname '.pdf']);


