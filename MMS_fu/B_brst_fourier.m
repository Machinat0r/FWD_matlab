clear;close all
clc

mms.db_init('local_file_db','/Volumes/SPART-WORK/Data/MMS/')
%%
ic=1:4;
iic = 4;
% Tsta='2018-08-22T15:34:32.60Z';   
% Tend='2018-08-22T15:34:33.80Z';
% Tint=irf.tint(Tsta,Tend);
% 
% Tsta1='2018-08-22T15:34:34.30Z';   
% Tend1='2018-08-22T15:34:35.20Z';
% Tint1=irf.tint(Tsta1,Tend1);
% Tint=irf.tint('2017-07-26T07:21:13.00Z/2017-07-26T07:38:42.00Z');
% Tint=irf.tint('2017-07-26T07:21:13.00Z/2017-07-26T07:28:00.00Z');
% % % Tint=irf.tint('2016-09-27T01:18:21.000Z/2016-09-27T01:18:21.700Z');
% Tint=irf.tint('2021-07-21T13:19:58.150Z/2021-07-21T13:19:58.200Z');

TT = '2017-07-23T19:49:36.114Z/2017-07-23T19:49:36.125Z';
Tint=irf.tint(TT);
% Tsta='2016-09-27T01:18:21.210Z';   
% Tend='2016-09-27T01:18:21.220Z';

%%
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tint);',ic);

c_eval(['Bts=Bxyz.abs;'],ic); 
c_eval(['Bt=irf.ts2mat(Bts);'],ic);

c_eval('Exyz_gse?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
c_eval('Exyz?=irf_gse2gsm(Exyz_gse?);');
c_eval('Bscm_gse?=mms.get_data(''B_gse_scm_brst_l2'',Tint,?);',ic);
c_eval('Bscm?=irf_gse2gsm(Bscm_gse?);');

c_eval('Exyz = Exyz?;',iic);
c_eval('Bscm = Bscm?;',iic);
%% filt for delta E / delta B
lf = 400; hf = 500;
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm.time.epochUnix));
% 
Exyzf = Exyz.filt(lf,hf,dfE,3);
Bscmf = Bscm.filt(lf,hf,dfB,3);

Bscmf = [double(irf_time(Bscmf.time,'epochtt>epoch')), Bscmf.data];
Exyzf = [double(irf_time(Exyzf.time,'epochtt>epoch')), Exyzf.data];

Bscmf = irf_abs(Bscmf); Exyzf = irf_abs(Exyzf);
%% filt for four SCs
c_eval('Exyzf? = Exyz?.filt(lf,hf,dfE,3);',ic);
c_eval('Bscmf? = Bscm?.filt(lf,hf,dfB,3);',ic);

c_eval("Bscmf? = [double(irf_time(Bscmf?.time,'epochtt>epoch')), Bscmf?.data];",ic);
c_eval("Exyzf? = [double(irf_time(Exyzf?.time,'epochtt>epoch')), Exyzf?.data];",ic);

c_eval("Bscmf? = irf_abs(Bscmf?); Exyzf? = irf_abs(Exyzf?);",ic);
%% Init figure 1
n_subplots=3;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% delta B
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bscmf(:,1) Bscmf(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Bscmf(:,1) Bscmf(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaB (nT)');
set(gca,'ylim',[-0.04 0.04]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{x}','\deltaB_{y}','\deltaB_{z}'},[0.1 0.12]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);

%% delta E
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Exyzf(:,1) Exyzf(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Exyzf(:,1) Exyzf(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Exyzf(:,1) Exyzf(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Exyzf(:,1) Exyzf(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Exyzf(:,1) Exyzf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaE (mV/m)');
set(gca,'Ylim',[-20 20]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaE_{x}','\deltaE_{y}','\deltaE_{z}'},[0.1 0.12]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);
%% delta E/B
V_wave = 1e6 * Exyzf(:,5)./Bscmf(:,5); % m/s
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
ylabel('Velocity (mV/m)');
irf_plot([Exyzf(:,1) V_wave], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Exyzf(:,1) Exyzf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
%%
colormap(jet);

% irf_plot_axis_align(h(1:n_subplots));
irf_zoom('x',TT,h(1:n_subplots));
irf_plot_axis_align(h(1:n_subplots))

set(h(1:n_subplots),'fontsize',12);
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gca,"XTickLabelRotation",0)


%% Init figure 2
n_subplots=3;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% delta B x
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bscmf1(:,1) Bscmf1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bscmf2(:,1) Bscmf2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bscmf3(:,1) Bscmf3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bscmf4(:,1) Bscmf4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaBx (nT)');
set(gca,'ylim',[-0.04 0.04]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.1 0.12]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);

%% delta B y
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bscmf1(:,1) Bscmf1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bscmf2(:,1) Bscmf2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bscmf3(:,1) Bscmf3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bscmf4(:,1) Bscmf4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaBy (nT)');
set(gca,'ylim',[-0.04 0.04]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);
%% delta B z
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bscmf1(:,1) Bscmf1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bscmf2(:,1) Bscmf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bscmf3(:,1) Bscmf3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bscmf4(:,1) Bscmf4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaBz (nT)');
set(gca,'ylim',[-0.04 0.04]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);
%%
colormap(jet);

% irf_plot_axis_align(h(1:n_subplots));
irf_zoom('x',TT,h(1:n_subplots));
irf_plot_axis_align(h(1:n_subplots))

set(h(1:n_subplots),'fontsize',12);
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gca,"XTickLabelRotation",0)



%% Init figure 3
n_subplots=3;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(3);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% delta E x
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Exyzf1(:,1) Exyzf1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Exyzf2(:,1) Exyzf2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Exyzf3(:,1) Exyzf3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Exyzf4(:,1) Exyzf4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaEx (mV/m)');
% set(gca,'ylim',[-0.04 0.04]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.1 0.12]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);

%% delta E y
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Exyzf1(:,1) Exyzf1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Exyzf2(:,1) Exyzf2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Exyzf3(:,1) Exyzf3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Exyzf4(:,1) Exyzf4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaEy (mV/m)');
% set(gca,'ylim',[-0.04 0.04]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);
%% delta E z
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Exyzf1(:,1) Exyzf1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Exyzf2(:,1) Exyzf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Exyzf3(:,1) Exyzf3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Exyzf4(:,1) Exyzf4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bscmf(:,1) Bscmf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaEz (mV/m)');
% set(gca,'ylim',[-0.04 0.04]);
irf_legend(gca,{'400-500 Hz'},[0.8 0.12]);
%%
colormap(jet);

% irf_plot_axis_align(h(1:n_subplots));
irf_zoom('x',TT,h(1:n_subplots));
irf_plot_axis_align(h(1:n_subplots))

set(h(1:n_subplots),'fontsize',12);
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gca,"XTickLabelRotation",0)
