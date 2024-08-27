clear;clc
close all
%% bug驱散令
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       南无电子阿弥陀佛驱散仿生bug
%                                _ooOoo_
%                               o8888888o
%                               88" . "88
%                               (| -_- |)
%                               O\  =  /O
%                            ____/`---'\____
%                          .'  \\|     |//  `.
%                         /  \\|||  :  |||//  \
%                        /  _||||| -:- |||||-  \
%                        |   | \\\  -  /// |   |
%                        | \_|  ''\-/''  |   |
%                        \  .-\__  `-`  ___/-. /
%                      ___`. .'  /-.-\  `. . __
%                   ."" '<  `.___\_<|>_/___.'  >'"".
%                  | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%                  \  \ `-.   \_ __\ /__ _/   .-` /  /
%             ======`-.____`-.___\_____/___.-`____.-'======
% 	                   `=-='
%                 天地玄宗，万气本根。广修亿劫，证吾神通。
%                 三界内外，惟道独尊。体有金光，覆映吾身。
%                 视之不见，听之不闻。包罗天地，养育群生。
%                 受持万遍，身有光明。三界侍卫，五帝司迎。
%                 万神朝礼，役使雷霆。鬼妖丧胆，精怪忘形。
%                 内有霹雳，雷神隐名。洞慧交彻，五炁腾腾。
%                金光速现，覆护真人。急急如律令，bug全去除！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/THEMIS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/THEMIS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TempDir = 'D:\THEMIS\tha\';mkdir(TempDir);
% TT = '2009-03-16T07:00:00Z\2009-03-16T09:00:00Z';
% TT = '2009-03-01T06:00:00Z\2009-03-01T08:00:00Z';
% TT = '2009-02-23T08:00:00Z\2009-02-23T010:00:00Z';
% TT = '2010-02-18T04:20:00Z\2010-02-18T05:40:00Z';
% TT = '2018-01-31T00:00:00Z\2018-01-31T23:59:59Z';
TT = '2024-05-10T00:00:00Z\2024-05-11T00:00:00Z';

Tsta = strsplit(TT,'\');
Tend = Tsta{2};Tsta = Tsta{1};
%% load data
tint = [iso2epoch(Tsta),iso2epoch(Tend)];
ic = {'a'};
c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','scm',TempDir)",ic);
c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','fgm',TempDir)",ic);
c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','efi',TempDir)",ic);
c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','esa',TempDir)",ic);
THEMISDataMove(TempDir,ParentDir)

mms.db_init('local_file_db',TempDir);
c_eval("B_? = th_read_l2_change_by_fwd('th?_fgl_gsm',tint);",ic);
c_eval("Vi_? = th_read_l2_change_by_fwd('th?_peir_velocity_gsm',tint);",ic);
c_eval("Ti_? = th_read_l2_change_by_fwd('th?_peir_t3',tint);",ic);
% c_eval("E_? = th_read_l2_change_by_fwd('th?_eff_dot0_gsm',tint);",ic);
c_eval("E_? = th_read_l2_change_by_fwd('th?_eff_dot0_gsm',tint);",ic);
% c_eval("Vires_? = irf_resamp(Vi_?,B_?);",ic);
% c_eval("Eres_? = irf_resamp(E_?,B_?);",ic);
% c_eval("Evxb_? = Eres_? + irf_cross(Vires_?,B_?)/1000;",ic);
% c_eval("Evxb_?(:,1) = Eres_?(:,1);",ic);

c_eval("Bres_? = irf_resamp(B_?,Vi_?);",ic);
c_eval("Eres_? = irf_resamp(E_?,Vi_?);",ic);
% c_eval("Evxb_? = Eres_? + irf_cross(Vi_?,Bres_?)/1000;",ic);
% c_eval("Evxb_?(:,1) = Eres_?(:,1);",ic);
%% lmn
% irf_minvar_gui(B_e);
% L=[-0.31 0.18 0.94];%最大变化方向
% N=[-0.95 -0.01 -0.31];%外法向
% M=cross(N,L);
% 
% c_eval('B_?=irf_newxyz(B_?,L,M,N);',ic);
% c_eval('Evxb_?=irf_newxyz(Evxb_?,L,M,N);',ic);
% c_eval('Vi_?=irf_newxyz(Vi_?,L,M,N);',ic);
%% Init figure
ic = 'b';
n=5;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Bx plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B_?(:,2))])-1 max([max(B_?(:,2))])+1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1]]);
% % % irf_legend(gca,{'B_x'},[0.97 0.92]);
% % % ylabel('B [nT]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Bx & Bz plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(B_?(:,2)) min(B_?(:,4))])-3 max([max(B_?(:,2)) max(B_?(:,4))])+3]);",ic);
set(gca,'Ylim',[-25 40], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% By plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B_?(:,1) B_?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min(B_?(:,3))-1 max(B_?(:,3))+1]);",ic);
    set(gca,'Ylim',[-2 6], 'ytick',[-8:4:8],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
irf_legend(gca,{'B_y'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%%  Bz plot
% % % h(i)=irf_subplot(n,1,-i);
% % % % c_eval("irf_plot([B_?(:,1) B_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B_?(:,1) B_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([B_?(:,1) 0*B_?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[0 max([max(B_?(:,4))])+1]);",ic);
% % % % c_eval("set(gca,'Ylim',[min([min(B_?(:,4))])-1 max([max(B_?(:,4))])+1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[1 0 0]]);
% % % irf_legend(gca,{'B_z'},[0.97 0.92]);
% % % ylabel('B [nT]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Vix plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Vi_?(:,1) Vi_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,2))])/10)*10-15 fix(max([max(Vi_?(:,2))])/10)*10+15],'fontsize',9);",ic);
set(gca,'Ylim',[-40 50], 'ytick',[-40:20:80],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1]]);
irf_legend(gca,{'Vi_x'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% Viy plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,3))])/10)*10-15 fix(min([max(Vi_?(:,3))])/10)*10+15],'fontsize',9);",ic);
set(gca,'Ylim',[-50 30], 'ytick',[-40:20:40],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 1 0]]);
irf_legend(gca,{'Vi_y'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% Viz plot
h(i)=irf_subplot(n,1,-i);
% c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([Vi_?(:,1) Vi_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Vi_?(:,4))])/10)*10-10 fix(max([max(Vi_?(:,4))])/10)*10+10],'fontsize',9);",ic);
set(gca,'Ylim',[-40 40], 'ytick',[-20 0 20],'fontsize',9);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[1 0 0]]);
irf_legend(gca,{'Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% E+vxb plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,3)], 'color','g', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,4)], 'color','r', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([Evxb_?(:,1) Evxb_?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[-100 100], 'ytick',[-100 -50 0 50 100],'fontsize',9);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(Evxb_?(:,2)) min(Evxb_?(:,3)) min(Evxb_?(:,4))])-2 max([max(Evxb_?(:,2)) max(Evxb_?(:,3)) max(Evxb_?(:,4))])+2]);",ic);
% % % set(gca,'Ylim',[-4 4], 'ytick',[-4:2:4],'fontsize',9);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('E [mV/m]','fontsize',12)
% % % i=i+1;
%%  出图保存部分

irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
Dir =   'C:\Matlab\bin\新建文件夹\fwd\Figures\';
figname = [Dir TT(1:10) '-' ic];
% print(gcf, '-dpdf', [figname '.pdf']);    