function Plot_TCS(tint,IC,Name,flagTime)
%------written by Wending Fu, Apr.2025 in Beijing------------
% see also SDCFilenames,SDCFilesDownload,SDCDataMove
%% 底下就是原来的overview程序
global OutputDir

ic=IC;
% load B
B1_ts=mms.get_data('B_gsm_brst',tint,1);B2_ts=mms.get_data('B_gsm_brst',tint,2);
B3_ts=mms.get_data('B_gsm_brst',tint,3);B4_ts=mms.get_data('B_gsm_brst',tint,4);
Bt1_ts = B1_ts.abs; Bt2_ts = B2_ts.abs;
Bt3_ts = B3_ts.abs; Bt4_ts = B4_ts.abs;
B1 = irf.ts2mat(B1_ts);B2 = irf.ts2mat(B2_ts);
B3 = irf.ts2mat(B3_ts);B4 = irf.ts2mat(B4_ts);
Bt1 = irf.ts2mat(Bt1_ts); Bt2 = irf.ts2mat(Bt2_ts); 
Bt3 = irf.ts2mat(Bt3_ts); Bt4 = irf.ts2mat(Bt4_ts); 

tStart = max([B1(1),B2(1),B3(1),B4(1)]);
B1(B1(:,1)<=tStart,:)=[];B2(B2(:,1)<=tStart,:)=[];
B3(B3(:,1)<=tStart,:)=[];B4(B4(:,1)<=tStart,:)=[];

tEnd = min([B1(end,1),B2(end,1),B3(end,1),B4(end,1)]);
B1(B1(:,1)>=tEnd,:)=[];B2(B2(:,1)>=tEnd,:)=[];
B3(B3(:,1)>=tEnd,:)=[];B4(B4(:,1)>=tEnd,:)=[];

Pos = mms.get_data('R_gsm',tint);
R1 = Pos.gsmR1;R2 = Pos.gsmR2;R3 = Pos.gsmR3;R4 = Pos.gsmR4;
R1 = [Pos.time.epochUnix R1(:,1:3)];R2 = [Pos.time.epochUnix R2(:,1:3)];
R3 = [Pos.time.epochUnix R3(:,1:3)];R4 = [Pos.time.epochUnix R4(:,1:3)];

% J
[J_B,~] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);
J_B = irf_abs(J_B);
J_B(:,2:5) = J_B(:,2:5)*1e9;
%% Init figure
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
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,2)) min(B2(:,2)) min(B3(:,2)) min(B4(:,2))])-0.1 ...
    max([max(B1(:,2)) max(B2(:,2)) max(B3(:,2)) max(B4(:,2))])+0.1]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bx [nT]','fontsize',8);
i=i+1;
%% By plot
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)],'k--', 'Linewidth',0.75);hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,3)) min(B2(:,3)) min(B3(:,3)) min(B4(:,3))])-0.1 ...
    max([max(B1(:,3)) max(B2(:,3)) max(B3(:,3)) max(B4(:,3))])+0.1]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('By [nT]','fontsize',8);
i=i+1;
%% Bz plot
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,4)) min(B2(:,4)) min(B3(:,4)) min(B4(:,4))])-0.1 ...
    max([max(B1(:,4)) max(B2(:,4)) max(B3(:,4)) max(B4(:,4))])+0.1]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bz [nT]','fontsize',8);
i=i+1;
%% Btotal plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bt2(:,1) Bt2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt3(:,1) Bt3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bt4(:,1) Bt4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(Bt1(:,2)) min(Bt2(:,2)) min(Bt3(:,2)) min(Bt4(:,2))])-0.1 ...
    max([max(Bt1(:,2)) max(Bt2(:,2)) max(Bt3(:,2)) max(Bt4(:,2))])+0.1]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Btotal [nT]','fontsize',8);
i=i+1;

%% J_B plot
h(i)=irf_subplot(n,1,-i);
irf_plot([J_B(:,1) J_B(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(J_B(:,2)) min(J_B(:,3)) min(J_B(:,4)) min(J_B(:,5))])-1 ...
    max([max(J_B(:,2)) max(J_B(:,3)) max(J_B(:,4)) max(J_B(:,5))])+1],'fontsize',9);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'J_x','J_y','J_z','|J|'},[0.97 0.92]);
ylabel('J [nA/m^2]','fontsize',8);
i=i+1;
%%
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)
irf_pl_mark(h(1:n),flagTime,'k')
set(gca,"XTickLabelRotation",0)
%%  出图保存部分
    set(gcf,'render','painters');
    set(gcf,'paperpositionmode','auto')
    % figname = [OutputDir,'OverviewFig\',Name(2:end-2)];    
    mkdir([OutputDir,'Figs/']);
    figname = [OutputDir,'Figs/',Name];  
    colormap(jet)
    print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
    close all
end