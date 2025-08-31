function ESWPlot(tint,IC,Name)
% see also SDCFilenames,SDCFilesDownload,SDCDataMove
%% 底下就是原来的overview程序
    global OutputDir

    for ic=IC
    % load B
    c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
    c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
    c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
    %  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
    c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);

    % load E
    c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
    %%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
    c_eval(['Et?_ts=E?_ts.abs;'],ic); 
    c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
    c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
    c_eval('Et? = irf_abs(E?);',ic);
    c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

    % load N
    Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
    Ne1=irf.ts2mat(Ne1_ts);
    end

%% Init figure
n=3;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% B plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) 0*Bt1(:,2)], 'k--', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
%irf_plot([Bbf1(:,1) Bbf1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))]))-1 fix(max(Bt1(:,2)))+1]);
%     set(gca,'Ylim',[-12 22], 'ytick',[-5 0 5 10 15 20 25],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B (nT)','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)

i=i+1;
%% N plot
h(i)=irf_subplot(n,1,-i);

c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
% c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2))])-0.02]) max([max(Ne?(:,2))])+0.02]);",ic)

  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
  set(gca,'xtick',[])
 irf_legend(gca,{'Ne'},[0.97 0.92]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
ylabel('N [cm^{-3}]','fontsize',8);
i=i+1;
%% Efac
h(i)=irf_subplot(n,1,-i);
irf_plot([Efac1(:,1) Efac1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Efac1(:,1) Efac1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Efac1(:,1) Efac1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Efac1(:,1) Efac1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Efac (mV/m)','fontsize',10)
%     set(gca,'Ylim',[-8 8], 'ytick',[-6 -4 -2 0 2 4 6 ]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
set(gca,'Ylim',[fix(min([min(Efac1(:,2)) min(Efac1(:,3)) min(Efac1(:,4))]))-1 fix(max(Efac1(:,2:4),[],'all'))+1]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.1 0.12]);

pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
i=i+1;
%%
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)
% for ii = 1:n  
%     tempTime = flagTime(1); %2s之内相近的符合判据的时间只画一条线
%     irf_pl_mark(h(ii),tempTime,'k')
%     id_flagTime = [1];
%     for jj = 2:length(flagTime)
%         if flagTime(jj) - tempTime >= 2
%             tempTime = flagTime(jj);
%             irf_pl_mark(h(ii),tempTime,'k');
%             id_flagTime(end+1) = jj;
%         end
%     end
% end
%%  出图保存部分
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
% figname = [OutputDir,'OverviewFig\',Name(2:end-2)];    
figname = [OutputDir,'OverviewFig\', Name];  
colormap(jet)
print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
close all
end