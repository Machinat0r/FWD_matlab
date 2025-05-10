 Tintr = irf.tint('2017-07-06T17:47:06.00Z/2017-07-06T17:47:14.00Z');
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tintr); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tintr); toc;',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',Tintr,?);',ic);
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',Tintr,?);',ic) 
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
load agyrotropy640
load energy_pad6


aa=agyrotropy40(:,2);
for ii = 1:size(aa,1) %数组行数
    tt(ii,1)=irf_time(agyrotropy40{ii,1},'epochtt>epoch');
    dd=aa{ii}.CData;
    dd3=sort(dd,2);  %各行元素升序排列
    c=nanmedian(dd3,2); %为该行从大到小排列的中间值
    for jj = 1:32
        Q_agyo(ii,jj)=(dd3(jj,26)-dd3(jj,6))/abs(c(jj));
    end   
end
n_subplots = 3;
i_subplot = 1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 5; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
Tintr = irf.tint('2017-07-06T17:47:06.00Z/2017-07-06T17:47:14.00Z');
%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
%irf_plot([B(:,1) B(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B(:,1) B(:,3)], 'color','g', 'Linewidth',0.75); hold on;

c_eval(['B4=irf.ts2mat(gsmB4);'],ic);
 irf_plot([B4(:,1) B4(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('B [nT]','fontsize',10);
%  set(gca,'Ylim',[fix(min([min(B(:,2)) min(B(:,3)) min(B(:,4))])/10)*10-10 fix(max(Bt(:,2))/10)*10+10]);
set(gca,'Ylim',[-5 15], 'ytick',[-25 -20 -10 -5  0 5 10 15 20 25 30 ]);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
 irf_legend(gca,{'B_x','B_y','B_z'},[0.05 0.97]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
%irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)


%% Te plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
c_eval(['Q4=irf.ts2mat(Q4);'],ic);

% irf_plot([Te_para(:,1) (Te_para(:,2)+2*Te_perp(:,2))/3], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Q4(:,1) Q4(:,2)], 'color','r', 'Linewidth',0.75); hold on;
%irf_plot([Te_perp1(:,1) Te_perp1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([Te(:,1) Te(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
ylabel({'Q';''},'fontsize',8);
%set(gca,'Ylim',[fix(min([min(Te_para1(:,2)) min(Te_para1(:,2)) min(Te_perp1(:,2))])/10)*10-10 ...
     %fix(max([max(Te_para1(:,2)) max(Te_para1(:,2)) max(Te_perp1(:,2))])/10)*10+10]);
% set(gca,'Ylim',[0 8000], 'ytick',[ 0  2000  4000 6000]);

% irf_legend(gca,'f',[0.01 0.98],'color','k','fontsize',12);
set(gca,'ColorOrder',[[1 0 0];[0 0 1];[0 0 0]]);
%irf_legend(gca,{'Te\{para}','Te\perp'},[0.05 0.92]);
%% Q_agyotropy

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i_subplot-1),jet);
for aa=1:size(Q_agyo,1)
    for bb=1:size(Q_agyo,2)
        if Q_agyo(aa,bb)>1
            Q_agyo(aa,bb)=NaN;
        end
    end
end

[ii,jj]=find(isnan(Q_agyo));
Q_agyo(ii(1),jj(1))=-max(max(Q_agyo));


Q_agy=struct('t',tt);
Q_agy.f= energy_pad(1,:);%energy levels,%1xM
Q_agy.p=Q_agyo;%data matrix% N*M
Q_agy.f_label='';
Q_agy.p_label={' ','(psd1-psd2)/psd0'};
[h(i_subplot-1), hcb]=irf_spectrogram(h(i_subplot-1),Q_agy);
ylabel('Q agyotropy','fontsize',10)
set(h(i_subplot-1),'yscale','log');
set(h(i_subplot-1),'ytick',[1e1 1e2 1e3 1e4]);
set(h(i_subplot-1),'Ylim',[1e2 1e4]);

caxis(h(i_subplot-1),[0 0.05]);
% irf_legend(h(i_subplot-1),'h',[0.99 0.98],'color','k','fontsize',12);
% poscbar=get(hcb,'pos');
% poscbar(3)=poscbar(3)*0.5;
% set(hcb,'pos',poscbar);
set(h(1:(i_subplot-1)),'fontsize',8);
  irf_adjust_panel_position;
  irf_plot_axis_align(h)
  irf_zoom(Tintr,'x',h(1:(i_subplot-1)));
  set(gcf,'color','w');
 set(gcf,'render','painters');
 set(gcf,'paperpositionmode','auto')
 figname=['20170712Q4(1111)'];
 print(gcf, '-dpdf', [figname '.pdf']);

