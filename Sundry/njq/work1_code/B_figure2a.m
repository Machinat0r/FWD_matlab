clear
clc
mms.db_init('local_file_db','E:\fu\data2')

%% load data
ic=1:4;
hh1= 20;
tintStr = '27-Jan-2017';
% tint=irf.tint('2017-01-27T12:05:42.20Z/2017-01-27T12:05:44.20Z');
tint=irf.tint('2017-01-27T12:05:42.20Z/2017-01-27T12:05:44.20Z');
% tint=irf.tint('2017-01-27T12:05:33.20Z/2017-01-27T12:05:33.60Z');

c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);

c_eval('B?(:,2)=smooth(B?(:,2),hh1);',ic);
c_eval('B?(:,3)=smooth(B?(:,3),hh1);',ic);
c_eval('B?(:,4)=smooth(B?(:,4),hh1);',ic);

n=length(B2(:,1)); 
for ii=1:n
    c_eval('B?(ii,5)=norm(B?(ii,2:4));',ic);
end

%% plot
h=irf_plot(7,'newfigure');
xSize=750; ySize=700;
set(gcf,'Position',[100 100 xSize ySize]);
mmscolors=[0 0 1; 0 1 0; 1 0 0; 0 0 0]; %蓝绿红黑
set(h,'ColorOrder',mmscolors)
lnwid = 1;
h(1)=irf_panel('Bt');
hold(h(1),'on');
irf_plot(h(1),[B2(:,1) B2(:,2)],'Linewidth',lnwid);
irf_plot(h(1),[B2(:,1) B2(:,3)],'Linewidth',lnwid);
irf_plot(h(1),[B2(:,1) B2(:,4)],'Linewidth',lnwid);
irf_plot(h(1),[B2(:,1) B2(:,5)],'Linewidth',lnwid);
% irf_plot(h(1),Bxyz2.x,'Linewidth',lnwid);
% irf_plot(h(1),Bxyz2.y,'Linewidth',lnwid);
% irf_plot(h(1),Bxyz2.z,'Linewidth',lnwid);
% irf_plot(h(1),Bt2,'Linewidth',lnwid);
hold(h(1),'off');
% ylabel(h(1),{'B_{tLMN}','(nT)'},'Interpreter','tex');
ylabel(h(1),{'B','(nT)'},'Interpreter','tex');
set(h(1),'ylim',[-25 33],'ytick',[-20:20:20]);
irf_legend(h(1),{'x   ','y   ','z   ','total'},[0.98 0.03],'fontsize',12)
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)
grid(h(1),'off');

xSize = 400; ySize=550;
set(gcf,'position',[200 50 xSize ySize])
set(gcf,'render','painters');
Paper_X = 18; Paper_Y = 20; 
coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
set(gcf,'PaperSize', [Paper_X Paper_Y]); 
set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])

irf_plot_ylabels_align(h)
irf_zoom(h,'x',tint);
% set(gcf,'paperpositionmode','auto')
set(h,'fontsize',12);
figname = ['MMS2_figure2a'];
% print(gcf, '-dpdf','-r600',[figname '.pdf']);


% c_eval('VeLMN?=irf.ts2mat(Ve_LMN?);',ic);
% Vout=abs(VeLMN1(1:545,2));
% Vin=abs(VeLMN1(589:1087,4));
% Vout1=mean(Vout);
%  Vin1=mean(Vin);
 
% squared_data = Vin.^2;
% mean_of_squares = mean(squared_data);
% % 计算平方平均数（均方根）
% Vin1 = sqrt(mean_of_squares);