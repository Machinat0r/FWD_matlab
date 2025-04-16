clear
clc
mms.db_init('local_file_db','E:\fu\data2')

%% load data
ic=1:4;
tintStr = '27-Jan-2017';
% tint=irf.tint('2017-01-27T12:05:42.20Z/2017-01-27T12:05:44.20Z');
tint=irf.tint('2017-01-27T12:05:26.40Z/2017-01-27T12:05:59.00Z');
% tint=irf.tint('2017-01-27T12:05:33.20Z/2017-01-27T12:05:33.60Z');

c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('Et?=Exyz?.abs;',ic);
c_eval('Et_brst?=irf.ts2mat(Et?);',ic);
c_eval('Ni?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval('Ne?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);

c_eval('E_brst?=irf.ts2mat(Exyz?);',ic);
%电子离子速度
c_eval('Vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('Ve?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval('V?=irf.ts2mat(Ve?);',ic);
c_eval('Vit?=irf.ts2mat(Vi?);',ic);

c_eval('B?=irf.ts2mat(Bxyz?);',ic);

c_eval('E?=irf_resamp(E_brst?,B?);',ic);
c_eval('V?=irf_resamp(V?,B?);',ic);
c_eval('Vit?=irf_resamp(Vit?,B?);',ic);


c_eval('dsp?(:,1)=B?(:,1);',ic);
c_eval('dsp?(:,2:4)=E?(:,2:4)+cross(V?(:,2:4),B?(:,2:4))/1000;',ic);
% c_eval('dsp?=irf_resamp(dsp?,Vit?);',ic);
c_eval('dspt?=irf_abs(dsp?);',ic);

c_eval('dspi?(:,1)=B?(:,1);',ic);
c_eval('dspi?(:,2:4)=E?(:,2:4)+cross(Vit?(:,2:4),B?(:,2:4))/1000;',ic);
% c_eval('dsp?=irf_resamp(dsp?,Vit?);',ic);
c_eval('dspit?=irf_abs(dspi?);',ic);

%% omni flux
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',tint);',ic)
c_eval('Omni_flux_e?(:,1)=irf_time(energy_e?.DEPEND_0.data,''ttns>epoch'');',ic)
c_eval('Omni_flux_e?(:,2:33)=energy_e?.data;;',ic)
c_eval('channel_e?=transpose(energy_e?.DEPEND_1.data(1,1:32));;',ic)

c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',tint);',ic)
c_eval('Omni_flux_i?(:,1)=irf_time(energy_i?.DEPEND_0.data,''ttns>epoch'');',ic)
c_eval('Omni_flux_i?(:,2:33)=energy_i?.data;;',ic)
c_eval('channel_i?=transpose(energy_i?.DEPEND_1.data(1,1:32));;',ic)

%% coordinates transform
% L=[0.03 0.24 0.97];% B
% M=[0.1 0.96 -0.25];
% N=[-0.99 0.11 0.01];

L=[0.05 -0.21 0.98]; % Vi
M=[-0.23 0.95 0.22];
N=[-0.97 -0.24 -0.01];
c_eval('B_LMN?=irf_newxyz(Bxyz?,L,M,N);',ic);
c_eval('E_LMN?=irf_newxyz(Exyz?,L,M,N);',ic);
c_eval('Vi_LMN?=irf_newxyz(Vi?,L,M,N);',ic);
c_eval('Ve_LMN?=irf_newxyz(Ve?,L,M,N);',ic);

[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
n=length(j(:,1)); 
j_total=zeros(n,1);
j_E=zeros(n,1);
for ii=1:n
    j_total(ii,1)= norm(j(ii,2:4));
end

for iii=1:n
    j_E(iii,1)= dot(E1(iii,2:4),j(iii,2:4));
end 

j_LMN=irf_newxyz(j(:,2:4),L,M,N);

%% plot
h=irf_plot(7,'newfigure');
xSize=750; ySize=700;
set(gcf,'Position',[100 100 xSize ySize]);
mmscolors=[0 0 1; 0 1 0; 1 0 0; 0 0 0]; %蓝绿红黑
set(h,'ColorOrder',mmscolors)
lnwid = 1;
h(1)=irf_panel('Bt');
hold(h(1),'on');
irf_plot(h(1),B_LMN2.x,'Linewidth',lnwid);
irf_plot(h(1),B_LMN2.y,'Linewidth',lnwid);
irf_plot(h(1),B_LMN2.z,'Linewidth',lnwid);
irf_plot(h(1),Bt2,'Linewidth',lnwid);
hold(h(1),'off');
% ylabel(h(1),{'B_{tLMN}','(nT)'},'Interpreter','tex');
ylabel(h(1),{'B','(nT)'},'Interpreter','tex');
% set(h(1),'ylim',[-5 90],'ytick',[0:20:80]);
irf_legend(h(1),{'L   ','M   ','N   ','total'},[0.98 0.03],'fontsize',12)
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)
grid(h(1),'off');

h(2)=irf_panel('Et');
hold(h(2),'on');
% irf_plot(h(1),Bt2,'Linewidth',lnwid);
% irf_plot(h(1),Bt3,'Linewidth',lnwid);
% irf_plot(h(1),Bt4,'Linewidth',lnwid);
irf_plot(h(2),E_LMN2.x,'Linewidth',lnwid);
irf_plot(h(2),E_LMN2.y,'Linewidth',lnwid);
irf_plot(h(2),E_LMN2.z,'Linewidth',lnwid);
irf_plot(h(2),Et2,'Linewidth',lnwid);
hold(h(2),'off');
% ylabel(h(1),{'B_{tLMN}','(nT)'},'Interpreter','tex');
ylabel(h(2),{'E','(mV/m)'},'Interpreter','tex');
set(h(2),'ylim',[-200 200],'ytick',[-150:150:150]);
irf_legend(h(2),{'L   ','M   ','N   ','total'},[0.98 0.03],'fontsize',12)
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)
grid(h(2),'off');
% h(2)=irf_panel('BL');
Blmn1=irf.ts2mat(B_LMN1);
% hold(h(2),'on');
% c_eval('irf_plot(h(2),B_LMN?.x,''Linewidth'',lnwid)',ic);
% irf_plot(h(2),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.5);
% hold(h(2),'off');
% ylabel(h(2),{'B_{L}','(nT)'},'Interpreter','tex');
% % set(h(2),'ylim',[-40 70],'ytick',[-20:20:40]);
% irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(2),'off');


h(3)=irf_panel('VeLMN');
hold(h(3),'on');
irf_plot(h(3),Ve_LMN2.x,'Linewidth',lnwid);
irf_plot(h(3),Ve_LMN2.y,'Linewidth',lnwid);
irf_plot(h(3),Ve_LMN2.z,'Linewidth',lnwid);
irf_plot(h(3),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.5);%%虚线
hold(h(3),'off');
% ylabel(h(2),{'V_{etLMN}','(km/s)'},'Interpreter','tex');
ylabel(h(3),{'V_{e}','(km/s)'},'Interpreter','tex');
set(h(3),'ylim',[-1200 1400],'ytick',[-1000:500:1000]);
irf_legend(h(3),{'L   ','M   ','N   '},[0.98 0.03],'fontsize',12)
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
grid(h(3),'off');


% h(4)=irf_panel('VeM');
% hold(h(4),'on');
% c_eval('irf_plot(h(4),Ve_LMN?.y,''Linewidth'',lnwid)',ic);
% irf_plot(h(4),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.5);
% hold(h(4),'off');
% ylabel(h(4),{'V_{eM}','(km/s)'},'Interpreter','tex');
% % set(h(4),'ylim',[-1000 1500],'ytick',[-500:500:1000]);
% irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(4),'off');

% h(5)=irf_panel('VeN');
% hold(h(5),'on');
% c_eval('irf_plot(h(5),Ve_LMN?.z,''Linewidth'',lnwid)',ic);
% irf_plot(h(5),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.5);
% hold(h(5),'off');
% ylabel(h(5),{'V_{eN}','(km/s)'},'Interpreter','tex');
% % set(h(5),'ylim',[-1000 1500],'ytick',[-500:500:1000]);
% irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(5),'off');


h(4)=irf_panel('ViL');
hold(h(4),'on');
% c_eval('irf_plot(h(6),Vi_LMN?.x,''Linewidth'',lnwid)',ic);
irf_plot(h(4),Vi_LMN2.x,'Linewidth',lnwid)
irf_plot(h(4),Vi_LMN2.y,'Linewidth',lnwid)
irf_plot(h(4),Vi_LMN2.z,'Linewidth',lnwid)
irf_plot(h(4),[Blmn1(:,1) Blmn1(:,2)*0+52],'k--', 'Linewidth',0.5);
hold(h(4),'off');
% ylabel(h(6),{'V_{itLMN}','(km/s)'},'Interpreter','tex');
ylabel(h(4),{'V_{i}','(km/s)'},'Interpreter','tex');
set(h(4),'ylim',[-200 400],'ytick',[-100:200:300]);
irf_legend(h(4),{'L   ','M   ','N   '},[0.98 0.03],'fontsize',12)
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)
grid(h(4),'off');

h(5)=irf_panel('density');
hold(h(5),'on');
irf_plot(h(5),Ni2,'Linewidth',lnwid,'color','b');
irf_plot(h(5),Ne2,'Linewidth',lnwid,'color','r');
hold(h(5),'off');
set(h(5),'ylim',[-4 15],'ytick',[0:5:15]);
ylabel(h(5),{'n','(cm^-3)'},'Interpreter','tex');
irf_legend(h(5),{'n_{i} '},[0.93 0.08],'color','b','fontsize',12)
irf_legend(h(5),{'n_{e} '},[0.99 0.08],'color','r','fontsize',12)
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)
grid(h(5),'off');

% h(6)=irf_panel('dissipation');
% hold(h(6),'on');
% irf_plot(h(6),[dspt1(:,1) dspt1(:,5)],'Linewidth',lnwid,'color','k');
% irf_plot(h(6),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'linewidth',0.5);
% hold(h(6),'off');
% % ylabel(h(6),{'E+V_{e} X B','mV/m'},'Interpreter','tex');
% ylabel(h(6),{'\color{black}|E+V_{e} X B|','\color{black}(mV/m)'},'Interpreter','tex','fontsize',12);
% % set(h(6),'ylim',[-2 60],'ytick',[0:25:50]);
% irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(6),'off');
% 
% h(7)=irf_panel('dissipationi');
% hold(h(7),'on');
% irf_plot(h(7),[dspit1(:,1) dspit1(:,5)],'Linewidth',lnwid,'color','k');
% irf_plot(h(7),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'linewidth',0.5);
% hold(h(7),'off');
% % ylabel(h(6),{'E+V_{e} X B','mV/m'},'Interpreter','tex');
% ylabel(h(7),{'\color{black}|E+V_{i} X B|','\color{black}(mV/m)'},'Interpreter','tex','fontsize',12);
% set(h(7),'ylim',[0 115],'ytick',[0:50:100]);
% irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(7),'off');


%% Omniflux_2
ic=2;
h(6)=irf_panel('Ve_s');

c_eval('specrec_p_e=struct(''t'',Omni_flux_e?(:,1));',ic)
c_eval('specrec_p_e.f=transpose(energy_e?.DEPEND_1.data(1,1:32))',ic);%energy levels
c_eval('specrec_p_e.p=Omni_flux_e?(:,2:33);',ic)%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i_plot), hcb1]=irf_spectrogram(h(i_plot),specrec_p_e,);
% [a5]=irf_spectrogram(h(6),specrec_p_e,'donotshowcolorbar');%'donotfitcolorbarlabel，替换句
[a5]=irf_spectrogram(h(6),specrec_p_e);%'donotfitcolorbarlabel%替换句
%[a5]=irf_spectrogram(h(4),specrec_p_e)
% [a5]=irf_spectrogram(h(4),specrec_p_e,'Location', 'east')
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',12);
grid off;
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),{'\color{black}E_{e}','\color{black}(eV)'},'Interpreter','tex','fontsize',12);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h(6),'position');
% set(a5,'position',[Pos_c(1) Pos_c(2)-0.004 Pos_s(3) Pos_c(4)]);
% set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2)-0.004 0.02 Pos_c(4)]);
% set(a5,'fontsize',12); set(b5,'fontsize',12);
colormap(h(6),jet)
caxis(a5,[5, 9])

%% 
% h(8)=irf_panel('ViN');
% hold(h(8),'on');
% c_eval('irf_plot(h(8),Vi_LMN?.z,''Linewidth'',lnwid)',ic);
% irf_plot(h(8),[Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.5);
% hold(h(8),'off');
% ylabel(h(8),{'V_{iN}','(km/s)'},'Interpreter','tex');
% set(h(8),'ylim',[-120 120],'ytick',[-300:200:100]);
% irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(8),'off');

h(7)=irf_panel('Vi_s');
c_eval('specrec_p_i=struct(''t'',Omni_flux_i?(:,1));',ic)
c_eval('specrec_p_i.f=transpose(energy_i?.DEPEND_1.data(1,1:32))',ic);%energy levels
c_eval('specrec_p_i.p=Omni_flux_i?(:,2:33);',ic)%data matrix
specrec_p_i.f_label='';
specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i_plot), hcb1]=irf_spectrogram(h(i_plot),specrec_p_e,);
% [a5]=irf_spectrogram(h(7),specrec_p_i,'donotshowcolorbar');%'donotfitcolorbarlabel，替换句
[a5]=irf_spectrogram(h(7),specrec_p_i);%替换句
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

irf_legend(h(7),'(g)',[0.99 0.98],'color','w','fontsize',12);
grid off;
set(h(7),'yscale','log');
set(h(7),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(7),{'\color{black}E_{i}','\color{black}(eV)'},'Interpreter','tex','fontsize',12);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h(7),'position');
% set(a5,'position',[Pos_c(1) Pos_c(2)-0.004 Pos_s(3) Pos_c(4)]);
% set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2)-0.004 0.02 Pos_c(4)]);
% set(a5,'fontsize',12); set(b5,'fontsize',12);
colormap(h(7),jet)
caxis(a5,[4, 8])



%% mark
tint2=irf.tint('2017-01-27T12:05:42.50Z/2017-01-27T12:05:43.80Z');
for ii = 1:7
    c_eval('irf_pl_mark(h(?),tint2,[0.7,0.7,0.7])',ii)
end

tint3=irf.tint('2017-01-27T12:05:30.40Z/2017-01-27T12:05:32.41Z');
for ii = 1:7
    c_eval('irf_pl_mark(h(?),tint3,[0.5,0.5,0.5])',ii)
end

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
figname = ['MMS2_colorbar'];
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