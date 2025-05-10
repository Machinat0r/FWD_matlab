clear;
clc;
mms.db_init('local_file_db','D:\Works\mms_db\data');
Units=irf_units;
me=Units.me;
tint1 = irf.tint('2017-07-17T07:53:00.00Z/2017-07-17T07:53:10.00Z');
tint2 = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
tint = irf.tint('2017-07-17T07:53:04.00Z/2017-07-17T07:53:07.00Z');
ic=1:4;
% load B

% load B
disp('Loading magnetic field...')
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);
c_eval(['B?=irf.ts2mat(gsmB?);'],ic);

% load E
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval(['gsmE?=irf_gse2gsm(gseE?);'],ic);
c_eval(['E?_fac = irf_convert_fac(gsmE?,gsmB?,[1 0 0]);'],ic);

c_eval('E?_fac = E?_fac.resample(gsmB1);',1:4);
c_eval(['E?_facmat=irf.ts2mat(E?_fac);'],ic);
c_eval(['E_perp? =sqrt(E?_facmat(:,2).^2+E?_facmat(:,3).^2);'],ic);

% 
R  = mms.get_data('R_gsm',tint1);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?);',1:4);
c_eval('R? = Rxyz?.resample(gsmB1);',1:4);
c_eval('R? = R?.tlim(tint);',1:4);

% caculate j with curl [Dunlop et al., 2005]
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','gsmB?');
j.data=1e9*j.data;
jmag=j.abs;
jmat=irf.ts2mat(j);
% gsmBav=gsmB1;
% gsmBav.data=(gsmB1.data+gsmB2.data+gsmB3.data+gsmB4.data)/4.0;
% j_fac = irf_convert_fac(j,gsmBav,[1 0 0]);
 c_eval(['j?_fac = irf_convert_fac(j,gsmB?,[1 0 0]);'],ic);
% c_eval(['j?_perp =sqrt(j?_fac(:,1).^2+j?_fac(:,2).^2);'],1:4);
c_eval('j?_fac = j?_fac.resample(gsmB1);',1:4);
c_eval(['j?_facmat=irf.ts2mat(j?_fac);'],ic);
 c_eval(['j_perp? =sqrt(j?_facmat(:,2).^2+j?_facmat(:,3).^2);'],ic);


%jdotE
c_eval('JdotE?=irf_dot([jmat(:,1) j?_fac.data(:,1) j?_fac.data(:,2) j?_fac.data(:,3)],[jmat(:,1) E?_fac.data(:,1) E?_fac.data(:,2) E?_fac.data(:,3)]);',1:4);
c_eval('JdotE?(:,2) = JdotE?(:,2)*1e-3;',1:4);%nW
%jdotE para
c_eval('JdotE?_para = [jmat(:,1) j?_fac.data(:,3).*E?_fac.data(:,3)];',1:4);
c_eval('JdotE?_para(:,2) = JdotE?_para(:,2)*1e-3;',1:4);%nW
%jdotE perp
c_eval('JdotE?_perp=irf_dot([jmat(:,1) j?_fac.data(:,1) j?_fac.data(:,2) 0*j?_fac.data(:,3)],[jmat(:,1) E?_fac.data(:,1) E?_fac.data(:,2) 0*E?_fac.data(:,3)]);',1:4);
c_eval('JdotE?_perp(:,2) = JdotE?_perp(:,2)*1e-3;',1:4);%nW

% c_eval('X?=JdotE?_para(:,2)./JdotE?(:,2);',1:4);%nW
% c_eval('X?=abs(X?)*1e2;',1:4);%nW
% c_eval('Y?=JdotE?_perp(:,2)./JdotE?(:,2);',1:4);%nW
% c_eval('Y?=abs(Y?)*1e2;',1:4);%nW
% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc
energy_e=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint);


%c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
%c_eval('ePDist?.data(ePDist?.data<ePDistErr?.data*1.1) = 0;',ic)

%% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval(['gsmTe?=irf_gse2gsm(gseTe?);'],ic);
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc
% c_eval(['gsmTi?=irf_gse2gsm(gseTi?);'],ic);
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gsmB?);',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gsmB?);',ic)
 %% Init figure
n_subplots=5;
i_subplot=1;
set(0,'DefaultAxesFontSize',5);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 10; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

if 1 % ne 
  h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',[[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{ne1},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',[[0 1 0];[0 0 1]])  
  %irf_legend(hca,{'n_e','n_i'},[0.1 0.2],'fontsize',12);
  grid(hca,'off');
end
if 0 % Te par perp
 h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  hca = irf_panel('Te');
  set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 0]])
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTe?.trace/3},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',[[0 1 0];[0 0 1];[0 0 0]])
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}','T_{e}'},[0.1 0.9],'fontsize',12);
  %hca.YLim = [0 4000];
  grid(hca,'off');
%     irf_zoom(hca,'y')
end

%Epara
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([E1_facmat(:,1) E1_facmat(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([E2_facmat(:,1) E2_facmat(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([E3_facmat(:,1) E3_facmat(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([E4_facmat(:,1) E4_facmat(:,4)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([E1_facmat(:,1) 0*E1_facmat(:,4)], 'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel({'E_{||}';'[mV/m]'},'fontsize',5);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
 irf_legend(gca,{'MMS1         ','MMS2         ','MMS3         ','MMS4'},[0.05 0.95],'fontsize',5); 


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([jmat(:,1) j1_facmat(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([jmat(:,1) j_perp1], 'color','b', 'Linewidth',0.75); hold off;
set(gca,'ColorOrder',[[0 1 0];[0 0 1]]);
irf_legend(gca,{'J||     ','J¡Í'},[0.05 0.95],'fontsize',5);
ylabel({'J';'[nA/m^2]'},'fontsize',5);

set(gca,'ylim',[-10 10]);
  
grid off;
%E.Jpara
 h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([JdotE1_para(:,1) JdotE1_para(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([JdotE2_para(:,1) JdotE2_para(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([JdotE3_para(:,1) JdotE3_para(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([JdotE4_para(:,1) JdotE4_para(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([JdotE1_para(:,1) 0*JdotE1_para(:,2)], 'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel({'(E\cdotJ)_{||}';'[nW/m^3]'},'fontsize',5); 
set(gca,'ylim',[-0.1 0.1]);
set(gca,'ytick',[-0.05 0 0.05]);   
%E.Jperp
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([JdotE1_perp(:,1) JdotE1_perp(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([JdotE2_perp(:,1) JdotE2_perp(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([JdotE3_perp(:,1) JdotE3_perp(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([JdotE4_perp(:,1) JdotE4_perp(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([JdotE1_para(:,1) 0*JdotE1_perp(:,2)], 'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel({'(E\cdotJ)_{\perp}';'[nW/m^3]'},'fontsize',5);
set(gca,'ylim',[-0.1 0.1]); 
set(gca,'ytick',[-0.05 0 0.05]);   
 % E.J
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([JdotE1(:,1) JdotE1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([JdotE2(:,1) JdotE2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([JdotE3(:,1) JdotE3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([JdotE4(:,1) JdotE4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([JdotE1(:,1) 0*JdotE1(:,2)], 'k--', 'Linewidth',0.75); hold off;
set(gca,'ylim',[-0.1 0.1]); 
grid off;
ylabel({'E\cdotJ';'[nW/m^3]'},'fontsize',5);
 grid off;
 set(gca,'ylim',[-0.1 0.1]); 
 set(gca,'ytick',[-0.05 0 0.05]);  
 
  % plot  Ee
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
  hca = irf_panel('e DEF omni 64');
specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
specrec_p_e.p=energy_e.data;%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
irf_spectrogram(h(i_subplot),specrec_p_e)
% ylabel({'Ee';'[ev]'},'fontsize',12);
hca.YLabel.String = {'Ee','[ev]'};   
set(h(i_subplot),'yscale','log');
set(h(i_subplot),'ytick',[1e1 1e2 1e3 1e4]);
set(h(i_subplot),'Ylim',[1e1 3e4]);
caxis(h(i_subplot),[4 8])

 
 set(h(1:(i_subplot-1)),'fontsize',5);
  irf_adjust_panel_position;
  irf_plot_axis_align(h)
  irf_zoom(tint2,'x',h(1:(i_subplot-1)));
  set(gcf,'color','w');
set(gcf,'render','painters');

figname=['EdotJ1'];
print(gcf, '-dpdf', [figname '.pdf']);