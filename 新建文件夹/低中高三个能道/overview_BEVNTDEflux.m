clear;
clc;
mms.db_init('local_file_db','D:\Matlab\xy-matlab\MMS\mms_db\data');
Units=irf_units;
me=Units.me;
% tintStr = '14-Nov-2015';
tint=irf.tint('2017-06-11T01:59:25Z/2017-06-11T01:59:55Z');
Tsta='2017-06-11T01:59:25Z';   
Tend='2017-06-11T01:59:55Z';
% tint=irf.tint('2017-08-12T05:23:55Z/2017-08-12T05:24:30Z');
% Tsta='2017-08-12T05:23:55Z';   
% Tend='2017-08-12T05:24:30Z';

for ic=1:1
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
% c_eval(['B?_gse=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
c_eval('E?_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?=irf.ts2mat(E?_ts);'],ic);
c_eval(['E?_gsm=irf_gse2gsm(E?);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);
c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);
c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
c_eval(['Ve?=irf_gse2gsm(Ve?);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);
c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['Vi?=irf_gse2gsm(Vi?);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);
c_eval('energy_low = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic);
c_eval('energy_mid = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic);
c_eval('energy_high = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic);
c_eval('energy_e = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',tint);',ic);
c_eval('energy_i = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',tint);',ic);
end
% irf_minvar_gui(B1);

% l=[-0.48 -0.88 0.00];
% m=[-0.87 0.48 0.11];
% n=[-0.09 0.05 -0.99];
% 9.3,7.5
% ic=1;
% c_eval('fpiFilei = dataobj(''D:\MATLAB\mms_db\data\mms?\fpi\brst\l2\dis-moms\2017\08\12\mms?_fpi_brst_l2_dis-moms_20170812052003_v3.2.0.cdf'');',ic);
% c_eval('energy_i = get_variable(fpiFilei,''mms?_dis_energyspectr_omni_brst'');',ic);
% ic=1;
% c_eval('fpiFilee = dataobj(''D:\MATLAB\mms_db\data\mms?\fpi\brst\l2\des-moms\2017\08\12\mms?_fpi_brst_l2_des-moms_20170812052003_v3.2.0.cdf'');',ic);
% c_eval('energy_low = get_variable(fpiFilee,''mms?_des_pitchangdist_lowen_brst'');',ic);
% c_eval('energy_mid = get_variable(fpiFilee,''mms?_des_pitchangdist_miden_brst'');',ic);
% c_eval('energy_high = get_variable(fpiFilee,''mms?_des_pitchangdist_highen_brst'');',ic);
% c_eval('energy_e = get_variable(fpiFilee,''mms?_des_energyspectr_omni_brst'');',ic);

% c_eval(['Vexb?=irf_cross(E?_resamp,B?);'],ic);
% c_eval(['Vexb?(:,2:4)=Vexb?(:,2:4)./[power(Bt?(:,2),2) power(Bt?(:,2),2) power(Bt?(:,2),2)]*1e-12*1e18;'],ic);
% c_eval(['Vexbt?=[Vexb?(:,1) sqrt(power(Vexb?(:,2),2)+power(Vexb?(:,3),2)+power(Vexb?(:,4),2))];'],ic);
% c_eval(['Energy_exb?=[Vexb?(:,1) 0.5*me*power(Vexbt?(:,2),2)/(1.602176565*1e-19)];'],ic);

% Vexb1=irf_cross(E1_resamp,B1_gse);
% Vexb1(:,2:4)=(Vexb1(:,2:4)./[power(Bt1(:,2),2) power(Bt1(:,2),2) power(Bt1(:,2),2)])*1e-12*1e18;
% Vexbt1=[Vexb1(:,1) sqrt(power(Vexb1(:,2),2)+power(Vexb1(:,3),2)+power(Vexb1(:,4),2))];
% Energy_exb1=[Vexbt1(:,1) 0.5*me*power(Vexbt1(:,2),2)/(1.602176565*1e-19)];
% 
% L=-[-0.562 0.719 -0.410];
% M=-[-0.555 -0.695 -0.458];
% N=[0.614 0.030 -0.789];


% for ic=1:1,
%     c_eval(['B?lmn=irf_newxyz(B?,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['E?lmn=irf_newxyz(E?_gsm,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['Ve?lmn=irf_newxyz(Ve?_gsm,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['Vi?lmn=irf_newxyz(Vi?_gsm,N,M,L);'],ic);
% end
 %% Init figure
n_subplots=12;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 18; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('B_{GSM} [nT]','fontsize',10)
set(h(1),'Ylim',[-20 20], 'ytick',[-10 0 10]);
% pos1=get(h(1),'pos');
 set(h(1),'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
 irf_legend(h(1),{'B_x','B_y','B_z','|B|'},[0.1 0.12])
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
irf_legend(h(1),'a',[0.99 0.98],'color','k','fontsize',12)
%% E plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([E1_gsm(:,1) E1_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([E1_gsm(:,1) E1_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E1_gsm(:,1) E1_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Et1(:,1) Et1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E1_gsm(:,1) E1_gsm(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('E_{GSM} [mV/m]','fontsize',10)
set(h(2),'Ylim',[-50 50], 'ytick',[-45 -30 -15 0 15 30 45]);
% pos1=get(h(1),'pos');
 set(h(2),'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
 irf_legend(h(2),{'E_x','E_y','E_z','|E|'},[0.1 0.12])
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
irf_legend(h(2),'b',[0.99 0.98],'color','k','fontsize',12)
%% Ve plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Ve1(:,1) Ve1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Ve1(:,1) Ve1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Ve1(:,1) Ve1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
irf_plot([Ve1(:,1) Ve1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Ve [km/s]','fontsize',10)
set(h(3),'Ylim',[-3500 3500]);
irf_legend(h(3),'c',[0.99 0.98],'color','k','fontsize',12);
set(h(3),'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(h(3),{'Ve_x','Ve_y','Ve_z','|Ve|'},[0.1 0.12]);
%% Vi plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Vi1(:,1) Vi1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Vi1(:,1) Vi1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Vi1(:,1) Vi1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
irf_plot([Vi1(:,1) Vi1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('Vi [km/s]','fontsize',10)
set(h(4),'Ylim',[-500 500]);
irf_legend(h(4),'d',[0.99 0.98],'color','k','fontsize',12);
set(h(4),'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(h(4),{'Vi_x','Vi_y','Vi_z','|Vi|'},[0.1 0.12]);

%% N plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Ni1(:,1) Ni1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Ne1(:,1) Ne1(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Ne1(:,1) Ne1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
grid off;
ylabel('N [cm^-3]','fontsize',10)
% set(h(3),'Ylim',[-500 500]);
irf_legend(h(5),'e',[0.99 0.98],'color','k','fontsize',12);
set(h(5),'ColorOrder',[[0 0 1];[0 1 0]]);
irf_legend(h(5),{'Ni','Ne'},[0.1 0.12]);
%% Ti plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Ti1(:,1) Ti1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Ti_para1(:,1) Ti_para1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Ti_perp1(:,1) Ti_perp1(:,2)], 'color','r', 'Linewidth',0.75); hold off;
grid off;
ylabel('Ti [eV]','fontsize',10)
% set(h(6),'Ylim',[20 80]);
irf_legend(h(6),'f',[0.99 0.98],'color','k','fontsize',12);
set(h(6),'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(h(6),{'Ti','Tipara','Tiperp'},[0.1 0.12]);
%% Te plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Te1(:,1) Te1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Te_para1(:,1) Te_para1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Te_perp1(:,1) Te_perp1(:,2)], 'color','r', 'Linewidth',0.75); hold off;
grid off;
ylabel('Te [eV]','fontsize',10)
% set(h(6),'Ylim',[20 80]);
irf_legend(h(7),'g',[0.99 0.98],'color','k','fontsize',12);
set(h(7),'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(h(7),{'Te','Tepara','Teperp'},[0.1 0.12]);
%% plot low e pad
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(8),jet)

specrec_p_elow=struct('t',irf_time(energy_low.time,'epochtt>epoch'));
specrec_p_elow.f=transpose([3 9 15 21 27 33 39 45 51 57 63 69 75 81 87 93 ...
    99 105 111 117 123 129 135 141 147 153 159 165 171 177]);%energy levels
% specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
specrec_p_elow.p=energy_low.data;%data matrix
specrec_p_elow.f_label='';
specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
[h(8) hcb8]=irf_spectrogram(h(8),specrec_p_elow);
hold off;
ylabel('PA','fontsize',10)
% set(h(7),'yscale','log');
set(h(8),'ytick',[0 90 180]);
%  caxis(h(7),[7.5 8.5])
irf_legend(h(8),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar8=get(hcb8,'pos');
poscbar8(3)=poscbar8(3)*0.5;
set(hcb8,'pos',poscbar8);
%% plot mid e pad
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(9),jet)

specrec_p_emid=struct('t',irf_time(energy_mid.time,'epochtt>epoch'));
specrec_p_emid.f=transpose([3 9 15 21 27 33 39 45 51 57 63 69 75 81 87 93 ...
    99 105 111 117 123 129 135 141 147 153 159 165 171 177]);%energy levels
specrec_p_emid.p=energy_mid.data;%data matrix
specrec_p_emid.f_label='';
specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
[h(9) hcb9]=irf_spectrogram(h(9),specrec_p_emid);
hold off;
ylabel('PA','fontsize',10)
% set(h(8),'yscale','log');
set(h(9),'ytick',[0 90 180]);
% caxis(h(2),[5.5 9.5])
irf_legend(h(9),'i',[0.99 0.98],'color','w','fontsize',12);
poscbar9=get(hcb9,'pos');
poscbar9(3)=poscbar9(3)*0.5;
set(hcb9,'pos',poscbar9);
%% plot high e pad
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(10),jet)

specrec_p_ehigh=struct('t',irf_time(energy_high.time,'epochtt>epoch'));
specrec_p_ehigh.f=transpose([3 9 15 21 27 33 39 45 51 57 63 69 75 81 87 93 ...
    99 105 111 117 123 129 135 141 147 153 159 165 171 177]);%energy levels
specrec_p_ehigh.p=energy_high.data;%data matrix
specrec_p_ehigh.f_label='';
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
[h(10) hcb10]=irf_spectrogram(h(10),specrec_p_ehigh);
hold off;
ylabel('PA','fontsize',10)
set(h(10),'ytick',[0 90 180]);
irf_legend(h(10),'j',[0.99 0.98],'color','w','fontsize',12);
poscbar10=get(hcb10,'pos');
poscbar10(3)=poscbar10(3)*0.5;
set(hcb10,'pos',poscbar10);
%% plot e energy spectrom
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(11),jet)

specrec_p_e=struct('t',irf_time(energy_e.time,'epochtt>epoch'));
specrec_p_e.f=transpose([6.5200000,8.5400000,11.170000,14.630000,19.150000,...
    25.070000,32.810001,42.950001,56.230000,73.599998,96.339996,126.12000,165.09000,...
    216.11000,282.89001,370.31000,484.73999,634.53998,830.63000,1087.3101,1423.3199,...
    1863.1600,2438.9199,3192.6101,4179.2002,5470.6802,7161.2500,9374.2500,12271.120,16063.200,21027.109,27525]);%energy levels
specrec_p_e.p=energy_e.data;%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
[h(11) hcb11]=irf_spectrogram(h(11),specrec_p_e);
grid off;
ylabel('Ee(ev)','fontsize',10)
set(h(11),'yscale','log');
set(h(11),'ytick',[1e1 1e2 1e3 1e4]);
% set(h(10),'Ylim',[1e3 3e4]);
% caxis(h(10),[3.5 7])
irf_legend(h(11),'k',[0.99 0.98],'color','w','fontsize',12);
poscbar11=get(hcb11,'pos');
poscbar11(3)=poscbar11(3)*0.5;
set(hcb11,'pos',poscbar11);
%% plot e energy spectrom
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(12),jet)

specrec_p_i=struct('t',irf_time(energy_i.time,'epochtt>epoch'));
specrec_p_i.f=transpose([2.1600001,3.9100001,7.0700002,10.930000,14.470000,19.160000,25.370001,33.590000,...
    44.480000,58.889999,77.980003,103.24000,136.70000,180.99001,239.63000,317.28000,...
    420.09000,556.21997,736.45001,975.08002,1291.0300,1709.3700,2263.2600,2996.6201,3967.6201,...
    5253.2402,6955.4600,9209.2402,12193.310,16144.310,21375.561,28301.891]);%energy levels
specrec_p_i.p=energy_i.data;%data matrix
specrec_p_i.f_label='';
specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
[h(12) hcb12]=irf_spectrogram(h(12),specrec_p_i);
grid off;
ylabel('Ei(ev)','fontsize',10)
set(h(12),'yscale','log');
set(h(12),'ytick',[1e1 1e2 1e3 1e4]);
% set(h(10),'Ylim',[1e3 3e4]);
% caxis(h(10),[3.5 7])
irf_legend(h(12),'l',[0.99 0.98],'color','w','fontsize',12);
poscbar12=get(hcb12,'pos');
poscbar12(3)=poscbar12(3)*0.5;
set(hcb12,'pos',poscbar12);

%% 

  set(h(1:12),'fontsize',11);
  irf_adjust_panel_position
  irf_zoom(tint,'x',h(1:12));
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:28.80Z')],'g');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.159Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.399Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.589Z')],'b');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.789Z')],'k'); 
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30Z')],'g');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30.209Z')],'k');
% %   add_position(h(3),gseR1), xlabel(h(3),'')
%   irf_zoom(tintlmn,'x',h(4:7))
 %% 
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['overview_BEVNTDEflux' '_' Tsta(1:4) Tsta(6:7) Tsta(9:10) '-' Tsta(12:13) Tsta(15:16) Tsta(18:19) '_' Tend(12:13) Tend(15:16) Tend(18:19)];
print(gcf, '-dpng', [figname '.png']);