clear;
clc;
mms.db_init('local_file_db','D:\MMS\');
Units=irf_units;
me=Units.me;
tint=irf.tint('2017-08-23T15:38:30.00Z/2017-08-23T15:39:15.00Z');
%
% tsta='2017-06-11T01:59:36.10Z';   
% tend='2017-06-11T01:59:44.10Z';

% tint1=irf.tint('2017-06-11T01:54:00Z/2017-06-11T02:00:00Z');
% R = mms.get_data('R_gsm',tint1);
% c_eval('R?_ts = irf.ts_vec_xyz(R.time,R.gsmR?);');

% for ic=1:4
% c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint1,?);'],ic);
% c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
% c_eval(['R?_ts=R?_ts.resample(B?_ts);'],ic); 
% c_eval(['R?_ts=R?_ts.tlim(tint);'],ic); 
% c_eval(['B?_ts=B?_ts.tlim(tint);'],ic); 
% c_eval(['B?=irf.ts2mat(B?_ts);'],ic);  
% c_eval(['R?=irf.ts2mat(R?_ts);'],ic);
% 
% c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
% c_eval('Bbf? = B?_ts.filt(0.8,1.2,dfB?,5);',ic);
% c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);
% end
% 
% irf_4_v_gui(Bbf1,Bbf2,Bbf3,Bbf4,R1,R2,R3,R4,'mms');

%zhuanhuan
N=[0.54,-0.82,0.17]
q=[0 0 1]
L=cross(N,q)
M=cross(N,L)



for ic=1
% load B


c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
% lvbo
c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,3);',ic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);

% c_eval('Bbff? = B?_ts.filt(0,0.8,dfB?,3);',ic);
% c_eval(['Bbff?=irf.ts2mat(Bbff?);'],ic);

c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);

c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

c_eval(['Vexb?=irf_cross(E?,B?);'],ic);
c_eval(['Vexb?(:,2:4)=1e3*Vexb?(:,2:4)./[Bt?_res(:,2).^2 Bt?_res(:,2).^2 Bt?_res(:,2).^2];'],ic);%km/s


c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
c_eval('Ebf? = E?_gsm.filt(0.8,1.1,dfE?,3);',ic);
c_eval(['Ebf?=irf.ts2mat(Ebf?);'],ic);



% c_eval('Ebff? = E?_gsm.filt(0,0.8,dfE?,3);',ic);
% c_eval(['Ebff?=irf.ts2mat(Ebf?);'],ic);

c_eval('Elmn?=irf_newxyz(Ebf1,L,M,N);',ic);


% load FPI
c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
c_eval('dfNe? = 1/median(diff(Ne?_ts.time.epochUnix));',ic);
c_eval('Nebf? = Ne?_ts.filt(0,1,dfNe?,5);',ic);
c_eval(['Nebf?=irf.ts2mat(Nebf?);'],ic);

c_eval('dfNi? = 1/median(diff(Ni?_ts.time.epochUnix));',ic);
c_eval('Nibf? = Ni?_ts.filt(0,1,dfNi?,5);',ic);
c_eval(['Nibf?=irf.ts2mat(Nibf?);'],ic);




c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);

% c_eval('dfTe_para? = 1/median(diff(Te_para?_ts.time.epochUnix));',ic);
% c_eval('Te_parabf? = Te_para?_ts.filt(0,1.5,dfE?,5);',ic);
% c_eval(['Te_parabf?=irf.ts2mat(Te_parabf?);'],ic);

c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);
c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
c_eval(['Ve?_gsm=irf_gse2gsm(Ve?);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);
c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['Vi?_gsm=irf_gse2gsm(Vi?);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);
ic=1;
% merge data/time from 2 cdf files
c_eval('energy_low1=mms.db_get_variable(''mms1_fpi_brst_l2_des-moms'',''mms1_des_pitchangdist_lowen_brst'',tint);',ic)
c_eval('energy_mid1=mms.db_get_variable(''mms1_fpi_brst_l2_des-moms'',''mms1_des_pitchangdist_miden_brst'',tint);',ic)
c_eval('energy_high1=mms.db_get_variable(''mms1_fpi_brst_l2_des-moms'',''mms1_des_pitchangdist_highen_brst'',tint);',ic)
c_eval('energy_e1=mms.db_get_variable(''mms1_fpi_brst_l2_des-moms'',''mms1_des_energyspectr_omni_brst'',tint);',ic)
c_eval('energy_i1=mms.db_get_variable(''mms1_fpi_brst_l2_dis-moms'',''mms1_dis_energyspectr_omni_brst'',tint);',ic)

%Press
%Pm
miu0 = 4*pi*10^(-7);
c_eval(['Pm? = 10^(-18)*Bt?(:,2).^2 / (2*miu0);'],ic); %10^(-18)Pa
%Pt
c_eval(["Pte? = 10^(-6)*Ne?(:,2).*Te?(:,2);"],ic);
c_eval(["Pti? = 10^(-6)*Ni?(:,2).*Ti?(:,2);"],ic);

end


energy_mid1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_pitchangdist_miden_brst',tint);
energy_high1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_pitchangdist_highen_brst',tint);
energy_e1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint);
energy_i1=mms.db_get_variable('mms1_fpi_brst_l2_dis-moms','mms1_dis_energyspectr_omni_brst',tint);

energy_low2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_lowen_brst',tint);
energy_mid2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_miden_brst',tint);
energy_high2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_highen_brst',tint);
energy_e2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_energyspectr_omni_brst',tint);
energy_i2=mms.db_get_variable('mms2_fpi_brst_l2_dis-moms','mms2_dis_energyspectr_omni_brst',tint);

energy_low3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_lowen_brst',tint);
energy_mid3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_miden_brst',tint);
energy_high3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_highen_brst',tint);
energy_e3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_energyspectr_omni_brst',tint);
energy_i3=mms.db_get_variable('mms3_fpi_brst_l2_dis-moms','mms3_dis_energyspectr_omni_brst',tint);

energy_low4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_lowen_brst',tint);
energy_mid4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_miden_brst',tint);
energy_high4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_highen_brst',tint);
energy_e4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_energyspectr_omni_brst',tint);
energy_i4=mms.db_get_variable('mms4_fpi_brst_l2_dis-moms','mms4_dis_energyspectr_omni_brst',tint);


% c_eval('dfE? = 1/median(diff(Exyz?.time.epochUnix));',ic);
% c_eval('dfB? = 1/median(diff(Bscm?.time.epochUnix));',ic);
% c_eval('Exyzfachf? = Exyzfac?.filt(9,12,dfE?,5);',ic);

c_eval('fpiFilee1 = dataobj(''D:\MMS\mms1\fpi\brst\l2\des-moms\2017\08\23\mms1_fpi_brst_l2_des-moms_20170823153713_v3.3.0.cdf'');',ic);
c_eval('energy_low1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_lowen_brst'');',ic);
c_eval('energy_mid1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_miden_brst'');',ic);
c_eval('energy_high1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_highen_brst'');',ic);
c_eval('energy_e1 = get_variable(fpiFilee1,''mms?_des_energyspectr_omni_brst'');',ic);

c_eval('fpiFilee2 = dataobj(''D:\MMS\mms1\fpi\brst\l2\des-moms\2017\08\23\mms1_fpi_brst_l2_des-moms_20170823153933_v3.3.0.cdf'');',ic);
c_eval('energy_low2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_lowen_brst'');',ic);
c_eval('energy_mid2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_miden_brst'');',ic);
c_eval('energy_high2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_highen_brst'');',ic);
c_eval('energy_e2 = get_variable(fpiFilee2,''mms?_des_energyspectr_omni_brst'');',ic);
% data merge
data1=energy_low1.data; data2=energy_low2.data; data=[data1;data2];energy_low=energy_low1; energy_low.data=data; energy_low.nrec=energy_low1.nrec+energy_low2.nrec;
data1=energy_mid1.data; data2=energy_mid2.data; data=[data1;data2];energy_mid=energy_mid1; energy_mid.data=data;  energy_mid.nrec=energy_mid1.nrec+energy_mid2.nrec;
data1=energy_high1.data; data2=energy_high2.data; data=[data1;data2];energy_high=energy_high1; energy_high.data=data;  energy_high.nrec=energy_high1.nrec+energy_high2.nrec;
data1=energy_e1.data; data2=energy_e2.data; data=[data1;data2];energy_e=energy_e1; energy_e.data=data;  energy_e.nrec=energy_e1.nrec+energy_e2.nrec;
% time merge
data1=energy_low1.DEPEND_0.data;data2=energy_low2.DEPEND_0.data; data=[data1;data2]; energy_low.DEPEND_0.data=data;
data1=energy_mid1.DEPEND_0.data;data2=energy_mid2.DEPEND_0.data; data=[data1;data2]; energy_mid.DEPEND_0.data=data;
data1=energy_high1.DEPEND_0.data;data2=energy_high2.DEPEND_0.data; data=[data1;data2]; energy_high.DEPEND_0.data=data;
data1=energy_e1.DEPEND_0.data;data2=energy_e2.DEPEND_0.data; data=[data1;data2]; energy_e.DEPEND_0.data=data;



% c_eval(['Vexb?=irf_cross(E?_resamp,B?);'],ic);
% c_eval(['Vexb?(:,2:4)=Vexb?(:,2:4)./[power(Bt?(:,2),2) power(Bt?(:,2),2) power(Bt?(:,2),2)]*1e-12*1e18;'],ic);
% c_eval(['Vexbt?=[Vexb?(:,1) sqrt(power(Vexb?(:,2),2)+power(Vexb?(:,3),2)+power(Vexb?(:,4),2))];'],ic);
% c_eval(['Energy_exb?=[Vexb?(:,1) 0.5*me*power(Vexbt?(:,2),2)/(1.602176565*1e-19)];'],ic);
% 
% Vexb1=irf_cross(E1_resamp,B1_gse);
% Vexb1(:,2:4)=(Vexb1(:,2:4)./[power(Bt1(:,2),2) power(Bt1(:,2),2) power(Bt1(:,2),2)])*1e-12*1e18;
% Vexbt1=[Vexb1(:,1) sqrt(power(Vexb1(:,2),2)+power(Vexb1(:,3),2)+power(Vexb1(:,4),2))];
% Energy_exb1=[Vexbt1(:,1) 0.5*me*power(Vexbt1(:,2),2)/(1.602176565*1e-19)];

% irf_minvar_gui(B1);
% L=[0.36 0.16 0.92];
% M=[0.40 -0.91 0.01];
% N=[0.84 0.37 -0.39];

% L=[0 0 1];
% M=[0 1 0];
% N=[1 0 0];
% for ic=1:1,
%     c_eval(['B?lmn=irf_newxyz(B?,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['E?lmn=irf_newxyz(E?,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['Ve?lmn=irf_newxyz(Ve?,N,M,L);'],ic);
% end
% for ic=1:1,
%     c_eval(['Vi?lmn=irf_newxyz(Vi?,N,M,L);'],ic);
% end
%% Init figure
n=4;
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
B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
%irf_plot([Bbf1(:,1) Bbf1(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
set(gca,'Ylim',[-12 22], 'ytick',[-5 0 5 10 15 20 25],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% N plot
h(i)=irf_subplot(n,1,-i);

%滤波
irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;

%非滤波
% irf_plot([Ne1(:,1) Ne1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
% irf_plot([Ni1(:,1) Ni1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
grid off;
%set(gca,'Ylim',[fix(min([min(Ne1(:,2)) min(Ni1(:,2))]))-2 fix(max([max(Ne1(:,2)) max(Ni1(:,2))]))+2]);
set(gca,'Ylim',[0.15 0.35], 'ytick',[0 0.2 0.3 0.4 0.6 0.8 1 1.2 1.6],'fontsize',9);
% pos1=get(h(1),'pos');
%  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
%  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
 irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
ylabel('N [cm^{-3}]','fontsize',12);
i=i+1;
%% plot mid e pad
%200-2000eV
h(i)=irf_subplot(n,1,-i);
%h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i),jet)

specrec_p_emid=struct('t',irf_time(energy_mid.DEPEND_0.data,'ttns>epoch'));
specrec_p_emid.f=transpose(energy_mid.DEPEND_1.data(1,1:30));%energy levels
specrec_p_emid.p=energy_mid.data;%data matrix
specrec_p_emid.f_label='';
specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb7]=irf_spectrogram(h(i),specrec_p_emid);
ylabel('PA mid','fontsize',12)
%set(gca,'yscale','log');
set(h(i),'ytick',[0 90 180]);
caxis(gca,[6.3 6.55]);
%irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar7=get(hcb7,'pos');
poscbar7(3)=poscbar7(3)*0.5;
set(hcb7,'pos',poscbar7);
i=i+1;
%% plot high e pad
%2k-30keV
h(i)=irf_subplot(n,1,-i);
%h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i),jet)

specrec_p_ehigh=struct('t',irf_time(energy_high.DEPEND_0.data,'ttns>epoch'));
specrec_p_ehigh.f=transpose(energy_high.DEPEND_1.data(1,1:30));%energy levels
specrec_p_ehigh.p=energy_high.data;%data matrix
specrec_p_ehigh.f_label='';
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb6]=irf_spectrogram(h(i),specrec_p_ehigh);
ylabel('PA high','fontsize',12)

set(h(i),'ytick',[0 90 180]);
caxis(gca,[6.95 7.25]);
%irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar6=get(hcb6,'pos');
poscbar6(3)=poscbar6(3)*0.5;
set(hcb6,'pos',poscbar6);
i=i+1;
 %% 
 irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
pause(2)
figname=['W1-Figure2'];
print(gcf, '-dpdf', [figname '.pdf']);

