function id_flagTime = SDCPlot(tint,tint2,IC,Name,flagTime,ip)
% see also SDCFilenames,SDCFilesDownload,SDCDataMove
%% 底下就是原来的overview程序
    global OutputDir
units = irf_units; 
mu0 = units.mu0; ep0 = units.eps0;
ic=IC;
    % load data
    c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
    c_eval('Bscm?=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
    c_eval('B_brst0?=irf.ts2mat(Bxyz?);',ic);
    try
        c_eval('B_scm?=irf.ts2mat(Bscm?);',ic);
    catch
        c_eval('B_scm?=irf.ts2mat(Bscm?{1});',ic);
    end
    c_eval('B_brst?=irf_resamp(B_brst0?,B_scm?);',ic);
    c_eval('B_fsm?(:,1)=B_brst?(:,1);',ic);
    c_eval('B_fsm?(:,2:4)=B_brst?(:,2:4)+B_scm?(:,2:4);',ic);
    % SCM——高频的磁场;FGM——低频磁场;FGM resample到SCM精度、然后加一起——低频+扰动【FSM】
    c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
    Pos = mms.get_data('R_gse',tint);
    c_eval('Pos? = Pos.gseR?;',ic) %% xyz + total 四列
    c_eval('R?(1,2:4)=Pos?(1,1:3);',ic); %修改
    c_eval('R?(1,1)=B_scm?(1,1);',ic); 
    c_eval('Bt?=Bxyz?.abs;',ic);
    c_eval('Bt_brst?=irf.ts2mat(Bt?);',ic);
    c_eval('Et?=Exyz?.abs;',ic);
    c_eval('Et_brst?=irf.ts2mat(Et?);',ic);
    c_eval('Ni?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
    c_eval('ni0?=irf.ts2mat(Ni?);',ic);
    c_eval('Ne?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
    c_eval('ne0?=irf.ts2mat(Ne?);',ic);
    c_eval('E_brst?=irf.ts2mat(Exyz?);',ic);
    %电子离子速度
    c_eval('vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
    c_eval('ve?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
    c_eval('Ve0?=irf.ts2mat(ve?);',ic);
    c_eval('Vi0?=irf.ts2mat(vi?);',ic);
    % E parallel d_dot_E
    c_eval('Efac?(:,:)=irf_convert_fac(E_brst?(:,1:4),B_brst?(:,1:4),R?(:,1:4));',ic)%z分量为平行分量
 %% 【dE/dT 】&【displacement current】 
c_eval('Efac?=irf_resamp(Efac?,Efac1);',ic)
c_eval('E_brst?=irf_resamp(E_brst?,Efac1);',ic)
n1=length(E_brst1(:,1)); 
for i = 1:n1-1
    c_eval('dE_dt?(i,1)=(E_brst?(i,1)+E_brst?(i+1,1))/2;',ic)
    c_eval('dE_dt?(i,2)=(E_brst?(i+1,2)-E_brst?(i,2))/(E_brst?(i+1,1)-E_brst?(i,1));',ic)%[(mV/m)*s^-1]/1e3=[(V/m)*s^-1]
    c_eval('dE_dt?(i,3)=(E_brst?(i+1,3)-E_brst?(i,3))/(E_brst?(i+1,1)-E_brst?(i,1));',ic)
    c_eval('dE_dt?(i,4)=(E_brst?(i+1,4)-E_brst?(i,4))/(E_brst?(i+1,1)-E_brst?(i,1));',ic)
    c_eval('dE_dt?(i,2:4)=dE_dt?(i,2:4)/1e3;',ic)
    c_eval('j_dis?(i,1)=dE_dt?(i,1);',ic)
    c_eval('j_dis?(i,2:4)=dE_dt?(i,2:4)*mu0*ep0*1e12;',ic)%[T/m]*1e12=[nT/km]
    j_dis_bary(i,1)=j_dis1(i,1);
    j_dis_bary(i,2)=(j_dis1(i,2)+j_dis2(i,2)+j_dis3(i,2)+j_dis4(i,2))/4;
    j_dis_bary(i,3)=(j_dis1(i,3)+j_dis2(i,3)+j_dis3(i,3)+j_dis4(i,3))/4;
    j_dis_bary(i,4)=(j_dis1(i,4)+j_dis2(i,4)+j_dis3(i,4)+j_dis4(i,4))/4;
end
%% 【delta X B】
d_cros_B=c_4_grad('R?','B_fsm?','curl'); % [nT/km]/1e9/1e3=[T/m]
c_eval('B_fsm?=irf_resamp(B_fsm?,B_fsm1);',ic)
%% j=nqv mean
c_eval('Ve?=irf_resamp(Ve0?,B_fsm?);',ic);
c_eval('Vi?=irf_resamp(Vi0?,B_fsm?);',ic);
c_eval('ni?=irf_resamp(ni0?,B_fsm?);',ic);
c_eval('ne?=irf_resamp(ne0?,B_fsm?);',ic);
% c_eval('j_v?(:,1)=ni?(:,1);',ic);
% n2=length(ni1(:,1));
% for i = 1:n2
%     c_eval('j_v?(i,2:4)=(ni?(i,2)*Vi?(i,2:4)-ne?(i,2).*Ve?(i,2:4))*1.6/1e10;',ic)
%     % [cm^(-3)*(km/s)*1.6/1e19]*1e6*1e3=[cm^(-3)*(km/s)*1.6/1e19]*1e9=[m^(-3)*m/s]*1.6/1e10
% end
% c_eval('j_v?(:,2:4)=j_v?(:,2:4)*mu0*1e12;',ic)%[T/m]*1e12=[nT/km]
%% 【rate=位移电流/dxB and 夹角】【4颗卫星】
c_eval('j_dis?=irf_resamp(j_dis?,d_cros_B);',ic)
c_eval('j_dis?(:,5)=vecnorm(j_dis?(:,2:4),2,2);',ic)
d_cros_B(:,5)=vecnorm(d_cros_B(:,2:4),2,2);
c_eval('rate?(:,1)=d_cros_B(:,1);',ic)
c_eval('angle?(:,1)=d_cros_B(:,1);',ic)
n3=length(rate1(:,1));
for i = 1:n3
    c_eval('rate?(i,2)=abs(j_dis?(i,5)/d_cros_B(i,5));',ic)
    c_eval('temp?(i,1)=dot(j_dis?(i,2:4),d_cros_B(i,2:4));',ic)
    c_eval('temp?(i,2)=acos(temp?(i,1)/(j_dis?(i,5)*d_cros_B(i,5)));',ic)
    c_eval('temp?(i,3)=rad2deg(temp?(i,2));',ic)
    c_eval('angle?(i,2)=temp?(i,3);',ic)
end

%% Init figure
   h=irf_plot(6,'newfigure');
xSize=1200; ySize=850;
set(gcf,'Position',[100 100 xSize ySize]);
mmscolors=[0 0 0; 1 0 0; 0 1 0; 0 0 1];%黑红绿蓝
set(h,'ColorOrder',mmscolors)
lnwid = 1;
i=1;

h(i)=irf_panel('Bfsm_GSE_x');
hold(h(i),'on');
irf_plot(h(i),[B_fsm1(:,1) B_fsm1(:,2)],'Linewidth',lnwid);
irf_plot(h(i),[B_fsm2(:,1) B_fsm2(:,2)],'Linewidth',lnwid);
irf_plot(h(i),[B_fsm3(:,1) B_fsm3(:,2)],'Linewidth',lnwid);
irf_plot(h(i),[B_fsm4(:,1) B_fsm4(:,2)],'Linewidth',lnwid);
hold(h(i),'off');
% ylabel(h(1),{'B_{tLMN}','(nT)'},'Interpreter','tex');
ylabel(h(i),{'B_{x}','(nT)'},'Interpreter','tex');
% set(h(i),'ylim',[-40 -35]);
irf_legend(h(i),{'MMS1   ','MMS2   ','MMS3   ','MMS4'},[0.98 0.03],'fontsize',12)
irf_legend(h(i),'(a)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');

i=i+1;
h(i)=irf_panel('Epara');
hold(h(i),'on');
irf_plot(h(i),[Efac1(:,1) Efac1(:,4)],'Linewidth',lnwid);
irf_plot(h(i),[Efac2(:,1) Efac2(:,4)],'Linewidth',lnwid);
irf_plot(h(i),[Efac3(:,1) Efac3(:,4)],'Linewidth',lnwid);
irf_plot(h(i),[Efac4(:,1) Efac4(:,4)],'Linewidth',lnwid);
hold(h(i),'off');
ylabel(h(i),{'E_{//}','(mV/m)'},'Interpreter','tex');
set(h(i),'ylim',[-200 200],'ytick',[-150:150:150]);
irf_legend(h(i),{'MMS1   ','MMS2   ','MMS3   ','MMS4'},[0.98 0.03],'fontsize',12)
irf_legend(h(i),'(b)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');
% 
% i=i+1;
% h(i)=irf_panel('rou');
% hold(h(i),'on');
% irf_plot(h(i),[rou(:,1) rou(:,2)],'Linewidth',0.5,'color','k');
% irf_plot(h(i),[rou(:,1) rou(:,2)*0],'r--', 'Linewidth',0.5);%%虚线
% hold(h(i),'off');
% ylabel(h(i),{'\rho ','(e/cm_{3})'},'Interpreter','tex');
% % set(h(i),'ylim',[-0.005 0.010],'ytick',[-0.005:0.005:0.010]);
% irf_legend(h(i),'(c)',[0.99 0.98],'color','k','fontsize',12);
% grid(h(i),'off');

i=i+1;
h(i)=irf_panel('d_cros_B_MMS1');
hold(h(i),'on');
% irf_plot(h(1),Bt2,'Linewidth',lnwid);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,2)],'Linewidth',lnwid);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,3)],'Linewidth',lnwid);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,4)],'Linewidth',lnwid);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,5)],'Linewidth',lnwid);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,2)*0],'k--', 'Linewidth',0.5);%%虚线
hold(h(i),'off');
ylabel(h(i),{'\nabla X B','(nT/km)'},'Interpreter','tex');
irf_legend(h(i),{'x   ','y   ','z   ','t   '},[0.98 0.03],'fontsize',12)
irf_legend(h(i),'(c)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');

% i=i+1;
% h(i)=irf_panel('x j=nqv');
% hold(h(i),'on');
% irf_plot(h(i),[j_v1(:,1) j_v1(:,2)], 'Linewidth',0.5);
% irf_plot(h(i),[j_v2(:,1) j_v2(:,2)], 'Linewidth',0.5);
% irf_plot(h(i),[j_v3(:,1) j_v3(:,2)], 'Linewidth',0.5);
% irf_plot(h(i),[j_v4(:,1) j_v4(:,2)], 'Linewidth',0.5); %event y
% irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,2)*0],'k--', 'Linewidth',0.5);%%虚线
% hold(h(i),'off');
% ylabel(h(i),{'(u_{0}j)_{x}','(nT/km)'},'Interpreter','tex');
% % set(h(i),'ylim',[0.03 0.11],'ytick',[0.04:0.02:0.1]);
% irf_legend(h(i),{'MMS1   ','MMS2   ','MMS3   ','MMS4'},[0.92 0.98],'fontsize',12)
% % irf_legend(h(i),{'MMS1   '},[0.95 0.98],'color','k','fontsize',12)
% irf_legend(h(i),'(e)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(i),'off');

i=i+1;
h(i)=irf_panel('displacement x');
hold(h(i),'on');
irf_plot(h(i),[j_dis1(:,1) j_dis1(:,2)],'Linewidth',0.5);
irf_plot(h(i),[j_dis2(:,1) j_dis2(:,2)],'Linewidth',0.5);
irf_plot(h(i),[j_dis3(:,1) j_dis3(:,2)],'Linewidth',0.5);
irf_plot(h(i),[j_dis4(:,1) j_dis4(:,2)],'Linewidth',0.5);
hold(h(i),'off');
% ylabel(h(2),{'V_{etLMN}','(km/s)'},'Interpreter','tex');
ylabel(h(i),{'u_{0}\epsilon_{0}(dE/dt)_{x} ','(nT/km)'},'Interpreter','tex');
% set(h(i),'ylim',[-1e-2 1e-2],'ytick',[-0.008:0.008:0.008]);
% irf_legend(h(i),{'MMS1   '},[0.95 0.98],'color','k','fontsize',12)
irf_legend(h(i),'(d)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');

i=i+1;
h(i)=irf_panel('rate');
hold(h(i),'on');
irf_plot(h(i),[rate1(:,1) rate1(:,2)],'Linewidth',0.5);
irf_plot(h(i),[rate2(:,1) rate2(:,2)],'Linewidth',0.5);
irf_plot(h(i),[rate3(:,1) rate3(:,2)],'Linewidth',0.5);
irf_plot(h(i),[rate4(:,1) rate4(:,2)],'Linewidth',0.5);
irf_plot(h(i),[rate1(:,1) rate1(:,2)*0+0.25],'k--', 'Linewidth',0.5);%%虚线
hold(h(i),'off');
ylabel(h(i),{'|rate|','   '},'Interpreter','tex');
set(h(i),'ylim',[0 0.3],'ytick',[0:0.2:0.3]);
% irf_legend(h(i),{'MMS1   ','MMS2   ','MMS3   ','MMS4'},[0.98 0.98],'fontsize',12)
% irf_legend(h(i),{'MMS1   '},[0.95 0.98],'color','k','fontsize',12)
irf_legend(h(i),'(e)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');

i=i+1;
h(i)=irf_panel('angle');
hold(h(i),'on');
irf_plot(h(i),[angle1(:,1) angle1(:,2)],'Linewidth',0.5);
irf_plot(h(i),[angle2(:,1) angle2(:,2)],'Linewidth',0.5);
irf_plot(h(i),[angle3(:,1) angle3(:,2)],'Linewidth',0.5);
irf_plot(h(i),[angle4(:,1) angle4(:,2)],'Linewidth',0.5);
irf_plot(h(i),[d_cros_B(:,1) d_cros_B(:,2)*0],'k--', 'Linewidth',0.5);%%虚线
hold(h(i),'off');
ylabel(h(i),{'angle','(deg)'},'Interpreter','tex');
% irf_legend(h(i),{'MMS1   ','MMS2   ','MMS3   ','MMS4'},[0.98 0.98],'fontsize',12)
% irf_legend(h(i),{'MMS1   '},[0.95 0.98],'color','k','fontsize',12)
set(h(i),'ylim',[0 180],'ytick',[0:90:180]);
irf_legend(h(i),'(f)',[0.99 0.98],'color','k','fontsize',12)
grid(h(i),'off');

 %% 标注电子洞时间
 for ii = 1:6 
        tempTime = flagTime; 
        irf_pl_mark(h(ii),tempTime,'k')  
 end
%%  出图保存部分
    set(gcf,'render','painters');
    set(gcf,'paperpositionmode','auto')
    % figname = [OutputDir,'OverviewFig\',Name(2:end-2)];
    outname = [Name(2:end-2),num2str(ip)];
    figname = [OutputDir,'OverviewFig\',outname];  
    colormap(jet)
    print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
    close all
end