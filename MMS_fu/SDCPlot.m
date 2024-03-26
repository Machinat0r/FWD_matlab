function id_flagTime = SDCPlot(tint,desmoms1,desmoms2,IC,Name,flagTime)
% see also SDCFilenames,SDCFilesDownload,SDCDataMove
%% 底下就是原来的overview程序
    %zhuanhuan
    N=[0.54,-0.82,0.17];
    q=[0 0 1];
    L=cross(N,q);
    M=cross(N,L);

    global OutputDir

    for ic=IC
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

    c_eval(['fpiFilee1 = dataobj(','''',desmoms1,'''',');'],ic);
    c_eval('energy_low1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_lowen_brst'');',ic);
    c_eval('energy_mid1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_miden_brst'');',ic);
    c_eval('energy_high1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_highen_brst'');',ic);
    c_eval('energy_e1 = get_variable(fpiFilee1,''mms?_des_energyspectr_omni_brst'');',ic);

    c_eval(['fpiFilee2 = dataobj(','''',desmoms2,'''',');'],ic);
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
    n=12;
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
    ylabel('B [nT]','fontsize',12);
    % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    
    i=i+1;

%% diffB plot
    % h(i)=irf_subplot(n,1,-i);
    % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % B1(:,5)=diff(B1(:,4));
    % 
    % irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,5)], 'color','g', 'Linewidth',0.75); hold on
    % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
    % % irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
    % % irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('B [nT]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % set(gca,'Ylim',[-5 15], 'ytick',[ -5 0 5 10]);
    % pos1=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;

%% Bsmall plot
    % h(i)=irf_subplot(n,1,-i);
    % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
    % irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('B [nT]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % set(gca,'Ylim',[5 10], 'ytick',[  0 5 6 7 8 9 10]);
    % pos1=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;

%% Bbf plot
    % h(i)=irf_subplot(n,1,-i);
    % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([Bbf1(:,1) Bbf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([Bbf1(:,1) Bbf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([Bbf1(:,1) Bbf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
    % % irf_plot([Bbft1(:,1) Bbft1], 'color','k', 'Linewidth',0.75); hold on;
    % % irf_plot([Bbf1(:,1) Bbf1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('Bbf [nT]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % set(gca,'Ylim',[-0.5 0.5], 'ytick',[0 5 6 7 8 9 10]);
    % pos1=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;


%% Blmn plot
    % h(i)=irf_subplot(n,1,-i);
    % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([Blmn1(:,1) Blmn1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([Blmn1(:,1) Blmn1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([Blmn1(:,1) Blmn1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
    % % irf_plot([Bbft1(:,1) Bbft1], 'color','k', 'Linewidth',0.75); hold on;
    % % irf_plot([Bbf1(:,1) Bbf1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('B2 [nT]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % set(gca,'Ylim',[-0.5 0.5], 'ytick',[-0.5 0 0.5]);
    % pos1=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'B_{\perp 1}','B_{\perp 2}','B_{||}'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;


%% dfN plot
    % h(i)=irf_subplot(n,1,-i);
    % irf_plot([Ne1(:,1) Ne1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
    % irf_plot([Ni1(:,1) Ni1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('N [cm^{-3}]','fontsize',10);
    % %set(gca,'Ylim',[fix(min([min(Ne1(:,2)) min(Ni1(:,2))]))-2 fix(max([max(Ne1(:,2)) max(Ni1(:,2))]))+2]);
    % set(gca,'Ylim',[0.3 1], 'ytick',[0.2 0.5 0.8]);
    % % pos1=get(h(1),'pos');
    %  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
    %  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;

%% N plot
    h(i)=irf_subplot(n,1,-i);

    %滤波
%     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
%     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;

    %非滤波
    irf_plot([Ne1(:,1) Ne1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
    irf_plot([Ni1(:,1) Ni1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
    grid off;
    set(gca,'Ylim',[min([min(Ne1(:,2)) min(Ni1(:,2))])-0.1 max([max(Ne1(:,2)) max(Ni1(:,2))])+0.1]);
%     set(gca,'Ylim',[0.15 0.35], 'ytick',[0 0.2 0.3 0.4 0.6 0.8 1 1.2 1.6],'fontsize',9);
    % pos1=get(h(1),'pos');
    %  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
    %  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
      set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
     irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
    % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
    ylabel('N [cm^{-3}]','fontsize',12);
    i=i+1;



%% deElectric field
    % h(i)=irf_subplot(n,1,-i);
    % irf_plot([Ebf1(:,1) Ebf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([Ebf1(:,1) Ebf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([Ebf1(:,1) Ebf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % irf_plot([Ebf1(:,1) Ebf1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('Ebf [mV/m]','fontsize',10)
    % % set(gca,'Ylim',[-1 3], 'ytick',[-0.5 0 1 2  ]);
    % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
    % %set(gca,'Ylim',[fix(min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))])/10)*10-10 fix(max(Et1(:,2))/10)*10+10]);
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'dfE_x','dfE_y','dfE_z'},[0.1 0.12]);
    % pos3=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 1 0]]);
    % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
    % i=i+1;
%% BandE plot
    % h(i)=irf_subplot(n,1,-i);
    % B1(:,5)=(Bbf1(:,2).^2+Bbf1(:,3).^2+Bbf1(:,4).^2).^0.5/9
    % B1(:,6)=(Ebf1(:,2).^2+Ebf1(:,3).^2+Ebf1(:,4).^2).^0.5
    % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,5)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([B1(:,1) B1(:,6)], 'color','g', 'Linewidth',0.75); hold on;
    % 
    % 
    % grid off;
    % ylabel('B [nT]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % % set(gca,'Ylim',[-5 15], 'ytick',[ -5 0 5 10]);
    % pos1=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;
%% Efac
    h(i)=irf_subplot(n,1,-i);
    irf_plot([Efac1(:,1) Efac1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    irf_plot([Efac1(:,1) Efac1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    irf_plot([Efac1(:,1) Efac1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    irf_plot([Efac1(:,1) Efac1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    grid off;
    ylabel('Efac [mV/m]','fontsize',10)
%     set(gca,'Ylim',[-8 8], 'ytick',[-6 -4 -2 0 2 4 6 ]);
    % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
    set(gca,'Ylim',[fix(min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))]))-1 fix(max(Et1(:,2)))+1]);
    set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
    irf_legend(gca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.1 0.12]);
    
    pos3=get(gca,'pos');
    set(gca,'ColorOrder',[[0 1 0]]);
    %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
    i=i+1;



%% E plot
    % h(i)=irf_subplot(n,1,-i);
    % E_fac = irf_convert_fac(E1,B1,[1 0 0]);
    % irf_plot([E_fac(:,1) E_fac(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([E_fac(:,1) E_fac(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([E_fac(:,1) E_fac(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % irf_plot([Et1(:,1) Et1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([E1(:,1) E1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('E [mV/m]','fontsize',10)
    % %set(gca,'Ylim',[-5 5], 'ytick',[-3 0 3]);
    % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
    % set(gca,'Ylim',[fix(min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))])/10)*10-10 fix(max(Et1(:,2))/10)*10+10]);
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.05 0.92]);
    % pos3=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 1 0]]);
    % % irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
    % i=i+1;



%% Ve plot
    h(i)=irf_subplot(n,1,-i);
    irf_plot([Ve1_gsm(:,1) Ve1_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    irf_plot([Ve1_gsm(:,1) Ve1_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    irf_plot([Ve1_gsm(:,1) Ve1_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    
    % irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
    irf_plot([Ve1_gsm(:,1) Ve1_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    
    grid off;
    ylabel('Ve [km/s]','fontsize',10);
    set(gca,'Ylim',[fix(min([min(Ve1_gsm(:,2)) min(Ve1_gsm(:,3)) min(Ve1_gsm(:,4))])/10)*10-10 fix(max(Vet1(:,2))/10)*10+10]);
    % set(gca,'Ylim',[-800 800], 'ytick',[-600 -400 -200 0 200 400 600]);
    %irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
    % irf_legend(gca,{'Ve_N','Ve_M','Ve_L','|Ve|','|Vexb|'},[0.1 0.12]);
    set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.05 0.92]);
    i=i+1;
%% Vexb plot
    % h(i)=irf_subplot(n,1,-i);
    % irf_plot([Vexb1(:,1) Vexb1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([Vexb1(:,1) Vexb1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    % irf_plot([Vexb1(:,1) Vexb1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
    % irf_plot([Vexb1(:,1) Vexb1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('Vexb [km/s]','fontsize',10);
    % set(gca,'Ylim',[fix(min([min(Vexb1(:,2)) min(Vexb1(:,3)) min(Vexb1(:,4))])/10)*10-10 fix(max(Vexb1(:,2))/10)*10+10]);
    % % set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
    % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
    % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
    % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'Vexb_x','Vexb_y','Vexb_z'},[0.05 0.92]);
    % i=i+1;

%% Electric field
% % %     h(i)=irf_subplot(n,1,-i);
% % %     irf_plot([E1(:,1) E1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % %     irf_plot([E1(:,1) E1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % %     irf_plot([E1(:,1) E1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % %     irf_plot([E1(:,1) E1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% % %     grid off;
% % %     % set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% % %     % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% % %     % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % %     set(gca,'Ylim',[min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))])-2 max([max(E1(:,2)) max(E1(:,3)) max(E1(:,4))])+2]);
% % %     set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % %     irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
% % %     pos3=get(gca,'pos');
% % %     set(gca,'ColorOrder',[[0 1 0]]);
% % %     %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % %     ylabel('E [mV/m]','fontsize',12)
% % %     i=i+1;
%% Press
    % h(i)=irf_subplot(n,1,-i);
    % irf_plot([Bt1(:,1) Pm1(:,1)], 'color','b', 'Linewidth',0.75); hold on;
    % irf_plot([Ti1(:,1) Pti1(:,1)], 'color','r', 'Linewidth',0.75); hold on;
    % % irf_plot([E1(:,1) E1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % % irf_plot([E1(:,1) E1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % set(gca,'Ylim',[-2*10^(-3) 2*10^(-3)],'fontsize',9);
    % % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
    % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
    % % set(gca,'Ylim',[fix(min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))])/10)*10-10 fix(max(Et1(:,2))/10)*10+10]);
    % set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
    % irf_legend(gca,{'Pm','Pthe','Ptotal'},[0.97 0.92]);
    % pos3=get(gca,'pos');
    % set(gca,'ColorOrder',[[0 1 0]]);
    % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
    % ylabel('P [N/m]','fontsize',12)
    % i=i+1; 
%% Vi plot
    h(i)=irf_subplot(n,1,-i);
    irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
    irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
    % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
    % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
    irf_plot([Vi1_gsm(:,1) Vi1_gsm(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    grid off;
    set(gca,'Ylim',[fix(min([min(Vi1_gsm(:,2)) min(Vi1_gsm(:,3)) min(Vi1_gsm(:,4))])/10)*10-10 fix(max(Vit1(:,2))/10)*10+10],'fontsize',9);
    %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
    % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
    % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
    % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
    set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
    ylabel('Vi [km/s]','fontsize',12);
    i=i+1;
%% Te plot
    h(i)=irf_subplot(n,1,-i);
    irf_plot([Te_para1(:,1) (Te_para1(:,2)+2*Te_perp1(:,2))/3], 'color','k', 'Linewidth',0.75); hold on;
    irf_plot([Te_para1(:,1) Te_para1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    irf_plot([Te_perp1(:,1) Te_perp1(:,2)], 'color','r', 'Linewidth',0.75); hold off;
    grid off;
    set(gca,'Ylim',[fix(min([min(Te_para1(:,2)) min(Te_para1(:,2)) min(Te_perp1(:,2))])/10)*10-10 ...
        fix(max([max(Te_para1(:,2)) max(Te_para1(:,2)) max(Te_perp1(:,2))])/10)*10+10],'fontsize',9);
    % set(gca,'Ylim',[500 2500]);
    % irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
    set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
    irf_legend(gca,{'Te','T_/_/','T_⊥'},[0.97 0.92]);
    ylabel('Te [eV]','fontsize',12);
    i=i+1;
%% Ti plot
    h(i)=irf_subplot(n,1,-i);
    irf_plot([Ti_para1(:,1) (Ti_para1(:,2)+2*Ti_perp1(:,2))/3], 'color','k', 'Linewidth',0.75); hold on;
    irf_plot([Ti_para1(:,1) Ti_para1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
    irf_plot([Ti_perp1(:,1) Ti_perp1(:,2)], 'color','r', 'Linewidth',0.75); hold off;
    grid off;
    ylabel('Ti [eV]','fontsize',10);
    set(gca,'Ylim',[fix(min([min(Ti_para1(:,2)) min(Ti_para1(:,2)) min(Ti_perp1(:,2))])/10)*10-10 ...
        fix(max([max(Ti_para1(:,2)) max(Ti_para1(:,2)) max(Ti_perp1(:,2))])/10)*10+10]);
    %set(gca,'Ylim',[-100 300]);
    % irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
    set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
    irf_legend(gca,{'Ti','Tipara','Tiperp'},[0.05 0.92]);
    i=i+1;
%% plot low e pad
%     %0-200eV
    h(i)=irf_subplot(n,1,-i);
    % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
    colormap(h(i),jet)
    specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
    specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
    specrec_p_elow.p=energy_low.data;%data matrix
    specrec_p_elow.f_label='';
    specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
    [h(i), hcb6]=irf_spectrogram(h(i),specrec_p_elow);
    ylabel('PA low','fontsize',10)
    % set(gca,'yscale','log');
    set(h(i),'ytick',[0 90 180]);
%     caxis(gca,[6.3 6.6]);
    %irf_legend(h(i),'g',[0.99 0.98],'color','w','fontsize',12);
    poscbar6=get(hcb6,'pos');
    poscbar6(3)=poscbar6(3)*0.5;
    set(hcb6,'pos',poscbar6);
    i=i+1;
%% plot mid e pad
%     %200-2000eV
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
%     caxis(gca,[6.2 6.6]);
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
%     caxis(gca,[7 7.3]);
    %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
    poscbar6=get(hcb6,'pos');
    poscbar6(3)=poscbar6(3)*0.5;
    set(hcb6,'pos',poscbar6);
    i=i+1;

%% plot high e pad2
    % h(i)=irf_subplot(n,1,-i);
    % %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
    % colormap(h(i),jet)
    % 
    % specrec_p_ehigh2=struct('t',irf_time(energy_high2.DEPEND_0.data,'ttns>epoch'));
    % specrec_p_ehigh2.f=transpose(energy_high2.DEPEND_1.data(1,1:30));%energy levels
    % specrec_p_ehigh2.p=energy_high2.data;%data matrix
    % specrec_p_ehigh2.f_label='';
    % specrec_p_ehigh2.p_label={' ','keV/(cm^2 s sr keV)'};
    % [h(i), hcb7]=irf_spectrogram(h(i),specrec_p_ehigh2);
    % ylabel('PA high','fontsize',10)
    % 
    % set(h(i),'ytick',[0 90 180]);
    %  caxis(gca,[6.3 7.1]);
    % %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
    % poscbar7=get(hcb7,'pos');
    % poscbar7(3)=poscbar7(3)*0.5;
    % set(hcb7,'pos',poscbar7);
    % i=i+1;

%% plot high e pad3
    % h(i)=irf_subplot(n,1,-i);
    % %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
    % colormap(h(i),jet)
    % 
    % specrec_p_ehigh3=struct('t',irf_time(energy_high3.DEPEND_0.data,'ttns>epoch'));
    % specrec_p_ehigh3.f=transpose(energy_high3.DEPEND_1.data(1,1:30));%energy levels
    % specrec_p_ehigh3.p=energy_high3.data;%data matrix
    % specrec_p_ehigh3.f_label='';
    % specrec_p_ehigh3.p_label={' ','keV/(cm^2 s sr keV)'};
    % [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_ehigh3);
    % ylabel('PA high','fontsize',10)
    % 
    % set(h(i),'ytick',[0 90 180]);
    %  caxis(gca,[6.3 7.1]);
    % %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
    % poscbar8=get(hcb8,'pos');
    % poscbar8(3)=poscbar8(3)*0.5;
    % set(hcb8,'pos',poscbar8);
    % i=i+1;

%% plot high e pad4
    % h(i)=irf_subplot(n,1,-i);
    % %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
    % colormap(h(i),jet)
    % 
    % specrec_p_ehigh4=struct('t',irf_time(energy_high4.DEPEND_0.data,'ttns>epoch'));
    % specrec_p_ehigh4.f=transpose(energy_high4.DEPEND_1.data(1,1:30));%energy levels
    % specrec_p_ehigh4.p=energy_high4.data;%data matrix
    % specrec_p_ehigh4.f_label='';
    % specrec_p_ehigh4.p_label={' ','keV/(cm^2 s sr keV)'};
    % [h(i), hcb9]=irf_spectrogram(h(i),specrec_p_ehigh4);
    % ylabel('PA high','fontsize',10)
    % 
    % set(h(i),'ytick',[0 90 180]);
    %  caxis(gca,[6.3 7.1]);
    % %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
    % poscbar9=get(hcb9,'pos');
    % poscbar9(3)=poscbar9(3)*0.5;
    % set(hcb9,'pos',poscbar9);
    % i=i+1;

%% pad pingjun plot
    % h(i)=irf_subplot(n,1,-i);
    % % specrec_p_ehigh1.p(:,32)=specrec_p_ehigh1.t(:,1)
    % % average_perp2_data=(specrec_p_ehigh2.p(:,14)+specrec_p_ehigh2.p(:,15)+specrec_p_ehigh2.p(:,16)+specrec_p_ehigh2.p(:,17))/4;
    % % specrec_p_ehigh3.p(:,31)=(specrec_p_ehigh3.p(:,14)+specrec_p_ehigh3.p(:,15)+specrec_p_ehigh3.p(:,16)+specrec_p_ehigh3.p(:,17))/4;
    % % specrec_p_ehigh4.p(:,31)=(specrec_p_ehigh4.p(:,14)+specrec_p_ehigh4.p(:,15)+specrec_p_ehigh4.p(:,16)+specrec_p_ehigh4.p(:,17))/4;
    % average_perp1_data=(specrec_p_ehigh1.p(:,14)+specrec_p_ehigh1.p(:,15)+specrec_p_ehigh1.p(:,16)+specrec_p_ehigh1.p(:,17))/4;
    % average_perp1=[irf_time(energy_high1.DEPEND_0.data,'ttns>epoch') average_perp1_data];
    % % average_perp1=irf_tlim(average_perp1,[iso2epoch('2017-06-11T01:59:34Z') iso2epoch('2017-06-11T01:59:46Z')]);
    % 
    % irf_plot([average_perp1(:,1) average_perp1(:,2)], 'color','k', 'Linewidth',0.75);
    % % irf_plot([specrec_p_ehigh1.t(:,1) specrec_p_ehigh1.p(:,31)], 'color','k', 'Linewidth',0.75); hold off;
    % % plot(specrec_p_ehigh1.p(:,32),specrec_p_ehigh1.p(:,31),'color','k','Linewidth',0.75);hold on;
    % % irf_plot([specrec_p_ehigh2.p(:,32) specrec_p_ehigh2.p(:,31)], 'color','b', 'Linewidth',0.75); hold on;
    % % irf_plot([specrec_p_ehigh3.t(:,1) specrec_p_ehigh3.p(:,31)], 'color','g', 'Linewidth',0.75); hold on;
    % % irf_plot([specrec_p_ehigh3.t(:,1) specrec_p_ehigh4.p(:,31)], 'color','r', 'Linewidth',0.75); hold on;
    % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
    % % irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
    % % irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
    % grid off;
    % ylabel('f [keV/(cm^2 s sr keV]','fontsize',10);
    % % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
    % % set(h(3),'Ylim',[3 8], 'ytick',[5 6]);
    % % set(h(9),'Yscale','log');
    % % pos1=get(gca,'pos');
    % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
    % % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
    % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
    % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
    % i=i+1;
%% plot e energy spectrom
    h(i)=irf_subplot(n,1,-i);
    colormap(h(i),jet)

    specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
    specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
    specrec_p_e.p=energy_e.data;%data matrix
    specrec_p_e.f_label='';
    specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
    [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);
    % hold on;
    % irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
    grid off;
    set(h(i),'yscale','log');
    set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
    ylabel('Ee(ev)','fontsize',12)
    set(gca,'Ylim',[1e1 3e4]);
    caxis(gca,[6.4 7.4])
    % irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
    poscbar8=get(hcb8,'pos');
    poscbar8(3)=poscbar8(3)*0.5;
    %poscbar6(1)=poscbar6(1)*0.5;
    set(hcb8,'pos',poscbar8);
    i=i+1;
%% plot ION energy spectrom
    h(i)=irf_subplot(n,1,-i);
    colormap(h(i),jet)

    specrec_p_i=struct('t',irf_time(energy_i1.DEPEND_0.data,'ttns>epoch'));
    specrec_p_i.f=transpose(energy_i1.DEPEND_1.data(1,1:32));%energy levels
    specrec_p_i.p=energy_i1.data;%data matrix
    specrec_p_i.f_label='';
    specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
    [h(i), hcb7]=irf_spectrogram(h(i),specrec_p_i);
    % hold on;
    % irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
    grid off;
    set(h(i),'yscale','log');
    set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
    ylabel('Ei(ev)','fontsize',12)
    set(gca,'Ylim',[1e1 3e4]);
    caxis(gca,[5 6])
    % irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
    poscbar7=get(hcb7,'pos');
    poscbar7(3)=poscbar7(3)*0.5;
    set(hcb7,'pos',poscbar7);
    i=i+1;
%% plot e energy spectrom
    % h(i)=irf_subplot(n,1,-i);
    % colormap(h(i),jet)
    % 
    % specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
    % specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
    % specrec_p_e.p=energy_e.data;%data matrix
    % specrec_p_e.f_label='';
    % specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
    % [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);
    % % hold on;
    % % irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
    % grid off;
    % set(h(i),'yscale','log');
    % set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
    % ylabel('Ee(ev)','fontsize',12)
    % set(gca,'Ylim',[1e1 3e4]);
    % caxis(gca,[6.4 7.4])
    % % irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
    % poscbar8=get(hcb8,'pos');
    % poscbar8(3)=poscbar8(3)*0.5;
    % %poscbar6(1)=poscbar6(1)*0.5;
    % set(hcb8,'pos',poscbar8);
    % i=i+1;

%% plot waves
    % c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
    % c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
    % c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
    % % Bscm=Bscm{1};            %Bscm??cell
    % c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
    % magB = Bxyz.abs;
    % 
    % %gse2gsm
    % c_eval(['Egse=irf.ts2mat(Exyz);'],ic);
    % c_eval(['Egsm=irf_gse2gsm(Egse);'],ic);
    % Exyz.data=Egsm(:,2:4);
    % 
    % % % c_eval(['Bscmgse=irf.ts2mat(Bscm);'],ic);
    % % % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
    % % % Bscm.data=Bscmgsm(:,2:4);
    % % % 
    % % % % Rotate E and B into field-aligned coordinates
    % % % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
    % % % Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);
    % % % % Bandpass filter E and B waveforms
    % % % dfE = 1/median(diff(Exyz.time.epochUnix));
    % % % dfB = 1/median(diff(Bscm.time.epochUnix));
    % % % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
    % % % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
    % % % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
    % % % 当Bscm发生bug时其会变为{1,2}的cell，点进去发现两部分是一样的，有时候重启matlab会好使有时候不好使就用下面这部分（到wave transforms之前）
    % c_eval(['Bscmgse=irf.ts2mat(Bscm{1,1});'],ic);
    % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
    % Bscm{1,1}.data=Bscmgsm(:,2:4);
    % 
    % % Rotate E and B into field-aligned coordinates
    % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
    % Bscmfac = irf_convert_fac(Bscm{1,1},Bxyz,[1 0 0]);
    % % Bandpass filter E and B waveforms
    % dfE = 1/median(diff(Exyz.time.epochUnix));
    % dfB = 1/median(diff(Bscm{1,1}.time.epochUnix));
    % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
    % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
    % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
    % 
    % 
    % % Wavelet transforms
    % nf = 100;
    % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 4000]);
    % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 50000]);
    % Bwavelet = irf_wavelet(Bscmfac,'nf',nf,'f',[5 4000]);
    % 
    % %compress wavelet transform data 10 point average
    % nc = 20;
    % idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
    % Ewavelettimes = Ewavelet.t(idx);
    % Ewaveletx = zeros(length(idx),nf);
    % Ewavelety = zeros(length(idx),nf);
    % Ewaveletz = zeros(length(idx),nf);
    % for ii = [1:length(idx)];
    %         Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    %         Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    %         Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    % end
    % specperpE=struct('t',Ewavelettimes);
    % specperpE.f=Ewavelet.f;
    % specperpE.p=Ewaveletx+Ewavelety;
    % specperpE.f_label='';
    % specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};
    % 
    % specparE=struct('t',Ewavelettimes);
    % specparE.f=Ewavelet.f;
    % specparE.p=Ewaveletz;
    % specparE.f_label='';
    % specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};
    % 
    % specE=struct('t',Ewavelettimes);
    % specE.f=Ewavelet.f;
    % specE.p=Ewaveletx+Ewavelety+Ewaveletz;
    % specE.f_label='';
    % specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};
    % 
    % 
    % idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
    % Bwavelettimes = Bwavelet.t(idx);
    % Bwaveletx = zeros(length(idx),nf);
    % Bwavelety = zeros(length(idx),nf);
    % Bwaveletz = zeros(length(idx),nf);
    % for ii = [1:length(idx)];
    %         Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    %         Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    %         Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
    % end
    % specB=struct('t',Bwavelettimes);
    % specB.f=Bwavelet.f;
    % specB.p=Bwaveletx+Bwavelety+Bwaveletz;
    % specB.f_label='';
    % specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
    % 
    % 
    % % Compute characteristic frequencies
    % Units=irf_units; % read in standard units
    % Me=Units.me;
    % Mp=Units.mp;
    % e=Units.e;
    % epso=Units.eps0;
    % mu0=Units.mu0;
    % Mp_Me = Mp/Me;
    % B_SI=magB.data*1e-9;
    % Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
    % Wce = e*B_SI/Me;
    % Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
    % Fce = Wce/2/pi;
    % Fce01=Fce*0.1;
    % Fce05=Fce*0.5;
    % Fpe = Wpe/2/pi;
    % Fcp = Fce/Mp_Me;
    % Fpp = Wpp/2/pi;
    % Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
    % Fpe = irf.ts_scalar(magB.time,Fpe);
    % Fce = irf.ts_scalar(magB.time,Fce);
    % Flh = irf.ts_scalar(magB.time,Flh);
    % Fpp = irf.ts_scalar(magB.time,Fpp);
    % Fce01=irf.ts_scalar(magB.time,Fce01);
    % Fce05=irf.ts_scalar(magB.time,Fce05);
    % 
    % h(i)=irf_subplot(n,1,-i);
    % colormap(h(i),jet)
    % [a8,b8]=irf_spectrogram(h(i),specE,'log');
    % 
    % hold(h(i),'on');
    % irf_plot(h(i),Fpe,'color','k','LineWidth',1)
    % irf_plot(h(i),Flh,'color','k','LineWidth',1)
    % irf_plot(h(i),Fce,'color','r','LineWidth',1)
    % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
    % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
    % hold(h(i),'off');
    % 
    % % irf_legend(h(i),'(h)',[0.99 0.97],'color','w','fontsize',12)
    % %%caxis(h(i),[-10 -2]);
    % set(h(i),'yscale','log');
    % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
    % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
    % ylabel(b8,{'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',10);
    % grid(h(i),'off');
    % poscbar8=get(b8,'pos');
    % poscbar8(3)=poscbar8(3)*0.5;
    % set(b8,'pos',poscbar8);
    % i=i+1;
    % 
    % h(i)=irf_subplot(n,1,-i);
    % colormap(h(i),jet)
    % [a9,b9]=irf_spectrogram(h(i),specB,'log');
    % %[h(i),b9]=irf_spectrogram(h(i),specB,'log');
    % 
    % hold(h(i),'on');
    % irf_plot(h(i),Flh,'color','k','LineWidth',1)
    % irf_plot(h(i),Fce,'color','r','LineWidth',1)
    % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
    % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
    % hold(h(i),'off');
    % 
    % % irf_legend(h(i),'(zhaomingjie)',[0.99 0.97],'color','w','fontsize',12)
    % caxis(h(i),[-15 -2]);
    % set(h(i),'yscale','log');
    % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
    % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
    % ylabel(b9,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',10);
    % grid(h(i),'off');
    % poscbar9=get(b9,'pos');
    % poscbar9(3)=poscbar9(3)*0.5;
    % set(b9,'pos',poscbar9);
    % i=i+1;
    % 
    %   set(h(1:n),'fontsize',8);
    % %   irf_zoom(tint,'x',h(1:n));
    irf_zoom(tint,'x',h(1:n));
    % irf_adjust_panel_position;
    % %   irf_plot_axis_align(h)
    irf_plot_axis_align(h)

    %   irf_pl_mark(h(1:n),[iso2epoch('2017-07-24T12:56:52.00Z')],'g');
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.159Z')],'k');
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.399Z')],'k');
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.589Z')],'b');
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.789Z')],'k'); 
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30Z')],'g');
    %   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30.209Z')],'k');
    %  add_position(gca,gseR1), xlabel(gca,'')
    %  irf_zoom(tintlmn,'x',h(4:7))
    for ii = 1:n  
        tempTime = flagTime(1); %2s之内相近的符合判据的时间只画一条线
        irf_pl_mark(h(ii),tempTime,'k')
        id_flagTime = [1];
        for jj = 2:length(flagTime)
            if flagTime(jj) - tempTime >= 2
                tempTime = flagTime(jj);
                irf_pl_mark(h(ii),tempTime,'k');
                id_flagTime(end+1) = jj;
            end
        end
    end
%%  出图保存部分
    set(gcf,'render','painters');
    set(gcf,'paperpositionmode','auto')
    % figname = [OutputDir,'OverviewFig\',Name(2:end-2)];    
    figname = [OutputDir,Name];  
    colormap(jet)
    print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
    close all
end