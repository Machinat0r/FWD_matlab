clear all



d=dataobj('D:\Matlab\xy-matlab\MMS\mms_db\data\mms1\fpi\brst\l2\des-dist\2017\06\11\mms1_fpi_brst_l2_des-dist_20170611015853_v3.2.0.cdf');
b=dataobj('D:\Matlab\xy-matlab\MMS\mms_db\data\mms1\fgm\brst\l2\2017\06\11\mms1_fgm_brst_l2_20170611015853_v5.92.0.cdf');
tint = irf.tint('2017-06-11T01:59:34Z/2017-06-11T01:59:45Z');



%%%%%%%%%%%%%%%%%%%%%   electron  %%%%%%%%%%%%%%%%%%%% 
diste = get_ts(d,'mms1_des_dist_brst');
energye0=get_variable(d,'mms1_des_energy0_brst');
energye1=get_variable(d,'mms1_des_energy1_brst');
phie=get_ts(d,'mms1_des_phi_brst');
thetae=get_variable(d,'mms1_des_theta_brst');
stepTablee=get_ts(d,'mms1_des_steptable_parity_brst');
Bxyz=get_ts(b,'mms1_fgm_b_gsm_brst_l2');


Bxyz1=TSeries(Bxyz.time,Bxyz.data(:,[1 2 3]),'TensorOrder',1,'TensorBasis','xyz');
% %   Bxyz1=get_ts(b,'mms1_fgm_b_gse_brst_l2');
%  

for ic=1:1
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
c_eval(['B?_gse=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
end

mirrorB=39;
for ib=1:length(Bt1(:,1))
    if Bt1(ib,2)==max(Bt1(:,2))
       mirrorB=35; 
    end
    if Bt1(ib,2)<=mirrorB
        alpha(ib,1:3)=[Bt1(ib,1) asind(Bt1(ib,2)/mirrorB) 180-asind(Bt1(ib,2)/mirrorB)];
    else
        alpha(ib,1:3)=[Bt1(ib,1) nan nan];
    end
end
  

energye0 = [10.934900,14.022900,17.982901,23.061100,29.573500,37.924801,48.634499,62.368599,79.981201,102.56700,131.53200,168.67599,216.30800,277.39301,355.72699,456.18201,585.00500,750.20599,962.06000,1233.7400,1582.1400,2028.9301,2601.8799,3336.6399,4278.8901,5487.2202,7036.7700,9023.9199,11572.200,14840.100,19030.900,24405.100];
energye1 = [12.383000,15.879900,20.364300,26.115101,33.489799,42.947102,55.075100,70.627998,90.572899,116.15000,148.95000,191.01300,244.95399,314.12701,402.83499,516.59302,662.47601,849.55499,1089.4600,1397.1200,1791.6600,2297.6101,2946.4500,3778.5000,4845.5298,6213.8799,7968.6401,10218.900,13104.700,16805.400,21551.100,27637];
thetae = [5.6250000,16.875000,28.125000,39.375000,50.625000,61.875000,73.125000,84.375000,95.625000,106.87500,118.12500,129.37500,140.62500,151.87500,163.12500,174.37500];


%% Compute PADS
[paddiste,thetae,energye,tinte] = get_pitchangledist_huang(diste,phie,thetae,stepTablee,energye0,energye1,Bxyz1,tint);
% [paddisti,thetai,energyi,tinti] = mms.get_pitchangledist(disti,phii,thetai,stepTablei,energyi0,energyi1,Bxyz,tint);

% Convert to DEflux
paddiste = paddiste/1e6/(5.486e-4)^2/0.53707;
% paddisti = paddisti/1e6/0.53707;

for ii = 1:length(paddiste.time),
    energytemp = energye(ii,:)'*ones(1,length(thetae));
    paddiste.data(ii,:,:) = squeeze(paddiste.data(ii,:,:)).*energytemp.^2;
end

% for ii = 1:length(paddisti.time),
%     energytemp = energyi(ii,:)'*ones(1,length(thetae));
%     paddisti.data(ii,:,:) = squeeze(paddisti.data(ii,:,:)).*energytemp.^2;
% end

%%

 %%%%%%%%  能量段 （eV）可以随意设置，如 E1 = [100 1000]，就会把100-1000内的数据都加起来画。

E1 = [10 13.305175];  
E2 = [13.305175,17.06250025];
E3 = [17.06250025,21.8808505];
E4 = [21.8808505,28.059874999999998];
E5 = [28.059874999999998,35.9838005];
E6 = [35.9838005,46.1453755];



% E7 = [1100,1500]
% E8 = [1500,1900]
% E9 = [1900,3000]
% E10 = [3000,4000]
% E11 = [4000,6000]
% E12 = [6000,8000]
% E13 = [8000,10000]
% E13 = [1700,7000]

E7 = [1500,1900]
E8 = [1900,2400]
E9 = [2400,3100]
E10 = [3100,4000]
E11 = [4000,5200]
E12 = [5200,6600]
E13 = [1700,7000]

% E21 = [1501.165,1925.08505]
% E22 = [1925.08505,2468.717525]
% E23 = [2468.717525,3165.86745];
% E24 = [3165.86745,4059.88995];
% 
% E25=[4059.88995,5206.38]
% E26=[5206.38,6676.62755]
% E27=[6676.62755,8562.0575]
% 
% E28=[8562.0575,10979.929975]
% E29=[10979.929975,14080.6]
% E30=[14080.6,18056.875]
% E31=[18056.875,26021.05]
% E32=[26021.05,30000]
idx1 = find(energye0 > E1(1) & energye0 < E1(2));
idx2 = find(energye0 > E2(1) & energye0 < E2(2));
idx3= find(energye0 > E3(1) & energye0 < E3(2));
idx4= find(energye0 > E4(1) & energye0 < E4(2));
idx5= find(energye0 > E5(1) & energye0 < E5(2));
idx6= find(energye0 > E6(1) & energye0 < E6(2));

idx7 = find(energye0 > E7(1) & energye0 < E7(2));
idx8 = find(energye0 > E8(1) & energye0 < E8(2));
idx9= find(energye0 > E9(1) & energye0 < E9(2));
idx10= find(energye0 > E10(1) & energye0 < E10(2));
idx11= find(energye0 > E11(1) & energye0 < E11(2));
idx12= find(energye0 > E12(1) & energye0 < E12(2));

idx13 = find(energye0 > E13(1) & energye0 < E13(2));



paddiste1 = squeeze(mean(paddiste.data(:,idx1,:),2));
paddiste2 = squeeze(mean(paddiste.data(:,idx2,:),2));
paddiste3 = squeeze(mean(paddiste.data(:,idx3,:),2));
paddiste4 = squeeze(mean(paddiste.data(:,idx4,:),2));
paddiste5 = squeeze(mean(paddiste.data(:,idx5,:),2));
paddiste6 = squeeze(mean(paddiste.data(:,idx6,:),2));

paddiste7 = squeeze(mean(paddiste.data(:,idx7,:),2));
paddiste8 = squeeze(mean(paddiste.data(:,idx8,:),2));
paddiste9 = squeeze(mean(paddiste.data(:,idx9,:),2));
paddiste10 = squeeze(mean(paddiste.data(:,idx10,:),2));
paddiste11 = squeeze(mean(paddiste.data(:,idx11,:),2));
paddiste12 = squeeze(mean(paddiste.data(:,idx12,:),2));

paddiste13 = squeeze(mean(paddiste.data(:,idx13,:),2));
% paddiste14 = squeeze(mean(paddiste.data(:,idx14,:),2));


a=cell(1,32);
a{1} = paddiste1();
a{2} = paddiste2();
a{3} = paddiste3();
a{4} = paddiste4();
a{5} = paddiste5();
a{6} = paddiste6();
a{7} = paddiste7();
a{8} = paddiste8();
a{9} = paddiste9();
a{10} = paddiste10();
a{11} = paddiste11();
a{12} = paddiste12();
a{13} = paddiste13();


length=length(paddiste1);

for k=1:13
    
for i=1:length;

    if i>=length
      for j=1:12
        aa{1,k}(i,j)=(a{1,k}(i,j)+a{1,k}(i,j))/2;      
      end
    else 
        for j=1:12
    aa{1,k}(i,j)=(a{1,k}(i,j)+a{1,k}(i+1,j))/2;
        end
    end
end

end






%%%%   组成结构体

specepad1=struct('t',paddiste.time.epochUnix);
specepad1.p = double(aa{1,1});
specepad1.p_label={'','keV/(cm^2 s sr keV)'};

specepad1.f_label={''};
specepad1.f = single(thetae);

specepad2=struct('t',paddiste.time.epochUnix);
specepad2.p = double(aa{1,2});
specepad2.f_label={''};
specepad2.f = single(thetae);
specepad2.p_label={'','keV/(cm^2 s sr keV)'};

specepad3=struct('t',paddiste.time.epochUnix);
specepad3.p = double(aa{1,3});
specepad3.f_label={''};
specepad3.f = single(thetae);
specepad2.p_label={'','keV/(cm^2 s sr keV)'};

specepad4=struct('t',paddiste.time.epochUnix);
specepad4.p = double(aa{1,4});
specepad4.f_label={''};
specepad4.f = single(thetae);
specepad4.p_label={'','keV/(cm^2 s sr keV)'};

specepad5=struct('t',paddiste.time.epochUnix);
specepad5.p = double(aa{1,5});
specepad5.f_label={''};
specepad5.f = single(thetae);
specepad5.p_label={'','keV/(cm^2 s sr keV)'};

specepad6=struct('t',paddiste.time.epochUnix);
specepad6.p = double(aa{1,6});
specepad6.f_label={''};
specepad6.f = single(thetae);
specepad6.p_label={'','keV/(cm^2 s sr keV)'};

specepad7=struct('t',paddiste.time.epochUnix);
specepad7.p = double(aa{1,7});
specepad7.p_label={' '};
specepad7.f_label={''};
specepad7.f = single(thetae);
specepad7.p_label={'','keV/(cm^2 s sr keV)'};

specepad8=struct('t',paddiste.time.epochUnix);
specepad8.p = double(aa{1,8});
specepad8.f_label={''};
specepad8.f = single(thetae);
specepad8.p_label={'','keV/(cm^2 s sr keV)'};

specepad9=struct('t',paddiste.time.epochUnix);
specepad9.p = double(aa{1,9});
specepad9.f_label={''};
specepad9.f = single(thetae);
specepad9.p_label={'','keV/(cm^2 s sr keV)'};

specepad10=struct('t',paddiste.time.epochUnix);
specepad10.p = double(aa{1,10});
specepad10.f_label={''};
specepad10.f = single(thetae);
specepad10.p_label={'','keV/(cm^2 s sr keV)'};

specepad11=struct('t',paddiste.time.epochUnix);
specepad11.p = double(aa{1,11});
specepad11.f_label={''};
specepad11.f = single(thetae);
specepad11.p_label={'','keV/(cm^2 s sr keV)'};

specepad12=struct('t',paddiste.time.epochUnix);
specepad12.p = double(aa{1,12});
specepad12.f_label={''};
specepad12.f = single(thetae);
specepad12.p_label={'','keV/(cm^2 s sr keV)'};

specepad13=struct('t',paddiste.time.epochUnix);
specepad13.p = double(aa{1,13});
specepad13.f_label={''};
specepad13.f = single(thetae);
specepad13.p_label={'','keV/(cm^2 s sr keV)'};


%% 
 %% Init figure
n_subplots=7;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 18; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])



% h=irf_plot(7,'newfigure');
% %h=irf_figure(540+ic,8);
% xSize=900; ySize=2400;
% set(gcf,'Position',[10 10 xSize ySize]);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(1)=irf_panel('Bxyz');
% set(h(1),'pos',[0.17 0.94 0.73 0.057]);
irf_plot(h(1),Bxyz1.abs,'linewidth',2);
set(h(1),'fontsize',13);
% ylabel(h(1),'B(nT)','Interpreter','tex','fontsize',17);
% irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'(1)',[0.99 0.98],'color','k','fontsize',13)
% irf_legend(h(1),'B(nT)',[-0.04 0.6],'color','w','fontsize',13)
ylabel(h(1),{'B[nT]'},'fontsize',12,'Interpreter','tex');

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(2)=irf_panel('epad1');
% set(h(2),'pos',[0.17 0.875 0.73 0.057]);
% irf_spectrogram(h(2),specepad1,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(2),'(2)',[0.99 0.98],'color','k','fontsize',13);
set(h(2),'yscale','lin');
set(h(2),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(2),{'12eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(2),'55eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(3)=irf_panel('epad2');
% set(h(3),'pos',[0.17 0.81 0.73 0.057]);
irf_spectrogram(h(3),specepad2,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(3),'(3)',[0.99 0.98],'color','k','fontsize',13);
set(h(3),'yscale','lin');
set(h(3),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(3),{'16eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(3),'66eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(3),'fontsize',13);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(4)=irf_panel('epad3');
% set(h(4),'pos',[0.17 0.745 0.73 0.057]);
irf_spectrogram(h(4),specepad3,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(4),'(4)',[0.99 0.98],'color','k','fontsize',13);
set(h(4),'yscale','lin');
set(h(4),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(4),{'20eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(4),'85eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(4),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(5)=irf_panel('epad4');
% set(h(5),'pos',[0.17 0.68 0.73 0.057]);
irf_spectrogram(h(5),specepad4,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(5),'(5)',[0.99 0.98],'color','k','fontsize',13);
set(h(5),'yscale','lin');
set(h(5),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(5),{'26eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(5),'109eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(5),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(6)=irf_panel('epad5');
% set(h(6),'pos',[0.17 0.615 0.73 0.057]);
irf_spectrogram(h(6),specepad5,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(6),'(6)',[0.99 0.98],'color','k','fontsize',13);
set(h(6),'yscale','lin');
set(h(6),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(6),{'33eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(6),'140eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(6),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(7)=irf_panel('epad6');
% set(h(7),'pos',[0.17 0.55 0.73 0.057]);
irf_spectrogram(h(7),specepad6,'log','donotfitcolorbarlabel');

hold on;
irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(7),'(7)',[0.99 0.98],'color','k','fontsize',13);
set(h(7),'yscale','lin');
set(h(7),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(7),{'1300eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(7),'179eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(7),'fontsize',13);

colormap(h(2),'jet');   %%%%%  调整颜色
colormap(h(3),'jet');
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(7),'jet');

irf_adjust_panel_position
irf_zoom(tint,'x',h(1:7))

set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['ePAD_step_1'];
print(gcf, '-dpdf', [figname '.pdf']);
%% 
%% Init figure
n_subplots=8;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 18; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


% h=irf_plot(7,'newfigure');
% %h=irf_figure(540+ic,8);
% xSize=900; ySize=2400;
% set(gcf,'Position',[10 10 xSize ySize]);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(1)=irf_panel('Bxyz');
% set(h(1),'pos',[0.17 0.94 0.73 0.057]);
irf_plot(h(1),Bxyz1.abs,'linewidth',2);
set(h(1),'fontsize',13);
% ylabel(h(1),'B(nT)','Interpreter','tex','fontsize',17);
% irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
irf_legend(h(1),'a',[0.99 0.98],'color','k','fontsize',13)
% irf_legend(h(1),'B(nT)',[-0.04 0.6],'color','w','fontsize',13)
ylabel(h(1),{'B[nT]'},'fontsize',12,'Interpreter','tex');

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(2)=irf_panel('epad7');
% set(h(2),'pos',[0.17 0.875 0.73 0.057]);
% [h(2),hcb]=irf_spectrogram(h(2),specepad7,'log','donotfitcolorbarlabel');

average_perp7(:,2)=(specepad7.p(:,6)+specepad7.p(:,7))/2;
smoothpad7(:,1)=smooth(average_perp7(:,2),5);
average_perp7(:,1)=specepad7.t(:,1)

irf_plot([average_perp7(:,1) smoothpad7(:,1)], 'color','r', 'Linewidth',0.75);hold on;
% set(hcb,'ytick',[-23 -22.5 -22]);
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};


irf_legend(h(2),'b',[0.99 0.98],'color','k','fontsize',13);
set(h(2),'yscale','lin');
set(h(2),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(2),{'1700eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(2),'55eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(3)=irf_panel('epad8');
% set(h(3),'pos',[0.17 0.81 0.73 0.057]);
% irf_spectrogram(h(3),specepad8,'log','donotfitcolorbarlabel');
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
average_perp8(:,2)=(specepad8.p(:,6)+specepad8.p(:,7))/2;
smoothpad8(:,1)=smooth(average_perp8(:,2),5);
average_perp8(:,1)=specepad8.t(:,1)

irf_plot([average_perp8(:,1) smoothpad8(:,1)], 'color','r', 'Linewidth',0.75);hold on;

irf_legend(h(3),'c',[0.99 0.98],'color','k','fontsize',13);
set(h(3),'yscale','lin');
set(h(3),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(3),{'2200eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(3),'66eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(3),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(4)=irf_panel('epad9');
% set(h(4),'pos',[0.17 0.745 0.73 0.057]);
% irf_spectrogram(h(4),specepad9,'log','donotfitcolorbarlabel');

average_perp9(:,2)=(specepad9.p(:,6)+specepad9.p(:,7))/2;
smoothpad9(:,1)=smooth(average_perp9(:,2),5);
average_perp9(:,1)=specepad9.t(:,1)
irf_plot([average_perp9(:,1) smoothpad9(:,1)], 'color','r', 'Linewidth',0.75);hold on;
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};


irf_legend(h(4),'d',[0.99 0.98],'color','k','fontsize',13);
set(h(4),'yscale','lin');
set(h(4),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(4),{'2700eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(4),'85eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(4),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(5)=irf_panel('epad10');
% set(h(5),'pos',[0.17 0.68 0.73 0.057]);
% irf_spectrogram(h(5),specepad10,'log','donotfitcolorbarlabel');
average_perp10(:,2)=(specepad10.p(:,6)+specepad10.p(:,7))/2;
smoothpad10(:,1)=smooth(average_perp10(:,2),5);
average_perp10(:,1)=specepad10.t(:,1)

irf_plot([average_perp10(:,1) smoothpad10(:,1)], 'color','r', 'Linewidth',0.75);hold on;


irf_legend(h(5),'e',[0.99 0.98],'color','k','fontsize',13);
set(h(5),'yscale','lin');
set(h(5),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(5),{'3600eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(5),'109eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(5),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(6)=irf_panel('epad11');
% set(h(6),'pos',[0.17 0.615 0.73 0.057]);
% irf_spectrogram(h(6),specepad11,'log','donotfitcolorbarlabel');

average_perp11(:,2)=(specepad11.p(:,6)+specepad11.p(:,7))/2;
smoothpad11(:,1)=smooth(average_perp11(:,2),5);
average_perp11(:,1)=specepad11.t(:,1)

irf_plot([average_perp11(:,1) smoothpad11(:,1)], 'color','r', 'Linewidth',0.75);hold on;

irf_legend(h(6),'f',[0.99 0.98],'color','k','fontsize',13);
set(h(6),'yscale','lin');
set(h(6),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(6),{'4600eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(6),'140eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(6),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(7)=irf_panel('epad12');
% set(h(7),'pos',[0.17 0.55 0.73 0.057]);
%  irf_spectrogram(h(7),specepad12,'log','donotfitcolorbarlabel');

average_perp12(:,2)=(specepad12.p(:,6)+specepad12.p(:,7))/2;
smoothpad12(:,1)=smooth(average_perp12(:,2),5);
average_perp12(:,1)=specepad12.t(:,1)
%  irf_plot([average_perp12(:,1) smoothpad12(:,1)], 'color','r', 'Linewidth',0.75);hold on;
irf_plot([average_perp8(:,1) smoothpad8(:,1)], 'color','r', 'Linewidth',0.75);hold on;
irf_plot([average_perp9(:,1) smoothpad9(:,1)], 'color','k', 'Linewidth',0.75);hold on;
irf_plot([average_perp10(:,1) smoothpad10(:,1)], 'color','g', 'Linewidth',0.75);hold on;
irf_plot([average_perp11(:,1) smoothpad11(:,1)], 'color','b', 'Linewidth',0.75);hold on;

 irf_legend(h(7),'g',[0.99 0.98],'color','k','fontsize',13);
set(h(7),'yscale','lin');
set(h(7),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(7),{'5900eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(7),'179eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(7),'fontsize',13);



% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % h(4)=irf_panel('epad9');
% % set(h(4),'pos',[0.17 0.745 0.73 0.057]);
% irf_spectrogram(h(8),specepad13,'log','donotfitcolorbarlabel');
% 
% % hold on;
% % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;
% 
% irf_legend(h(8),'h',[0.99 0.98],'color','k','fontsize',13);
% set(h(8),'yscale','lin');
% set(h(8),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
% ylabel(h(8),{'9000eV'},'fontsize',12,'Interpreter','tex');
% % irf_legend(h(4),'85eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
% set(h(8),'fontsize',13);
% specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};



h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_spectrogram(h(8),specepad13,'log','donotfitcolorbarlabel');



irf_legend(h(8),'h',[0.99 0.98],'color','k','fontsize',13);
set(h(8),'yscale','lin');
set(h(8),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
ylabel(h(8),{'7500eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(4),'85eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(8),'fontsize',13);
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};



average_perp13(:,2)=(specepad13.p(:,6)+specepad13.p(:,7))/2;
smoothpad13(:,1)=smooth(average_perp13(:,2),5);
average_perp13(:,1)=specepad13.t(:,1)

irf_plot([average_perp13(:,1) smoothpad13(:,1)], 'color','r', 'Linewidth',0.75);hold on;
% irf_plot([average_perp9(:,1) smoothpad9(:,1)], 'color','k', 'Linewidth',0.75);hold on;
% irf_plot([average_perp10(:,1) smoothpad10(:,1)], 'color','g', 'Linewidth',0.75);hold on;




colormap(h(2),'jet');   %%%%%  调整颜色
colormap(h(3),'jet');
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(7),'jet');
colormap(h(8),'jet');


irf_adjust_panel_position
irf_zoom(tint,'x',h(1:8))


set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['ePAD_step_2'];
print(gcf, '-dpdf', [figname '.pdf']);
%% 


