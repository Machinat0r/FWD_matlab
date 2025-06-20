close all
clear,clc

ParentDir = 'D:/MMS/'; 
mms.db_init('local_file_db',ParentDir);

%%%%  按照路径读取好文件，设置好tint的时间  %%%%

% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20151023145734_v2.1.0.cdf');
% b=dataobj('F:\mms\mms1\fgm\brst\l2\2015\10\mms1_fgm_brst_l2_20151023145734_v*.**.*.cdf');
% tint = irf.tint('2015-10-23T14:59:33.20Z/2015-10-23T14:59:35.20Z');

% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20150911110144_v2.1.0.cdf');
% b=dataobj('F:\mms\mms1\fgm\brst\l2\2016\01\mms1_fgm_brst_l2_20150911110144_v*.**.*.cdf');
% tint = irf.tint('2015-10-25T12:47:43.50Z/2015-10-25T12:47:45.00Z');

%%%  3  %%%

% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20150911110144_v2.1.0.cdf');
% b=dataobj('F:\mms\mms1\fgm\brst\l2\2016\01\mms1_fgm_brst_l2_20150911110144_v4.18.1.cdf');
% tint = irf.tint('2015-10-19T14:07:54.00Z/2015-10-19T14:07:59.00Z');

%%%%  1  %%%%
% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20151211235324_v2.1.0.cdf');
% b=dataobj('F:\mms\mms1\fgm\brst\l2\2015\12\mms1_fgm_brst_l2_20151211235324_v4.**.*.cdf');
% tint = irf.tint('2015-12-11T23:54:48.00Z/2015-12-11T23:54:59.00Z');


%%%%  2  %%%%
% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20150925092754_v2.1.0.cdf');
% b=dataobj('D:\MMS_Data\mms1\afg\brst\2016\01\25\mms1_afg_brst_ql_20150925092754_v*.**.*.cdf');
% tint = irf.tint('2015-09-25T09:29:34.00Z/2015-09-25T09:29:36.00Z');
%%%% 

% d=dataobj('F:\mms\mms1\fpi\brst\l2\des-dist\mms1_fpi_brst_l2_des-dist_20151016130524_v2.1.0.cdf');
% b=dataobj('F:\mms\mms1\fgm\brst\l2\2016\01\mms1_fgm_brst_l2_20151016130524_v4.**.*.cdf');
% tint = irf.tint('2015-10-16T13:05:40.00Z/2015-10-16T13:06:10.00Z');
% d=dataobj('D:\MATLAB\mms_db\data\mms1\fpi\fast\l2\des-dist\2015\12\mms1_fpi_fast_l2_des-dist_20151231000000_v2.1.0.cdf');
% energye=get_variable(d,'mms1_des_energy_fast');
% % % 
d=dataobj('D:\MMS\mms1\fpi\brst\l2\des-dist\2019\08\05\mms1_fpi_brst_l2_des-dist_20190805162313_v3.4.0.cdf');
b=dataobj('D:/MMS/mms1/fgm/brst/l2/2019/08/05/mms1_fgm_brst_l2_20190805162313_v5.202.0.cdf');
tint = irf.tint('2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z');

% % % c_eval("d=dataobj('/Volumes/172.17.190.41/Data/MMS/mms1/fpi/brst/l2/des-dist/2018/07/03/mms1_fpi_brst_l2_des-dist_20180703154943_v3.3.0.cdf');");
% % % c_eval("b=dataobj('/Volumes/172.17.190.41/Data/MMS/mms1/fgm/brst/l2/2018/07/03/mms1_fgm_brst_l2_20180703154943_v5.146.0.cdf');");
% % % tint = irf.tint('2018-07-03T15:50:10.00Z/2018-07-03T15:50:25.00Z');

% % % d=dataobj('/Users/fwd/Documents/MATLAB/MMS/mms1/fpi/brst/l2/des-dist/2021/07/14/mms1_fpi_brst_l2_des-dist_20210714180703_v3.4.0.cdf');
% % % b=dataobj('/Users/fwd/Documents/MATLAB/MMS/mms1/fgm/brst/l2/2021/07/14/mms1_fgm_brst_l2_20210714180703_v5.304.0.cdf');
% % % tint = irf.tint('2021-07-14T18:08:00.00Z/2021-07-14T18:08:15.00Z');
% % d=dataobj('D:\MMS\mms1\fpi\brst\l2\des-dist\2021\07\14\mms1_fpi_brst_l2_des-dist_20210714172943_v3.4.0.cdf');
% % b=dataobj('D:\MMS\mms1\fgm\brst\l2\2021\07\14\mms1_fgm_brst_l2_20210714172943_v5.304.0.cdf');
% % tint = irf.tint('2021-07-14T17:30:00.00Z/2021-07-14T17:31:00.00Z');

%%%%%%%%%%%%%%%%%%   ion   %%%%%%%%%%%%%%%%%%%%%%%
% diste = get_ts(d,'mms1_dis_dist_brst');
% energye0=get_variable(d,'mms1_dis_energy0_brst');
% energye1=get_variable(d,'mms1_dis_energy1_brst');
% phie=get_ts(d,'mms1_dis_phi_brst');
% thetae=get_variable(d,'mms1_dis_theta_brst');
% stepTablee=get_ts(d,'mms1_dis_steptable_parity_brst');
% Bxyz=get_ts(b,'mms1_fgm_b_gse_brst_l2');

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
c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
end

% % % mirrorB=39;
% % % for ib=1:length(Bt1(:,1))
% % %     if Bt1(ib,2)==max(Bt1(:,2))
% % %        mirrorB=35; 
% % %     end
% % %     if Bt1(ib,2)<=mirrorB
% % %         alpha(ib,1:3)=[Bt1(ib,1) asind(Bt1(ib,2)/mirrorB) 180-asind(Bt1(ib,2)/mirrorB)];
% % %     else
% % %         alpha(ib,1:3)=[Bt1(ib,1) nan nan];
% % %     end
% % % end
  

energye0 = [10.934900,14.022900,17.982901,23.061100,29.573500,37.924801,48.634499,62.368599,79.981201,102.56700,131.53200,168.67599,216.30800,277.39301,355.72699,456.18201,585.00500,750.20599,962.06000,1233.7400,1582.1400,2028.9301,2601.8799,3336.6399,4278.8901,5487.2202,7036.7700,9023.9199,11572.200,14840.100,19030.900,24405.100];
energye1 = [12.383000,15.879900,20.364300,26.115101,33.489799,42.947102,55.075100,70.627998,90.572899,116.15000,148.95000,191.01300,244.95399,314.12701,402.83499,516.59302,662.47601,849.55499,1089.4600,1397.1200,1791.6600,2297.6101,2946.4500,3778.5000,4845.5298,6213.8799,7968.6401,10218.900,13104.700,16805.400,21551.100,27637];
% % % thetae = [5.6250000,16.875000,28.125000,39.375000,50.625000,61.875000,73.125000,84.375000,95.625000,106.87500,118.12500,129.37500,140.62500,151.87500,163.12500,174.37500];
% energye_fast=[11.636500,14.922600,19.136600,24.540600,31.470800,40.357899,51.754799,66.370003,85.112503,109.14800,139.97000,179.49699,230.18600,295.18900,378.54901,485.44800,622.53601,798.33698,1023.7800,1312.8900,1683.6400,2159.0901,2768.8101,3550.7100,4553.3999,5839.2598,7488.2300,9602.8496,12314.600,15792.200,20251.801,25970.801];
% for iii=1:32
%    energyestep(1,iii)=energye0(1,iii)+energye1(1,iii);  
% end
% for iii=1:32
%     if (iii==1)||(iii==32)
%     energyenew(1,iii)=energyestep(1,iii)/2;
%     else
%         energyenew(1,iii)=(energyestep(1,iii-1)+energyestep(1,iii))/4;
%     end
%     
% end

% for iii=1:32
%     if (iii==1)||(iii==32)
%     energye_fast_new(1,iii)=energye_fast(1,iii);
%     else
%         energye_fast_new(1,iii)=(energye_fast(1,iii-1)+energye_fast(1,iii))/2;
%     end
%     
% end


%% Compute moments
% Units = irf_units; % Use IAU and CODATA values for fundamental constants.
% qe = Units.e;
% mp = Units.mp;
% imoments = mms.psd_moments(disti,phii,thetai,stepTablei,energyi0,energyi1,SCpot,'ion');
% ni = imoments.n_psd;
% Vi = imoments.V_psd;
% Ti = imoments.T_psd;
% Tipp = mms.rotate_tensor(Ti,'fac',Bxyz,'pp'); 
% Tiparperp = TSeries(Ti.time,[Tipp.data(:,1,1) Tipp.data(:,2,2)]);
% 
% emoments = mms.psd_moments(diste,phie,thetae,stepTablee,energye0,energye1,SCpot,'electron');
% ne = emoments.n_psd;
% Ve = emoments.V_psd;
% Te = emoments.T_psd;
% Tepp = mms.rotate_tensor(Te,'fac',Bxyz,'pp'); 
% Teparperp = TSeries(Te.time,[Tepp.data(:,1,1) Tepp.data(:,2,2)]);
% 
diste.data = diste.data*1e30; % Unit conversion
% disti.data = disti.data*1e30;
% 
% energyspec = ones(length(diste.time),1)*energye0;
% for ii = 1:length(diste.time);
%     if stepTablee.data(ii),
%         energyspec(ii,:) = energye1;
%     end
% end
% 
% energyspeci = ones(length(disti.time),1)*energyi0;
% for ii = 1:length(disti.time);
%     if stepTablei.data(ii),
%         energyspeci(ii,:) = energyi1;
%     end
% end
% 
% % define angles
% dangle = pi/16;
% lengthphi = 32;
% 
% z2 = ones(lengthphi,1)*sind(thetae);
% solida = dangle*dangle*z2;
% allsolidi = zeros(size(disti.data));
% allsolide = zeros(size(diste.data));
% 
% for ii = 1:length(disti.time);
%     for jj=1:length(energyi0);
%         allsolidi(ii,jj,:,:) = solida;
%     end
% end
% 
% for ii = 1:length(diste.time);
%     for jj=1:length(energye0);
%         allsolide(ii,jj,:,:) = solida;
%     end
% end
% 
% distis = disti.data.*allsolidi;
% distes = diste.data.*allsolide;
% 
% % Electron analysis - OMNI
% for ii = 1:length(diste.time);
%     disttemp = squeeze(distes(ii,:,:,:));
%     PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
% end
% 
% % Ion analysis - OMNI
% PSDiomni = zeros(length(disti.time),length(energyi0));
% for ii = 1:length(disti.time);
%     disttemp = squeeze(distis(ii,:,:,:));
%     PSDiomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
% end
% 
% efluxomni = PSDomni.*energyspec.^2;
% efluxomni = efluxomni; %convert to normal units
% 
% ifluxomni = PSDiomni.*energyspeci.^2;
% ifluxomni = ifluxomni/1e6/0.53707; %convert to normal units

%% Compute PADS
[paddiste,thetae,energye,tinte] = mms_get_pitchangledist_my_change(diste,phie,thetae,stepTablee,energye0,energye1,Bxyz1,tint);
% % [paddiste,thetae,energye,tinte] = mms.get_pitchangledist(diste,phie,thetae,stepTablee,energye0,energye1,Bxyz1,tint);
% [paddisti,thetai,energyi,tinti] = mms.get_pitchangledist(disti,phii,thetai,stepTablei,energyi0,energyi1,Bxyz,tint);

% Convert to DEflux
paddiste = paddiste/1e6/(5.486e-4)^2/0.53707;
% paddiste = paddiste/(5.486e-4)^2/0.53707;
% paddisti = paddisti/1e6/0.53707;

for ii = 1:length(paddiste.time)
    energytemp = energye(ii,:)'*ones(1,length(thetae));
    paddiste.data(ii,:,:) = squeeze(paddiste.data(ii,:,:)).*energytemp.^2;
end

% for ii = 1:length(paddisti.time),
%     energytemp = energyi(ii,:)'*ones(1,length(thetae));
%     paddisti.data(ii,:,:) = squeeze(paddisti.data(ii,:,:)).*energytemp.^2;
% end

%%

 %%%%%%%%  能量段 （eV）可以随意设置，如 E1 = [100 1000]，就会把100-1000内的数据都加起来画。

% E1 = [30 300];  
% E2 = [300 20000];
% E3 = [20000 30000];
E1 = [30 200];  
E2 = [200 20000];
E3 = [20000 30000];
% E4 = [21.8808505,28.059874999999998];
% E5 = [28.059874999999998,35.9838005];
% E6 = [35.9838005,46.1453755];

% % E1 = [1000,1500];  
% % E2 = [1500,1900];
% % E3 = [3000,5000];
% % E4 = [14000,18000];
% % E5 = [19000,22000];
% % E6 = [35.9838005,46.1453755];


idx1 = find(energye0 > E1(1) & energye0 < E1(2));
idx2 = find(energye0 > E2(1) & energye0 < E2(2));
idx3= find(energye0 > E3(1) & energye0 < E3(2));
% idx4= find(energye0 > E4(1) & energye0 < E4(2));
% idx5= find(energye0 > E5(1) & energye0 < E5(2));
% idx6= find(energye0 > E6(1) & energye0 < E6(2));

paddiste1 = squeeze(mean(paddiste.data(:,idx1,:),2));
paddiste2 = squeeze(mean(paddiste.data(:,idx2,:),2));
paddiste3 = squeeze(mean(paddiste.data(:,idx3,:),2));

a=cell(1,3);
a{1} = paddiste1();
a{2} = paddiste2();
a{3} = paddiste3();
% a{4} = paddiste4();
% a{5} = paddiste5();
% a{6} = paddiste6();


length=length(paddiste1);

for k=1:3
    
for i=1:length

    if i>=length
      for j=1:16
        aa{1,k}(i,j)=(a{1,k}(i,j)+a{1,k}(i,j))/2;      
      end
    else 
        for j=1:16
        aa{1,k}(i,j)=(a{1,k}(i,j)+a{1,k}(i+1,j))/2;
        end
    end
end

end




% paddistilow = squeeze(mean(paddisti.data(:,1:idxlow(end),:),2));
% paddistimid = squeeze(mean(paddisti.data(:,idxmid,:),2));
% paddistihigh = squeeze(mean(paddisti.data(:,idxhigh,:),2));

% Make structures for plotting

% speceomni=struct('t',diste.time.epochUnix);
% speceomni.p = double(efluxomni);
% speceomni.p_label={'e log(dEF)','keV/(cm^2 s sr keV)'};
% speceomni.f_label={''};
% speceomni.f = single(energyspec);

%%%%   组成结构体

specepad1=struct('t',paddiste.time.epochUnix);
specepad1.p = double(aa{1,1});
specepad1.p_label={'keV/(cm^2 s sr keV)'};

specepad1.f_label={''};
specepad1.f = single(thetae);

specepad2=struct('t',paddiste.time.epochUnix);
specepad2.p = double(aa{1,2});
specepad2.f_label={''};
specepad2.f = single(thetae);
specepad2.p_label={'keV/(cm^2 s sr keV)'};

specepad3=struct('t',paddiste.time.epochUnix);
specepad3.p = double(aa{1,3});
specepad3.f_label={''};
specepad3.f = single(thetae);
specepad3.p_label={'keV/(cm^2 s sr keV)'};

% specepad4=struct('t',paddiste.time.epochUnix);
% specepad4.p = double(aa{1,4});
% specepad4.f_label={''};
% specepad4.f = single(thetae);
% specepad4.p_label={'','keV/(cm^2 s sr keV)'};
% 
% specepad5=struct('t',paddiste.time.epochUnix);
% specepad5.p = double(aa{1,5});
% specepad5.f_label={''};
% specepad5.f = single(thetae);
% specepad5.p_label={'','keV/(cm^2 s sr keV)'};
% 
% specepad6=struct('t',paddiste.time.epochUnix);
% specepad6.p = double(aa{1,6});
% specepad6.f_label={''};
% specepad6.f = single(thetae);
% specepad6.p_label={'','keV/(cm^2 s sr keV)'};



% % % % % speciomni=struct('t',disti.time.epochUnix);
% % % % % speciomni.p = double(ifluxomni);
% % % % % speciomni.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
% % % % % speciomni.f_label={''};
% % % % % speciomni.f = single(energyspeci);
% % % % % 
% % % % % specipadl=struct('t',paddisti.time.epochUnix);
% % % % % specipadl.p = double(paddistilow);
% % % % % specipadl.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
% % % % % specipadl.f_label={''};
% % % % % specipadl.f = single(thetai);
% % % % % 
% % % % % specipadm=struct('t',paddisti.time.epochUnix);
% % % % % specipadm.p = double(paddistimid);
% % % % % specipadm.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
% % % % % specipadm.f_label={''};
% % % % % specipadm.f = single(thetai);
% % % % % 
% % % % % specipadh=struct('t',paddisti.time.epochUnix);
% % % % % specipadh.p = double(paddistihigh);
% % % % % specipadh.p_label={'i log(dEF)','keV/(cm^2 s sr keV)'};
% % % % % specipadh.f_label={''};
% % % % % specipadh.f = single(thetai);
%% Init figure
n_subplots=4;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 60; ySize = 40; coef=floor(min(600/xSize,600/ySize));
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
irf_legend(h(1),'(a)',[0.99 0.9],'color','k','fontsize',13)
set(gca,'Ylim',[0,25])
% irf_legend(h(1),'B(nT)',[-0.04 0.6],'color','w','fontsize',13)
ylabel(h(1),{'|B| [nT]'},'fontsize',12,'Interpreter','tex','color','k');
grid off
%%
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(2)=irf_panel('epad1');
% set(h(2),'pos',[0.17 0.875 0.73 0.057]);
    irf_spectrogram(h(2),specepad1,'log','donotfitcolorbarlabel');

% % % hold on;
% % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(2),'(b)',[0.99 0.9],'color','k','fontsize',13);
set(h(2),'yscale','lin');
set(h(2),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(h(2),[6.2 6.7]);
ylabel(h(2),{'PAD(deg)'},'fontsize',13,'Interpreter','tex');
irf_legend(h(2),'0.03-0.6 keV',[0.04 0.85],'fontsize',15,'color','k');
%%
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(3)=irf_panel('epad2');
% set(h(3),'pos',[0.17 0.81 0.73 0.057]);
irf_spectrogram(h(3),specepad2,'log','donotfitcolorbarlabel');

hold on;
% B0 = 12.1863 * ones(7680,1); %16:24:23-16:24:30
% B0 = 11.7607 * ones(7680,1); %16:24:23-16:24:28
B01=4.7;%16:24:05-16:24:15
% B01 = 18;
% B01 = 5;
B02 = 5.7;
% B0=11;
% % % B0=10.7;
alpha1 = real(asind(sqrt(B01./Bt1(:,2))));
irf_plot([Bt1(:,1),alpha1],'w--', 'Linewidth',1); hold on;
irf_plot([Bt1(:,1),180-alpha1],'w--', 'Linewidth',1); hold on;

alpha2 = real(asind(sqrt(B02./Bt1(:,2))));
irf_plot([Bt1(:,1),alpha2],'k--', 'Linewidth',1); hold on;
irf_plot([Bt1(:,1),180-alpha2],'k--', 'Linewidth',1); hold on;

irf_legend(h(3),'(c)',[0.99 0.9],'color','k','fontsize',13);
set(h(3),'yscale','lin');
set(h(3),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(h(3),[6.95 7.4]);
% % % caxis(h(3),[6.825 7.35]);
% caxis(h(3),[7.3 7.48]);
ylabel(h(3),{'PAD(deg)'},'fontsize',13,'Interpreter','tex');
irf_legend(h(3),'0.6-16 keV',[0.04 0.85],'fontsize',15,'color','k');
%set(h(3),'fontsize',13);
%%
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(4)=irf_panel('epad3');
% set(h(4),'pos',[0.17 0.745 0.73 0.057]);
irf_spectrogram(h(4),specepad3,'log','donotfitcolorbarlabel');

% % hold on;
% % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',13);
set(h(4),'yscale','lin');
set(h(4),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(h(4),[5 7.2]);
ylabel(h(4),{'PAD(deg)'},'fontsize',13,'Interpreter','tex');
irf_legend(h(4),'16+ keV',[0.04 0.85],'color','k','fontsize',15);
%set(hh,'Box','on')
%set(h(4),'fontsize',13);
% 
% 
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % h(5)=irf_panel('epad4');
% % set(h(5),'pos',[0.17 0.68 0.73 0.057]);
% irf_spectrogram(h(5),specepad4,'log','donotfitcolorbarlabel');
% 
% % % % hold on;
% % % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;
% 
% irf_legend(h(5),'(5)',[0.99 0.98],'color','k','fontsize',13);
% set(h(5),'yscale','lin');
% set(h(5),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
% ylabel(h(5),{'25eV'},'fontsize',12,'Interpreter','tex');
% % irf_legend(h(5),'109eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
% set(h(5),'fontsize',13);
% 
% 
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % h(6)=irf_panel('epad5');
% % set(h(6),'pos',[0.17 0.615 0.73 0.057]);
% irf_spectrogram(h(6),specepad5,'log','donotfitcolorbarlabel');
% 
% % % % hold on;
% % % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;
% 
% irf_legend(h(6),'(6)',[0.99 0.98],'color','k','fontsize',13);
% set(h(6),'yscale','lin');
% set(h(6),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
% ylabel(h(6),{'32eV'},'fontsize',12,'Interpreter','tex');
% % irf_legend(h(6),'140eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
% set(h(6),'fontsize',13);
% 
% 
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % h(7)=irf_panel('epad6');
% % set(h(7),'pos',[0.17 0.55 0.73 0.057]);
% irf_spectrogram(h(7),specepad6,'log','donotfitcolorbarlabel');
% 
% % % % hold on;
% % % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;


colormap(h(2),'jet');   %%%%%  调整颜色
colormap(h(3),'jet');
colormap(h(4),'jet');
% colormap(h(5),'jet');
% colormap(h(6),'jet');

irf_adjust_panel_position
irf_zoom(tint,'x',h(1:4))

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['W2F2'];
% print(gcf, '-dpdf', [figname '.pdf']);
