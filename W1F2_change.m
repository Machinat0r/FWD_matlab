close all
clear,clc


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


d=dataobj('D:\MMS\20170823\mms1_fpi_brst_l2_des-dist_20170823153713_v3.3.0.cdf');
b=dataobj('D:\MMS\20170823\mms1_fgm_brst_l2_20170823153713_v5.100.0.cdf');
tint = irf.tint('2017-08-23T15:38:52.00Z/2017-08-23T15:39:02.00Z');

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
% % % thetae = [11.2500   22.5000   33.7500   45.0000   56.2500   67.5000   78.7500   90.0000  101.2500  112.5000  123.7500  135.0000  146.2500  157.5000  168.7500  180.0000];
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


% %% Compute moments
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
[paddiste,thetae,energye,tinte] = get_pitchangledist(diste,phie,thetae,stepTablee,energye0,energye1,Bxyz1,tint);
% [paddiste,thetae,energye,tinte] = get_pitchangledist(diste,Bxyz1,tint);
% [paddisti,thetai,energyi,tinti] = mms.get_pitchangledist(disti,phii,thetai,stepTablei,energyi0,energyi1,Bxyz,tint);

% Convert to DEflux
paddiste = paddiste/1e6/(5.486e-4)^2/0.53707;
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

E1 = [10 13.305175];  
E2 = [13.305175,17.06250025];
E3 = [17.06250025,21.8808505];
E4 = [21.8808505,28.059874999999998];
E5 = [28.059874999999998,35.9838005];
E6 = [35.9838005,46.1453755];

% % E1 = [1000,1500];  
% % E2 = [1500,1900];
% % E3 = [3000,5000];
% % E4 = [14000,18000];
% % E5 = [19000,22000];
% % E6 = [35.9838005,46.1453755];
% 

E7 = [46.1453755,59.176549]
E8 = [59.176549,75.88767425]
E9 = [75.88767425,97.317775]
E10 = [97.317775,124.79975]
E11 = [124.79975,160.0427475]
E12 = [160.0427475,205.237745]
% % 
% % E7 = [1100,1500]
% % E8 = [1500,1900]
% % E9 = [1900,3000]
% % E10 = [3000,4000]
% % E11 = [4000,6000]
% % E12 = [6000,8000]
% % % E13 = [8000,10000]
% % E13 = [1700,7000]
% 
% 
E13 = [205.237745,263.1955025]
E14 = [263.1955025,337.5205]
E15 = [337.5205,432.8342525]
E16 = [432.8342525,555.06401]
E17 = [555.06401,711.8104975]
E18 = [711.8104975,912.820245]

%
E19 = [912.820245,1170.595]
E20 = [1170.595,1501.165]
E21 = [1501.165,1925.08505]
E22 = [1925.08505,2468.717525]
E23 = [2468.717525,3165.86745];
E24 = [3165.86745,4059.88995];

E25=[4059.88995,5206.38]
E26=[5206.38,6676.62755]
E27=[6676.62755,8562.0575]

E28=[8562.0575,10979.929975]
E29=[10979.929975,14080.6]
E30=[14080.6,18056.875]
E31=[18056.875,22021.05]
E32=[22021.05,30000]

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
idx14 = find(energye0 > E14(1) & energye0 < E14(2));
idx15= find(energye0 > E15(1) & energye0 < E15(2));
idx16= find(energye0 > E16(1) & energye0 < E16(2));
idx17= find(energye0 > E17(1) & energye0 < E17(2));
idx18= find(energye0 > E18(1) & energye0 < E18(2));

idx19 = find(energye0 > E19(1) & energye0 < E19(2));
idx20 = find(energye0 > E20(1) & energye0 < E20(2));
idx21= find(energye0 > E21(1) & energye0 < E21(2));
idx22= find(energye0 > E22(1) & energye0 < E22(2));
idx23= find(energye0 > E23(1) & energye0 < E23(2));
idx24= find(energye0 > E24(1) & energye0 < E24(2));

idx25= find(energye0 > E25(1) & energye0 < E25(2));
idx26= find(energye0 > E26(1) & energye0 < E26(2));
idx27= find(energye0 > E27(1) & energye0 < E27(2));
idx28= find(energye0 > E28(1) & energye0 < E28(2));
idx29= find(energye0 > E29(1) & energye0 < E29(2));
idx30= find(energye0 > E30(1) & energye0 < E30(2));
idx31= find(energye0 > E31(1) & energye0 < E31(2));
idx32= find(energye0 > E32(1) & energye0 < E32(2));


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
paddiste14 = squeeze(mean(paddiste.data(:,idx14,:),2));
paddiste15 = squeeze(mean(paddiste.data(:,idx15,:),2));
paddiste16 = squeeze(mean(paddiste.data(:,idx16,:),2));
paddiste17 = squeeze(mean(paddiste.data(:,idx17,:),2));
paddiste18 = squeeze(mean(paddiste.data(:,idx18,:),2));

paddiste19 = squeeze(mean(paddiste.data(:,idx19,:),2));
paddiste20 = squeeze(mean(paddiste.data(:,idx20,:),2));
paddiste21 = squeeze(mean(paddiste.data(:,idx21,:),2));
paddiste22 = squeeze(mean(paddiste.data(:,idx22,:),2));
paddiste23 = squeeze(mean(paddiste.data(:,idx23,:),2));
paddiste24 = squeeze(mean(paddiste.data(:,idx24,:),2));


paddiste25 = squeeze(mean(paddiste.data(:,idx25,:),2));
paddiste26 = squeeze(mean(paddiste.data(:,idx26,:),2));
paddiste27 = squeeze(mean(paddiste.data(:,idx27,:),2));
paddiste28 = squeeze(mean(paddiste.data(:,idx28,:),2));
paddiste29 = squeeze(mean(paddiste.data(:,idx29,:),2));
paddiste30 = squeeze(mean(paddiste.data(:,idx30,:),2));
paddiste31 = squeeze(mean(paddiste.data(:,idx31,:),2));
paddiste32 = squeeze(mean(paddiste.data(:,idx32,:),2));


a=cell(1,6);
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
a{14} = paddiste14();
a{15} = paddiste15();
a{16} = paddiste16();
a{17} = paddiste17();
a{18} = paddiste18();
a{19} = paddiste19();
a{20} = paddiste20();
a{21} = paddiste21();
a{22} = paddiste22();
a{23} = paddiste23();
a{24} = paddiste24();

a{25} = paddiste25();
a{26} = paddiste26();
a{27} = paddiste27();
a{28} = paddiste28();
a{29} = paddiste29();
a{30} = paddiste30();
a{31} = paddiste31();
a{32} = paddiste32();

length=length(paddiste1);

for k=1:32
    
for i=1:length;

    if i>=length
      for j=1:30
        aa{1,k}(i,j)=(a{1,k}(i,j)+a{1,k}(i,j))/2;      
      end
    else 
        for j=1:30
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
specepad3.p_label={'','keV/(cm^2 s sr keV)'};

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

specepad14=struct('t',paddiste.time.epochUnix);
specepad14.p = double(aa{1,14});
specepad14.p_label={' '};
specepad14.f_label={''};
specepad14.f = single(thetae);

specepad15=struct('t',paddiste.time.epochUnix);
specepad15.p = double(aa{1,15});
specepad15.p_label={' '};
specepad15.f_label={''};
specepad15.f = single(thetae);

specepad16=struct('t',paddiste.time.epochUnix);
specepad16.p = double(aa{1,16});
specepad16.p_label={' '};
specepad16.f_label={''};
specepad16.f = single(thetae);

specepad17=struct('t',paddiste.time.epochUnix);
specepad17.p = double(aa{1,17});
specepad17.p_label={' '};
specepad17.f_label={''};
specepad17.f = single(thetae);

specepad18=struct('t',paddiste.time.epochUnix);
specepad18.p = double(aa{1,18});
specepad18.p_label={' '};
specepad18.f_label={''};
specepad18.f = single(thetae);

specepad19=struct('t',paddiste.time.epochUnix);
specepad19.p = double(aa{1,19});
specepad19.p_label={' '};
specepad19.f_label={''};
specepad19.f = single(thetae);

specepad20=struct('t',paddiste.time.epochUnix);
specepad20.p = double(aa{1,20});
specepad20.p_label={' '};
specepad20.f_label={''};
specepad20.f = single(thetae);

specepad21=struct('t',paddiste.time.epochUnix);
specepad21.p = double(aa{1,21});
specepad21.p_label={' '};
specepad21.f_label={''};
specepad21.f = single(thetae);

specepad22=struct('t',paddiste.time.epochUnix);
specepad22.p = double(aa{1,22});
specepad22.p_label={' '};
specepad22.f_label={''};
specepad22.f = single(thetae);

specepad23=struct('t',paddiste.time.epochUnix);
specepad23.p = double(aa{1,23});
specepad23.p_label={' '};
specepad23.f_label={''};
specepad23.f = single(thetae);

specepad24=struct('t',paddiste.time.epochUnix);
specepad24.p = double(aa{1,24});
specepad24.p_label={' '};
specepad24.f_label={''};
specepad24.f = single(thetae);

specepad25=struct('t',paddiste.time.epochUnix);
specepad25.p = double(aa{1,25});
specepad25.p_label={' '};
specepad25.f_label={''};
specepad25.f = single(thetae);

specepad26=struct('t',paddiste.time.epochUnix);
specepad26.p = double(aa{1,26});
specepad26.p_label={' '};
specepad26.f_label={''};
specepad26.f = single(thetae);

specepad27=struct('t',paddiste.time.epochUnix);
specepad27.p = double(aa{1,27});
specepad27.p_label={' '};
specepad27.f_label={''};
specepad27.f = single(thetae);

specepad28=struct('t',paddiste.time.epochUnix);
specepad28.p = double(aa{1,28});
specepad28.p_label={' '};
specepad28.f_label={''};
specepad28.f = single(thetae);

specepad29=struct('t',paddiste.time.epochUnix);
specepad29.p = double(aa{1,29});
specepad29.p_label={' '};
specepad29.f_label={''};
specepad29.f = single(thetae);

specepad30=struct('t',paddiste.time.epochUnix);
specepad30.p = double(aa{1,30});
specepad30.p_label={' '};
specepad30.f_label={''};
specepad30.f = single(thetae);


specepad31=struct('t',paddiste.time.epochUnix);
specepad31.p = double(aa{1,31});
specepad31.p_label={' '};
specepad31.f_label={''};
specepad31.f = single(thetae);


specepad32=struct('t',paddiste.time.epochUnix);
specepad32.p = double(aa{1,32});
specepad32.p_label={' '};
specepad32.f_label={''};
specepad32.f = single(thetae);

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



%% 
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
% h=irf_plot(7,'newfigure');
% %h=irf_figure(540+ic,8);
% xSize=900; ySize=2400;
% set(gcf,'Position',[10 10 xSize ySize]);
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(1)=irf_panel('Bxyz');
% set(h(1),'pos',[0.17 0.94 0.73 0.057]);
irf_plot(h(1),Bxyz1.abs,'linewidth',2);
set(h(1),'fontsize',8);
% ylabel(h(1),'B(nT)','Interpreter','tex','fontsize',17);
% irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.1 0.12])
set(h(1),'Ylim',[8 13], 'ytick',[-5 0 5 10 12 15 20 25],'fontsize',9);
irf_legend(h(1),'(1)',[0.99 0.98],'color','k','fontsize',13)
% irf_legend(h(1),'B(nT)',[-0.04 0.6],'color','w','fontsize',13)
ylabel(h(1),{'B[nT]'},'fontsize',8,'Interpreter','tex');

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(5)=irf_panel('epad22');
% set(h(5),'pos',[0.17 0.68 0.73 0.057]);
irf_spectrogram(h(2),specepad22,'log','donotfitcolorbarlabel');

% % % hold on;
% % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(1),'(2)',[0.99 0.98],'color','k','fontsize',13);
set(h(2),'yscale','lin');
set(h(2),'ytick',[45 90 135],'fontsize',8,'TickDir','in');
caxis(gca,[7.05 7.35])
ylabel(h(2),{'2197eV'},'fontsize',8,'Interpreter','tex');
% irf_legend(h(5),'109eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(2),'fontsize',8);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(6)=irf_panel('epad123');
% set(h(6),'pos',[0.17 0.615 0.73 0.057]);
irf_spectrogram(h(3),specepad23,'log','donotfitcolorbarlabel');

% % % hold on;
% % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(3),'(3)',[0.99 0.98],'color','k','fontsize',13);
set(h(3),'yscale','lin');
set(h(3),'ytick',[45 90 135],'fontsize',8,'TickDir','in');
caxis(gca,[7.2 7.4])
ylabel(h(3),{'2817eV'},'fontsize',8,'Interpreter','tex');
% irf_legend(h(6),'140eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(3),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(7)=irf_panel('epad24');
% set(h(7),'pos',[0.17 0.55 0.73 0.057]);
irf_spectrogram(h(4),specepad24,'log','donotfitcolorbarlabel');

% % % hold on;
% % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(4),'(4)',[0.99 0.98],'color','k','fontsize',13);
set(h(4),'yscale','lin');
set(h(4),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[7.2 7.5])
ylabel(h(4),{'3612eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(7),'179eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(4),'fontsize',13);




h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(2)=irf_panel('epad19');
% set(h(2),'pos',[0.17 0.875 0.73 0.057]);
irf_spectrogram(h(5),specepad25,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(5),'(5)',[0.99 0.98],'color','k','fontsize',13);
set(h(5),'yscale','lin');
set(h(5),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[7.3 7.5])
ylabel(h(5),{'4633eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(2),'55eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(3)=irf_panel('epad20');
% set(h(3),'pos',[0.17 0.81 0.73 0.057]);
irf_spectrogram(h(6),specepad26,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(6),'(6)',[0.99 0.98],'color','k','fontsize',13);
set(h(6),'yscale','lin');
set(h(6),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[7.25 7.5])
ylabel(h(6),{'5941eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(3),'66eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(6),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(4)=irf_panel('epad21');
% set(h(4),'pos',[0.17 0.745 0.73 0.057]);
irf_spectrogram(h(7),specepad27,'log','donotfitcolorbarlabel');

% % % hold on;
% % % irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% % % irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(7),'(7)',[0.99 0.98],'color','k','fontsize',13);
set(h(7),'yscale','lin');
set(h(7),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[7.2 7.45])
ylabel(h(7),{'7619eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(4),'85eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(7),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(5)=irf_panel('epad22');
% set(h(5),'pos',[0.17 0.68 0.73 0.057]);
irf_spectrogram(h(8),specepad28,'log','donotfitcolorbarlabel');
% 
% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(8),'(8)',[0.99 0.98],'color','k','fontsize',13);
set(h(8),'yscale','lin');
set(h(8),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[7 7.4])
ylabel(h(8),{'9771eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(5),'109eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(8),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(6)=irf_panel('epad123');
% set(h(6),'pos',[0.17 0.615 0.73 0.057]);
irf_spectrogram(h(9),specepad29,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(9),'(9)',[0.99 0.98],'color','k','fontsize',13);
set(h(9),'yscale','lin');
set(h(9),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[6.7 7.2])
ylabel(h(9),{'12531eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(6),'140eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(9),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(7)=irf_panel('epad24');
% set(h(7),'pos',[0.17 0.55 0.73 0.057]);
irf_spectrogram(h(10),specepad30,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(10),'(10)',[0.99 0.98],'color','k','fontsize',13);
set(h(10),'yscale','lin');
set(h(10),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[6.4 7])
ylabel(h(10),{'16069eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(7),'179eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
set(h(10),'fontsize',13);


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(2)=irf_panel('epad19');
% set(h(2),'pos',[0.17 0.875 0.73 0.057]);
irf_spectrogram(h(11),specepad31,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(11),'(11)',[0.99 0.98],'color','k','fontsize',13);
set(h(11),'yscale','lin');
set(h(11),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[6 7])
ylabel(h(11),{'20039eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(2),'55eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% h(3)=irf_panel('epad20');
% set(h(3),'pos',[0.17 0.81 0.73 0.057]);
irf_spectrogram(h(12),specepad32,'log','donotfitcolorbarlabel');

% hold on;
% irf_plot([alpha(:,1) alpha(:,2)],'k--', 'Linewidth',0.75); hold on;
% irf_plot([alpha(:,1) alpha(:,3)],'k--', 'Linewidth',0.75); hold off;

irf_legend(h(12),'(12)',[0.99 0.98],'color','k','fontsize',13);
set(h(12),'yscale','lin');
set(h(12),'ytick',[45 90 135],'fontsize',12,'TickDir','in');
caxis(gca,[5.7 6.8])
ylabel(h(12),{'26010eV'},'fontsize',12,'Interpreter','tex');
% irf_legend(h(3),'66eV',[-0.04 0.7],'Rotation',90,'color','k','fontsize',13);
% % % set(h(3),'fontsize',13);


colormap(h(2),'jet');   %%%%%  调整颜色
colormap(h(3),'jet');
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(7),'jet');
colormap(h(8),'jet');
colormap(h(9),'jet');
colormap(h(10),'jet');
colormap(h(11),'jet');
colormap(h(12),'jet');

irf_adjust_panel_position
irf_zoom(tint,'x',h(1:7))

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['C:\Users\fwd\Desktop\Ti~mor~\M\Formation of the rolling-pin distribution of suprathermal electrons behind dipolarization fronts\obs\4633-16069'];
print(gcf, '-dpdf', [figname '.pdf']);

