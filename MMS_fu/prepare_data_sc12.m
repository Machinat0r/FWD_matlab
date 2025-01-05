ic = 1:4;
units = irf_units;

%% SC4 distance
% function d=distance(x,y)
r1=(gseR1.x-gseR2.x)^2+(gseR1.y-gseR2.y)^2+(gseR1.z-gseR2.z)^2;
r2=(gseR3.x-gseR2.x)^2+(gseR3.y-gseR2.y)^2+(gseR3.z-gseR2.z)^2;
r3=(gseR3.x-gseR4.x)^2+(gseR3.y-gseR4.y)^2+(gseR3.z-gseR4.z)^2;
r4=(gseR1.x-gseR4.x)^2+(gseR1.y-gseR4.y)^2+(gseR1.z-gseR4.z)^2;
r5=(gseR4.x-gseR2.x)^2+(gseR4.y-gseR2.y)^2+(gseR4.z-gseR2.z)^2;
r6=(gseR1.x-gseR3.x)^2+(gseR1.y-gseR3.y)^2+(gseR1.z-gseR3.z)^2;
d12=sqrt(r1.data);
d23=sqrt(r2.data);
d34=sqrt(r3.data);
d14=sqrt(r4.data);
d42=sqrt(r5.data);
d13=sqrt(r6.data);
A=[d12 d23 d34 d14 d42 d13];
d=mean(A,2);


%% Current
% Current from magnetic field (curlometer) 
if all(ic==[1:4])
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
c_eval(['gsmJcurl=irf_gse2gsm(gseJcurl);'],ic);
gsmJcurlfac = irf_convert_fac(gsmJcurl,gsmB1,[1 0 0]); 
end

ic = 1:4;
ic_B = 1:4;
% Currents from moments, use ne also for Ji 
% 


% c_eval('gsmJe? = -units.e*ne?*gsmVe?*1e3*1e6*1e9; gsmJe?.units = ''nA/m^2''; gsmJe?.coordinateSystem = ''gsm'';',ic);
% c_eval('gsmJi? = units.e*ne?*gsmVi?.resample(ne?.time)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);

c_eval('gsmJe? = -units.e*ne?.resample(gsmB?)*gsmVe?.resample(gsmB?)*1e3*1e6*1e9; gsmJe?.units = ''nA/m^2''; gsmJe?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJi? = units.e*ne?.resample(gsmB?)*gsmVi?.resample(gsmB?)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);
% c_eval('gsmJi? = units.e*ni?*gsmVi?.resample(ni?.time)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);

% c_eval('gsmJe? = -units.e*ne?*gsmVe?*1e3*1e6*1e9; gsmJe?.units = ''nA/m^2''; gsmJe?.coordinateSystem = ''gsm'';',ic);
% c_eval('gsmJi? = units.e*ne?*gsmVi?.resample(gsmB?.time)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJ? = (gsmJe?+gsmJi?);',ic);
    
 
%% Perpendicular and parallel decomposition
% Celocity and current

c_eval('[gsmVe?par,gsmVe?perp] = irf_dec_parperp(gsmB?,gsmVe?); gsmVe?par.name = ''Ve par''; gsmVe?perp.name = ''Ve perp'';',ic)
c_eval('[gsmVe?fac] = irf_convert_fac(gsmVe?,gsmB?,[1 0 0]); ',ic)

c_eval('[gsmVi?par,gsmVi?perp] = irf_dec_parperp(gsmB?,gsmVi?); gsmVi?par.name = ''Vi par''; gsmVi?perp.name = ''Vi perp'';',ic)
c_eval('[gsmJ?par,gsmJ?perp] = irf_dec_parperp(gsmB?,gsmJ?); gsmJ?par.name = ''J par''; gsmJ?perp.name = ''J perp'';',ic)


c_eval('[gsmVi?fac] = irf_convert_fac(gsmVi?,gsmB?,[1 0 0]); ',ic)
c_eval('[gsmJ?fac] = irf_convert_fac(gsmJ?,gsmB?,[1 0 0]); ',ic)

c_eval('[gsmJi?fac] = irf_convert_fac(gsmJi?,gsmB?,[1 0 0]); ',ic)
c_eval('[gsmJe?fac] = irf_convert_fac(gsmJe?,gsmB?,[1 0 0]); ',ic)

% Electric fields
c_eval('[gsmE?par,gsmE?perp] = irf_dec_parperp(gsmB?,gsmE?); gsmE?par.name = ''E par''; gsmE?perp.name = ''E perp'';',ic_B)
c_eval('[gsmE?fac] = irf_convert_fac(gsmE?,gsmB?,[1 0 0]); ',ic_B)
% Wave magnetic field
%gsmB2scm = gsmB2scm{2};
% c_eval('[gsmB?scmpar,gsmB?scmperp] = irf_dec_parperp(gsmB?,gsmB?scm); gsmB?scmpar.name = ''B par scm''; gsmB?scmperp.name = ''B perp scm'';',ic_B)
try
c_eval('[gsmE?hmfepar,gsmE?hmfeperp] = irf_dec_parperp(gsmB?,gsmE?hmfe); gsmE?hmfepar.name = ''E par hmfe''; gsmE?hmfeperp.name = ''E perp hmfe'';',ic_B)
end


%% LMN

% irf_minvar_gui(gsmB1.tlim(tint));
L=[1,0,0];
M=[0,1,0];
N=[0,0,1];

% c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')
c_eval('mvaB? = irf.ts_vec_xyz(gsmB?.time,[gsmB?.dot(L).data gsmB?.dot(M).data gsmB?.dot(N).data]); mvaB?.name = ''B LMN'';',ic_B)
c_eval('mvaE? = irf.ts_vec_xyz(gsmE?.time,[gsmE?.dot(L).data gsmE?.dot(M).data gsmE?.dot(N).data]); mvaE?.name = ''E LMN'';',ic)
c_eval('mvaE?l= mvaE?.x',ic)
c_eval('mvaE?m= mvaE?.y',ic)
c_eval('mvaE?n= mvaE?.z',ic)
% smooth
c_eval('mvaE?lsmoothdata= smooth(mvaE?l.data, 200);',ic)
c_eval('mvaE?lsmooth=irf.ts_scalar(mvaE?l.time, mvaE?lsmoothdata); ',ic)
c_eval('mvaE?msmoothdata= smooth(mvaE?m.data,200);',ic)
c_eval('mvaE?msmooth=irf.ts_scalar(mvaE?m.time, mvaE?msmoothdata); ',ic)
c_eval('mvaE?nsmoothdata= smooth(mvaE?n.data,200);',ic)
c_eval('mvaE?nsmooth=irf.ts_scalar(mvaE?n.time, mvaE?nsmoothdata); ',ic)
% resample
c_eval('mvaE?resample = mvaE?.resample(gsmVe?); mvaE?resample.units = ''mV/m'';',ic)

c_eval('mvaVe? = irf.ts_vec_xyz(gsmVe?.time,[gsmVe?.dot(L).data gsmVe?.dot(M).data gsmVe?.dot(N).data]); mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = irf.ts_vec_xyz(gsmVi?.time,[gsmVi?.dot(L).data gsmVi?.dot(M).data gsmVi?.dot(N).data]); mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaJ? = irf.ts_vec_xyz(gsmJ?.time,[gsmJ?.dot(L).data gsmJ?.dot(M).data gsmJ?.dot(N).data]); mvaJ?.units = gsmJ?.units; mvaJ?.name = ''J LMN'';',ic)
c_eval('mvaJ?l= mvaJ?.x',ic)
c_eval('mvaJ?m= mvaJ?.y',ic)
c_eval('mvaJ?n= mvaJ?.z',ic)



c_eval('mvaJe? = irf.ts_vec_xyz(gsmJe?.time,[gsmJe?.dot(L).data gsmJe?.dot(M).data gsmJe?.dot(N).data]); mvaJe?.units = gsmJe?.units; mvaJe?.name = ''Je LMN'';',ic)
c_eval('mvaJi? = irf.ts_vec_xyz(gsmJi?.time,[gsmJi?.dot(L).data gsmJi?.dot(M).data gsmJi?.dot(N).data]); mvaJi?.units = gsmJi?.units; mvaJi?.name = ''Ji LMN'';',ic)
c_eval('mvaJcurl = irf.ts_vec_xyz(gsmJcurl.time,[gsmJcurl.dot(L).data gsmJcurl.dot(M).data gsmJcurl.dot(N).data]); mvaJcurl.units = gsmJcurl.units; mvaJcurl.name = ''Jcurl LMN'';')
c_eval('mvaJcurl_l= mvaJcurl.x',ic)
c_eval('mvaJcurl_m= mvaJcurl.y',ic)
c_eval('mvaJcurl_n= mvaJcurl.z',ic)
c_eval('[mvaJ?fac] = irf_convert_fac(mvaJ?,mvaB?,[1 0 0]); ',ic)


c_eval('mvaPe? = irf.ts_vec_xyz(gsePe?.time,[gsePe?.dot(L).data gsePe?.dot(M).data gsePe?.dot(N).data]); mvaPe?.units = gsePe?.units; mvaPe?.name = ''Pe LMN'';',ic)
c_eval(['mvaPe?=irf_gse2gsm(mvaPe?);'],ic);
c_eval('mvaPe?l= mvaPe?.x',ic)
c_eval('mvaPe?m= mvaPe?.y',ic)
c_eval('mvaPe?n= mvaPe?.z',ic)


c_eval('mvaR? = irf.ts_vec_xyz(gsmR?.time,[gsmR?.dot(L).data gsmR?.dot(M).data gsmR?.dot(N).data]); mvaR?.units = gsmR?.units; mvaR?.name = ''R LMN'';',ic)

%% Cross products
% ExB drift
c_eval('gsmVExB? = cross(gsmE?.resample(gsmB?.time),gsmB?)/gsmB?.abs/gsmB?.abs*1e3; gsmVExB?.units = '''';',ic) % km/s

% Convective electric fields
c_eval('gsmVexB? = gsmVe?.cross(gsmB?.resample(gsmVe?))*1e-3; gsmVexB?.units = ''mV/m'';',ic)
c_eval('gsmVixB? = gsmVi?.cross(gsmB?.resample(gsmVi?))*1e-3; gsmVixB?.units = ''mV/m'';',ic)

c_eval('mvaVexB? = irf.ts_vec_xyz(gsmVexB?.time,[gsmVexB?.dot(L).data gsmVexB?.dot(M).data gsmVexB?.dot(N).data]); mvaVexB?.name = ''VexB LMN'';',ic)
c_eval('mvaVixB? = irf.ts_vec_xyz(gsmVixB?.time,[gsmVixB?.dot(L).data gsmVixB?.dot(M).data gsmVixB?.dot(N).data]); mvaVixB?.name = ''VixB LMN'';',ic)
% c_eval('mvaVixB?l= mvaVixB?.x',ic)
% c_eval('mvaVixB?m= mvaVixB?.y',ic)
% c_eval('mvaVixB?n= mvaVixB?.z',ic)
% Non-ideal electric field, E+VexB
c_eval('gsmEVexB? = gsmE?.resample(gsmVexB?.time)+gsmVexB?; gsmEVexB?.name = ''E+VexB'';',ic)
% c_eval('testmvaEVexB? = irf.ts_vec_xyz(gsmEVexB?.time,[gsmEVexB?.dot(L).data gsmEVexB?.dot(M).data gsmEVexB?.dot(N).data]); testmvaEVexB?.name = ''testE+VexB LMN'';',ic)

c_eval('mvaEVexB? = mvaE?.resample(mvaVexB?.time)+mvaVexB?; mvaEVexB?.name = ''E+VexB'';',ic)


% JxB
c_eval('gsmJxB? = gsmJ?.resample(gsmE?).cross(gsmB?.resample(gsmE?));',ic)
c_eval('mvaJxB? = irf.ts_vec_xyz(gsmJxB?.time,[gsmJxB?.dot(L).data gsmJxB?.dot(M).data gsmJxB?.dot(N).data]); mvaJxB?.name = ''JxB LMN'';',ic)

%% ohmslaw
% JxB/ne
c_eval('gsmOhmJxB? = gsmJxB?/ne?/units.e*1e-9*1e-9*1e-6*1e3; gseOhmJxB?.units = ''mV/m'';',ic)
c_eval('mvaOhmJxB? = mvaJxB?/ne?/units.e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB?.units = ''mV/m'';',ic)
% c_eval('mvaOhmJxB?l= mvaJxB?.x',ic)
% c_eval('mvaOhmJxB?m= mvaJxB?.y',ic)
% c_eval('mvaOhmJxB?n= mvaJxB?.z',ic)

% ?Pe/ne
% sc1
% deltaPe=diff(Pe);
% deltaPe_deltaN=deltaPe/?t*V

% % delta Pe
GPeXX = c_4_grad('gseR?','gsePe?.xx','grad');
GPeXY = c_4_grad('gseR?','gsePe?.xy','grad');
GPeXZ = c_4_grad('gseR?','gsePe?.xz','grad');
GPeYY = c_4_grad('gseR?','gsePe?.yy','grad');
GPeYZ = c_4_grad('gseR?','gsePe?.yz','grad');
GPeZZ = c_4_grad('gseR?','gsePe?.zz','grad');
GPeX = -(GPeXX.data(:,1)+GPeXY.data(:,2)+GPeXZ.data(:,3))*1e-3;
GPeY = -(GPeXY.data(:,1)+GPeYY.data(:,2)+GPeYZ.data(:,3))*1e-3;
GPeZ = -(GPeXZ.data(:,1)+GPeYZ.data(:,2)+GPeZZ.data(:,3))*1e-3;%nPa/m
divPe = irf.ts_vec_xyz(GPeXX.time,[GPeX GPeY GPeZ]);
divPe = irf_gse2gsm(divPe);
c_eval('divPe?ne = divPe/ne?/units.e*1e-12; divPe?ne.units = ''mV/m'';', ic);
% divPe=irf.ts2mat(divPe);
c_eval('divPe?nel = divPe?ne.x;',ic);
c_eval('divPe?nem = divPe?ne.y;',ic);
c_eval('divPe?nen = divPe?ne.z;',ic);
% c_eval('divPe?ne= irf.ts_scalar(divPe?t, {divPe?nel, divPe?nem, divPe?nen}); divPe?ne.units = ''mV/m'';divPe?ne.tensorOrder=0;',ic)

% E total
c_eval('mvaE?t_data=-mvaVixB?.data+mvaOhmJxB?.resample(mvaVixB?.time).data-divPe?nen.resample(mvaVixB?.time).data;',ic)
c_eval('mvaE?t=irf.ts_vec_xyz(mvaVixB?.time, mvaE?t_data);mvaE?t.units = ''mV/m'';mvaE?t.name = ''E ohm LMN''; ',ic)
% sc2
% DmvaPe12=mvaPe1.resample(mvaPe2.time)-mvaPe2;
% DmvaPe23=mvaPe2.resample(mvaPe3.time)-mvaPe3;
DmvaPe21=mvaPe2.resample(mvaPe1.time)-mvaPe1;
Rmva21=mvaR2-mvaR1;
divPeS2=DmvaPe21/11.1736/ne1/units.e*1e-15;



%% Magnetic field curvature 
% % if all(ic==[1:4])
% c_eval('R? = gsmR?.resample(gsmB4);',1:4)
% c_eval('B? = gsmB?.resample(gsmB4);',1:4)
% [gsmCurvB,avB]=c_4_grad('R?','B?','curvature'); gsmCurvB.name = 'curv B'; gsmCurvB.coordinateSystem = 'gsm';
% curvBradius = 1/gsmCurvB.abs; curvBradius.name = 'R_c';
% % end
%% fpi moments
if 0
    species='ion';% 'ion' or 'electron'

filepath_and_filename = mms.get_filepath(['mms',num2str(ic),'_fpi_brst_l2_d',species(1),'s-dist'],tint);
[PDist,PDistError] = mms.make_pdist(filepath_and_filename);
% c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
% c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

e_phi = PDist.depend{1,2};
e_theta = PDist.depend{1,3};
e_energy0 =PDist.ancillary.energy0;
e_energy1 =PDist.ancillary.energy1;
esteptable=PDist.ancillary.esteptable;
e_phi=irf.ts_scalar(PDist.time,e_phi);
%e_theta=irf.ts_scalar(PDist.time,e_theta);
esteptable=irf.ts_scalar(PDist.time,esteptable);

energyrange1=[2000,29000]% 1
%  energyrange1=[3000,29000];%2
% energyrange1=[3500,29000]% 0722
% energyrange1=[2000,29000];%set energy range
% energyrange=[3000,29000];%set energy range
pl_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange1);
c_eval('Tfac= mms.rotate_tensor(pl_moments.T_psd,''fac'',dmpaB?);',ic);
T_scalar=1/3*(Tfac.xx.data+Tfac.yy.data+Tfac.zz.data);
T_scalar_psd=irf.ts_scalar(Tfac.time,T_scalar);
T_para_psd=irf.ts_scalar(Tfac.time,Tfac.xx.data);
T_perp_psd=irf.ts_scalar(Tfac.time,1/2*(Tfac.yy.data+Tfac.zz.data));
% c_eval('T_scalar_psd?=irfd.Nan_remove(T_scalar_psd);',ic)
% c_eval('T_para_psd?=irfd.Nan_remove(T_para_psd);',ic)
% c_eval('T_perp_psd?=irfd.Nan_remove(T_perp_psd);',ic)
% c_eval('N_psd?=irfd.Nan_remove(pl_moments.n_psd);',ic)
% c_eval('V_psd?=irfd.Nan_remove(pl_moments.V_psd);',ic)
c_eval('T1_scalar_psd?=T_scalar_psd;',ic)
c_eval('T1_para_psd?=T_para_psd;',ic)
c_eval('T1_perp_psd?=T_perp_psd;',ic)
c_eval('N1_psd?=pl_moments.n_psd;',ic)
c_eval('V1_psd0?=pl_moments.V_psd;',ic)
c_eval('V1_psd?=irf_gse2gsm(V1_psd0?);',ic);
c_eval('[gsmV1_psd?par,gsmV1_psd?perp] = irf_dec_parperp(gsmB?,V1_psd?); V1_psd?par.name = ''V1_psd par''; V1_psd?perp.name = ''V1_psd perp'';',ic)
c_eval('[gsmV1_psd?fac] = irf_convert_fac(V1_psd?,gsmB?,[1 0 0]); ',ic)
% Units=irf_units;
% pressure=N1_psd1.data.*Units.kB.*T1_scalar_psd1.data;
c_eval('Pi1_psd?=1e6*11604.5*1e9*1.38e-23*N1_psd?.data.*T1_scalar_psd?.data; ',ic);%nPa
c_eval('Pi1_psd?=irf.ts_scalar(T1_scalar_psd1.time,Pi1_psd?);',ic)
c_eval('Pd1_psd?= irf.ts_scalar(N1_psd?.time,2e-6*double(N1_psd?.data).*V1_psd?.abs2.data);',ic)

c_eval('PB? = gsmB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('Ptmp?=irf_add(1,PB?,1,Pi1_psd?);',ic);
c_eval('Pt1_psd?=irf_add(1,Ptmp?,1,Pd1_psd?);',ic);
%
energyrange2=[10,900]% 1 two cold ion
% energyrange2=[400,2900];% 2
% energyrange2=[300,1000]% 0722
% energyrange2=[400,1200];%set energy range
% energyrange=[400,2700];%set energy range
pl_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange2);
c_eval('Tfac= mms.rotate_tensor(pl_moments.T_psd,''fac'',dmpaB?);',ic);
T_scalar=1/3*(Tfac.xx.data+Tfac.yy.data+Tfac.zz.data);
T_scalar_psd=irf.ts_scalar(Tfac.time,T_scalar);
T_para_psd=irf.ts_scalar(Tfac.time,Tfac.xx.data);
T_perp_psd=irf.ts_scalar(Tfac.time,1/2*(Tfac.yy.data+Tfac.zz.data));
c_eval('T2_scalar_psd?=T_scalar_psd;',ic)
c_eval('T2_para_psd?=T_para_psd;',ic)
c_eval('T2_perp_psd?=T_perp_psd;',ic)
c_eval('N2_psd?=pl_moments.n_psd;',ic)
c_eval('V2_psd0?=pl_moments.V_psd;',ic)
c_eval('V2_psd?=irf_gse2gsm(V2_psd0?);',ic);
c_eval('[gsmV2_psd?par,gsmV2_psd?perp] = irf_dec_parperp(gsmB?,V2_psd?); V2_psd?par.name = ''V2_psd par''; V2_psd?perp.name = ''V2_psd perp'';',ic)
c_eval('[gsmV2_psd?fac] = irf_convert_fac(V2_psd?,gsmB?,[1 0 0]); ',ic)

c_eval('Pi2_psd?=1e6*11604.5*1e9*1.38e-23*N2_psd?.data.*T2_scalar_psd?.data; ',ic);%nPa
c_eval('Pi2_psd?=irf.ts_scalar(T2_scalar_psd?.time,Pi2_psd?);',ic)
c_eval('Pd2_psd?= irf.ts_scalar(N2_psd?.time,2e-6*double(N2_psd?.data).*V2_psd?.abs2.data);',ic)


%
energyrange0=[10,28301.8];%set energy range
% energyrange=[400,2700];%set energy range
plt_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange0);
c_eval('TTfac= mms.rotate_tensor(plt_moments.T_psd,''fac'',dmpaB?);',ic);
TT_scalar=1/3*(TTfac.xx.data+TTfac.yy.data+TTfac.zz.data);
TT_scalar_psd=irf.ts_scalar(TTfac.time,TT_scalar);
TT_para_psd=irf.ts_scalar(TTfac.time,TTfac.xx.data);
TT_perp_psd=irf.ts_scalar(TTfac.time,1/2*(TTfac.yy.data+TTfac.zz.data));
c_eval('TT_scalar_psd?=TT_scalar_psd;',ic)
c_eval('NT_psd?=plt_moments.n_psd;',ic)
c_eval('VT_psd0?=pl_moments.V_psd;',ic)
c_eval('VT_psd?=irf_gse2gsm(VT_psd0?);',ic);
c_eval('[gsmVT_psd?par,gsmVT_psd?perp] = irf_dec_parperp(gsmB?,VT_psd?); VT_psd?par.name = ''VT_psd par''; VT_psd?perp.name = ''VT_psd perp'';',ic)

% c_eval('PiT_psd?=1e6*11604.5*1e9*1.38e-23*NT_psd?.data.*TT_scalar_psd?.data; ',ic);%nPa
% c_eval('PiT_psd?=irf.ts_scalar(TT_scalar_psd?.time,PiT_psd?);',ic)
c_eval('PdT_psd?= irf.ts_scalar(NT_psd?.time,2e-6*double(NT_psd?.data).*VT_psd?.abs2.data);',ic)
c_eval('PiT_psd?=irf_add(1,Pi1_psd?,1,Pi2_psd?); ',ic);

c_eval('PTtmp?=irf_add(1,PB?,1,PiT_psd?);',ic);
c_eval('PtT_psd?=irf_add(1,PTtmp?,1,PdT_psd?);',ic);
end
%% Remove low counts
if 0
eint = [000 40000];
iDist = iPDist1.tlim(tint).elim(eint);
% first try: just remove everything below certain energy
iDist_rmE10 = iDist;
iDist_rmE10.data(:,1:10,:,:) = 0;
iDist_rmE12 = iDist;
iDist_rmE12.data(:,1:12,:,:) = 0;

% Make reduced distributions

% 1D: vipar
% 2D: vix, vit, viz
eint = [000 40000];
vint = [-Inf Inf];
par = dmpaB1.resample(iDist).norm;

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% vg = [-2200:50:2200];
vg = linspace(-3e4,3e4,100); % km/s 


% 
% tic; if1D = iDist.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
% tic; if2D_xy = iDist.reduce('2D',x,y,'vg',vg,'base','cart'); toc
% tic; if2D_xz = iDist.reduce('2D',x,z,'vg',vg,'base','cart'); toc
% tic; if2D_yz = iDist.reduce('2D',y,z,'vg',vg,'base','cart'); toc
% 
% tic; if1D_rmE10 = iDist_rmE10.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
% tic; if2D_xy_rmE10 = iDist_rmE10.reduce('2D',x,y,'vg',vg,'base','cart'); toc
% tic; if2D_xz_rmE10 = iDist_rmE10.reduce('2D',x,z,'vg',vg,'base','cart'); toc
% tic; if2D_yz_rmE10 = iDist_rmE10.reduce('2D',y,z,'vg',vg,'base','cart'); toc

tic; if1D_rmE12 = iDist_rmE12.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
tic; if2D_xy_rmE12 = iDist_rmE12.reduce('2D',x,y,'vg',vg,'base','cart'); toc
tic; if2D_xz_rmE12 = iDist_rmE12.reduce('2D',x,z,'vg',vg,'base','cart'); toc
tic; if2D_yz_rmE12 = iDist_rmE12.reduce('2D',y,z,'vg',vg,'base','cart'); toc
end

%% Pitchangle distributions
% c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
%   c_eval('iPitch? = iPDist?.pitchangles(dmpaB?,15);',ic)
if 0
  load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
  %load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
elseif 1
  %%
  ictmp=ic;
%   ic = 1;
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
  c_eval('iPitch? = iPDist?.pitchangles(dmpaB?,15);',ic)
  %c_eval('ePitch?par = ePDist?.pitchangles(dmpaB?,[0 15]);',ic)
  %c_eval('ePitch?perp = ePDist?.pitchangles(dmpaB?,[75 105]);',ic)
  %c_eval('ePitch?apar = ePDist?.pitchangles(dmpaB?,[165 180]);',ic)
  ic = ictmp;
end

%% Calculate some additional parameters, irf_plasma_calc
if 1
% Speeds
c_eval('matB? = gsmB?.abs.data;',ic)
c_eval('matParTe? = facTe?.xx.resample(gsmB?.time).data;',ic)
c_eval('matParTi? = facTi?.xx.resample(gsmB?.time).data;',ic)
c_eval('matPerTe? = (facTe?.yy.resample(gsmB?.time).data + facTe?.zz.resample(gsmB?.time).data)/2;',ic)
c_eval('matPerTi? = (facTi?.yy.resample(gsmB?.time).data + facTi?.zz.resample(gsmB?.time).data)/2;',ic)
c_eval('matTe? = facTe?.trace.resample(gsmB?.time).data/3;',ic)
c_eval('matTi? = facTi?.trace.resample(gsmB?.time).data/3;',ic)
c_eval('matNe? = ne?.resample(gsmB?.time).data;',ic)

c_eval('vte?perp = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Vte''); vte?perp = irf.ts_scalar(gsmB?.time,vte?perp)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte?par = irf_plasma_calc(matB?,matNe?,0,matParTe?,matTi?,''Vte''); vte?par = irf.ts_scalar(gsmB?.time,vte?par)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vte? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Vte''); vte? = irf.ts_scalar(gsmB?.time,vte?)*1e-3; vte?.units = ''km/s'';',ic)
c_eval('vtp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matPerTi?,''Vtp''); vtp? = irf.ts_scalar(gsmB?.time,vtp?)*1e-3; vtp?.units = ''km/s'';',ic)
c_eval('vA? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Va''); vA? = irf.ts_scalar(gsmB?.time,vA?)*1e-3; vA?.units = ''km/s'';',ic)

% Frequencies
c_eval('flh? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Flh''); flh? = irf.ts_scalar(gsmB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fce''); fce? = irf.ts_scalar(gsmB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fcp''); fcp? = irf.ts_scalar(gsmB?.time,fcp?);',ic)
c_eval('fpe? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpe''); fpe? = irf.ts_scalar(gsmB?.time,fpe?);',ic)
c_eval('fpp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpp''); fpp? = irf.ts_scalar(gsmB?.time,fpp?);',ic)

% Length scales
c_eval('Lp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Li''); Lp? = irf.ts_scalar(gsmB?.time,Lp?)*1e-3; Lp?.units = ''km''; Lp?.name=''p inertial length'';',ic)
c_eval('Le? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Le''); Le? = irf.ts_scalar(gsmB?.time,Le?)*1e-3; Le?.units = ''km''; Le?.name=''e inertial length'';',ic)
c_eval('Ld? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Ld''); Ld? = irf.ts_scalar(gsmB?.time,Ld?)*1e-3; Ld?.units = ''km''; Ld?.name=''Debye length'';',ic)
c_eval('re? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Roe''); re? = irf.ts_scalar(gsmB?.time,re?)*1e-3; re?.units = ''km''; re?.name=''e gyroradius'';',ic)
c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gsmB?.time,rp?)*1e-3; rp?.units = ''km''; rp?.name=''p gyroradius'';',ic)

%c_eval('beta? = (re?/Le?).^2;',ic)

% Pressure 
c_eval('PB? = gsmB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('Pi?para = irf.ts_scalar(facPi?.time,facPi?.xx.data);',ic)
% c_eval('Pi_perp?=1e6*11604.5*1e9*1.38e-23*Ni?.data.*Ti_perp?.data; ',ic);%nPa
% c_eval('Pi_perp?=irf.ts_scalar(Ti_perp1.time,Pi_perp?);',ic)
% c_eval('Pic?para = irf.ts_scalar(ni?.time,1e6*11604.5*1e9*1.38e-23*ni?.data.*Ti_para(:,2));',ic)
c_eval('Pi?perp = irf.ts_scalar(facPi?.time,(facPi?.yy.data+facPi?.zz.data)/2);',ic)
c_eval('Pd?= irf.ts_scalar(ni?.time,2e-6*double(ni?.data).*gsmVi?.abs2.data);',ic)
%    Ptmp=irf_add(1,Pcis,1,PB);Pt=irf_add(1,Ptmp,1,Pd);
c_eval('Pi?=irf_add(1/3,Pi?para,2/3,Pi?perp);',ic);
c_eval('Ptmp?=irf_add(1,PB?,1,Pi?);',ic);
c_eval('Pt?=irf_add(1,Ptmp?,1,Pd?);',ic);

% % Pressure hot
% c_eval('PiH?para = irf.ts_scalar(N1_psd?.time,1e6*11604.5*1e9*1.38e-23*N1_psd?.data.*T1_para_psd?.data);',ic)
% c_eval('PiH?perp = irf.ts_scalar(N1_psd?.time,1e6*11604.5*1e9*1.38e-23*N1_psd?.data.*T1_perp_psd?.data);',ic)
% c_eval('PdH?= irf.ts_scalar(N1_psd?.time,2e-6*double(N1_psd?.data).*V1_psd?.abs2.data);',ic)
% c_eval('PiH?=irf_add(1/3,PiH?para,2/3,PiH?perp);',ic);
% c_eval('PtmpH?=irf_add(1,PB?,1,PiH?);',ic);
% c_eval('PtH?=irf_add(1,PtmpH?,1,PdH?);',ic);
% % Pressure cold
% c_eval('PiC?para = irf.ts_scalar(N2_psd?.time,1e6*11604.5*1e9*1.38e-23*N2_psd?.data.*T2_para_psd?.data);',ic)
% c_eval('PiC?perp = irf.ts_scalar(N2_psd?.time,1e6*11604.5*1e9*1.38e-23*N2_psd?.data.*T2_perp_psd?.data);',ic)
% c_eval('PdC?= irf.ts_scalar(N2_psd?.time,2e-6*double(N2_psd?.data).*V2_psd?.abs2.data);',ic)
% c_eval('PiC?=irf_add(1/3,PiC?para,2/3,PiC?perp);',ic);
% c_eval('PtmpC?=irf_add(1,PB?,1,PiC?);',ic);
% c_eval('PtC?=irf_add(1,PtmpC?,1,PdC?);',ic);

% Magnetic moment
c_eval('mag_mom? = 0.5*units.me*vte?perp.^2*10^6/(gsmB?.abs*1e-9)*1e9;  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)
end
%% EDR signatures
if 1
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

% Compute electron Mach number
c_eval('Me?perp = gsmVe?perp.abs/vte?perp;',ic);
c_eval('Me?par = gsmVe?par.abs/vte?par;',ic);

% Compute current density and J.E

c_eval('EdotJi? = gsmE?.resample(gsmJi?).dot(gsmJi?); EdotJi?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJe? = gsmE?.resample(gsmJe?).dot(gsmJe?); EdotJe?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)


c_eval('EdotJ?curl = mvaE?.resample(mvaJcurl).dot(mvaJcurl); mvaEdotJ?curl.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ? = gsmE?.resample(gsmJ?).dot(gsmJ?); EdotJ?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% c_eval('EdotJ? = gsmE?.dot(gsmJ?.resample(gsmE?))/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
%% VixB dot J
c_eval('VixBdotJ?x = gsmVixB?.x.resample(gsmJ?)*gsmJ?.x;', ic);
c_eval('VixBdotJ?y = gsmVixB?.y.resample(gsmJ?)*gsmJ?.y;', ic);
c_eval('VixBdotJ?z = gsmVixB?.z.resample(gsmJ?)*gsmJ?.z;', ic);
c_eval('VixBdotJ? = VixBdotJ?x + VixBdotJ?y + VixBdotJ?z;', ic);

%% delta Pe / ne dot J
c_eval('divPe?nedotJx = divPe?ne.x.resample(gsmJ?)*gsmJ?.x;', ic);
c_eval('divPe?nedotJy = divPe?ne.y.resample(gsmJ?)*gsmJ?.y;', ic);
c_eval('divPe?nedotJz = divPe?ne.z.resample(gsmJ?)*gsmJ?.z;', ic);
c_eval('divPe?nedotJ = divPe?nedotJx + divPe?nedotJy + divPe?nedotJz;', ic);

%% check JxB/ne . J
% c_eval('JxBdotJ? = gsmOhmJxB?.dot(gsmJ?.resample(gsmOhmJxB?))/1000; JxBdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('JxBdotJ?x = gsmOhmJxB?.x*(gsmJ?.x.resample(gsmOhmJxB?.x)); JxBdotJ?x.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('JxBdotJ?y = gsmOhmJxB?.y*(gsmJ?.y.resample(gsmOhmJxB?.y)); JxBdotJ?x.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('JxBdotJ?z = gsmOhmJxB?.z*(gsmJ?.z.resample(gsmOhmJxB?.z)); JxBdotJ?x.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('JxBdotJ? = JxBdotJ?x + JxBdotJ?y +JxBdotJ?z; JxBdotJ?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)


% c_eval('EdotJ? = gsmE?.resample(gsmB?).dot(gsmJ?.resample(gsmB?))/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)


c_eval('gsmEdotJ?x = gsmE?.x.resample(gsmJ?.x)*(gsmJ?.x); gsmEdotJ?x.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('gsmEdotJ?y = gsmE?.y.resample(gsmJ?.y)*(gsmJ?.y); gsmEdotJ?y.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('gsmEdotJ?z = gsmE?.z.resample(gsmJ?.z)*(gsmJ?.z); gsmEdotJ?z.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('gsmEdotJ?par = gsmE?fac.z.resample(gsmJ?fac.z)*(gsmJ?fac.z); gamdotJ?par.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('gsmEdotJ?perp = gsmE?fac.x.resample(gsmJ?fac.x)*(gsmJ?fac.x)+gsmE?fac.y.resample(gsmJ?fac.y)*(gsmJ?fac.y); gamdotJ?perp.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

% % % c_eval('gsmEdotJ?x = gsmE?.x*(gsmJ?.x.resample(gsmE?.x))/1000; gsmEdotJ?x.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% % % c_eval('gsmEdotJ?y = gsmE?.y*(gsmJ?.y.resample(gsmE?.y))/1000; gsmEdotJ?y.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% % % c_eval('gsmEdotJ?z = gsmE?.z*(gsmJ?.z.resample(gsmE?.z))/1000; gsmEdotJ?z.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% % % c_eval('gsmEdotJ?par = gsmE?fac.z*(gsmJ?fac.z.resample(gsmE?fac.z))/1000; gamdotJ?par.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% % % c_eval('gsmEdotJ?perp = gsmE?fac.x*(gsmJ?fac.x.resample(gsmE?fac.x))/1000+gsmE?fac.y*(gsmJ?fac.y.resample(gsmE?fac.y))/1000; gamdotJ?perp.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

c_eval('EdotJ?par = gsmE?par*gsmJ?par.resample(gsmE?par); EdotJ?par.units = ''pW/m^3''; EdotJ?par.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% c_eval('EdotJ?perp = gsmE?perp.dot(gsmJ?perp.resample(gsmE?perp))/1000; EdotJ?perp.units = ''nW/m^3''; EdotJ?perp.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)


c_eval('RedotJ? = gsmEVexB?.resample(gsmJ?).dot(gsmJ?); RedotJ?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% E?J lmn
c_eval('mvaEdotJ?l = mvaE?l.resample(mvaJ?l)*(mvaJ?l); mvaEdotJ?l.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('mvaEdotJ?m = mvaE?m.resample(mvaJ?m)*(mvaJ?m); mvaEdotJ?m.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('mvaEdotJ?n = mvaE?n.resample(mvaJ?n)*(mvaJ?n); mvaEdotJ?n.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% c_eval('mvaEdotJ?par = mvaE?par.resample(mvaJ?par)*mvaJ?par/1000; mvaEdotJ?par.units = ''nW/m^3''; mvaEdotJ?par.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% c_eval('mvaEdotJ?perp = mvaE?perp.resample(mvaJ?perp).dot(mvaJ?perp)/1000; mvaEdotJ?perp.units = ''nW/m^3''; mvaEdotJ?perp.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

c_eval('mvaEdotJ?curl_l = mvaE?l.resample(mvaJcurl_l)*(mvaJcurl_l); mvaEdotJ?curl_l.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('mvaEdotJ?curl_m = mvaE?m.resample(mvaJcurl_m)*(mvaJcurl_m); mvaEdotJ?curl_m.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('mvaEdotJ?curl_n = mvaE?l.resample(mvaJcurl_n)*(mvaJcurl_n); mvaEdotJ?curl_n.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
% Calculate epsilon and delta parameters
c_eval('oce? = fce?*2*pi;',ic)
c_eval('EdotVe? = gsmE?.resample(gsmVe?).dot(gsmVe?);',ic);
c_eval('epsilone? = abs(6*pi*EdotVe?/(oce?.resample(gsmVe?)*(facTe?.trace)));',ic);
%% Calculate E - JxB/ne & (E - JxBohm).J

c_eval('Residual? = gsmE?.resample(gsmJ?) - gsmOhmJxB?.resample(gsmJ?) + gsmVixB?.resample(gsmJ?) - divPe?ne.resample(gsmJ?); E_JxB?.units = ''mV/m'';',ic)
c_eval('ResidualdotJ?x = Residual?.x*(gsmJ?.x); gsmEdotJ?x.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('ResidualdotJ?y = Residual?.y*(gsmJ?.y); gsmEdotJ?y.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('ResidualdotJ?z = Residual?.z*(gsmJ?.z); gsmEdotJ?z.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('ResidualdotJ? = Residual?.dot(gsmJ?); E_JxBdotJ?.units = ''pW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

% c_eval('gsmE_JxBdotJ? = gsmE_JxBdotJ?x + gsmE_JxBdotJ?y + gsmE_JxBdotJ?z; gsmE_JxBdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
%% Calculate E - Hall - vixB - divPe
c_eval('elec_inertia_0? = gsmE? - gsmOhmJxB? + gsmVixB?.resample(gsmE?); E_JxB?.units = ''mV/m'';',ic)
% % % c_eval('elec_inertia_0x? = elec_inertia_0?.x; E_JxBx?.units = ''mV/m'';',ic)
% % % c_eval('elec_inertia_0y? = elec_inertia_0?.y; E_JxBy?.units = ''mV/m'';',ic)
% % % c_eval('elec_inertia_0z? = elec_inertia_0?.z; E_JxBz?.units = ''mV/m'';',ic)
% % % 
% % % c_eval('elec_inertia_x? = elec_inertia_0x? + divPe?nel.resample(elec_inertia_0x?).data;',ic)
% % % c_eval('elec_inertia_y? = elec_inertia_0y? + divPe?nem.resample(elec_inertia_0y?).data;',ic)
% % % c_eval('elec_inertia_z? = elec_inertia_0z? + divPe?nen.resample(elec_inertia_0z?).data;',ic)




c_eval('deltae? = gsmVexB?/(gsmVe?perp.abs*gsmB?.resample(gsmVe?).abs*1e-9);',ic);
%c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

% Plasma beta and magnetic pressure
%c_eval('beta?_ = (re?/Le?).^2;',ic) % this is beta_e

% c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

end
% c_eval('wavVe?par = irf_wavelet(gsmVe?par.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
% c_eval('wavVe?par.f_units = ''Hz''; wavVe?par.f_label = ''f [Hz]''; wavVe?par.p_label = {''log_{10} v_{e,||}^2'',''(km/s)^2/Hz''};',ic)
% c_eval('wavVe?perp = irf_wavelet(gsmVe?perp.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
% c_eval('wavVe?perp.f_units = ''Hz''; wavVe?perp.f_label = ''f [Hz]''; wavVe?perp.p_label = {''log_{10} v_{e,\perp}^2'',''(km/s)^2/Hz''};',ic)

%% Assume normal electric fields are directly proportional to ion pressure gradient at boundary
if 0
c_eval('gsmGradPi?_fromE = gsmE?.resample(ne?)*ne?*units.e*1e-3*1e6*1e9*1e3; gsmGradPi?_fromE.units = ''nPa/km'';',ic)
% c_eval('Pi?perp = irf.ts_scalar(facPi?.time,(facPi?.yy.data+facPi?.zz.data)/2);',ic)
c_eval('LP?= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp.abs/units.e*1e-9*1e-6*1e3;',ic)
c_eval('gsmE?perp_filt = gsmE?perp.filt(0,3,[],3);',ic)
c_eval('LP?_filt= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp_filt.resample(Pi?perp).abs/units.e*1e-9*1e-6*1e3;',ic)
end
%%
disp('Done preparing data. Not MVA system.')
