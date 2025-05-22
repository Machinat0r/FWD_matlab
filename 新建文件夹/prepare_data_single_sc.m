ic = 1;
units = irf_units;

%% Current
% Currents from moments, use ne also for Ji 
c_eval('gsmJe? = -units.e*ne?*gsmVe?*1e3*1e6*1e9; gsmJe?.units = ''nA/m^2''; gsmJe?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJi? = units.e*ne?*gsmVi?.resample(ne?.time)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJ? = (gsmJe?+gsmJi?);',ic);
 
%% Perpendicular and parallel decomposition
% Celocity and current
c_eval('[gsmVe?par,gsmVe?perp] = irf_dec_parperp(gsmB?,gsmVe?); gsmVe?par.name = ''Ve par''; gsmVe?perp.name = ''Ve perp'';',ic)
c_eval('[gsmVe?fac] = irf_convert_fac(gsmVe?,gsmB?,[1 0 0]); ',ic)

c_eval('[gsmVi?par,gsmVi?perp] = irf_dec_parperp(gsmB?,gsmVi?); gsmVi?par.name = ''Vi par''; gsmVi?perp.name = ''Vi perp'';',ic)

c_eval('[gsmJ?par,gsmJ?perp] = irf_dec_parperp(gsmB?,gsmJ?); gsmJ?par.name = ''J par''; gsmJ?perp.name = ''J perp'';',ic)
c_eval('[gsmJ?fac] = irf_convert_fac(gsmJ?,gsmB?,[1 0 0]); ',ic)
% Electric fields
c_eval('[gsmE?par,gsmE?perp] = irf_dec_parperp(gsmB?,gsmE?); gsmE?par.name = ''E par''; gsmE?perp.name = ''E perp'';',ic)
c_eval('[gsmE?fac] = irf_convert_fac(gsmE?,gsmB?,[1 0 0]); ',ic)
% Wave magnetic field
%gsmB2scm = gsmB2scm{2};
c_eval('[gsmB?scmpar,gsmB?scmperp] = irf_dec_parperp(gsmB?,gsmB?scm); gsmB?scmpar.name = ''B par scm''; gsmB?scmperp.name = ''B perp scm'';',ic)
try
c_eval('[gsmE?hmfepar,gsmE?hmfeperp] = irf_dec_parperp(gsmB?,gsmE?hmfe); gsmE?hmfepar.name = ''E par hmfe''; gsmE?hmfeperp.name = ''E perp hmfe'';',ic)
end

%% Cross products
% ExB drift
c_eval('gsmVExB? = cross(gsmE?.resample(gsmB?.time),gsmB?)/gsmB?.abs/gsmB?.abs*1e3; gsmVExB?.units = '''';',ic) % km/s

% Convective electric fields
c_eval('gsmVexB? = gsmVe?.cross(gsmB?.resample(gsmVe?))*1e-3; gsmVexB?.units = ''mV/m'';',ic)
c_eval('gsmVixB? = gsmVi?.cross(gsmB?.resample(gsmVi?))*1e-3; gsmVixB?.units = ''mV/m'';',ic)

% Non-ideal electric field, E+VexB
c_eval('gsmEVexB? = gsmE?.resample(gsmVexB?.time)+gsmVexB?; gsmEVexB?.name = ''E+VexB'';',ic)

% JxB
c_eval('gsmJxB? = gsmJ?.cross(gsmB?.resample(gsmJ?));',ic)

% Magnetic field curvature 
% if all(ic==[1:4])
% c_eval('R? = gsmR?.resample(gsmB1);',1:4)
% c_eval('B? = gsmB?.resample(gsmB1);',1:4)
% [gsmCurvB,avB]=c_4_grad('R?','B?','curvature'); gsmCurvB.name = 'curv B'; gsmCurvB.coordinateSystem = 'gsm';
% curvBradius = 1/gsmCurvB.abs; curvBradius.name = 'R_c';
% end
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

energyrange=[3000,29000];%set energy range
pl_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange);
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
c_eval(['V1_psd?=irf_gse2gsm(V1_psd0?);'],ic);
c_eval('[gsmV1_psd?par,gsmV1_psd?perp] = irf_dec_parperp(gsmB?,V1_psd?); V1_psd?par.name = ''V1_psd par''; V1_psd?perp.name = ''V1_psd perp'';',ic)
%
energyrange=[300,1300];%set energy range
pl_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange);
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
c_eval('V2_psd0?=pl_moments.V_psd;',ic)
c_eval(['V2_psd?=irf_gse2gsm(V2_psd0?);'],ic);
c_eval('[gsmV2_psd??par,gsmV2_psd?perp] = irf_dec_parperp(gsmB?,V2_psd?); V2_psd?par.name = ''V2_psd par''; V2_psd?perp.name = ''V2_psd perp'';',ic)
if 1
energyrange=[1600,3000];%set energy range
pl_moments= mms.psd_moments(PDist,e_phi,e_theta,esteptable,e_energy0,e_energy1,scPot1,species,'energyrange',energyrange);
c_eval('Tfac= mms.rotate_tensor(pl_moments.T_psd,''fac'',dmpaB?);',ic);
T_scalar=1/3*(Tfac.xx.data+Tfac.yy.data+Tfac.zz.data);
T_scalar_psd=irf.ts_scalar(Tfac.time,T_scalar);
T_para_psd=irf.ts_scalar(Tfac.time,Tfac.xx.data);
T_perp_psd=irf.ts_scalar(Tfac.time,1/2*(Tfac.yy.data+Tfac.zz.data));
c_eval('T3_scalar_psd?=T_scalar_psd;',ic)
c_eval('T3_para_psd?=T_para_psd;',ic)
c_eval('T3_perp_psd?=T_perp_psd;',ic)
c_eval('N3_psd?=pl_moments.n_psd;',ic)
c_eval('V3_psd?=pl_moments.V_psd;',ic)
end
end

%% Pitchangle distributions
if 0
  load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
  %load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
elseif 1
  %%
  ictmp=ic;
  ic = 1;
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
if 0
% Pressure hot
c_eval('PiH?para = irf.ts_scalar(N1_psd?.time,1e6*11604.5*1e9*1.38e-23*N1_psd?.data.*T1_para_psd?.data);',ic)
c_eval('PiH?perp = irf.ts_scalar(N1_psd?.time,1e6*11604.5*1e9*1.38e-23*N1_psd?.data.*T1_perp_psd?.data);',ic)
c_eval('PdH?= irf.ts_scalar(N1_psd?.time,2e-6*double(N1_psd?.data).*V1_psd?.abs2.data);',ic)
c_eval('PiH?=irf_add(1/3,PiH?para,2/3,PiH?perp);',ic);
c_eval('PtmpH?=irf_add(1,PB?,1,PiH?);',ic);
c_eval('PtH?=irf_add(1,PtmpH?,1,PdH?);',ic);
% Pressure cold
c_eval('PiC?para = irf.ts_scalar(N2_psd?.time,1e6*11604.5*1e9*1.38e-23*N2_psd?.data.*T2_para_psd?.data);',ic)
c_eval('PiC?perp = irf.ts_scalar(N2_psd?.time,1e6*11604.5*1e9*1.38e-23*N2_psd?.data.*T2_perp_psd?.data);',ic)
c_eval('PdC?= irf.ts_scalar(N2_psd?.time,2e-6*double(N2_psd?.data).*V2_psd?.abs2.data);',ic)
c_eval('PiC?=irf_add(1/3,PiC?para,2/3,PiC?perp);',ic);
c_eval('PtmpC?=irf_add(1,PB?,1,PiC?);',ic);
c_eval('PtC?=irf_add(1,PtmpC?,1,PdC?);',ic);
end
% Magnetic moment
c_eval('mag_mom? = 0.5*units.me*vte?perp.^2*10^6/(gsmB?.abs*1e-9)*1e9;  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)
end
%% EDR signatures
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
c_eval('EdotJ? = gsmE?.resample(gsmJ?).dot(gsmJ?)/1000; EdotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?par = gsmE?par.resample(gsmJ?par)*gsmJ?par/1000; EdotJ?par.units = ''nW/m^3''; EdotJ?par.name = ''E*J par'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ?perp = gsmE?perp.resample(gsmJ?perp).dot(gsmJ?perp)/1000; EdotJ?perp.units = ''nW/m^3''; EdotJ?perp.name = ''E*J perp'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('RedotJ? = gsmEVexB?.resample(gsmJ?).dot(gsmJ?)/1000; RedotJ?.units = ''nW/m^3'';',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)

% Calculate epsilon and delta parameters
c_eval('oce? = fce?*2*pi;',ic)
c_eval('EdotVe? = gsmE?.resample(gsmVe?).dot(gsmVe?);',ic);
c_eval('epsilone? = abs(6*pi*EdotVe?/(oce?.resample(gsmVe?)*(facTe?.trace)));',ic);

c_eval('deltae? = gsmVexB?/(gsmVe?perp.abs*gsmB?.resample(gsmVe?).abs*1e-9);',ic);
%c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

% Plasma beta and magnetic pressure
%c_eval('beta?_ = (re?/Le?).^2;',ic) % this is beta_e

% c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)


% c_eval('wavVe?par = irf_wavelet(gsmVe?par.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
% c_eval('wavVe?par.f_units = ''Hz''; wavVe?par.f_label = ''f [Hz]''; wavVe?par.p_label = {''log_{10} v_{e,||}^2'',''(km/s)^2/Hz''};',ic)
% c_eval('wavVe?perp = irf_wavelet(gsmVe?perp.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
% c_eval('wavVe?perp.f_units = ''Hz''; wavVe?perp.f_label = ''f [Hz]''; wavVe?perp.p_label = {''log_{10} v_{e,\perp}^2'',''(km/s)^2/Hz''};',ic)

%% Assume normal electric fields are directly proportional to ion pressure gradient at boundary
c_eval('gsmGradPi?_fromE = gsmE?.resample(ne?)*ne?*units.e*1e-3*1e6*1e9*1e3; gsmGradPi?_fromE.units = ''nPa/km'';',ic)
% c_eval('Pi?perp = irf.ts_scalar(facPi?.time,(facPi?.yy.data+facPi?.zz.data)/2);',ic)
c_eval('LP?= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp.abs/units.e*1e-9*1e-6*1e3;',ic)
c_eval('gsmE?perp_filt = gsmE?perp.filt(0,3,[],3);',ic)
c_eval('LP?_filt= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp_filt.resample(Pi?perp).abs/units.e*1e-9*1e-6*1e3;',ic)

%%
disp('Done preparing data. Not MVA system.')
