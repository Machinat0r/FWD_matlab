% L=[0.20 0.18 0.96];
% M=[-0.02 -0.98 0.19];
% N=[0.98 -0.06 -0.19];
% for ic=1:4
%     c_eval(['gsmB?=irf_newxyz(gsmB?,M,N,L);'],ic);
% 
%     c_eval(['gsmE?=irf_newxyz(gsmE?,M,N,L);'],ic);
% 
% end

ic = 1:4;
units = irf_units;

%% Current
% Currents from moments, use ne also for Ji 
c_eval('gsmJe? = -units.e*ne?*gsmVe?*1e3*1e6*1e9; gsmJe?.units = ''nA/m^2''; gsmJe?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJi? = units.e*ne?*gsmVi?.resample(ne?.time)*1e3*1e6*1e9; gsmJi?.units = ''nA/m^2''; gsmJi?.coordinateSystem = ''gsm'';',ic);
c_eval('gsmJ? = (gsmJe?+gsmJi?);',ic);

%% Perpendicular and parallel decomposition
% Celocity and current
c_eval('[gsmVe?par,gsmVe?perp] = irf_dec_parperp(gsmB?,gsmVe?); gsmVe?par.name = ''Ve par''; gsmVe?perp.name = ''Ve perp'';',ic)
c_eval('[gsmVi?par,gsmVi?perp] = irf_dec_parperp(gsmB?,gsmVi?); gsmVi?par.name = ''Vi par''; gsmVi?perp.name = ''Vi perp'';',ic)
c_eval('[gsmJ?par,gsmJ?perp] = irf_dec_parperp(gsmB?,gsmJ?); gsmJ?par.name = ''J par''; gsmJ?perp.name = ''J perp'';',ic)

% Electric fields
c_eval('[gsmE?par,gsmE?perp] = irf_dec_parperp(gsmB?,gsmE?); gsmE?par.name = ''E par''; gsmE?perp.name = ''E perp'';',ic)
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
if all(ic==[1:4])
c_eval('R? = gsmR?.resample(gsmB1);',1:4)
c_eval('B? = gsmB?.resample(gsmB1);',1:4)
[gsmCurvB,avB]=c_4_grad('R?','B?','curvature'); gsmCurvB.name = 'curv B'; gsmCurvB.coordinateSystem = 'gsm';
curvBradius = 1/gsmCurvB.abs; curvBradius.name = 'R_c';
end
%% Pitchangle distributions
if 0
%   load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
  %load /Users/Cecilia/Data/MMS/20151112071854_2017-03-11_ePitch15.mat
elseif 1
  %%
  ictmp=ic;
  ic = 1:4;
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
  %c_eval('ePitch?par = ePDist?.pitchangles(dmpaB?,[0 15]);',ic)
  %c_eval('ePitch?perp = ePDist?.pitchangles(dmpaB?,[75 105]);',ic)
  %c_eval('ePitch?apar = ePDist?.pitchangles(dmpaB?,[165 180]);',ic)
  ic = ictmp;
end

%% Calculate some additional parameters, irf_plasma_calc
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
c_eval('PB? = gsmB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

% Magnetic moment
c_eval('mag_mom? = 0.5*units.me*vte?perp.^2*10^6/(gsmB?.abs*1e-9)*1e9;  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)

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

c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)


c_eval('wavVe?par = irf_wavelet(gsmVe?par.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?par.f_units = ''Hz''; wavVe?par.f_label = ''f [Hz]''; wavVe?par.p_label = {''log_{10} v_{e,||}^2'',''(km/s)^2/Hz''};',ic)
c_eval('wavVe?perp = irf_wavelet(gsmVe?perp.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 15],''nf'',100);',ic)
c_eval('wavVe?perp.f_units = ''Hz''; wavVe?perp.f_label = ''f [Hz]''; wavVe?perp.p_label = {''log_{10} v_{e,\perp}^2'',''(km/s)^2/Hz''};',ic)

%% Assume normal electric fields are directly proportional to ion pressure gradient at boundary
c_eval('gsmGradPi?_fromE = gsmE?.resample(ne?)*ne?*units.e*1e-3*1e6*1e9*1e3; gsmGradPi?_fromE.units = ''nPa/km'';',ic)
c_eval('Pi?perp = irf.ts_scalar(facPi?.time,(facPi?.yy.data+facPi?.zz.data)/2);',ic)
c_eval('LP?= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp.abs/units.e*1e-9*1e-6*1e3;',ic)
c_eval('gsmE?perp_filt = gsmE?perp.filt(0,3,[],3);',ic)
c_eval('LP?_filt= Pi?perp/ne?.resample(Pi?perp)/gsmE?perp_filt.resample(Pi?perp).abs/units.e*1e-9*1e-6*1e3;',ic)

%%
disp('Done preparing data. Not MVA system.')
