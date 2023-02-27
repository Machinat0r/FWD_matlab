function ang_tilt= dipole_tilt_get(T_iso)

% the result unit is degeree

mo=str2num(substr(T_iso,6,2)); dy=str2num(substr(T_iso,9,2));
hr=str2num(substr(T_iso,12,2)); mi=str2num(substr(T_iso,15,2));  

dayspermonth = [31 28 31 30 31 30 31 31 30 31 30 31];
doy = sum(dayspermonth(1:mo-1))+dy;
ut=hr+mi/60;


% Reference: Nowada et al. (2009), Effects of dipole tilt angle on geomagnetic activity
Phi_yr=23.4*cos((doy-172)*(2*pi/365.25));
Phi_dy=11.2*cos((ut-16.72)*(2*pi/24));

Phi_tilt=Phi_yr+Phi_dy;


ang_tilt=Phi_tilt;