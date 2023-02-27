units = irf_units;
B = 0.05*1e9;%nT
Energy_e = 5000; %eV
v_e = sqrt(2*Energy_e*units.e/units.me); %m/s
R_e = units.me * v_e /(1e3*units.e*B*1e-9) %km
Energy_p = 1000;
v_p = sqrt(2*Energy_p*units.e/units.mp); %m/s
R_p = units.mp * v_p /(1e3*units.e*B*1e-9) %km