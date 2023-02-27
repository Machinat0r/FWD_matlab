clear;clc;
V = spdfcdfread('E:\OMNI\omni2_h0_mrg1hr_20210101_v01.cdf','Variables', 'V');
V_lon = spdfcdfread('E:\OMNI\omni2_h0_mrg1hr_20210101_v01.cdf','Variables', 'PHI-V');
V_lat = spdfcdfread('E:\OMNI\omni2_h0_mrg1hr_20210101_v01.cdf','Variables','THETA-V');
B_lon = spdfcdfread('E:\OMNI\omni2_h0_mrg1hr_20210101_v01.cdf','Variables', 'PHI_AV');
B_lat = spdfcdfread('E:\OMNI\omni2_h0_mrg1hr_20210101_v01.cdf','Variables','THETA_AV');

V_lat(V_lat>180)=0;
V_lon(V_lon>180)=0;
B_lat(B_lat>180)=0;
B_lon(B_lon>180)=0;

subplot(5,1,1)
plot(V)
subplot(5,1,2)
plot(V_lat)
subplot(5,1,3)
plot(V_lon)
subplot(5,1,4)
plot(B_lat)
subplot(5,1,5)
plot(B_lon)
