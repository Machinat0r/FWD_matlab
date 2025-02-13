clear;clc;close all
IMFdata = importdata('/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/Sundry/IMF_20151018.dat');
IMFdata = IMFdata.data;
%% Download Bow Shock Location
TT = '2015-10-18T00:00:00.00Z/2015-10-18T12:11:00.00Z';
tint = irf.tint(TT);
BSNx= irf_get_data_omni_modified(tint,'bsnx','omni_min');
% BSNy= irf_get_data_omni_modified(tint,'bsny','omni_min');
% BSNz= irf_get_data_omni_modified(tint,'bsnz','omni_min');
BSNx = inpaint_nans(BSNx);
%% Time Shifting (only X direction)
units = irf_units;
dx = BSNx(:,2) - 32; % shift to 32Re
dx = dx * units.RE / 1000;
dt = dx ./ IMFdata(:,11);
Tshifted = BSNx(:,1) - dt;
Tvector = irf_time(Tshifted,'epoch>vector9');
IMFdata(:,1:7) = Tvector(:,1:7);
%% Save Data
save('./IMF_20151018_shifted.dat', 'IMFdata', '-mat')