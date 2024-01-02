clc; close all

% 连续小波变换
figure
fs = 8192;
wave = Bt1(:,2);

wavename='amor';
[cfs, frq] = cwt(wave,wavename,fs,FrequencyLimits=[500 3e3]);
tms = (0:numel(wave)-1)/fs;
surface(tms,frq,log10(abs(cfs)))
shading flat
set(gca,'yscale','log');
colormap jet
colorbar
clim([-4,-1.5])