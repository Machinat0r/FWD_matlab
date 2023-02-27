clear;clc;close all
Eion = [762.5 1561.9 2957 5290 7240 9560 12060 14580 22540 25290 28000 31920 34830 ...
    37840 44100 47206 122200 131000 140500 152600 163000 173600 188100 195200 851800 895161];
units = irf_units;
Eioni = Eion * 1000 * 1e-7 / (6.02*1e23);%单位：erg
T = 1e6;
n = 1.3*1e22;% 假定气态铁被束缚在一个标准大气压下
%  由于精度不够，取自然对数进行计算
for i = 1:length(Eioni)
    temp =  2.4 * 1e21 * T^(1.5) * exp(-Eioni(i)/(units.kB*T));
    n(end+1) = 0.5*(-temp+sqrt(temp^2+4*temp*n(end)));
end
rate = n(2:end)/n(1);
plot(rate)
ylabel('Ionization Rate','fontsize',12);
xlabel('Ionization state','fontsize',12);

E = sum(Eion.* rate);