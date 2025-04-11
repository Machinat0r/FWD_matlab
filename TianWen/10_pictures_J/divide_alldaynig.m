clear; clc;
load('E:\Martian\programs\2_Maven_Download\5_2_T_BRsphere_za_year/T_BRsphere_za_year_15_23.mat');

% for x=1:length(T_BRsphere_za_year)
%     if T_BRsphere_za_year(x,8) >= 120
%         T_BRsphere_za_year_day = [T_BRsphere_za_year_day;T_BRsphere_za_year(x,:)];
%     end
%     if T_BRsphere_za_year(x,8) <= 80
%         T_BRsphere_za_year_nig = [T_BRsphere_za_year_nig;T_BRsphere_za_year(x,:)];
%     end
% end
x = 1:length(T_BRsphere_za_year);
nig = find(T_BRsphere_za_year(x,8) >= 120);
T_BRsphere_za_year_nig = T_BRsphere_za_year(nig, :);

day = find(T_BRsphere_za_year(x,8) <= 80);
T_BRsphere_za_year_day = T_BRsphere_za_year(day, :);

savename_nig='E:\Martian\programs\2_Maven_Download\10_pictures_J\T_BRsphere_za_year_15_23_nig.mat';
save(savename_nig,'T_BRsphere_za_year_nig');

savename_day='E:\Martian\programs\2_Maven_Download\10_pictures_J\T_BRsphere_za_year_15_23_day.mat';
save(savename_day,'T_BRsphere_za_year_day');

savename_all='E:\Martian\programs\2_Maven_Download\10_pictures_J\T_BRsphere_za_year_15_23.mat';
save(savename_all,'T_BRsphere_za_year');