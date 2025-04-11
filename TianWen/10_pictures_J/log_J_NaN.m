clear; clc; close all;

gridNuma = 360; gridNume = 180; gridNumh = 5; 

load(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '.mat']);
load(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_day', '.mat']);
load(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_nig', '.mat']);

%% NaN
inda = find((Ja < 1 & Ja >= 0 & ~isnan(Ja)) | (Ja > -1 & Ja < 0 & ~isnan(Ja)));
inda_day = find((Ja_day < 1 & Ja_day >= 0 & ~isnan(Ja_day)) | (Ja_day > -1 & Ja_day < 0 & ~isnan(Ja_day)));
inda_nig = find((Ja_nig < 1 & Ja_nig >= 0 & ~isnan(Ja_nig)) | (Ja_nig > -1 & Ja_nig < 0 & ~isnan(Ja_nig)));

inde = find((Je < 1 & Je >= 0 & ~isnan(Je)) | (Je > -1 & Je < 0 & ~isnan(Je)));
inde_day = find((Je_day < 1 & Je_day >= 0 & ~isnan(Je_day)) | (Je_day > -1 & Je_day < 0 & ~isnan(Je_day)));
inde_nig = find((Je_nig < 1 & Je_nig >= 0 & ~isnan(Je_nig)) | (Je_nig > -1 & Je_nig < 0 & ~isnan(Je_nig)));

indh = find((Jh < 1 & Jh >= 0 & ~isnan(Jh)) | (Jh > -1 & Jh < 0 & ~isnan(Jh)));
indh_day = find((Jh_day < 1 & Jh_day >= 0 & ~isnan(Jh_day)) | (Jh_day > -1 & Jh_day < 0 & ~isnan(Jh_day)));
indh_nig = find((Jh_nig < 1 & Jh_nig >= 0 & ~isnan(Jh_nig)) | (Jh_nig > -1 & Jh_nig < 0 & ~isnan(Jh_nig)));

Ja(inda) = NaN; Ja_day(inda_day) = NaN; Ja_nig(inda_nig) = NaN;
Je(inde) = NaN; Je_day(inde_day) = NaN; Je_nig(inde_nig) = NaN;
Jh(indh) = NaN; Jh_day(indh_day) = NaN; Jh_nig(indh_nig) = NaN;

%% log
loga = log10(abs(Ja)); loga_day = log10(abs(Ja_day)); loga_nig = log10(abs(Ja_nig));
loge = log10(abs(Je)); loge_day = log10(abs(Je_day)); loge_nig = log10(abs(Je_nig));
logh = log10(abs(Jh)); logh_day = log10(abs(Jh_day)); logh_nig = log10(abs(Jh_nig));

logJa = sign(Ja).*loga; logJa_day = sign(Ja_day).*loga_day; logJa_nig = sign(Ja_nig).*loga_nig;
logJe = sign(Je).*loge; logJe_day = sign(Je_day).*loge_day; logJe_nig = sign(Je_nig).*loge_nig;
logJh = sign(Jh).*logh; logJh_day = sign(Jh_day).*logh_day; logJh_nig = sign(Jh_nig).*logh_nig;

save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J_log_NaN\Ja' ...
   , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '.mat'], ...
   "logJa", "logJe", "logJh");
save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J_log_NaN\Ja' ...
   , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_day', '.mat'], ...
   "logJa_day", "logJe_day", "logJh_day");
save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\J_log_NaN\Ja' ...
   , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_nig', '.mat'], ...
   "logJa_nig", "logJe_nig", "logJh_nig");