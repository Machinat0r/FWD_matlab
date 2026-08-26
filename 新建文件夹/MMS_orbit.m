clear
clc

tint=irf.tint('2021-07-22T11:16:00.000Z/2021-07-22T11:24:00.000Z');

mms.db_init('local_file_db','/Volumes/100.114.210.8/SPART-WORK/Data/MMS/');
h = mms.mms4_pl_conf(tint);
% h = mms.mms4_pl_conf('config3D')
h = mms.mms4_pl_conf('gsm')
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['Relative_location'];
% print(gcf, '-dpdf', [figname '.pdf']);