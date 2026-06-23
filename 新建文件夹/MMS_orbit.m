clear
clc

tint=irf.tint('2020-08-03T01:45:23.000Z/2020-08-03T01:45:38.000Z');

mms.db_init('local_file_db','/Volumes/SPART-WORK/Data/MMS/');
h = mms.mms4_pl_conf(tint);
% h = mms.mms4_pl_conf('config3D')
h = mms.mms4_pl_conf('gsm')
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['Relative_location'];
% print(gcf, '-dpdf', [figname '.pdf']);