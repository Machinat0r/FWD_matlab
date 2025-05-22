clear
clc

tint=irf.tint('2018-07-05T09:15:00.00Z/2015-09-08T11:35:00.00Z');

mms.db_init('local_file_db','/Volumes/172.17.190.41/Data/MMS/');
h = mms.mms4_pl_conf(tint);
% h = mms.mms4_pl_conf('config3D')
h = mms.mms4_pl_conf('gsm')
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['Relative_location'];
% print(gcf, '-dpdf', [figname '.pdf']);