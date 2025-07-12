clear
clc
%tint=irf.tint('2017-07-12T11:54:33.30Z/2017-07-12T11:54:34.50Z');
% tint=irf.tint('2017-07-17T07:53:06.00Z/2017-07-17T07:53:11.00Z');
% tint=irf.tint('2017-06-11T17:54:00.000Z/2017-06-11T17:57:00.000Z');
% tint=irf.tint('2019-07-29T16:05:50.00Z/2019-07-29T16:07:00.00Z');

tint=irf.tint('2016-01-27T03:03:26.000Z/2016-01-27T03:03:28.000Z');

mms.db_init('local_file_db','/Volumes/SPART-NAS/Data/MMS/'); 
h = mms.mms4_pl_conf(tint);
h = mms.mms4_pl_conf('gse');
% set(h(1),'Xlim',[-30 30]);
% set(h(3),'Xlim',[-30 30]);
set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['locat0624'];
 %print(gcf, '-dpdf', [figname '.pdf']);