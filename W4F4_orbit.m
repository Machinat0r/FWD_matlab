clear;clc

global ParentDir 
ParentDir = 'Z:/Data/MMS/'; 
% TempDir = 'E:\MMS\temp\';mkdir(TempDir);
% TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% TT = '2017-08-07T16:31:30.00Z/2017-08-07T16:33:00.00Z';
% TT = '2017-07-06T17:31:30.000Z/2017-07-06T17:32:30.000Z';
% TT = '2017-08-20T02:01:00.00Z/2017-08-20T02:03:00.00Z';
TT = '2019-01-16T04:09:00.000Z/2019-01-16T04:10:00.000Z';
tint=irf.tint(TT);

%%
mms.db_init('local_file_db',ParentDir);
% h = irf_gse2gsm(mms.mms4_pl_conf(tint));
h = mms.mms4_pl_conf(tint);
h = mms.mms4_pl_conf('gse');
% h = mms.mms4_pl_conf('config3D');
% h = mms.mms4_pl_eb(tint)
%% Neutral Sheet
% x = -25:0.1:0;
% tilt = 25.474825;%07
% % tilt = 5.4347405;%20
% h0 = 12.6/pi; 
% y=0;
% dz2NS = -h0 * sind(tilt) * atan(x/5) * (2 * cos(y/6));
% plot(h(1),x,dz2NS,'k-.','linewidth',0.5)
%% 
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

figname=['Relative_location'];
% print(gcf, '-dpdf', [figname '.pdf']);