close all
clear;clc

global ParentDir 
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

TT = '2015-11-19T13:40:00.000Z/2015-11-19T13:41:00.000Z';
% TT = '2015-11-19T14:20:00.000Z/2015-11-19T14:21:00.000Z';
% TT = '2015-11-19T13:10:00.000Z/2015-11-19T13:12:00.000Z';
% TT = '2015-11-19T14:07:30.000Z/2015-11-19T14:09:30.000Z';
% TT = '2015-11-19T14:02:30.000Z/2015-11-19T14:07:30.000Z';
% TT = '2019-03-29T13:23:00.000Z/2019-03-29T13:23:30.000Z';
% TT = '2017-07-18T13:04:50.000Z/2017-07-18T13:05:00.000Z';
% TT = '2015-11-19T13:20:00.000Z/2015-11-19T13:25:00.000Z';

tint=irf.tint(TT);
mms.db_init('local_file_db',ParentDir);
% h = irf_gse2gsm(mms.mms4_pl_conf(tint));
h = mms.mms4_pl_conf(tint);
h = mms.mms4_pl_conf('gse');

% Bx_omni= irf_get_data_omni_modified(tint,'Bx','omni_min');
% By_omni= irf_get_data_omni_modified(tint,'ByGSE','omni_min');
% Bz_omni= irf_get_data_omni_modified(tint,'BzGSE','omni_min');
% B_omni = [Bx_omni, By_omni(:,2), Bz_omni(:,2)];