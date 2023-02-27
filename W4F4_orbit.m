clear;clc

global ParentDir 
ParentDir = 'E:\MMS\'; 
TempDir = 'E:\MMS\temp\';mkdir(TempDir);
% TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% TT = '2017-08-07T16:31:30.00Z/2017-08-07T16:33:00.00Z';
TT = '2019-01-16T04:09:55.00Z/2019-01-16T04:09:56.00Z';
% TT = '2017-08-20T02:01:00.00Z/2017-08-20T02:03:00.00Z';
tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
if Datelist{1} == Datelist{2} %防止有跨天的情况出现
    Date = [Datelist{1},'/',Datelist{1}];
else
    Date = [Datelist{1},'/',Datelist{2}];
end
filenames_srvy = SDCFilenames(Date,1:4,'inst','fgm','drm','srvy'); %坐标存在fgm的srvy数据里
SDCFilesDownload(filenames_srvy,TempDir)
SDCDataMove(TempDir,ParentDir)
cd  C:\Matlab\bin\新建文件夹\fwd\  %改一个任意文件夹都可，否则会出现该文件夹正在被占用无法删除的情况
rmdir(TempDir,'s');
%%
mms.db_init('local_file_db','E:\MMS\');
% h = irf_gse2gsm(mms.mms4_pl_conf(tint));
h = mms.mms4_pl_conf(tint);
h = mms.mms4_pl_conf('gsm');
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