%% MATLAB/IRFU 真实原 CDF 接口测试；输出只写临时日志
root='C:\Users\Administrator\Documents\irfu-matlab-master';
addpath(root,fullfile(root,'irf'),fullfile(root,'plots'),fullfile(root,'mission','mms'), ...
    fullfile(root,'mission','cluster'),fullfile(root,'contrib','nasa_cdf_patch'));
setenv('CDF_LEAPSECONDSTABLE',fullfile(root,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
global MMS_DB; MMS_DB=mms_db; MMS_DB.add_db(mms_local_file_db('Z:\SPART-WORK\Data\MMS\'));
tint=irf.tint('2026-07-20T04:00:00Z/2026-07-20T06:00:00Z');
v=mms.get_data('Vi_gse_fpi_fast_l2',tint,1);disp(v);disp(v.data(1:3,:));
q=mms.db_get_variable('mms1_fpi_fast_l2_dis-moms','mms1_dis_energyspectr_omni_fast',tint);
disp(fieldnames(q));disp(size(q.data));disp(q.UNITS);disp(fieldnames(q.DEPEND_0));disp(size(q.DEPEND_1.data));disp(q.DEPEND_1.data(1,:));
t=mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_temppara_fast',tint);disp(t);
disp('CDF_PROBE_COMPLETE');
