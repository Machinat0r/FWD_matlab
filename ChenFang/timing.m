clc
clear
ic = 1:4;
mms.db_init('local_file_db','D:\Works\mms_db\data');
db_info = datastore('mms_db');   

% tint = irf.tint('2017-06-11T17:55:17.40Z/2017-06-11T17:55:18.40Z');
tint1 = irf.tint('2017-07-17T07:53:00.00Z/2017-07-17T07:53:20.00Z');
%  tint = irf.tint('2017-07-17T07:53:05.00Z/2017-07-17T07:53:07.00Z');
% tint = irf.tint('2017-07-17T07:53:07.00Z/2017-07-17T07:53:08.50Z');
%tint = irf.tint('2017-07-17T07:53:13.60Z/2017-07-17T07:53:14.60Z');
% tint = irf.tint('2017-07-17T07:53:14.70Z/2017-07-17T07:53:15.00Z');
tint = irf.tint('2017-07-17T07:53:21.00Z/2017-07-17T07:53:23.00Z');

disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval(['gsmE?=irf_gse2gsm(gseE?);'],ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);


R  = mms.get_data('R_gsm',tint1);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?);',1:4);
c_eval('R? = Rxyz?.resample(gsmE1);',1:4);
c_eval('R? = R?.tlim(tint);',1:4);
figure('name','timing');
np=1;
zoomy = [];
h=irf_plot(np);
isub=1;
hca=h(isub); isub=isub+1;
    zoomy = [zoomy isub];
irf_4_v_gui(gsmE1,gsmE2,gsmE3,gsmE4,R1,R2,R3,R4,'mms');