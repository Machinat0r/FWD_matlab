clear;
clc;
mms.db_init('local_file_db','D:\matlab2019a\MMS\');
Units=irf_units;
me=Units.me;
tint=irf.tint('2018-07-06T12:37:40Z/2018-07-06T12:39:00Z');

tint1=irf.tint('2018-07-06T12:37:40Z/2018-07-06T12:38:00Z');
R = mms.get_data('R_gsm',tint1);
c_eval('R?_ts = irf.ts_vec_xyz(R.time,R.gsmR?);');

for ic=1:4
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint1,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['R?_ts=R?_ts.resample(B?_ts);'],ic); 
c_eval(['R?_ts=R?_ts.tlim(tint);'],ic); 
c_eval(['B?_ts=B?_ts.tlim(tint);'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);  
c_eval(['R?=irf.ts2mat(R?_ts);'],ic);

c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
c_eval('Bbf? = B?_ts.filt(0.8,1.2,dfB?,5);',ic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);
end

irf_4_v_gui(Bbf1,Bbf2,Bbf3,Bbf4,R1,R2,R3,R4,'mms');