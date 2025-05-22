clc
clear
ic = 1;
% tint = irf.tint('2020-07-16T18:08:33.00Z/2020-07-16T18:12:53.00Z');

% tint = irf.tint('2021-07-10T12:41:23.00Z/2021-07-10T12:41:50.00Z');
% tint = irf.tint('2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z');
% tint = irf.tint('2017-06-21T23:27:00.00Z/2017-06-21T23:28:00.00Z');
% tint = irf.tint('2017-06-22T00:30:40.00Z/2017-06-22T00:31:30.00Z');
% tint = irf.tint('2017-06-22T02:33:05.00Z/2017-06-22T02:33:23.00Z');
% tint = irf.tint('2017-06-22T03:01:20.00Z/2017-06-22T03:01:43.00Z');
tint = irf.tint('2021-07-22T12:44:30.00Z/2021-07-22T12:45:30.00');
% tint = irf.tint('2017-06-17T04:53:40.00Z/2017-06-17T04:55:00.00Z');
% tint = irf.tint('2017-07-25T22:10:05.50Z/2017-07-25T22:10:06.50Z');
% tint = irf.tint('2017-07-18T13:04:40.00Z/2017-07-18T13:05:20.00Z');
%  tint=irf.tint('2017-07-20T11:43:33.00Z/2017-07-20T11:45:03.00Z');

%% Load datastore
% mms.db_init('local_file_db','G:\data\mms_db\data');
mms.db_init('local_file_db','D:\MMS\');
% mms.db_init('local_file_db','/Users/xuyin/Documents/data/mms_db/data');
db_info = datastore('mms_db');   

%% Load defatt, for coordinate tranformation
if 0 % not nessecary unless E needs to be recalibrated
  disp('Loading defatt...')
  %load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat
  c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
  c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
  c_eval('defatt? = mms_removerepeatpnts(defatt?);',ic)
end

%% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB?s = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint); toc;',ic);

%c_eval('tic; gseB?scm = mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint); toc',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval(['gsmB?scm=irf_gse2gsm(gseB?scm);'],ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval(['gsmE?=irf_gse2gsm(gseE?);'],ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

%% Load spacecraft position
% disp('Loading spacecraft position...')
% 
% c_eval('R?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_r_gse_brst_l2'',tint);',ic);% gse
% c_eval('R?=irf.ts_vec_xyz(R?.time,R?.data(:,1:3));',ic);
% % % R = mms.get_data('R_gse',tint);
% % % 
% % % if size(R.gseR1,2) == 4
% % %   c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,1:3));',ic); % dfg_srvy_l2pre
% % %   c_eval(['gsmR?=irf_gse2gsm(gseR?);'],ic);
% % % else
% % %   c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',ic); % mec
% % %   c_eval(['gsmR?=irf_gse2gsm(gseR?);'],ic);
% % % end

%% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

%% Particle moments
% Skymap distributions
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
energy_e=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint);
energy_i=mms.db_get_variable('mms1_fpi_brst_l2_dis-moms','mms1_dis_energyspectr_omni_brst',tint);

% c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
% c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
% c_eval('tic; iPDist?=iPDist?.tlim(tint); toc',ic)
% c_eval('tic; ePDist?=ePDist?.tlim(tint); toc',ic)

% ePDistErr1=ePDistErr1.tlim(tint);
% Remove all one-count "noise"
%c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
%c_eval('ePDist?.data(ePDist?.data<ePDistErr?.data*1.1) = 0;',ic)
%


%% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
% c_eval(['gsmPe?=irf_gse2gsm(gsePe?);'],ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval(['gsmTe?=irf_gse2gsm(gseTe?);'],ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
% c_eval(['gsmPi?=irf_gse2gsm(gsePi?);'],ic);
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc
% c_eval(['gsmTi?=irf_gse2gsm(gseTi?);'],ic);
%
c_eval('Ti?_para_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
c_eval(['Ti?_para=irf.ts2mat(Ti?_para_ts);'],ic);
c_eval('Ti?_perp_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
c_eval(['Ti?_perp=irf.ts2mat(Ti?_perp_ts);'],ic);
c_eval(['Ti?=[Ti?_para(:,1),(Ti?_para(:,2)+2*Ti?_perp(:,2))/3.0];'],ic);
%
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gsmB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gsmB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gsmB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gsmB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Ion species
disp('Loading HPCA data...')
c_eval('nH? = mms.get_data(''Nhplus_hpca_brst_l2'',tint,?);',ic);
c_eval('nO? = mms.get_data(''Noplus_hpca_brst_l2'',tint,?);',ic);
c_eval('nHe? = mms.get_data(''Nheplus_hpca_brst_l2'',tint,?);',ic);

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval(['gsmVe?=irf_gse2gsm(gseVe?);'],ic);
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval(['gsmVi?=irf_gse2gsm(gseVi?);'],ic);
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

%c_eval('tic; gseVe?fast = mms.get_data(''Ve_gse_fpi_fast_l2'',fastTint,?); toc;',ic)
%c_eval('tic; gseVi?fast = mms.get_data(''Vi_gse_fpi_fast_l2'',fastTint,?); toc;',ic)

disp('Done loading data.');

