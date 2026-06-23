% electron_vorticity_main
close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


TT = '2020-08-03T01:45:23.000Z/2020-08-03T01:45:38.000Z'; % case 16-long
% TT = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z'; 


tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:3;
iic = 1:3;
% % % try
% % % filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% % % filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% % % filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% % % filenames4 = SDCFilenames(Date,iic,'inst','edp','drm','brst','dpt','dce');
% % % % filenames5 = SDCFilenames(Date,ic,'inst','feeps','drm','brst','dpt','electron');
% % % % filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% % % % filenames_fast1 = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms');
% % % % filenames_fast2 = SDCFilenames(Date,ic,'inst','edp','drm','fast');
% % % filenames = [filenames1, filenames2, filenames3, filenames4];
% % % % filenames_fast = [filenames_fast1, filenames_fast2];
% % % % % % 
% % % [filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% % % % [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% % % % [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);
% % % 
% % % SDCFilesDownload_NAS(filenames,TempDir, 'CheckSize', 0, 'Threads', 16)
% % % % SDCFilesDownload(filenames,TempDir)
% % % % % % 
% % % % SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 64, 'CheckSize', 0)
% % % % SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 64, 'CheckSize', 0)
% % % % % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
% % % catch
% % %     warning('no files have been downloaded')
% % % end
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
c_eval(['Bt?_ts=B?_ts.abs;'],iic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],iic);
 % c_eval(['B?=irf_gse2gsm(B?);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],iic);


% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);

c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

c_eval('E?_err_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);
c_eval('E?_err=irf.ts2mat(E?_err_ts);',ic);

c_eval(['Vexb?=irf_cross(E?,B?);'],ic);
c_eval(['Vexb?(:,2:4)=1e3*Vexb?(:,2:4)./[Bt?_res(:,2).^2 Bt?_res(:,2).^2 Bt?_res(:,2).^2];'],ic);%km/s


% load FPI
c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);

c_eval('Ve?_ts = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
c_eval(['gsmVe?_ts=irf_gse2gsm(Ve?_ts);'],ic);
c_eval(['gsmVe?=irf.ts2mat(gsmVe?_ts);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);


c_eval('Pe?_ts = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic)
c_eval(["Pe?=[double(irf_time(Pe1_ts.time.epoch,'ttns>epoch')) double(Pe?_ts.data(:,1,1))" + ...
    " double(Pe?_ts.data(:,2,2)) double(Pe?_ts.data(:,3,3))];"],ic);
c_eval(['gsmPe?=irf_gse2gsm(Pe?);'],ic);


% R
units = irf_units;
Pos = mms.get_data('R_gsm',tint);
c_eval('R? = Pos.gsmR?;',iic)
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];',iic)
c_eval('R? = irf_resamp(R?,Ne1);',iic)
%%
Lvec = [0.23, 0.96, 0.14];
Mvec = [0.37, 0.04, -0.93];
Nvec = [-0.90, 0.26, -0.35];

out = calc_electron_vorticity_phi_3sc( ...
    R1, R2, R3, ...
    Ne1, Ne2, Ne3, ...
    gsmPe1, gsmPe2, gsmPe3, ...
    gsmVe1, gsmVe2, gsmVe3, ...
    B1, B2, B3, ...
    Lvec, Mvec, Nvec);
%%
figure('Color','w');
hold on;

t = out.t(:);

baro = abs(out.baro_phi(:));
em   = abs(out.eMag_phi(:));

% plot(out.t, out.baro_phi, 'LineWidth', 1.5);
% plot(out.t, out.eMag_phi, 'LineWidth', 1.5);

t_fine = linspace(min(t), max(t), size(B1,1)).';
baro_fine = 10.^interp1(t, log10(baro), t_fine, 'makima');
em_fine   = 10.^interp1(t, log10(em), t_fine, 'makima');

plot(t_fine, baro_fine, 'LineWidth', 1.8);
plot(t_fine, em_fine,   'LineWidth', 1.8);


legend('Baroclinic term', 'Electron EM term');
% set(gca, 'YScale', 'log');
xlabel('Time');
ylabel('Angular acceleration, s^{-2}');
irf_timeaxis(gca);
grid off;