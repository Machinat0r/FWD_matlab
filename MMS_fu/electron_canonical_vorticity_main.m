%% ================================================================
% Electron canonical vorticity diagnostic using 3 MMS spacecraft
% Local coordinates:
%   er   = N
%   ephi = L
%   ez   = M
%
% Assumption:
%   d/dphi = 0
%
% Main quantities:
%   omega_e = curl(Ve)
%   W_e     = omega_e + (qe/me)*B
%
% Canonical vorticity equation:
%   dW_e/dt = curl(Ve x W_e) + (grad ne x grad pe)/(me ne^2)
%%
close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


% TT = '2020-08-03T01:45:23.000Z/2020-08-03T01:45:38.000Z'; % case 16-long
TT = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z'; 


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
%% ================================================================

Lvec = [0.23, 0.96, 0.14];
Mvec = [0.37, 0.04, -0.93];
Nvec = [-0.90, 0.26, -0.35];

%% ================================================================
% Estimate Omega0 from electron azimuthal flow
% using FOTE estimated distance to the structure center
%
% Omega0(t) = Ve_phi(t) / r
%
% Local coordinates:
%   er   = N
%   ephi = L
%   ez   = M
%
% Units:
%   Ve_phi : km/s
%   r      : km
%   Omega  : s^-1
%% ================================================================

r_center_km = 235.257762572608;

%% construct local orthonormal basis, same as the function
er = Nvec(:) ./ norm(Nvec);

ephi = Lvec(:) - dot(Lvec(:), er).*er;
ephi = ephi ./ norm(ephi);

ez_tmp = cross(er, ephi);
ez_tmp = ez_tmp ./ norm(ez_tmp);

Munit = Mvec(:) ./ norm(Mvec);
if dot(ez_tmp, Munit) < 0
    ez_tmp = -ez_tmp;
end
ez = ez_tmp;

Q = [er, ephi, ez];

%% choose reference time
% Use MMS1 electron velocity time as the reference
t_omega = gsmVe1(:,1);

%% resample electron velocities to the same time
gsmVe1_omega = gsmVe1;

if exist('gsmVe2','var')
    gsmVe2_omega = irf_resamp(gsmVe2, gsmVe1_omega);
end

if exist('gsmVe3','var')
    gsmVe3_omega = irf_resamp(gsmVe3, gsmVe1_omega);
end

if exist('gsmVe4','var')
    gsmVe4_omega = irf_resamp(gsmVe4, gsmVe1_omega);
end

%% project GSM electron velocity into local coordinates
Ve1_local = gsmVe1_omega(:,2:4) * Q;
Vephi1 = Ve1_local(:,2);

Omega1 = Vephi1 ./ r_center_km;

Omega_all = Omega1;

if exist('gsmVe2_omega','var')
    Ve2_local = gsmVe2_omega(:,2:4) * Q;
    Vephi2 = Ve2_local(:,2);
    Omega2 = Vephi2 ./ r_center_km;
    Omega_all = [Omega_all, Omega2];
end

if exist('gsmVe3_omega','var')
    Ve3_local = gsmVe3_omega(:,2:4) * Q;
    Vephi3 = Ve3_local(:,2);
    Omega3 = Vephi3 ./ r_center_km;
    Omega_all = [Omega_all, Omega3];
end

if exist('gsmVe4_omega','var')
    Ve4_local = gsmVe4_omega(:,2:4) * Q;
    Vephi4 = Ve4_local(:,2);
    Omega4 = Vephi4 ./ r_center_km;
    Omega_all = [Omega_all, Omega4];
end

%% remove bad values and extreme spikes
Omega_abs_all = abs(Omega_all);

idx_good = isfinite(Omega_abs_all);

Omega_vec_signed = Omega_all(idx_good);
Omega_vec_abs    = Omega_abs_all(idx_good);

% optional robust spike removal
% Keep central 98 percent of the absolute Omega distribution
p_low  = prctile(Omega_vec_abs, 1);
p_high = prctile(Omega_vec_abs, 99);

idx_robust = Omega_vec_abs >= p_low & Omega_vec_abs <= p_high;

Omega_vec_signed_robust = Omega_vec_signed(idx_robust);
Omega_vec_abs_robust    = Omega_vec_abs(idx_robust);

%% characteristic Omega0
% signed value preserves rotation direction
Omega0_signed_median = median(Omega_vec_signed_robust, 'omitnan');
Omega0_signed_mean   = mean(Omega_vec_signed_robust, 'omitnan');

% absolute value gives rotation strength
Omega0_abs_median = median(Omega_vec_abs_robust, 'omitnan');
Omega0_abs_mean   = mean(Omega_vec_abs_robust, 'omitnan');
Omega0_rms        = sqrt(mean(Omega_vec_signed_robust.^2, 'omitnan'));

%% recommended choice
% If you want the Coriolis term to preserve rotation direction, use signed median.
% If sign changes are mostly noise, use absolute median.
Omega0 = Omega0_signed_median;

% Alternative:
% Omega0 = Omega0_abs_median;

fprintf('\n===== Omega0 estimated from electron flow =====\n');
fprintf('r_center = %.6f km\n', r_center_km);
fprintf('Omega0 signed median = %.6g s^-1\n', Omega0_signed_median);
fprintf('Omega0 signed mean   = %.6g s^-1\n', Omega0_signed_mean);
fprintf('Omega0 abs median    = %.6g s^-1\n', Omega0_abs_median);
fprintf('Omega0 abs mean      = %.6g s^-1\n', Omega0_abs_mean);
fprintf('Omega0 rms           = %.6g s^-1\n', Omega0_rms);
fprintf('Chosen Omega0        = %.6g s^-1\n', Omega0);

%% save time series for checking
Omega0_ts = [t_omega, Omega_all];

%% plot diagnostic
figure('Color','w','Position',[100 100 900 450]);
hold on;

plot(t_omega, Omega_all, 'LineWidth', 1.0);
yline(Omega0_signed_median, 'k--', 'LineWidth', 1.2);
yline(Omega0_abs_median, 'k:', 'LineWidth', 1.2);

ylabel('\Omega_e = V_{e\phi}/r  (s^{-1})');
xlabel('Time');
title('Estimate of characteristic electron angular velocity');
grid off;

legend_entries = {};
legend_entries{end+1} = '\Omega_e MMS1';

if exist('Omega2','var')
    legend_entries{end+1} = '\Omega_e MMS2';
end

if exist('Omega3','var')
    legend_entries{end+1} = '\Omega_e MMS3';
end

if exist('Omega4','var')
    legend_entries{end+1} = '\Omega_e MMS4';
end

legend_entries{end+1} = 'signed median';
legend_entries{end+1} = 'abs median';

legend(legend_entries, 'Location','best');

irf_timeaxis(gca);
%%
out = calc_electron_canonical_vorticity_phi_3sc( ...
    R1, R2, R3, ...
    Ne1, Ne2, Ne3, ...
    gsmPe1, gsmPe2, gsmPe3, ...
    gsmVe1, gsmVe2, gsmVe3, ...
    B1, B2, B3, ...
    Lvec, Mvec, Nvec, Omega0);

%% ================================================================
% Plot all canonical vorticity panels in one figure
%% ================================================================

t = out.t(:);

figure('Color','w','Position',[100 50 900 1050]);

np = 6;
h = gobjects(np,1);

%% Panel 1: axial shear of azimuthal electron flow
h(1) = subplot(np,1,1);
hold on;

plot(t, out.dVephi_dz, 'LineWidth', 1.4);
plot(t, out.omega_r,   'LineWidth', 1.4);

ylabel('s^{-1}');
legend( ...
    '\partial V_{e\phi}/\partial z', ...
    '\omega_{e,r}', ...
    'Location','best');
title('Axial shear of azimuthal electron flow');
grid off;

%% Panel 2: r component of canonical vorticity
h(2) = subplot(np,1,2);
hold on;

plot(t, out.omega_r, 'LineWidth', 1.4);
plot(t, out.mag_r,   'LineWidth', 1.4);
plot(t, out.W_r,     'LineWidth', 1.6);

ylabel('s^{-1}');
legend( ...
    '\omega_{e,r}', ...
    '(q_e/m_e)B_r', ...
    'W_{e,r}', ...
    'Location','best');
title('r component of electron canonical vorticity');
grid off;

%% Panel 3: phi component of canonical vorticity
h(3) = subplot(np,1,3);
hold on;

plot(t, out.omega_phi, 'LineWidth', 1.4);
plot(t, out.mag_phi,   'LineWidth', 1.4);
plot(t, out.W_phi,     'LineWidth', 1.6);

ylabel('s^{-1}');
legend( ...
    '\omega_{e,\phi}', ...
    '(q_e/m_e)B_\phi', ...
    'W_{e,\phi}', ...
    'Location','best');
title('\phi component of electron canonical vorticity');
grid off;

%% Panel 4: magnetic dominance and compensation diagnostics
h(4) = subplot(np,1,4);
hold on;

plot(t, out.mag_over_omega_phi, 'LineWidth', 1.4);
plot(t, out.comp_phi,           'LineWidth', 1.4);

yline(1,'--','LineWidth',1.0);

ylabel('dimensionless');
legend( ...
    '|(q_e/m_e)B_\phi| / |\omega_{e,\phi}|', ...
    '|W_{e,\phi}| / (|\omega_{e,\phi}|+|(q_e/m_e)B_\phi|)', ...
    '1', ...
    'Location','best');
title('Magnetic dominance and compensation in W_{e,\phi}');
grid off;

%% Panel 5: canonical vorticity equation, phi component
% smoothed in the same way as electron_rotating_vorticity_main
% Note:
%   This function uses the rotating-frame equation:
%   dW_phi/dt = conv_phi + rot_phi + baro_phi + res_phi
%   where
%   conv_phi       = [curl(Ve x W_e)]_phi
%   rot_phi        = [curl(Ve x 2Omega0)]_phi
%   Vexomegae_phi  = [curl(Ve x omega_e)]_phi
h(5) = subplot(np,1,5);
hold on;

y1 = abs(out.dWphi_dt(:));
y2 = abs(out.conv_phi(:));
y3 = abs(out.rot_phi(:));
y4 = abs(out.baro_phi(:));
y5 = abs(out.res_phi(:));

if isfield(out,'Vexomegae_phi')
    y6 = abs(out.Vexomegae_phi(:));
else
    error(['out.Vexomegae_phi does not exist. ', ...
           'Please add Vexomegae_phi output in calc_electron_canonical_vorticity_phi_3sc.']);
end

% avoid log interpolation problems
y1(y1<=0) = NaN;
y2(y2<=0) = NaN;
y3(y3<=0) = NaN;
y4(y4<=0) = NaN;
y5(y5<=0) = NaN;
y6(y6<=0) = NaN;

% use B1 sampling number to make smooth curves, same style as rotating code
t_fine = linspace(min(t), max(t), size(B1,1)).';

% interpolate in log-space
valid1 = isfinite(t) & isfinite(y1);
valid2 = isfinite(t) & isfinite(y2);
valid3 = isfinite(t) & isfinite(y3);
valid4 = isfinite(t) & isfinite(y4);
valid5 = isfinite(t) & isfinite(y5);
valid6 = isfinite(t) & isfinite(y6);

y1_fine = NaN(size(t_fine));
y2_fine = NaN(size(t_fine));
y3_fine = NaN(size(t_fine));
y4_fine = NaN(size(t_fine));
y5_fine = NaN(size(t_fine));
y6_fine = NaN(size(t_fine));

if sum(valid1) >= 2
    y1_fine = 10.^interp1(t(valid1), log10(y1(valid1)), t_fine, 'makima', NaN);
end

if sum(valid2) >= 2
    y2_fine = 10.^interp1(t(valid2), log10(y2(valid2)), t_fine, 'makima', NaN);
end

if sum(valid3) >= 2
    y3_fine = 10.^interp1(t(valid3), log10(y3(valid3)), t_fine, 'makima', NaN);
end

if sum(valid4) >= 2
    y4_fine = 10.^interp1(t(valid4), log10(y4(valid4)), t_fine, 'makima', NaN);
end

if sum(valid5) >= 2
    y5_fine = 10.^interp1(t(valid5), log10(y5(valid5)), t_fine, 'makima', NaN);
end

if sum(valid6) >= 2
    y6_fine = 10.^interp1(t(valid6), log10(y6(valid6)), t_fine, 'makima', NaN);
end

plot(t_fine, y1_fine, 'LineWidth', 1.8);
plot(t_fine, y2_fine, 'LineWidth', 1.8);
plot(t_fine, y3_fine, 'LineWidth', 1.5);
plot(t_fine, y4_fine, 'LineWidth', 1.8);
plot(t_fine, y5_fine, 'LineWidth', 1.2);
plot(t_fine, y6_fine, 'LineWidth', 1.8);

set(gca,'YScale','log');

ylabel('s^{-2}');
legend( ...
    '|\partial W_{e,\phi}/\partial t|', ...
    '|[\nabla\times(V_e\times W_e)]_\phi|', ...
    '|[\nabla\times(V_e\times 2\Omega_0)]_\phi|', ...
    '|Baroclinic_\phi|', ...
    '|Residual|', ...
    '|[\nabla\times(V_e\times\omega_e)]_\phi|', ...
    'Location','best');

title('Electron canonical vorticity equation, \phi component');
grid off;

%% Panel 6: three-spacecraft gradient quality
h(6) = subplot(np,1,6);
hold on;

plot(t, out.condG, 'LineWidth', 1.4);

ylabel('cond(G)');
xlabel('Time');
title('Three-spacecraft 2D gradient quality');
grid off;

%% Format time axis
linkaxes(h,'x');
xlim([min(t) max(t)]);

for ii = 1:np-1
    set(h(ii),'XTickLabel',[]);
end

irf_timeaxis(h(end));

%% Optional: mark event core interval
% coreTint = irf.tint('2020-08-03T01:45:27.800Z/2020-08-03T01:45:28.500Z');
% for ii = 1:np
%     irf_pl_mark(h(ii), coreTint, [0.85 0.85 0.85]);
% end

%% Optional: export figure
% exportgraphics(gcf,'canonical_vorticity_all_panels.png','Resolution',300);