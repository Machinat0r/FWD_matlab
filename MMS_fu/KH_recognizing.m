%% ================= KH recognizing (LMN Ve, based on GSM data) =================
clear; clc; close all

%% ---------------- Paths ----------------
global ParentDir
ParentDir   = '/Volumes/SPART-WORK/Data/MMS/';
DownloadDir = '/Volumes/SPART-WORK/Data/MMS/';
TempDir     = [DownloadDir,'temp\'];
if ~exist(TempDir,'dir'); mkdir(TempDir); end

SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

%% ---------------- Time ----------------
% TT   = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z';
TT = '2018-07-03T17:07:20.000Z/2018-07-03T17:07:35.000Z';
tint = irf.tint(TT);
ic = 1;   % MMS1

%% ================= User-defined LMN =================
% L_hat = [ 0.23 0.96 0.14];   % 最大变化方向
% N_hat = [-0.90 0.26 -0.35];   % 外法向
% M_hat = cross(N_hat,L_hat);
L_hat = [ 0.28 0.7 0.66];   % 最大变化方向
N_hat = [ -0.77 0.57 -0.29];   % 外法向
M_hat = cross(N_hat,L_hat);

% Transformation matrix:
% GSM vector V_gsm -> LMN vector V_lmn = [dot(V,L), dot(V,M), dot(V,N)]
T_gsm2lmn = [L_hat; M_hat; N_hat];

disp('LMN basis vectors:')
disp(['L = [',num2str(L_hat,'% .6f '),']'])
disp(['M = [',num2str(M_hat,'% .6f '),']'])
disp(['N = [',num2str(N_hat,'% .6f '),']'])

%% ================= B (GSM) =================
% !!!!! DO NOT MODIFY ORIGINAL DATA READING PART !!!!!
B_ts  = mms.get_data('B_gsm_brst', tint, ic);
Bt_ts = B_ts.abs;

B  = irf.ts2mat(B_ts);    % [t Bx By Bz] GSM, nT
Bt = irf.ts2mat(Bt_ts);   % [t |B|] nT

%% ================= Electron velocity (GSE -> GSM) =================
Ve_gse_ts = mms.get_data('Ve_gse_fpi_brst_l2', tint, ic);
Vi_gse_ts = mms.get_data('Vi_gse_fpi_brst_l2', tint, ic);

% coordinate transform ONLY: GSE -> GSM
Ve_gsm_ts = irf_gse2gsm(Ve_gse_ts);
Ve = irf.ts2mat(Ve_gsm_ts);   % [t Vx Vy Vz] km/s (GSM)

Vi_gsm_ts = irf_gse2gsm(Vi_gse_ts);
Vi = irf.ts2mat(Vi_gsm_ts);   % [t Vx Vy Vz] km/s (GSM)

%% ================= Density =================
Ne_ts = mms.get_data('Ne_fpi_brst_l2', tint, ic);
Ne = irf.ts2mat(Ne_ts);       % [t Ne] cm^-3

%% ================= Resample to B timeline =================
Ve = irf_resamp(Ve,B);
Ne = irf_resamp(Ne,B);

ne = Ne(:,2) * 1e6;   % cm^-3 -> m^-3

%% ================= Convert GSM -> LMN =================
% B in LMN
B_gsm_vec = B(:,2:4);              % nT
B_lmn_vec = (T_gsm2lmn * B_gsm_vec')';   % nT

% Ve in LMN
Ve_gsm_vec = Ve(:,2:4);                 % km/s
Ve_lmn_vec = (T_gsm2lmn * Ve_gsm_vec')'; % km/s

Vi_gsm_vec = Vi(:,2:4);                 % km/s
Vi_lmn_vec = (T_gsm2lmn * Vi_gsm_vec')'; % km/s

% Rebuild time series matrices
B_lmn  = [B(:,1)  B_lmn_vec];
Ve_lmn = [Ve(:,1) Ve_lmn_vec];
Vi_lmn = [Vi(:,1) Vi_lmn_vec];

Vi_lmn = irf_resamp(Vi_lmn, Ve_lmn);
%% ================= Physical constants =================
mu0 = 4*pi*1e-7;
me  = 9.10938356e-31;

%% ================= Electron Alfven speed (LMN components) =================
% component form: C_Ae,i = |B_i| / sqrt(mu0 * me * ne)
CAe_lmn = abs(B_lmn(:,2:4))*1e-9 ./ sqrt(mu0*me.*ne);   % m/s
CAe_lmn = CAe_lmn / 1e3;                                % km/s

CAe_L = CAe_lmn(:,1);
CAe_M = CAe_lmn(:,2);
CAe_N = CAe_lmn(:,3);

%% ================= Extract LMN components =================
B_L = B_lmn(:,2);
B_M = B_lmn(:,3);
B_N = B_lmn(:,4);

Ve_L = Ve_lmn(:,2);
Ve_M = Ve_lmn(:,3);
Ve_N = Ve_lmn(:,4);

Vi_L = Vi_lmn(:,2);
Vi_M = Vi_lmn(:,3);
Vi_N = Vi_lmn(:,4);

%% ================= Plot =================
figure('color','w','position',[200 80 900 950])
n = 4; i = 1;

%% ---- (a) B (LMN) ----
h(i) = irf_subplot(n,1,-i);
irf_plot([Bt(:,1) Bt(:,2)], 'k','LineWidth',0.9); hold on
irf_plot([B_lmn(:,1) B_L], 'b');
irf_plot([B_lmn(:,1) B_M], 'g');
irf_plot([B_lmn(:,1) B_N], 'r');
hold off
ylabel('B [nT]')
irf_legend(gca,{'|B|','B_L','B_M','B_N'},[0.97 0.92])
title(['2020-08-03   MMS',num2str(ic),'   LMN'])

i = i + 1;

%% ---- (b) Ve_L  CAe_L ----
h(i) = irf_subplot(n,1,-i);
irf_plot([Ve_lmn(:,1) Ve_L - Vi_L], 'b'); hold on
irf_plot([Ve_lmn(:,1) CAe_L/2], 'k','LineWidth',1.0);
hold off
ylabel('V_{e,L} (km/s)')
irf_legend(gca,{'V_{e,L}','|C_{Ae,L}|/2'},[0.05 0.92]);

i = i + 1;

%% ---- (c) Ve_M  CAe_M ----
h(i) = irf_subplot(n,1,-i);
irf_plot([Ve_lmn(:,1) Ve_M - Vi_M], 'g'); hold on
irf_plot([Ve_lmn(:,1) CAe_M/2], 'k','LineWidth',1.0);
hold off
ylabel('V_{e,M} (km/s)')
irf_legend(gca,{'V_{e,M}','|C_{Ae,M}|/2'},[0.05 0.92]);

i = i + 1;

%% ---- (d) Ve_N  CAe_N ----
h(i) = irf_subplot(n,1,-i);
irf_plot([Ve_lmn(:,1) Ve_N - Vi_N], 'r'); hold on
irf_plot([Ve_lmn(:,1) CAe_N/2], 'k','LineWidth',1.0);
hold off
ylabel('V_{e,N} (km/s)')
xlabel('UTC')
irf_legend(gca,{'V_{e,N}','|C_{Ae,N}|/2'},[0.05 0.92]);

%% ---- Align ----
irf_zoom(h,'x',tint);
irf_plot_axis_align(h);

disp('KH recognizing (LMN Ve) finished.');