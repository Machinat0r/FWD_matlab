close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

TT = '2024-05-10T00:00:00.000Z/2024-05-11T00:00:00.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
try
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','hpca','drm','brst');
filenames4 = SDCFilenames(Date,iic,'inst','edp','drm','brst','dpt','dce');
% filenames5 = SDCFilenames(Date,ic,'inst','feeps','drm','brst','dpt','electron');
filenames_srvy = SDCFilenames(Date,iic,'inst','hpca','drm','srvy'); 
% filenames_fast1 = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms');
% filenames_fast2 = SDCFilenames(Date,ic,'inst','edp','drm','fast');
filenames = [filenames1, filenames2, filenames3, filenames4];
% filenames_fast = [filenames_fast1, filenames_fast2];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
[filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

SDCFilesDownload_NAS(filenames,TempDir, 'CheckSize', 0, 'Threads', 16)
% SDCFilesDownload(filenames,TempDir)
% % % 
% SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 64, 'CheckSize', 0)
SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 16, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
catch
    warning('no files have been downloaded')
end
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

%% plot_mms_hpca_tof_spectrogram_irfu
% Purpose:
%   Plot MMS HPCA TOF-COUNTS as a "mass-spectrogram-like" plot:
%   x = time, y = q/m proxy (~1/TOF^2) or TOF, color = counts (log).
%
% Requirements:
%   - IRFU-MATLAB (irfu-matlab) on MATLAB path
%   - Local MMS database initialized with mms.db_init('local_file_db', rootPath)
%     (rootPath should contain the MMS SDC-like directory tree)
%
% References:
%   - irfu-matlab usage suggests running "irf" each session. :contentReference[oaicite:1]{index=1}
%   - Example db_init syntax appears in irfu-matlab issues. :contentReference[oaicite:2]{index=2}
%   - irf_spectrogram supports specrec.plot_type='log' etc. :contentReference[oaicite:3]{index=3}

%% 1) Point IRFU MMS database to your local MMS data root
% Example:
%   D:\mms\mms2\hpca\srvy\l2\tof-counts\2024\05\mms2_hpca_srvy_l2_tof-counts_20240510....
mms.db_init('local_file_db', ParentDir);  % :contentReference[oaicite:5]{index=5}

%% 2) Settings
ic  = 1; % spacecraft: 1..4
tint = irf.tint('2024-05-10T00:00:00Z/2024-05-11T00:00:00Z');

% Dataset / variable names (HPCA survey L2 TOF-COUNTS)
% NOTE: Some installations use "tof-counts" (hyphen) as in file names.
% If your local DB uses underscores, change "tof-counts" -> "tof_counts".
prod = sprintf('mms%d_hpca_srvy_l2_tof-counts', ic);
var  = sprintf('mms%d_hpca_tof_counts_allA_atE', ic);   % preferred spectrogram-ready variable

use_qm_proxy = true;   % true: y = 1/TOF^2 proxy; false: y = TOF bins (or TOF axis)
sum_over_energy_dim = true; % if data is [t x tof x E], sum over E (else pick one index)
energy_index = 1;      % used only if sum_over_energy_dim = false

%% 3) Load data via IRFU MMS DB

% Works if variable can be represented as TSeries (time + numeric array)
hpca_file = dataobj('/Volumes/SPART-WORK/Data/MMS/mms1/hpca/srvy/l2/tof-counts/2024/05/mms1_hpca_srvy_l2_tof-counts_20240510000000_v4.3.0.cdf');
ts = get_variable(hpca_file, 'mms1_hpca_tof_counts', tint);
Epoch = get_variable(hpca_file, 'Epoch', tint);
Energy = get_variable(hpca_file, 'mms1_hpca_tof_energy', tint);
Bin_index = get_variable(hpca_file, 'mms1_hpca_tof_bin_index', tint);
t = Epoch.data;
P = double(ts.data);

specrec_p_e=struct('t',t);
specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
specrec_p_e.p=energy_e.data;%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','log10(keV/(cm^2 s sr keV))'};
[h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);




%% 4) Reduce dimensions to [Nt x Ny] if needed
sz = size(P);
if numel(sz) == 2
  % [Nt x Ny] OK
  P2 = P;
elseif numel(sz) == 3
  % Assume [Nt x Ny x Ne] OR [Nt x Ne x Ny] depending on file.
  % Heuristic: choose the dimension that looks like "time".
  Nt = numel(t);
  if sz(1) ~= Nt
    error('Time dimension mismatch: size(P,1)=%d, Nt=%d. Please inspect variable dimensions.', sz(1), Nt);
  end

  % Now P is [Nt x A x B]. We assume one of A/B is TOF bins.
  if sum_over_energy_dim
    P2 = squeeze(nansum(P, 3));  % sum over 3rd dim
  else
    P2 = squeeze(P(:,:,energy_index));
  end
else
  error('Unexpected dimensionality for %s: %s', var, mat2str(sz));
end

% Replace non-physical values for log plotting
P2(P2 <= 0) = NaN;

%% 5) Build y-axis (TOF or q/m proxy)
Ny = size(P2,2);

% Try to read TOF axis from metadata (DEPEND_1) if available
tof_axis = (1:Ny)';
try
  v1 = mms.db_get_variable(prod, var, tint);
  if isfield(v1,'DEPEND_1')
    if isstruct(v1.DEPEND_1) && isfield(v1.DEPEND_1,'data')
      cand = double(v1.DEPEND_1.data(:));
    else
      cand = double(v1.DEPEND_1(:));
    end
    if numel(cand) == Ny
      tof_axis = cand;
    end
  end
catch
  % keep default bin index
end

if use_qm_proxy
  y = 1 ./ (tof_axis.^2);
  ylab = 'q/m proxy ~ 1/TOF^2 (arb.)';
else
  y = tof_axis;
  ylab = 'TOF (bin or axis units)';
end

%% 6) Prepare IRFU spectrogram record and plot
specrec = struct();
specrec.t = t;            % EpochTT time array
specrec.f = y;            % y-axis
specrec.p = P2.';         % irf_spectrogram commonly expects [Ny x Nt]
specrec.f_label = ylab;
specrec.p_label = 'counts';
specrec.plot_type = 'log'; % :contentReference[oaicite:6]{index=6}

h = irf_plot(1, 'newfigure');
irf_spectrogram(h(1), specrec);       % :contentReference[oaicite:7]{index=7}
irf_zoom(h(1), 'x', tint);
irf_timeaxis(h(1), 'utc');
grid(h(1), 'on');
% title(h(1), sprintf('MMS%d HPCA TOF-COUNTS (%s)', ic, char(tint)));
