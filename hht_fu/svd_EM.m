function [n,theta,t_wavelet,freq_wavelet] = svd_EM(B0, E0, B, E, freq_int)
%------written by Wending Fu, Mar.2024 in Beijing------------
% %  Input:
%
%    B0       - low-frequency background magnetic field columns (t bx by bz)
%    E        - high-frequency wave electric field, columns (t ex ey ez)
%    B        - high-frequency wave magnetic field, columns (t bx by bz)
%    freq_int - frequency interval: [fmin fmax]
%
%    Citation: Santolı  ́k, O., M. Parrot, and F. Lefeuvre, Singular value 
%    decomposition methods for wave propagation analysis, Radio Sci., 38(1), 
%    1010, doi:10.1029/2000RS002523, 2003.
%% convert TS
warning('off','all')
% Begin  format to older format (Must include spacecraft position)
if isa(E,'TSeries')
  ttemp = E.time.epochUnix;
  datatemp = double(E.data);
  E = [ttemp, double(datatemp)];
end
if isa(E0,'TSeries')
  ttemp = E0.time.epochUnix;
  datatemp = double(E0.data);
  E0 = [ttemp, double(datatemp)];
end
if isa(B,'TSeries')
  ttemp = B.time.epochUnix;
  datatemp = double(B.data);
  B = [ttemp, datatemp];
end
if isa(B0,'TSeries')
  ttemp = B0.time.epochUnix;
  datatemp = double(B0.data);
  B0 = [ttemp, datatemp];
end
%% resamp
SampleE=1/(E(2,1)-E(1,1));
SampleB=1/(B(2,1)-B(1,1));

if SampleB > 1.5*SampleE
    E=irf_resamp(E,B); B0=irf_resamp(B0,B); E0=irf_resamp(E0,B);
    disp('irf_pl_ebs: interpolating e to b');
elseif SampleE > 1.5*SampleB 
    B=irf_resamp(B,E); B0=irf_resamp(B0,E); E0=irf_resamp(E0,E);
    disp('irf_pl_ebs: interpolating b to e');
elseif SampleE == SampleB && E(1) ~= B(1)
    E=irf_resamp(E,B); B0=irf_resamp(B0,B); E0=irf_resamp(E0,B);
    disp(['b and e have the same sampling, but with different start time,' ...
        'irf_pl_ebs: interpolating e to b']);
elseif SampleE == SampleB && E(1) == B(1)
    disp(['b and e have the same sampling, no need for interpolating']);
else   
    inSampling=SampleE;
    t=max(E(1,1),B(1,1)):1/inSampling:min(E(end,1),B(end,1)); 
    t=t';
    E=irf_resamp(E,t); B=irf_resamp(B,t); B0=irf_resamp(B0,t); E0=irf_resamp(E0,t);
    irf_log('proc','interpolating B and E to E sampling');
end

B(:,2:4) = B(:,2:4) - B0(:,2:4);
E(:,2:4) = E(:,2:4) - E0(:,2:4);
B = irf_convert_fac(B,B0,[1,0,0]);
E = irf_convert_fac(E,B0,[1,0,0]);
%% MODWPT IMF
B = MODWPT_HHT(B,4);
E = MODWPT_HHT(E,4);
%% calculate wave spectrom by wavelet
% nT^2/Hz
WaveB = irf_wavelet(B,'returnpower',0,'cutedge',1,'f',[freq_int(1) freq_int(2)],'fs',8196);
% (mV/m)^2/Hz
WaveE = irf_wavelet(E,'returnpower',0,'cutedge',1,'f',[freq_int(1) freq_int(2)],'fs',8196);

freq_wavelet = WaveB.f;
t_wavelet = WaveB.t;
WaveB_wavelet_cell = WaveB.p;
WaveE_wavelet_cell = WaveE.p;

% convert to SI
units = irf_units; %#ok
% c*B (c*T^2/Hz)
c_eval('WaveB_wavelet_cell{?} = units.c * WaveB_wavelet_cell{?} * (1e-9)^2 * 1e6;',1:3);
% E (V/m)^2/Hz
c_eval('WaveE_wavelet_cell{?} = WaveE_wavelet_cell{?} * (1e-3)^2 * 1e6;',1:3);

WaveB_wavelet_cell1 = WaveB_wavelet_cell{1};WaveB_wavelet_cell2 = WaveB_wavelet_cell{2};WaveB_wavelet_cell3 = WaveB_wavelet_cell{3};
WaveE_wavelet_cell1 = WaveE_wavelet_cell{1};WaveE_wavelet_cell2 = WaveE_wavelet_cell{2};WaveE_wavelet_cell3 = WaveE_wavelet_cell{3};
%% parfor
if isempty(gcp('nocreate'))
    parpool;
end

n = zeros(length(t_wavelet),length(freq_wavelet),3);
theta = zeros(length(t_wavelet),length(freq_wavelet));

parfor_progress(length(freq_wavelet));
parfor freq = 1:length(freq_wavelet)
WaveB = [t_wavelet, WaveB_wavelet_cell1(:,freq), WaveB_wavelet_cell2(:,freq), WaveB_wavelet_cell3(:,freq)];
WaveE = [t_wavelet, WaveE_wavelet_cell1(:,freq), WaveE_wavelet_cell2(:,freq), WaveE_wavelet_cell3(:,freq)];
%% calculate A
%equation (15)-(18)

S = [WaveB(:,2:4), WaveE(:,2:4)];
Q = zeros(size(S,1),6,6);

% % % AE = zeros(size(S,1),6,3);
% % % b = zeros(size(S,1),6,1);
% % % UE = zeros(size(S,1),6,3);
AE = zeros(size(S,1),36,3);
b = zeros(size(S,1),36,1);
UE = zeros(size(S,1),36,3);

WE = zeros(size(S,1),3,3);
VE = zeros(size(S,1),3,3);

n_tt = zeros(size(S,1),3);
theta_tt = zeros(size(S,1),1);

for row = 1:size(Q,1)
    tempQ = S(row,:)' * conj(S(row,:));
    Q(row,:,:) = tempQ;

    tempAE = [0, real(tempQ(6,1)), -real(tempQ(5,1)); 0, imag(tempQ(6,1)), -imag(tempQ(5,1));
              0, real(tempQ(6,2)), -real(tempQ(5,2)); 0, imag(tempQ(6,2)), -imag(tempQ(5,2));
              0, real(tempQ(6,3)), -real(tempQ(5,3)); 0, imag(tempQ(6,3)), -imag(tempQ(5,3));
              0, real(tempQ(6,4)), -real(tempQ(5,4)); 0, imag(tempQ(6,4)), -imag(tempQ(5,4));
              0, real(tempQ(6,5)), -real(tempQ(5,5)); 0, imag(tempQ(6,5)), -imag(tempQ(5,5));
              0, real(tempQ(6,6)), -real(tempQ(5,6)); 0, imag(tempQ(6,6)), -imag(tempQ(5,6));
              -real(tempQ(6,1)), 0, real(tempQ(4,1)); -imag(tempQ(6,1)), 0, imag(tempQ(4,1));
              -real(tempQ(6,2)), 0, real(tempQ(4,2)); -imag(tempQ(6,2)), 0, imag(tempQ(4,2));
              -real(tempQ(6,3)), 0, real(tempQ(4,3)); -imag(tempQ(6,3)), 0, imag(tempQ(4,3));
              -real(tempQ(6,4)), 0, real(tempQ(4,4)); -imag(tempQ(6,4)), 0, imag(tempQ(4,4));
              -real(tempQ(6,5)), 0, real(tempQ(4,5)); -imag(tempQ(6,5)), 0, imag(tempQ(4,5));
              -real(tempQ(6,6)), 0, real(tempQ(4,6)); -imag(tempQ(6,6)), 0, imag(tempQ(4,6));
              real(tempQ(5,1)), -real(tempQ(4,1)), 0; imag(tempQ(5,1)), -imag(tempQ(4,1)), 0;
              real(tempQ(5,2)), -real(tempQ(4,2)), 0; imag(tempQ(5,2)), -imag(tempQ(4,2)), 0;
              real(tempQ(5,3)), -real(tempQ(4,3)), 0; imag(tempQ(5,3)), -imag(tempQ(4,3)), 0;
              real(tempQ(5,4)), -real(tempQ(4,4)), 0; imag(tempQ(5,4)), -imag(tempQ(4,4)), 0;
              real(tempQ(5,5)), -real(tempQ(4,5)), 0; imag(tempQ(5,5)), -imag(tempQ(4,5)), 0;
              real(tempQ(5,6)), -real(tempQ(4,6)), 0; imag(tempQ(5,6)), -imag(tempQ(4,6)), 0;];
        % % % tempAE = [0, real(tempQ(6,2)), -real(tempQ(5,2)); 0, imag(tempQ(6,2)), -imag(tempQ(5,2));
        % % %       -real(tempQ(6,3)), 0, real(tempQ(4,3)); -imag(tempQ(6,3)), 0, imag(tempQ(4,3));
        % % %       real(tempQ(5,1)), -real(tempQ(4,1)), 0; imag(tempQ(5,1)), -imag(tempQ(4,1)), 0;];
    AE(row,:,:) = tempAE;

    % % % tempb = [real(tempQ(1,2)); imag(tempQ(1,2));
    % % %          real(tempQ(2,3)); imag(tempQ(2,3));
    % % %          real(tempQ(3,1)); imag(tempQ(3,1));];
    tempb = [real(tempQ(1,1)); imag(tempQ(1,1));
             real(tempQ(1,2)); imag(tempQ(1,2));
             real(tempQ(1,3)); imag(tempQ(1,3));
             real(tempQ(1,4)); imag(tempQ(1,4));
             real(tempQ(1,5)); imag(tempQ(1,5));
             real(tempQ(1,6)); imag(tempQ(1,6));
             real(tempQ(2,1)); imag(tempQ(2,1));
             real(tempQ(2,2)); imag(tempQ(2,2));
             real(tempQ(2,3)); imag(tempQ(2,3));
             real(tempQ(2,4)); imag(tempQ(2,4));
             real(tempQ(2,5)); imag(tempQ(2,5));
             real(tempQ(2,6)); imag(tempQ(2,6));
             real(tempQ(3,1)); imag(tempQ(3,1));
             real(tempQ(3,2)); imag(tempQ(3,2));
             real(tempQ(3,3)); imag(tempQ(3,3));
             real(tempQ(3,4)); imag(tempQ(3,4));
             real(tempQ(3,5)); imag(tempQ(3,5));
             real(tempQ(3,6)); imag(tempQ(3,6));];
    b(row,:,:) = tempb;
    
    if any(any(isnan(tempAE)))
        % % % tempU = NaN*zeros(6,3); tempW = NaN*zeros(3,3); tempV = NaN*zeros(3,3);
        tempU = NaN*zeros(36,3); tempW = NaN*zeros(3,3); tempV = NaN*zeros(3,3);
    else
        [tempU, tempW, tempV] = svd(tempAE, 'econ');
        % [tempU, tempW, tempV] = svd(tempAE, 0);
    end
    UE(row,:,:) = tempU;
    WE(row,:,:) = tempW;
    VE(row,:,:) = tempV;
    
    tempn = transpose(tempV' / inv(tempW) * tempU' * tempb); 
    n_tt(row,:) = tempn;

    temptheta = acosd(dot(tempn, B0(row,2:4))/(norm(tempn)*norm(B0(row,2:4))));
    theta_tt(row) = temptheta;
end
    n(:,freq,:) = n_tt;
    theta(:,freq) = theta_tt;
    parfor_progress;
end
parfor_progress(0);
end

%% MODWPT improved HHT
function B_MOD = MODWPT_HHT(Bfac,lev)
%% MODWPT
% lev = floor(log2(numel(Bfac1(:,1)))); % default, about 16 for 10 sec
% 2^lev frequency bands
[wpt_x, ~, ~] = modwpt(Bfac(:,2), lev, TimeAlign = true);
[wpt_y, ~, ~] = modwpt(Bfac(:,3), lev, TimeAlign = true);
[wpt_z, ~, ~] = modwpt(Bfac(:,4), lev, TimeAlign = true);
%% emd and choose IMFs
for wpt_i = 1:2^lev
[imf_x{wpt_i}, res_x{wpt_i}, info_x{wpt_i}] = emd(wpt_x(wpt_i,:)', 'MaxNumIMF', 20);
[imf_y{wpt_i}, res_y{wpt_i}, info_y{wpt_i}] = emd(wpt_y(wpt_i,:)', 'MaxNumIMF', 20);
[imf_z{wpt_i}, res_z{wpt_i}, info_z{wpt_i}] = emd(wpt_z(wpt_i,:)', 'MaxNumIMF', 20);

% c_eval('[imf_x?, res_x?, info_x?] = vmd(Bfac?(:,2), ''NumIMF'', 10);',ic);
% c_eval('[imf_y?, res_y?, info_y?] = vmd(Bfac?(:,3), ''NumIMF'', 10);',ic);
% c_eval('[imf_z?, res_z?, info_z?] = vmd(Bfac?(:,4), ''NumIMF'', 10);',ic);
% c_eval('[imf_t?, res_t?, info_t?] = vmd(Bt?(:,2), ''NumIMF'', 10);',ic);

corr_x{wpt_i} = corr(imf_x{wpt_i},Bfac(:,2));
corr_y{wpt_i} = corr(imf_y{wpt_i},Bfac(:,3));
corr_z{wpt_i} = corr(imf_z{wpt_i},Bfac(:,4));

u_threshold_x(wpt_i) = threshold_cal(corr_x{wpt_i});
u_threshold_y(wpt_i) = threshold_cal(corr_y{wpt_i});
u_threshold_z(wpt_i) = threshold_cal(corr_z{wpt_i});

imf_x_phy{wpt_i} = imf_x{wpt_i}(:,corr_x{wpt_i} >= u_threshold_x(wpt_i));
imf_y_phy{wpt_i} = imf_y{wpt_i}(:,corr_y{wpt_i} >= u_threshold_y(wpt_i));
imf_z_phy{wpt_i} = imf_z{wpt_i}(:,corr_z{wpt_i} >= u_threshold_z(wpt_i));
end
%% 
imfx = zeros(size(wpt_x));imfy = zeros(size(wpt_y));imfz = zeros(size(wpt_z));
for wpt_i = 1:2^lev
    if ~isempty(imf_x_phy{wpt_i}), imfx(wpt_i,:) = sum(imf_x_phy{wpt_i},2);end
    if ~isempty(imf_y_phy{wpt_i}), imfy(wpt_i,:) = sum(imf_y_phy{wpt_i},2);end
    if ~isempty(imf_z_phy{wpt_i}), imfz(wpt_i,:) = sum(imf_z_phy{wpt_i},2);end
end
BMOD_x = imodwpt(imfx);
BMOD_y = imodwpt(imfy);
BMOD_z = imodwpt(imfz);

B_MOD = [Bfac(:,1),BMOD_x', BMOD_y', BMOD_z'];
end

%% Threshold Calculation
function threshold =  threshold_cal(corr)
if max(corr) >= 0
    threshold = max(corr)/(10 * max(corr) - 3);
else
    threshold = max(corr) / 3;
end
end