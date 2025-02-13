function [k,theta,t_wavelet,freq_wavelet, B, E] = svd_EM(B0, E0, B, E, freq_int)
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
    disp('b and e have the same sampling, no need for interpolating');
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
% [B,E] = ortho_denoise(B,E,4);
%% calculate wave spectrom by wavelet
% convert to SI
units = irf_units; 

% T^2/Hz
WaveB = irf_wavelet(units.c*B*1e-9,'returnpower',0,'cutedge',1,'f',[freq_int(1) freq_int(2)],'fs',8196);
% (V/m)^2/Hz
WaveE = irf_wavelet(E*1e-3,'returnpower',0,'cutedge',1,'f',[freq_int(1) freq_int(2)],'fs',8196);

freq_wavelet = WaveB.f;
t_wavelet = WaveB.t;
WaveB_wavelet_cell = WaveB.p;
WaveE_wavelet_cell = WaveE.p;

% c*B (c*T^2/Hz)
c_eval('WaveB_wavelet_cell{?} = WaveB_wavelet_cell{?};',1:3);
% E (V/m)^2/Hz
c_eval('WaveE_wavelet_cell{?} = WaveE_wavelet_cell{?};',1:3);

WaveB_wavelet_cell1 = WaveB_wavelet_cell{1};WaveB_wavelet_cell2 = WaveB_wavelet_cell{2};WaveB_wavelet_cell3 = WaveB_wavelet_cell{3};
WaveE_wavelet_cell1 = WaveE_wavelet_cell{1};WaveE_wavelet_cell2 = WaveE_wavelet_cell{2};WaveE_wavelet_cell3 = WaveE_wavelet_cell{3};
%% parfor
if isempty(gcp('nocreate'))
    parpool;
end

n = zeros(length(t_wavelet),length(freq_wavelet),3);
theta = zeros(length(t_wavelet),length(freq_wavelet));

% % % parfor_progress(length(freq_wavelet));
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
    
    tempn = transpose(tempV' / tempW * tempU' * tempb); 
    n_tt(row,:) = tempn;

    % temptheta = acosd(dot(tempn, B0(row,2:4))/(norm(tempn)*norm(B0(row,2:4))));
    temptheta = acosd(dot(tempn, [1, 0, 0])/(norm(tempn)*norm([1, 0, 0])));
    theta_tt(row) = temptheta;
end
    n(:,freq,:) = n_tt;
    theta(:,freq) = theta_tt;
    % % % parfor_progress;
end
k = nan*zeros(size(n));
freq_time_mat = repmat(freq_wavelet',[size(n,1),1]);
for comp = 1:3
    k(:,:,comp) = n(:,:,comp) * 2 * pi .* freq_time_mat / units.c;
end
% % % parfor_progress(0);
end

%% Orthogonal Noise Reduction
function [B_MOD, E_MOD] = ortho_denoise(Bfac,Efac, lev)
%% MODWPT
% lev = floor(log2(numel(Bfac1(:,1)))); % default, about 16 for 10 sec
% 2^lev frequency bands
[wpt_Bx, ~, ~] = modwpt(Bfac(:,2), lev, TimeAlign = true);
[wpt_By, ~, ~] = modwpt(Bfac(:,3), lev, TimeAlign = true);
[wpt_Bz, ~, ~] = modwpt(Bfac(:,4), lev, TimeAlign = true);

[wpt_Ex, ~, ~] = modwpt(Efac(:,2), lev, TimeAlign = true);
[wpt_Ey, ~, ~] = modwpt(Efac(:,3), lev, TimeAlign = true);
[wpt_Ez, ~, ~] = modwpt(Efac(:,4), lev, TimeAlign = true);

corr_wpt_x = ortho(wpt_Bx,wpt_Ex); %#ok
corr_wpt_y = ortho(wpt_By,wpt_Ey); %#ok
corr_wpt_z = ortho(wpt_Bz,wpt_Ez); %#ok
%% emd and choose IMFs
for wpt_i = 1:2^lev
% B and E should be orthogonalized
threshold_orth = 0.3;

% x
[tempB_imf_x, ~, ~] = emd(wpt_Bx(wpt_i,:)', 'MaxNumIMF', 20);
[tempE_imf_x, ~, ~] = emd(wpt_Ex(wpt_i,:)', 'MaxNumIMF', 20);

NumIMF_x = min(size(tempB_imf_x,2),size(tempE_imf_x,2));
corr_imf_x = ones(max(size(tempB_imf_x,2),size(tempE_imf_x,2)),1);
corr_imf_x(1:NumIMF_x) = ortho(tempB_imf_x(:,1:NumIMF_x),tempE_imf_x(:,1:NumIMF_x));

corr_wpt_Bx = corr(tempB_imf_x,Bfac(:,2));
threshold_Bx = threshold_cal(corr_wpt_Bx);
corr_wpt_Ex = corr(tempE_imf_x,Efac(:,2));
threshold_Ex = threshold_cal(corr_wpt_Ex);

tempB_imf_x = tempB_imf_x(:,corr_imf_x(1:size(tempB_imf_x,2)) < threshold_orth & corr_wpt_Bx >= threshold_Bx); 
tempE_imf_x = tempE_imf_x(:,corr_imf_x(1:size(tempE_imf_x,2)) < threshold_orth & corr_wpt_Ex >= threshold_Ex);
% % % tempB_imf_x = tempB_imf_x(:,corr_wpt_Bx >= threshold_Bx); 
% % % tempE_imf_x = tempE_imf_x(:,corr_wpt_Ex >= threshold_Ex);

wpt_Bx(wpt_i,:) = sum(tempB_imf_x,2); 
wpt_Ex(wpt_i,:) = sum(tempE_imf_x,2);

% y
[tempB_imf_y, ~, ~] = emd(wpt_By(wpt_i,:)', 'MaxNumIMF', 20);
[tempE_imf_y, ~, ~] = emd(wpt_Ey(wpt_i,:)', 'MaxNumIMF', 20);

NumIMF_y = min(size(tempB_imf_y,2),size(tempE_imf_y,2));
corr_imf_y = ones(max(size(tempB_imf_y,2),size(tempE_imf_y,2)),1);
corr_imf_y(1:NumIMF_y) = ortho(tempB_imf_y(:,1:NumIMF_y),tempE_imf_y(:,1:NumIMF_y));

corr_wpt_By = corr(tempB_imf_y,Bfac(:,3));
threshold_By = threshold_cal(corr_wpt_By);
corr_wpt_Ey = corr(tempE_imf_y,Efac(:,3));
threshold_Ey = threshold_cal(corr_wpt_Ey);

tempB_imf_y = tempB_imf_y(:,corr_imf_y(1:size(tempB_imf_y,2)) < threshold_orth & corr_wpt_By >= threshold_By); 
tempE_imf_y = tempE_imf_y(:,corr_imf_y(1:size(tempE_imf_y,2)) < threshold_orth & corr_wpt_Ey >= threshold_Ey);
% % % tempB_imf_y = tempB_imf_y(:,corr_wpt_By >= threshold_By); 
% % % tempE_imf_y = tempE_imf_y(:,corr_wpt_Ey >= threshold_Ey);

wpt_By(wpt_i,:) = sum(tempB_imf_y,2); 
wpt_Ey(wpt_i,:) = sum(tempE_imf_y,2);

% z
[tempB_imf_z, ~, ~] = emd(wpt_Bz(wpt_i,:)', 'MaxNumIMF', 20);
[tempE_imf_z, ~, ~] = emd(wpt_Ez(wpt_i,:)', 'MaxNumIMF', 20);

NumIMF_z = min(size(tempB_imf_z,2),size(tempE_imf_z,2));
corr_imf_z = ones(max(size(tempB_imf_z,2),size(tempE_imf_z,2)),1);
corr_imf_z(1:NumIMF_z) = ortho(tempB_imf_z(:,1:NumIMF_z),tempE_imf_z(:,1:NumIMF_z));

corr_wpt_Bz = corr(tempB_imf_z,Bfac(:,4));
threshold_Bz = threshold_cal(corr_wpt_Bz);
corr_wpt_Ez = corr(tempE_imf_z,Efac(:,4));
threshold_Ez = threshold_cal(corr_wpt_Ez);

tempB_imf_z = tempB_imf_z(:,corr_imf_z(1:size(tempB_imf_z,2)) < threshold_orth & corr_wpt_Bz >= threshold_Bz); 
tempE_imf_z = tempE_imf_z(:,corr_imf_z(1:size(tempE_imf_z,2)) < threshold_orth & corr_wpt_Ez >= threshold_Ez);
% % % tempB_imf_z = tempB_imf_z(:,corr_wpt_Bz >= threshold_Bz); 
% % % tempE_imf_z = tempE_imf_z(:,corr_wpt_Ez >= threshold_Ez);

wpt_Bz(wpt_i,:) = sum(tempB_imf_z,2); 
wpt_Ez(wpt_i,:) = sum(tempE_imf_z,2);
end
%% inverse MODWPT
BMOD_x = imodwpt(wpt_Bx);
BMOD_y = imodwpt(wpt_By);
BMOD_z = imodwpt(wpt_Bz);
B_MOD = [Bfac(:,1),BMOD_x', BMOD_y', BMOD_z'];

EMOD_x = imodwpt(wpt_Ex);
EMOD_y = imodwpt(wpt_Ey);
EMOD_z = imodwpt(wpt_Ez);
E_MOD = [Efac(:,1),EMOD_x', EMOD_y', EMOD_z'];
end
%% Threshold Calculation
function threshold =  threshold_cal(corr)
if max(corr) >= 0
    threshold = max(corr)/(10 * max(corr) - 3);
else
    threshold = max(corr) / 3;
end
end
%% Orthogonalization Examination
function corr_coef = ortho(wpt_B,wpt_E)
if any(size(wpt_B) ~= size(wpt_E))
    warning('Different Levels of MODWPT to B and E!')
end

if size(wpt_B,1) < size(wpt_B,2)
    wpt_B = transpose(wpt_B);
    wpt_E = transpose(wpt_E);
end

corr_coef = zeros(size(wpt_B,2),1);
for cfre = 1:size(wpt_B,2)
    corr_coef(cfre) = abs(corr(wpt_B(:,cfre), wpt_E(:,cfre))); % correlative coeffecient is the same as the orthogonality
    % corr_coef(cfre) = dot(wpt_B(:,cfre), wpt_E(:,cfre))/(norm(wpt_B(:,cfre))*norm(wpt_E(:,cfre)));
end
end