close all
clear;clc;
%% import data
OMNI_Data = importdata("C:\Users\pc\Documents\Timor\Recovery-Work-DF-MP\4-Differential Rotation\数理统计\OMNI2_H0_MRG1HR_1426744.txt");
Dst_Data = importdata("C:\Users\pc\Documents\Timor\Recovery-Work-DF-MP\4-Differential Rotation\数理统计\WWW_dstae00456773.dat.txt");
AE_Data = importdata("C:\Users\pc\Documents\Timor\Recovery-Work-DF-MP\4-Differential Rotation\数理统计\OMNI2_H0_MRG1HR_465521.txt");
lines = readlines("C:\Users\pc\Documents\Timor\Recovery-Work-DF-MP\4-Differential Rotation\数理统计\Kp_ap.txt");
time_length = 1:365*24*10;
time_length2 = 365*24*10:365*24*11;
%% Kp & ap Data
lines = lines(2:end);

date_strs = extractBetween(lines, 1, 8);
base_dates = datetime(date_strs, 'InputFormat', 'yyyyMMdd');

kp_tokens = regexp(lines, '[0-9][+-]?', 'match'); 
kp_tokens = cellfun(@(x) x(1:8), kp_tokens, 'UniformOutput', false); 

kp_matrix_str = vertcat(kp_tokens{:}); 
kp_matrix_str = string(kp_matrix_str);

kp_base = double(extract(kp_matrix_str, digitsPattern));
kp_sign = zeros(size(kp_matrix_str));
kp_sign(endsWith(kp_matrix_str, "+")) = 1/3;
kp_sign(endsWith(kp_matrix_str, "-")) = -1/3;
kp_val = kp_base + kp_sign;

kp_hours = hours(0:3:21);                      % 1 × 8
time_matrix = base_dates + kp_hours;          % n × 8 datetime
time_flat = reshape(time_matrix', [], 1);     % m × 1 datetime
kp_flat   = reshape(kp_val', [], 1);          % m × 1 double
kp = [kp_flat];

ap_tokens = regexp(lines, '\d+', 'match');    % 每行提取所有数字
ap_vals = cellfun(@(x) str2double(x(end-9:end-2)), ap_tokens, 'UniformOutput', false); % 提取 ap1~ap8
ap_matrix = vertcat(ap_vals{:});              % n × 8 double
ap_flat = reshape(ap_matrix', [], 1);         % m × 1
ap = [ap_flat];

ap(ap>=1000) = nan;
kp = repelem(kp,3);
ap = repelem(ap,3);
kp1 = kp(time_length);
ap1 = ap(time_length);
kp2 = kp(time_length2);
ap2 = ap(time_length2);
%% Dst Data
textdata = Dst_Data.textdata; 
data = Dst_Data.data(:, 2:25);
Dst_mean = Dst_Data.data(:, 26);
Dst_mean = repelem(Dst_mean,24);

str = string(textdata(:,1));
yy = 2000 + str2double(extractBetween(str, 4, 5));
mm = str2double(extractBetween(str, 6, 7));
dd = str2double(extractBetween(str, 9, 10));
YY = repelem(yy, 24);
MM = repelem(mm, 24);
DD = repelem(dd, 24);
HH = repmat(transpose(0:23), size(data, 1), 1);
time = datetime(YY, MM, DD, HH, 30, 0);

Dst_flat = reshape(data', [], 1);
Dst = [irf_time(datenum(time),'date>epoch'), Dst_flat];
Dst = Dst(1:end-24*31,:); %delete the first month in 2025
Dst = Dst(:,2);
Dst_var = Dst - Dst_mean(1:end-24*31);

Dst1 = Dst(time_length,:);
Dst2 = Dst(time_length2,:);
Dst_var1 = Dst_var(time_length,:);
Dst_var2 = Dst_var(time_length2,:);
%% OMNI Data
n = numel(OMNI_Data)-5;
joined = strjoin(OMNI_Data(3:end-3), newline);
fmt = '%s %s %f %f %f %f %f %f %f';
parsed = textscan(joined, fmt, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);

time_str = strcat(parsed{1}, {' '}, parsed{2});
time_dt = datetime(time_str, 'InputFormat', 'dd-MM-yyyy HH:mm:ss.SSS');
time_num = irf_time(datenum(time_dt),'date>epoch'); 

% Bx = [time_num, parsed{3}];
% By = [time_num, parsed{4}];
% Bz = [time_num, parsed{5}];
% T = [time_num, parsed{6}];
% n = [time_num, parsed{7}];
% V = [time_num, parsed{8}];
% Pdyn = [time_num, parsed{9}];

Bx = [parsed{3}];
By = [parsed{4}];
Bz = [parsed{5}];
T = [parsed{6}];
n = [parsed{7}];
V = [parsed{8}];
Pdyn = [parsed{9}];

Pdyn(Pdyn>99) = nan;
V(V>=9990) = nan;
Bz(Bz>=999) = nan;
Bx(Bx>=999) = nan;
By(By>=999) = nan;
n(n>=999) = nan;

Pdyn1 = Pdyn(time_length,:); n1 = n(time_length,:);
Bz1 = Bz(time_length,:); V1 = V(time_length,:);
Bx1 = Bx(time_length,:); By1 = By(time_length,:);
Pdyn2 = Pdyn(time_length2,:); n2 = n(time_length2,:);
Bz2 = Bz(time_length2,:); V2 = V(time_length2,:);
Bx2 = Bx(time_length2,:); By2 = By(time_length2,:);
%% AE Data
n = numel(AE_Data);
joined = strjoin(AE_Data, newline);
fmt = '%s %s %f %f %f %f %f %f %f';
parsed = textscan(joined, fmt, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);

SunSpot = [parsed{3}];
AE = [parsed{4}];
SunSpot1 = SunSpot(time_length);
AE1 = AE(time_length);
SunSpot2 = SunSpot(time_length2);
AE2 = AE(time_length2);
%% Linear Regression
dt = 24*30;
% AE = smooth(AE, dt);
% Dst = smooth(Dst, dt);
% SunSpot = smooth(SunSpot, dt);
% Bx = smooth(Bx, dt);
% Bz = smooth(Bz, dt);
% By = smooth(By, dt);
% kp = smooth(kp, dt);
% ap = smooth(ap, dt);

% idx = find(abs(Dst)<=50);
% Dst = Dst(idx);
% Dst_var = Dst_var(idx);
% Pdyn = Pdyn(idx);
% n = n(idx);
% V = V(idx);
% Bz = Bz(idx);
% Bx = Bx(idx);
% By = By(idx);
% V2 = V.^2;