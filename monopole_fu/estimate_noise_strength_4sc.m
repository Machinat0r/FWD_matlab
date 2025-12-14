function out = estimate_noise_strength_4sc(B1, B2, B3, B4, opts)
% out = estimate_noise_strength_4sc(B1,B2,B3,B4, opts)
% ------------------------------------------------------------
% Estimate noise strength from a "no-event" time interval for 4 spacecraft.
%
% Inputs:
%   B1..B4 : Nx3 or Nx4 arrays.
%            If Nx4, first column is time, cols 2:4 are [Bx By Bz] in nT.
%            If Nx3, interpreted as [Bx By Bz] in nT (same timeline assumed).
%   opts (optional struct):
%       .removeCommonMode (default true)
%           If true: subtract 4-s/c mean field at each time: b_i <- b_i - mean(b_all)
%           This removes common background fluctuations shared by all spacecraft.
%       .detrendMode (default "linear")
%           "none" | "constant" | "linear"
%           Applied per satellite and per component after common-mode removal.
%       .hpWindowSec (default [])
%           If set, remove slow variation by subtracting moving average with this window (seconds).
%           Requires time column in B1..B4. Does NOT need toolboxes.
%       .alignTime (default false)
%           If true and time columns exist but not identical, interpolates to common time grid (B1 time).
%       .robust (default false)
%           If true, use MAD-based sigma (more robust to outliers).
%
% Outputs (units):
%   out.sigma_sc_xyz   : 4x3 std (nT) per spacecraft & component
%   out.sigma_scalar   : scalar std (nT) pooled (useful for isotropic assumption)
%   out.SigmaVar_scalar: scalar variance (nT^2) = out.sigma_scalar^2
%   out.SigmaCov_sc    : 3x3x4 covariance matrices (nT^2) per spacecraft
%   out.SigmaCov_all   : 3x3 covariance from concatenating all residuals (nT^2)
%   out.dt             : estimated sampling interval (sec) if time provided
%   out.N              : number of samples used

arguments
    B1 double
    B2 double
    B3 double
    B4 double
    opts.removeCommonMode (1,1) logical = true
    opts.detrendMode (1,1) string = "linear"
    opts.hpWindowSec double = []
    opts.alignTime (1,1) logical = false
    opts.robust (1,1) logical = false
end

% ---- Extract time and XYZ
[t1, X1] = take_t_xyz(B1);
[t2, X2] = take_t_xyz(B2);
[t3, X3] = take_t_xyz(B3);
[t4, X4] = take_t_xyz(B4);

% ---- Optional time alignment (interp to t1)
if opts.alignTime && ~isempty(t1) && ~isempty(t2) && ~isempty(t3) && ~isempty(t4)
    X2 = interp1(t2, X2, t1, 'linear', 'extrap');
    X3 = interp1(t3, X3, t1, 'linear', 'extrap');
    X4 = interp1(t4, X4, t1, 'linear', 'extrap');
    t = t1;
else
    % assume same sampling & length
    t = t1; % may be empty
end

% ---- Basic checks
N = size(X1,1);
if any([size(X2,1), size(X3,1), size(X4,1)] ~= N)
    error('B1..B4 must have the same number of samples (or set opts.alignTime=true with time columns).');
end

% ---- Remove common mode: b_i <- b_i - mean(b_all)
if opts.removeCommonMode
    Xm = (X1 + X2 + X3 + X4) / 4;
    X1 = X1 - Xm;
    X2 = X2 - Xm;
    X3 = X3 - Xm;
    X4 = X4 - Xm;
end

% ---- Optional high-pass via moving average subtraction
dt = [];
if ~isempty(opts.hpWindowSec)
    if isempty(t)
        error('opts.hpWindowSec requires time column in inputs (Nx4).');
    end
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('Invalid time vector for high-pass windowing.');
    end
    win = max(3, round(opts.hpWindowSec / dt));
    X1 = X1 - movmean(X1, win, 1, 'omitnan');
    X2 = X2 - movmean(X2, win, 1, 'omitnan');
    X3 = X3 - movmean(X3, win, 1, 'omitnan');
    X4 = X4 - movmean(X4, win, 1, 'omitnan');
end

% ---- Detrend
X1 = apply_detrend(X1, opts.detrendMode);
X2 = apply_detrend(X2, opts.detrendMode);
X3 = apply_detrend(X3, opts.detrendMode);
X4 = apply_detrend(X4, opts.detrendMode);

% ---- Robust or standard sigma per sat per component
sigma_sc_xyz = zeros(4,3);
SigmaCov_sc  = zeros(3,3,4);

Xcell = {X1,X2,X3,X4};
for i = 1:4
    Xi = Xcell{i};

    if opts.robust
        % MAD -> sigma: sigma ~= 1.4826 * median(|x - median(x)|)
        sigma_sc_xyz(i,:) = 1.4826 * median(abs(Xi - median(Xi,1,'omitnan')), 1, 'omitnan');
    else
        sigma_sc_xyz(i,:) = std(Xi, 0, 1, 'omitnan');
    end

    % covariance (nT^2)
    SigmaCov_sc(:,:,i) = cov(Xi, 'omitrows');
end

% ---- Pooled / concatenated covariance
Xall = [X1; X2; X3; X4];   % (4N)x3
SigmaCov_all = cov(Xall, 'omitrows');

% ---- Scalar sigma (isotropic-ish): average variance across sats & components
vars = sigma_sc_xyz.^2; % 4x3
sigma_scalar = sqrt(mean(vars(:), 'omitnan')); % nT
SigmaVar_scalar = sigma_scalar^2;              % nT^2

out.sigma_sc_xyz = sigma_sc_xyz;
out.sigma_scalar = sigma_scalar;
out.SigmaVar_scalar = SigmaVar_scalar;
out.SigmaCov_sc = SigmaCov_sc;
out.SigmaCov_all = SigmaCov_all;
out.dt = dt;
out.N = N;

end

% ---------------- helpers ----------------
function [t, X] = take_t_xyz(B)
if size(B,2) == 4
    t = B(:,1);
    X = B(:,2:4);
elseif size(B,2) == 3
    t = [];
    X = B(:,1:3);
else
    error('Each B must be Nx3 or Nx4 (t + Bx By Bz).');
end
end

function Xo = apply_detrend(X, mode)
switch lower(mode)
    case "none"
        Xo = X;
    case "constant"
        Xo = X - mean(X,1,'omitnan');
    case "linear"
        % detrend each column
        Xo = zeros(size(X));
        for k=1:3
            Xo(:,k) = detrend(X(:,k));
        end
    otherwise
        error('opts.detrendMode must be "none" | "constant" | "linear".');
end
end
