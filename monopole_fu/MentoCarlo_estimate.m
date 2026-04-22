function res = MentoCarlo_estimate(Rin, Bin, Sigma, K, opts)
% res = mm_event_confidence_monopole(Rin, Bin, Sigma, K, opts)
% ------------------------------------------------------------
% Event confidence for monopole detection using MC-calibrated LRT under H0: Q=0
% Matching your units/formula:
%   B in nT, R in km, Bmono = Q/(4*pi) * (r-r0)/|r-r0|^3  (=> Q has unit nT*km^2)
%
% Inputs
%   Rin   : 4x3 positions (km) or 4x4 with time column (t, x, y, z)
%   Bin   : 4x3 B-field  (nT) or 4x4 with time column (t, Bx, By, Bz)
%   Sigma : noise model in nT. Supported:
%           - scalar variance: Sigma = sigma^2
%           - 3x3 covariance shared by all sats
%           - 4x3 std per sat per axis: Sigma(i,:)=[sx sy sz]
%           - 3x3x4 covariance per sat
%   K     : #MC replicates (e.g., 5000~20000)
%   opts  : struct with optional fields:
%       .seed        (default 1)
%       .Qsign       (default 0)  0 no constraint, +1 force Q>=0, -1 force Q<=0
%       .r0_init     (default mean(R))
%       .Q_init      (default 1e4) used only if fitMethod='lsqnonlin4'
%       .fitMethod   (default 'profile_fminsearch')  'profile_fminsearch' or 'lsqnonlin4'
%       .maxIter     (default 800) iterations for optimizer
%       .tolX        (default 1e-10)
%       .tolFun      (default 1e-10)
%       .r_eps_km    (default 1e-6 * typical scale) softening to avoid singularity
%       .saveT       (default false) save MC T distribution
%
% Output res fields
%   res.Qhat, res.r0hat
%   res.Tobs, res.pvalue, res.confidence (=1-p), res.Z (one-sided sigma)
%   res.chi2_null, res.chi2_alt
%   res.Tmc (optional)

arguments
    Rin double
    Bin double
    Sigma
    K (1,1) double {mustBeInteger,mustBePositive} = 5000
    opts.seed (1,1) double = 1
    opts.Qsign (1,1) double = 0
    opts.r0_init double = []
    opts.Q_init (1,1) double = 1e4
    opts.fitMethod (1,1) string = "profile_fminsearch"
    opts.maxIter (1,1) double = 800
    opts.tolX (1,1) double = 1e-10
    opts.tolFun (1,1) double = 1e-10
    opts.r_eps_km double = []
    opts.saveT (1,1) logical = false
end

rng(opts.seed);

% ---- Extract (km, nT) vectors, remove time column if exists
R = take_xyz(Rin);
Bobs = take_xyz(Bin);

if size(R,1) ~= 4 || size(R,2) ~= 3
    error('R must be 4x3 (or 4x4 with time column).');
end
if size(Bobs,1) ~= 4 || size(Bobs,2) ~= 3
    error('B must be 4x3 (or 4x4 with time column).');
end

% softening (km)
if isempty(opts.r_eps_km)
    % typical scale ~ median distance to centroid
    scale = median(vecnorm(R - mean(R,1), 2, 2));
    if ~isfinite(scale) || scale<=0, scale = 1; end
    r_eps = 1e-6 * scale;
else
    r_eps = opts.r_eps_km;
end

% init
if isempty(opts.r0_init), r0_init = mean(R,1);
else, r0_init = reshape(opts.r0_init,1,3);
end

% ---- Fit H1 on observed data
fitObs = fit_alt(R, Bobs, Sigma, r0_init, opts, r_eps);
chi2_alt = fitObs.chi2;
chi2_null = chi2_under_null(Bobs, Sigma);

Tobs = chi2_null - chi2_alt;

% ---- MC under H0: B = noise
Tmc = zeros(K,1);

% For speed: use observed r0hat as init for all MC fits
r0_init_mc = fitObs.r0hat;

% parfor_progress(K);
for k = 1:K
    E = sample_noise(4, Sigma);  % 4x3 in nT
    Bk = E;

    fitK = fit_alt(R, Bk, Sigma, r0_init_mc, opts, r_eps);
    chi2_alt_k = fitK.chi2;
    chi2_null_k = chi2_under_null(Bk, Sigma);

    Tmc(k) = chi2_null_k - chi2_alt_k;
    
    clc;k
    % parfor_progress;
end
parfor_progress(0);

% one-sided p-value with +1 smoothing
p = (sum(Tmc >= Tobs) + 1) / (K + 1);
conf = 1 - p;

% one-sided Z without requiring Statistics Toolbox:
% Z = norminv(1-p) = sqrt(2)*erfcinv(2*p)
Z = sqrt(2) * erfcinv(2*p);

res.Qhat = fitObs.Qhat;
res.r0hat = fitObs.r0hat;
res.Tobs = Tobs;
res.pvalue = p;
res.confidence = conf;
res.Z = Z;
res.chi2_null = chi2_null;
res.chi2_alt = chi2_alt;
res.r_eps_km = r_eps;

if opts.saveT
    res.Tmc = Tmc;
end

end

%% ----------------- helper: fit under H1 -----------------
function out = fit_alt(R, Bobs, Sigma, r0_init, opts, r_eps)
switch lower(opts.fitMethod)
    case "profile_fminsearch"
        % Optimize r0 only; Q is solved analytically (because model is linear in Q)
        obj = @(x) chi2_profile_r0(x(:).', R, Bobs, Sigma, opts.Qsign, r_eps);

        o = optimset('MaxIter', opts.maxIter, 'TolX', opts.tolX, 'TolFun', opts.tolFun, 'Display', 'off');
        [r0hat, chi2min] = fminsearch(obj, r0_init, o);

        % recover Qhat at r0hat
        [~, Qhat] = chi2_profile_r0(r0hat, R, Bobs, Sigma, opts.Qsign, r_eps);

        out.Qhat = Qhat;
        out.r0hat = reshape(r0hat,1,3);
        out.chi2 = chi2min;

    case "lsqnonlin4"
        % Direct 4-parameter lsqnonlin on t=[Q, x0, y0, z0]
        if ~exist('lsqnonlin','file')
            error('fitMethod=lsqnonlin4 requires Optimization Toolbox (lsqnonlin).');
        end

        t0 = [opts.Q_init, r0_init];
        options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', ...
            'OptimalityTolerance',opts.tolFun,'FunctionTolerance',opts.tolFun, ...
            'MaxFunctionEvaluations',2e4,'MaxIterations',opts.maxIter,'Display','off');

        % optional Q sign bounds
        lb = [-inf,-inf,-inf,-inf]; ub = [inf,inf,inf,inf];
        if opts.Qsign > 0
            lb(1) = 0;
        elseif opts.Qsign < 0
            ub(1) = 0;
        end

        [that, ~, resvec] = lsqnonlin(@(t) resid_whitened_4par(t, R, Bobs, Sigma, r_eps), t0, lb, ub, options);

        out.Qhat = that(1);
        out.r0hat = that(2:4);
        out.chi2 = sum(resvec(:).^2);

    otherwise
        error('Unknown opts.fitMethod: %s', opts.fitMethod);
end
end

%% ---- chi2 profile: given r0, solve Q analytically, compute chi2 ----
function [chi2, Qhat] = chi2_profile_r0(r0, R, Bobs, Sigma, Qsign, r_eps)
% Build stacked (whitened) vectors for Bobs and g
bw = whiten_stack(Bobs, Sigma);      % 12x1
gw = whiten_stack(g_stack(R, r0, r_eps), Sigma); % 12x1 for unit of (dr/|dr|^3)

den = gw' * gw;
if den <= 0 || ~isfinite(den)
    Qhat = 0;
else
    % Model: B = (Q/(4*pi)) * g
    Qhat = 4*pi * (gw' * bw) / den;
end

% apply optional sign constraint consistent with your idx_flag usage
if Qsign > 0
    Qhat = max(Qhat, 0);
elseif Qsign < 0
    Qhat = min(Qhat, 0);
end

rw = bw - (Qhat/(4*pi)) * gw;
chi2 = rw' * rw;
end

%% ---- residual for lsqnonlin4 ----
function rw = resid_whitened_4par(t, R, Bobs, Sigma, r_eps)
Q = t(1);
r0 = t(2:4);
g  = g_stack(R, r0, r_eps); % 4x3
Bpred = (Q/(4*pi)) * g;     % 4x3, in nT (Q carries nT*km^2)
res = Bobs - Bpred;         % 4x3
rw = whiten_stack(res, Sigma); % 12x1
end

%% ---- g = (r-r0)/|r-r0|^3  (km^-2) ----
function g = g_stack(R, r0, r_eps)
dr = R - reshape(r0,1,3);                 % 4x3 (km)
d2 = sum(dr.^2, 2) + r_eps^2;             % 4x1 (km^2), softened
inv_r3 = 1 ./ (d2 .* sqrt(d2));           % 1/km^3
g = dr .* inv_r3;                         % 4x3 (1/km^2)
end

%% ---- chi2 under null: Q=0 => residual=Bobs ----
function chi2 = chi2_under_null(Bobs, Sigma)
bw = whiten_stack(Bobs, Sigma);
chi2 = bw' * bw;
end

%% ---- whiten: stack 4x3 to 12x1 and apply Sigma^-1/2 ----
function v = whiten_stack(X, Sigma)
% X can be 4x3; returns 12x1 whitened vector
if isequal(size(X),[4 3])
    % ok
elseif isequal(size(X),[4 4])
    X = X(:,2:4);
else
    error('whiten_stack expects 4x3 (or 4x4 with time column).');
end

v = zeros(12,1);
for i = 1:4
    xi = X(i,:).'; % 3x1
    wi = whiten_vec3(xi, Sigma, i);
    v( (i-1)*3 + (1:3) ) = wi;
end
end

function w = whiten_vec3(x, Sigma, iSat)
% returns Sigma^{-1/2} x such that ||w||^2 = x' Sigma^{-1} x
if isscalar(Sigma)
    % Sigma is variance (nT^2)
    w = x / sqrt(Sigma);
elseif isequal(size(Sigma),[3 3])
    % shared covariance
    L = chol(Sigma, 'lower');
    w = L \ x;
elseif ndims(Sigma)==2 && isequal(size(Sigma),[4 3])
    % per sat per axis std (nT)
    s = Sigma(iSat,:).';
    w = x ./ s;
elseif ndims(Sigma)==3 && isequal(size(Sigma),[3 3 4])
    L = chol(Sigma(:,:,iSat),'lower');
    w = L \ x;
else
    error('Unsupported Sigma format.');
end
end

%% ---- noise sampler under H0: N(0, Sigma) in nT ----
function E = sample_noise(nSat, Sigma)
E = zeros(nSat, 3);
if isscalar(Sigma)
    s = sqrt(Sigma);
    E = s * randn(nSat,3);
elseif isequal(size(Sigma),[3 3])
    L = chol(Sigma,'lower');
    for i=1:nSat
        E(i,:) = (L * randn(3,1)).';
    end
elseif ndims(Sigma)==2 && isequal(size(Sigma),[4 3])
    E = randn(nSat,3) .* Sigma;
elseif ndims(Sigma)==3 && isequal(size(Sigma),[3 3 4])
    for i=1:nSat
        L = chol(Sigma(:,:,i),'lower');
        E(i,:) = (L * randn(3,1)).';
    end
else
    error('Unsupported Sigma format.');
end
end

%% ---- strip time column if present ----
function X = take_xyz(Xin)
if size(Xin,2) == 4
    X = Xin(:,2:4);
elseif size(Xin,2) == 3
    X = Xin;
else
    error('Input must be Nx3 or Nx4 (time + xyz).');
end
end
