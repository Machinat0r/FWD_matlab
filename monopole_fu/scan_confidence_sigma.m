function out = scan_confidence_sigma(Rin, Bin, K, nSigma, BGnoise, seed)
% out = scan_confidence_sigma(Rin, Bin, K, nSigma, seed)
% --------------------------------------------------------
% 1) Input: one event R,B (units: R in km, B in nT)
% 2) Scan noise std sigma from 0.1 to 10 nT
% 3) Plot:
%    (a) observed B RMS vs noise B RMS across sigma
%    (b) confidence vs sigma
%
% Requires:
%   - Your main function: MentoCarlo_estimate(R,B,Sigma,K,opts)
%     where scalar Sigma is interpreted as variance (nT^2).
%
% Parameters:
%   K      : MC replicates per sigma (e.g., 2000, 5000, 10000)
%   nSigma : number of sigma points between 0.1 and 10 (default 25~40)
%   seed   : base RNG seed for reproducibility

if nargin < 3 || isempty(K),      K = 5000; end
if nargin < 4 || isempty(nSigma), nSigma = 30; end
if nargin < 5 || isempty(nSigma), BGnoise = 4; end
if nargin < 6 || isempty(seed),   seed = 1; end

R = take_xyz(Rin);   % 4x3
B = take_xyz(Bin);   % 4x3

if ~isequal(size(R),[4 3]) || ~isequal(size(B),[4 3])
    error('R and B must be 4x3 (or 4x4 with time column).');
end

% ---- sigma grid: log-spaced (recommended for 0.1~10)
sigma_list = logspace(log10(3.93), log10(3.93), nSigma);  % nT

% ---- observed field amplitude (constant across sigma)
Bmag = vecnorm(B,2,2);          % |B_i| for each satellite
% % Bmag_rms = sqrt(mean(Bmag.^2)); % RMS of |B| across 4 sats
% sigma_eq = Bmag_rms / sqrt(3);   % nT, equivalent per-component noise std
sigma_eq = BGnoise;
Bmag_rms = sigma_eq * sqrt(3);

% ---- arrays to store results
pvals = nan(nSigma,1);
confs = nan(nSigma,1);
Zs    = nan(nSigma,1);

% For plot 1: expected noise RMS(|e|) under independent components
% If e_x,e_y,e_z ~ N(0, sigma^2), then E[|e|^2] = 3*sigma^2
noise_mag_rms_theory = sqrt(3) .* sigma_list(:); % nT

% Optional: also show one random noise realization RMS(|e|) per sigma
noise_mag_rms_sample = nan(nSigma,1);

fprintf('Scanning sigma in [0.1, 10] nT with %d points, K=%d MC each...\n', nSigma, K);

parfor i = 1:nSigma
    sigma = sigma_list(i);     % nT (std)
    SigmaVar = sigma^2;        % nT^2 (variance) -- IMPORTANT

    % generate one sample noise field just for visualization
    rng(seed + 100000 + i);
    En = sigma * randn(4,3);
    noise_mag_rms_sample(i) = sqrt(mean(vecnorm(En,2,2).^2));

    % call your main routine to compute confidence
    opts = struct();
    opts.seed = seed + i;          % different seed per sigma
    opts.saveT = false;
    % opts.fitMethod = "profile_fminsearch"; % if your main supports it
    % opts.Qsign = 0;               % or +1/-1 if you constrain Q sign

    res = MentoCarlo_estimate(R, B, SigmaVar, K);

    pvals(i) = res.pvalue;
    confs(i) = res.confidence;
    Zs(i)    = res.Z;

    fprintf('sigma=%.4g nT: p=%.3e, conf=%.6f, Z=%.3f\n', sigma, pvals(i), confs(i), Zs(i));
end

% ---------------- Plot 1: Observed vs Noise amplitude ----------------
figure('Name','Observed B vs Noise B across sigma');
hold on; grid on; box on;

% Observed is constant line
semilogx(sigma_list, Bmag_rms * ones(size(sigma_list)), 'LineWidth', 2);

% Noise RMS magnitude: theory and one sample
semilogx(sigma_list, noise_mag_rms_theory, 'LineWidth', 2);
semilogx(sigma_list, noise_mag_rms_sample, '--', 'LineWidth', 1.5);

xlabel('\sigma (nT)  [noise std per component]');
ylabel('RMS magnitude (nT)');
legend('Observed RMS(|B|)', 'Noise RMS(|e|) theory (\surd3\sigma)', 'Noise RMS(|e|) one sample', ...
       'Location','best');
title('Observed field strength vs assumed noise level');
set(gca,'XScale','log');  

% ---------------- Plot 2: Confidence curve ----------------
figure('Name','Confidence vs sigma');
hold on; grid off; box on;

semilogx(sigma_list, confs, 'LineWidth', 2);
ylim([0 1]);
xlabel('\sigma (nT)  [noise std per component]');
ylabel('Confidence = 1 - p-value');
% title(sprintf('Event confidence vs noise std (K=%d per point)', K));

set(gca,'XScale','log');

% 竖虚线：sigma_eq = observed RMS(|B|)/sqrt(3)
xline(sigma_eq, '--', 'LineWidth', 1.5);
text(sigma_eq, 0.05, sprintf('  \\sigma_{eq}=%.3g nT', sigma_eq), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom');


% package outputs
out.sigma_list = sigma_list(:);
out.pvalue = pvals;
out.confidence = confs;
out.Z = Zs;
out.Bmag_rms = Bmag_rms;
out.noise_mag_rms_theory = noise_mag_rms_theory;
out.noise_mag_rms_sample = noise_mag_rms_sample;
end

% ---------- helper: strip time column ----------
function X = take_xyz(Xin)
if size(Xin,2) == 4
    X = Xin(:,2:4);
elseif size(Xin,2) == 3
    X = Xin;
else
    error('Input must be 4x3 or 4x4 (time + xyz).');
end
end
