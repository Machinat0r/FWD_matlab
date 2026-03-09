function total = combine_noise_folder(folderPath, outFile)
% total = combine_noise_folder(folderPath, outFile)
% ------------------------------------------------------------
% Combine many daily "noise" structs saved in xxx.mat (xxx numeric).
% Each .mat contains a single variable: noise
%
% Pooling rule:
%   If cov computed with MATLAB cov() on M samples, it uses 1/(M-1).
%   So pooled covariance = sum_i (cov_i * (M_i-1)) / sum_i (M_i-1)
%
% This function combines:
%   - total.SigmaCov_all  (3x3)  pooled across all spacecraft & days
%   - total.SigmaCov_sc   (3x3x4) pooled per spacecraft (if available)
%   - total.sigma_scalar  (nT)    sqrt(mean(diag(SigmaCov_all)))
%   - total.sigma_xyz     (1x3)   std per component from SigmaCov_all
%   - total.sigma_sc_xyz  (4x3)   std per spacecraft/component (if available)
%
% Inputs:
%   folderPath : path containing *.mat files
%   outFile    : (optional) save combined result to this .mat file
%
% Example:
%   total = combine_noise_folder('D:\noise_daily\', 'total_noise.mat');

if nargin < 1 || isempty(folderPath)
    folderPath = pwd;
end
if nargin < 2
    outFile = '';
end

files = dir(fullfile(folderPath, '*.mat'));
if isempty(files)
    error('No .mat files found in folder: %s', folderPath);
end

% --- sort by numeric filename if possible (xxx.mat where xxx is number)
names = {files.name};
nums = nan(size(names));
for i = 1:numel(names)
    [~, base, ~] = fileparts(names{i});
    v = str2double(base);
    if isfinite(v), nums(i) = v; end
end
if any(isfinite(nums))
    [~, idx] = sort(nums);
    files = files(idx);
end

% --- pooled accumulators
S_all = zeros(3,3);      dof_all = 0;   % for SigmaCov_all
S_sc  = zeros(3,3,4);    dof_sc  = zeros(4,1); % for SigmaCov_sc

nUsed = 0;
badFiles = strings(0);

for k = 1:numel(files)
    fp = fullfile(folderPath, files(k).name);
    D = load(fp);

    if ~isfield(D, 'noise')
        badFiles(end+1) = string(files(k).name); %#ok<AGROW>
        continue;
    end
    noise = D.noise;

    % ---- combine SigmaCov_all (preferred)
    if isfield(noise, 'SigmaCov_all') && isequal(size(noise.SigmaCov_all), [3 3])
        C = noise.SigmaCov_all;

        % infer sample count used in cov
        % In estimate_noise_strength_4sc: Xall is (4N)x3, so M = 4N
        if isfield(noise, 'N') && ~isempty(noise.N) && isfinite(noise.N) && noise.N > 1
            M = 4 * noise.N;
            dof = max(M - 1, 1);
        else
            % fallback if N missing: treat as equal weight
            dof = 1;
        end

        S_all = S_all + C * dof;
        dof_all = dof_all + dof;
        nUsed = nUsed + 1;
    end

    % ---- combine SigmaCov_sc (if available)
    if isfield(noise, 'SigmaCov_sc') && isequal(size(noise.SigmaCov_sc), [3 3 4]) ...
            && isfield(noise, 'N') && ~isempty(noise.N) && isfinite(noise.N) && noise.N > 1
        % per spacecraft cov computed on N samples => dof = N-1
        dof_day_sc = max(noise.N - 1, 1);
        for s = 1:4
            S_sc(:,:,s) = S_sc(:,:,s) + noise.SigmaCov_sc(:,:,s) * dof_day_sc;
            dof_sc(s)   = dof_sc(s)   + dof_day_sc;
        end
    end
end

if dof_all <= 0
    error('No valid "noise.SigmaCov_all" found in the folder.');
end

% --- final pooled covariance
SigmaCov_all = S_all / dof_all;

total = struct();
total.folderPath = folderPath;
total.nFilesFound = numel(files);
total.nFilesUsed = nUsed;
total.badFiles = badFiles;

total.SigmaCov_all = SigmaCov_all;
total.sigma_xyz = sqrt(max(diag(SigmaCov_all), 0)).';           % 1x3, nT
total.sigma_scalar = sqrt(mean(max(diag(SigmaCov_all), 0)));    % nT
total.SigmaVar_scalar = total.sigma_scalar^2;                   % nT^2
total.dof_all = dof_all;

% --- per spacecraft results if available
if all(dof_sc > 0)
    SigmaCov_sc = zeros(3,3,4);
    sigma_sc_xyz = zeros(4,3);
    for s = 1:4
        SigmaCov_sc(:,:,s) = S_sc(:,:,s) / dof_sc(s);
        sigma_sc_xyz(s,:)  = sqrt(max(diag(SigmaCov_sc(:,:,s)), 0)).';
    end
    total.SigmaCov_sc = SigmaCov_sc;
    total.sigma_sc_xyz = sigma_sc_xyz;
    total.dof_sc = dof_sc;
end

% --- optionally save
if ~isempty(outFile)
    save(outFile, 'total');
end

end
