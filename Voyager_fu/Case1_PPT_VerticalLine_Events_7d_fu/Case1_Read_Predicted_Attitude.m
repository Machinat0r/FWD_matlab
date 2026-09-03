function attitude = Case1_Read_Predicted_Attitude(epoch, cfg)
%Case1_Read_Predicted_Attitude Voyager 1 predicted bus orientation from CK.
%   EPOCH is an explicit UTC datetime vector. Pass the original LECP Epoch
%   under the current approved policy. The result is an instantaneous prediction,
%   not a daily average or telemetry-reconstructed attitude.
%
%   C_SC_to_RTN maps SC_BUS column vectors into the COHO RTN basis.
%   No LECP sector direction is assumed or produced by this function.
%   Requires an empty SPICE kernel pool; never clears the caller's pool.

%% input and existing IRFU/MICE routines
if nargin < 2, cfg = Case1_Config; end
validateattributes(epoch, {'datetime'}, {'vector', 'nonempty'});
if any(isnat(epoch)), error('VoyagerCase1:AttitudeEpoch', 'NaT epoch.'); end
if ~cfg.PredictedAttitudeApproved
    error('VoyagerCase1:AttitudeApproval', ...
        'Official predicted attitude has not been approved.');
end
time = epoch(:);
time.TimeZone = 'UTC';
addpath(fullfile(cfg.IRFURoot, 'contrib', 'mice'));
if exist('mice', 'file') ~= 3
    error('VoyagerCase1:MICEMissing', 'IRFU MICE binary is unavailable.');
end
[~, poolHasVariables] = cspice_gnpool('*', 1, 1);
if cspice_ktotal('ALL') ~= 0 || poolHasVariables
    error('VoyagerCase1:KernelPoolNotEmpty', ...
        'Use an isolated MATLAB run with no preloaded SPICE kernels.');
end
files = Case1_Attitude_Files(cfg);
kernelFiles = files(files.Role ~= "COMMENT", :);
if any(~isfile(kernelFiles.SourceFile))
    error('VoyagerCase1:AttitudeKernelMissing', ...
        'Run Case1_Download_Attitude before reading predicted attitude.');
end
files.SHA256 = strings(height(files), 1);
for ii = 1:height(files)
    if isfile(files.SourceFile(ii))
        files.SHA256(ii) = Case1_File_SHA256(char(files.SourceFile(ii)));
    end
end

%% load exact kernels; cleanup only the files loaded by this call
cleanup = onCleanup(@() unloadFiles(kernelFiles.SourceFile));
for ii = 1:height(kernelFiles)
    cspice_furnsh(char(kernelFiles.SourceFile(ii)));
end
time.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSSSSS';
utcText = cellstr(char(time));
utcText = cellfun(@(x) [x, ' UTC'], utcText, 'UniformOutput', false);
et = cspice_str2et(utcText);
clock = cspice_sce2c(-31, et);
ckFile = char(files.SourceFile(files.Role == "CK"));
coverage = cspice_ckcov(ckFile, -31000, false, ...
    'INTERVAL', 0, 'TDB', 100000);

%% allocate missing as NaN; no extrapolation or gap filling
nTime = numel(time);
C_SC_to_J2000 = nan(3, 3, nTime);
C_J2000_to_RTN = nan(3, 3, nTime);
C_SC_to_RTN = nan(3, 3, nTime);
found = false(nTime, 1);
clockReturned = nan(nTime, 1);
positionJ2000 = nan(nTime, 3);
status = repmat("outside predicted CK coverage", nTime, 1);
predictionStart = datetime(1990, 4, 16, 'TimeZone', 'UTC');

%% official type-3 rotation interpolation, evaluated at requested epochs
for ii = 1:nTime
    if time(ii) < predictionStart || ~cspice_wnelmd(et(ii), coverage)
        continue
    end
    [C, clockOut, hasAttitude] = cspice_ckgp( ...
        -31000, clock(ii), 0, 'J2000');
    if ~hasAttitude, continue; end
    % CK C maps J2000 -> SC_BUS. Its transpose maps SC_BUS -> J2000.
    C_SC_to_J2000(:, :, ii) = C.';
    clockReturned(ii) = clockOut;
    try
        [position, ~] = cspice_spkpos('-31', et(ii), ...
            'J2000', 'NONE', 'SUN');
    catch ME
        if contains(ME.message, 'SPICE(SPKINSUFFDATA)')
            status(ii) = "missing SPK coverage; RTN unavailable";
            continue
        end
        rethrow(ME)
    end
    positionJ2000(ii, :) = position.';
    R = position / norm(position);
    solarRotation = cspice_pxform('IAU_SUN', 'J2000', et(ii));
    sunNorth = solarRotation(:, 3);
    T = cross(sunNorth, R);
    if norm(T) < 1e-12
        error('VoyagerCase1:RTNDegenerate', 'RTN basis is undefined.');
    end
    T = T / norm(T);
    N = cross(R, T);
    Q = [R.'; T.'; N.'];
    C_J2000_to_RTN(:, :, ii) = Q;
    C_SC_to_RTN(:, :, ii) = Q * C.';
    % Numerical validation only; this is not a scientific data threshold.
    if norm(C*C.'-eye(3), 'fro') > 1e-10 || abs(det(C)-1) > 1e-10 ...
            || norm(Q*Q.'-eye(3), 'fro') > 1e-10 || det(Q) < 0
        error('VoyagerCase1:AttitudeRotation', 'Invalid rotation matrix.');
    end
    found(ii) = true;
    status(ii) = "official predicted attitude";
end

%% provenance and explicit scientific meaning
attitude = struct;
attitude.TimeUTC = time;
attitude.EphemerisTimeTDB = et(:);
attitude.SpacecraftClock = clock(:);
attitude.ReturnedSpacecraftClock = clockReturned;
attitude.Found = found;
attitude.Status = status;
attitude.Predicted = found;
attitude.C_SC_to_J2000 = C_SC_to_J2000;
attitude.C_J2000_to_RTN = C_J2000_to_RTN;
attitude.C_SC_to_RTN = C_SC_to_RTN;
attitude.PositionJ2000_km = positionJ2000;
attitude.CoverageTDB = coverage;
attitude.Sources = files;
attitude.Toolkit = cspice_tkvrsn('TOOLKIT');
attitude.Method = 'Official predicted CK; native type-3 rotation interpolation; tolerance=0';
attitude.RTNDefinition = 'R=Sun to V1; T=unit(SunNorth cross R); N=R cross T';
attitude.AberrationCorrection = 'NONE';
attitude.SectorGeometryApplied = false;
end

function unloadFiles(files)
for ii = numel(files):-1:1
    try
        cspice_unload(char(files(ii)));
    catch
        % Do not mask the original load/evaluation exception.
    end
end
end
