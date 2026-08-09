function state = Voyager_SPICE_Position(epoch, spacecraft, varargin)
% Voyager_SPICE_Position Evaluate a Voyager state from JPL SPICE kernels.
%
%   state = Voyager_SPICE_Position(epoch, spacecraft)
%   state = Voyager_SPICE_Position(..., 'KernelDirectory', folder)
%
% Inputs
%   epoch       Nonempty datetime array. Zoned datetimes are converted to
%               UTC. An unzoned datetime is interpreted as UTC.
%   spacecraft  Numeric scalar equal to 1 or 2.
%
% Name-value input
%   KernelDirectory
%       Directory containing the downloaded naif0012.tls, DE440,
%       broad/x2100 trajectory, and encounter-reconstruction BSP files for
%       the selected spacecraft. The default is:
%         Z:\SPART-WORK\Data\Voyager\voyagerN\ephemeris\spice\kernels
%
% Output fields
%   time_utc                 N-by-1 datetime array in UTC
%   ephemeris_time_tdb_s     N-by-1 TDB seconds past J2000
%   position_j2000_km        N-by-3 geometric position relative to the Sun
%   velocity_j2000_kms       N-by-3 geometric velocity relative to the Sun
%   one_way_light_time_s     N-by-1 geometric range divided by light speed
%   spacecraft, naif_id, observer, reference_frame, kernel_files
%
% Load priority is broad trajectory, DE440 planets, 2022 x2100 trajectory,
% then encounter reconstructions. Later kernels therefore supply the
% highest-quality available flyby solution in overlapping intervals.
% No light-time or stellar-aberration correction is applied.
%
% MICE example
%   t = datetime(2026, 7, 24, 0, 0, 0, 'TimeZone', 'UTC');
%   s = Voyager_SPICE_Position(t, 1);
%   disp(s.position_j2000_km)
%
% This function requires NASA NAIF MICE for MATLAB:
%   https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html

validateattributes(epoch, {'datetime'}, {'nonempty'}, mfilename, 'epoch', 1);
if any(isnat(epoch(:)))
    error('Voyager:InvalidEpoch', ...
        'epoch must not contain NaT values.');
end

validateattributes(spacecraft, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer'}, ...
    mfilename, 'spacecraft', 2);
if ~ismember(double(spacecraft), [1, 2])
    error('Voyager:InvalidSpacecraft', ...
        'spacecraft must be the numeric value 1 or 2.');
end
spacecraft = double(spacecraft);

defaultKernelDirectory = fullfile( ...
    'Z:\SPART-WORK\Data\Voyager', sprintf('voyager%d', spacecraft), ...
    'ephemeris', 'spice', 'kernels');

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'KernelDirectory', defaultKernelDirectory, ...
    @isTextScalar);
parse(parser, varargin{:});
kernelDirectory = char(parser.Results.KernelDirectory);
if isempty(strtrim(kernelDirectory))
    error('Voyager:InvalidKernelDirectory', ...
        'KernelDirectory must be a nonempty folder path.');
end

requiredMiceFunctions = { ...
    'cspice_furnsh', 'cspice_unload', ...
    'cspice_str2et', 'cspice_spkezr'};
missingMiceFunctions = requiredMiceFunctions(cellfun( ...
    @(name) exist(name, 'file') == 0, requiredMiceFunctions));
if ~isempty(missingMiceFunctions)
    error('Voyager:MICENotInstalled', [ ...
        'NASA NAIF MICE for MATLAB is not installed or is not on the ' ...
        'MATLAB path. Missing function(s): %s.\nDownload MICE from ' ...
        'https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html and add its ' ...
        'mice/src/mice and mice/lib directories to the MATLAB path.'], ...
        strjoin(missingMiceFunctions, ', '));
end

if spacecraft == 1
    spacecraftName = 'Voyager 1';
    naifId = -31;
    mergedFilename = 'Voyager_1.a54206u_V0.2_merged.bsp';
    extendedFilename = 'vgr1.x2100.bsp';
    encounterFilenames = {'vgr1_jup230.bsp', 'vgr1.sat427.bsp'};
else
    spacecraftName = 'Voyager 2';
    naifId = -32;
    mergedFilename = 'Voyager_2.m05016u.merged.bsp';
    extendedFilename = 'vgr2.x2100.bsp';
    encounterFilenames = { ...
        'vgr2_jup230.bsp', 'vgr2.sat427.bsp', ...
        'vgr2.ura182.bsp', 'vgr2_nep097.bsp'};
end

kernelFiles = { ...
    fullfile(kernelDirectory, 'naif0012.tls'), ...
    fullfile(kernelDirectory, mergedFilename), ...
    fullfile(kernelDirectory, 'de440s.bsp'), ...
    fullfile(kernelDirectory, extendedFilename)};
kernelFiles = [kernelFiles, cellfun( ...
    @(name) fullfile(kernelDirectory, name), encounterFilenames, ...
    'UniformOutput', false)];
missingKernel = ~cellfun(@isfile, kernelFiles);
if any(missingKernel)
    missingText = sprintf('  %s\n', kernelFiles{missingKernel});
    error('Voyager:SPICEKernelMissing', [ ...
        'Required SPICE kernel file(s) are missing:\n%s' ...
        'Download the spice_spk product, or pass the folder containing ' ...
        'the kernels with the KernelDirectory name-value input.'], ...
        missingText);
end

timeUtc = epoch(:);
% Assigning UTC to an unzoned datetime preserves its displayed clock time
% and explicitly records the UTC assumption. For a zoned datetime it
% converts the same instants to UTC.
timeUtc.TimeZone = 'UTC';
timeUtc.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS';
utcText = cellstr(char(timeUtc));
utcText = cellfun(@(value) [strtrim(value), ' UTC'], ...
    utcText, 'UniformOutput', false);

loadedFiles = {};
try
    % Kernel load order controls SPK segment priority in overlapping spans.
    for ii = 1:numel(kernelFiles)
        cspice_furnsh(kernelFiles{ii});
        loadedFiles{end + 1} = kernelFiles{ii}; %#ok<AGROW>
    end
catch cause
    unloadKernels(loadedFiles);
    exception = MException('Voyager:SPICEKernelLoadFailed', ...
        'MICE could not load the required kernels from %s.', ...
        kernelDirectory);
    exception = addCause(exception, cause);
    throwAsCaller(exception);
end

try
    ephemerisTime = cspice_str2et(utcText);
    [stateVector, lightTime] = cspice_spkezr( ...
        sprintf('%d', naifId), ephemerisTime, ...
        'J2000', 'NONE', 'SUN');
catch cause
    unloadKernels(loadedFiles);
    exception = MException('Voyager:SPICEEvaluationFailed', [ ...
        'MICE could not evaluate the %s state. Check that all epochs lie ' ...
        'within the supplied BSP coverage and inspect the attached MICE ' ...
        'error for details.'], spacecraftName);
    exception = addCause(exception, cause);
    throwAsCaller(exception);
end
unloadKernels(loadedFiles);

state = struct;
state.time_utc = timeUtc;
state.ephemeris_time_tdb_s = ephemerisTime(:);
state.position_j2000_km = stateVector(1:3, :).';
state.velocity_j2000_kms = stateVector(4:6, :).';
state.one_way_light_time_s = lightTime(:);
state.spacecraft = spacecraftName;
state.naif_id = naifId;
state.observer = 'SUN';
state.reference_frame = 'J2000';
state.aberration_correction = 'NONE';
state.kernel_files = kernelFiles(:);
end

function tf = isTextScalar(value)
tf = (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value));
end

function unloadKernels(kernelFiles)
% Reverse order mirrors loading and restores any earlier kernel priority.
for ii = numel(kernelFiles):-1:1
    try
        cspice_unload(kernelFiles{ii});
    catch
        % Preserve the primary load/evaluation result if cleanup fails.
    end
end
end
