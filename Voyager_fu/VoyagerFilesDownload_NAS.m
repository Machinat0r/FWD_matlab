function report = VoyagerFilesDownload_NAS(Date, OutputFilesDir, varargin)
% VoyagerFilesDownload_NAS Download Voyager science files from NASA archives.
%
% This function mirrors the ACEFilesDownload_NAS MATLAB interface while its
% Python backend performs discovery, retries, staging and atomic publishing.
%
% Example
%   report = VoyagerFilesDownload_NAS( ...
%       '2018-01-01/2020-12-31', ...
%       'Z:\SPART-WORK\Data\Voyager', ...
%       'Spacecraft', 2, ...
%       'Products', 'coho1hr,position1day', ...
%       'Threads', 5);

parser = inputParser;
parser.CaseSensitive = false;

addRequired(parser, 'Date', @(x) ischar(x) || isstring(x));
addRequired(parser, 'OutputFilesDir', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'Products', 'coho1hr,position1day', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(parser, 'StageDir', fullfile(tempdir, 'Voyager_staging'), ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'ManifestName', 'voyager_download_manifest.json', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'Threads', 5, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && x <= 5 && fix(x) == x);
addParameter(parser, 'CheckSize', true, @isBooleanScalar);
addParameter(parser, 'Force', false, @isBooleanScalar);
addParameter(parser, 'ListOnly', false, @isBooleanScalar);
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
parse(parser, Date, OutputFilesDir, varargin{:});

opts = parser.Results;
dateText = char(opts.Date);
outputDir = char(opts.OutputFilesDir);
stageDir = char(opts.StageDir);
pythonScript = char(opts.PythonScript);
manifestName = validateManifestName(opts.ManifestName);

if isempty(pythonScript)
    pythonScript = fullfile(fileparts(mfilename('fullpath')), ...
        'download_voyager_files.py');
end
if ~isfile(pythonScript)
    error('Voyager:PythonBackendMissing', ...
        'Cannot find Python backend: %s', pythonScript);
end

spacecraft = unique(double(opts.Spacecraft(:).'), 'stable');
if isempty(spacecraft) || any(~ismember(spacecraft, [1 2]))
    error('Voyager:InvalidSpacecraft', ...
        'Spacecraft must contain only 1 and/or 2.');
end
spacecraftText = strjoin(arrayfun(@num2str, spacecraft, ...
    'UniformOutput', false), ',');

if iscell(opts.Products)
    products = strjoin(cellfun(@char, opts.Products, ...
        'UniformOutput', false), ',');
else
    products = char(opts.Products);
end

pythonCmd = resolvePythonCommand(char(opts.PythonExe));
cmd = [pythonCmd, ' ', shellquote(pythonScript), ...
    ' --date ', shellquote(dateText), ...
    ' --spacecraft ', shellquote(spacecraftText), ...
    ' --product ', shellquote(products), ...
    ' --out ', shellquote(outputDir), ...
    ' --stage ', shellquote(stageDir), ...
    ' --manifest-name ', shellquote(manifestName), ...
    ' --threads ', num2str(round(opts.Threads)), ...
    ' --check-size ', num2str(logical(opts.CheckSize)), ...
    ' --force ', num2str(logical(opts.Force))];

if logical(opts.ListOnly)
    cmd = [cmd, ' --list-only --json'];
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('Voyager:FileDiscoveryFailed', ...
            'Voyager file discovery failed (exit %d).\n%s', status, cmdout);
    end
    try
        report = jsondecode(strtrim(cmdout));
    catch exception
        error('Voyager:InvalidBackendJSON', ...
            'Cannot decode the downloader file list: %s\nBackend output:\n%s', ...
            exception.message, cmdout);
    end
    return
end

fprintf('========== Voyager download started ==========\n');
fprintf('Date:       %s\n', dateText);
fprintf('Spacecraft: %s\n', spacecraftText);
fprintf('Products:   %s\n', products);
fprintf('Data root:  %s\n', outputDir);
fprintf('Staging:    %s\n', stageDir);
fprintf('Manifest:   %s\n', manifestName);

status = system(cmd);
if status ~= 0
    error('Voyager:DownloadFailed', ...
        'Voyager download failed. Python backend exit status: %d', status);
end
fprintf('========== Voyager download finished ==========\n');

manifestFile = fullfile(outputDir, manifestName);
if isfile(manifestFile)
    try
        report = jsondecode(fileread(manifestFile));
    catch exception
        warning('Voyager:ManifestReadFailed', ...
            'Download completed, but the manifest could not be read: %s', ...
            exception.message);
        report = struct('status', 'completed', ...
            'manifest', manifestFile);
    end
else
    report = struct('status', 'completed', 'manifest', '');
end
end

function value = validateManifestName(input)
value = strtrim(char(input));
[folder, base, extension] = fileparts(value);
invalidCharacters = ['\', '/', ':', '*', '?', '"', '<', '>', '|'];
if isempty(value) || ~isempty(folder) || isempty([base, extension]) || ...
        any(ismember(value, invalidCharacters))
    error('Voyager:InvalidManifestName', ...
        'ManifestName must be a filename without directory components.');
end
end

function tf = isBooleanScalar(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isfinite(double(value)) && ismember(double(value), [0 1]);
end

function pythonCmd = resolvePythonCommand(pythonExe)
pythonExe = strtrim(char(pythonExe));
if isempty(pythonExe)
    pythonExe = strtrim(getenv('VOYAGER_PYTHON'));
end

if ~isempty(pythonExe)
    pythonCmd = formatPythonCommand(pythonExe);
    if pythonCommandWorks(pythonCmd)
        return
    end
    error('Voyager:PythonInvalid', ...
        'PythonExe is not a working Python 3 interpreter: %s', pythonExe);
end

try
    environment = pyenv;
    if ~isempty(environment.Executable)
        pythonCmd = shellquote(char(environment.Executable));
        if pythonCommandWorks(pythonCmd)
            return
        end
    end
catch
end

knownFiles = knownPythonFiles();
for ii = 1:numel(knownFiles)
    if isfile(knownFiles{ii})
        pythonCmd = shellquote(knownFiles{ii});
        if pythonCommandWorks(pythonCmd)
            return
        end
    end
end

if ispc
    candidates = {'python', 'python3', 'py -3', 'py'};
else
    candidates = {'python3', 'python'};
end
for ii = 1:numel(candidates)
    pythonCmd = candidates{ii};
    if pythonCommandWorks(pythonCmd)
        return
    end
end

error('Voyager:PythonNotFound', [ ...
    'Cannot find Python 3. Set MATLAB pyenv("Version", executable), ', ...
    'pass ''PythonExe'', or set the VOYAGER_PYTHON environment variable.']);
end

function files = knownPythonFiles()
files = {};
if ispc
    userProfile = getenv('USERPROFILE');
    if ~isempty(userProfile)
        files{end + 1} = fullfile(userProfile, '.cache', 'codex-runtimes', ...
            'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
    end
end
end

function ok = pythonCommandWorks(pythonCmd)
testCode = 'import sys; raise SystemExit(0 if sys.version_info[0] >= 3 else 1)';
[status, ~] = system([pythonCmd, ' -c ', shellquote(testCode)]);
ok = (status == 0);
end

function output = formatPythonCommand(pythonExe)
if isfile(pythonExe) || contains(pythonExe, '\') || contains(pythonExe, '/')
    output = shellquote(pythonExe);
else
    output = pythonExe;
end
end

function output = shellquote(input)
input = char(input);
if ispc
    dq = char(34);
    input = strrep(input, dq, [dq, dq]);
    output = [dq, input, dq];
else
    q = char(39);
    dq = char(34);
    input = strrep(input, q, [q, dq, q, dq, q]);
    output = [q, input, q];
end
end
