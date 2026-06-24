function [filenames, fileInfo] = ACEFilenames(Date, varargin)
% ACEFilenames Crawl ACE CDAWeb directory listings and return matched files.
%
% Example:
%   files = ACEFilenames('2026-01-01/2026-01-03', 'product', 'mfi_h0')
%   [files, info] = ACEFilenames('2024-05-10/2024-05-12', ...
%       'product', 'mfi_h0,swe_h0')

parser = inputParser;
parser.CaseSensitive = false;

addRequired(parser, 'Date', @(x) ischar(x) || isstring(x));
addParameter(parser, 'product', 'mfi_h0', @(x) ischar(x) || isstring(x));
addParameter(parser, 'inst', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'drm', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
parse(parser, Date, varargin{:});

Date = char(parser.Results.Date);
Product = char(parser.Results.product);
Instrument = char(parser.Results.inst);
Mode = char(parser.Results.drm);
PythonScript = char(parser.Results.PythonScript);
PythonExe = char(parser.Results.PythonExe);

if ~isempty(Instrument)
    if isempty(Mode)
        Product = Instrument;
    else
        Product = [Instrument, '_', Mode];
    end
end

if isempty(PythonScript)
    func_path = mfilename('fullpath');
    func_dir = fileparts(func_path);
    PythonScript = fullfile(func_dir, 'download_ace_files_new.py');
end

if ~isfile(PythonScript)
    error('Cannot find Python backend: %s', PythonScript);
end

pythonCmd = resolvePythonCommand(PythonExe);
cmd = [pythonCmd, ' ', shellquote(PythonScript), ...
    ' --date ', shellquote(Date), ...
    ' --product ', shellquote(Product), ...
    ' --list-only --json'];

[status, cmdout] = system(cmd);
if status ~= 0
    error('ACE filename crawl failed. Python backend exit status: %d\n%s', status, cmdout);
end

fileInfo = jsondecode(cmdout);
if isempty(fileInfo)
    filenames = {};
elseif isstruct(fileInfo)
    filenames = {fileInfo.filename};
    filenames = filenames(:);
else
    filenames = {};
end
end

function pythonCmd = resolvePythonCommand(pythonExe)
pythonExe = strtrim(char(pythonExe));
if isempty(pythonExe)
    pythonExe = strtrim(getenv('ACE_PYTHON'));
end

if ~isempty(pythonExe)
    pythonCmd = formatPythonCommand(pythonExe);
    if pythonCommandWorks(pythonCmd)
        return
    end
    error('PythonExe is not a working Python 3 interpreter: %s', pythonExe);
end

try
    pe = pyenv;
    if ~isempty(pe.Executable)
        pythonCmd = shellquote(char(pe.Executable));
        if pythonCommandWorks(pythonCmd)
            return
        end
    end
catch
end

candidateFiles = knownPythonFiles();
for ii = 1:numel(candidateFiles)
    if isfile(candidateFiles{ii})
        pythonCmd = shellquote(candidateFiles{ii});
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

error(['Cannot find a working Python 3 interpreter. Set MATLAB pyenv("Version", ' ...
    '"C:\path\to\python.exe"), pass ''PythonExe'', ''C:\path\to\python.exe'', ' ...
    'or set the ACE_PYTHON environment variable.']);
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

function out = formatPythonCommand(pythonExe)
if isfile(pythonExe) || contains(pythonExe, '\') || contains(pythonExe, '/')
    out = shellquote(pythonExe);
else
    out = pythonExe;
end
end

function out = shellquote(in)
in = char(in);
if ispc
    dq = char(34);
    in = strrep(in, dq, [dq, dq]);
    out = [dq, in, dq];
else
    q = char(39);
    dq = char(34);
    in = strrep(in, q, [q, dq, q, dq, q]);
    out = [q, in, q];
end
end
