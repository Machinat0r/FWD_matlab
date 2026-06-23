function ACEFilesDownload_NAS(Date, OutputFiles_dir, varargin)
% ACEFilesDownload_NAS Download ACE CDF files from NASA CDAWeb.
%
% Example:
%   ACEFilesDownload_NAS('2026-01-01/2026-01-03', './ACE/', ...
%       'product', 'mfi_h0', 'Threads', 8, 'CheckSize', 1)
%
% Multiple products can be requested with comma-separated names:
%   ACEFilesDownload_NAS('2024-05-10/2024-05-12', './ACE/', ...
%       'product', 'mfi_h0,swe_h0')
%
% MMS-style aliases are also accepted:
%   ACEFilesDownload_NAS('2026-01-01/2026-01-01', './ACE/', ...
%       'inst', 'mfi', 'drm', 'h0')

parser = inputParser;
parser.CaseSensitive = false;

addRequired(parser, 'Date', @(x) ischar(x) || isstring(x));
addRequired(parser, 'OutputFiles_dir', @(x) ischar(x) || isstring(x));
addParameter(parser, 'product', 'mfi_h0', @(x) ischar(x) || isstring(x));
addParameter(parser, 'inst', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'drm', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'CheckSize', 1, @isnumeric);
addParameter(parser, 'Threads', 8, @isnumeric);
addParameter(parser, 'KeepTree', 0, @isnumeric);
addParameter(parser, 'ListOnly', 0, @isnumeric);
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
parse(parser, Date, OutputFiles_dir, varargin{:});

Date = char(parser.Results.Date);
OutputFiles_dir = char(parser.Results.OutputFiles_dir);
Product = char(parser.Results.product);
Instrument = char(parser.Results.inst);
Mode = char(parser.Results.drm);
CheckSize = parser.Results.CheckSize;
Threads = parser.Results.Threads;
KeepTree = parser.Results.KeepTree;
ListOnly = parser.Results.ListOnly;
PythonScript = char(parser.Results.PythonScript);

if ~isempty(Instrument)
    if isempty(Mode)
        Product = Instrument;
    else
        Product = [Instrument, '_', Mode];
    end
end

if ~isfolder(OutputFiles_dir)
    mkdir(OutputFiles_dir);
end

if isempty(PythonScript)
    func_path = mfilename('fullpath');
    func_dir = fileparts(func_path);
    PythonScript = fullfile(func_dir, 'download_ace_files_new.py');
end

if ~isfile(PythonScript)
    error('Cannot find Python backend: %s', PythonScript);
end

disp('========== ACE download started ==========')
cmd = ['python3 ', shellquote(PythonScript), ...
    ' --date ', shellquote(Date), ...
    ' --product ', shellquote(Product), ...
    ' --out ', shellquote(OutputFiles_dir), ...
    ' --threads ', num2str(Threads), ...
    ' --check-size ', num2str(CheckSize), ...
    ' --keep-tree ', num2str(KeepTree)];

if ListOnly
    cmd = [cmd, ' --list-only'];
end

status = system(cmd);
if status ~= 0
    error('ACE download failed. Python backend exit status: %d', status);
end
disp('========== ACE download finished ==========')
end

function out = shellquote(in)
in = char(in);
q = char(39);
dq = char(34);
in = strrep(in, q, [q, dq, q, dq, q]);
out = [q, in, q];
end
