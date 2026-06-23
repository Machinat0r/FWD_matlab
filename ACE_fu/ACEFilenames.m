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
parse(parser, Date, varargin{:});

Date = char(parser.Results.Date);
Product = char(parser.Results.product);
Instrument = char(parser.Results.inst);
Mode = char(parser.Results.drm);
PythonScript = char(parser.Results.PythonScript);

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

cmd = ['python3 ', shellquote(PythonScript), ...
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

function out = shellquote(in)
in = char(in);
q = char(39);
dq = char(34);
in = strrep(in, q, [q, dq, q, dq, q]);
out = [q, in, q];
end
