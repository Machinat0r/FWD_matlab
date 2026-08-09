function [filenames, fileInfo] = VoyagerFilenames(Date, varargin)
% VoyagerFilenames Discover NASA/SPDF Voyager CDF files without downloading.
%
% Example
%   [names, info] = VoyagerFilenames('2018-01-01/2018-02-28', ...
%       'Spacecraft', [1 2], 'Products', 'highres');

parser = inputParser;
parser.CaseSensitive = false;

addRequired(parser, 'Date', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'Products', 'coho1hr,position1day', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(parser, 'PythonScript', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
parse(parser, Date, varargin{:});

scratchRoot = fullfile(tempdir, 'Voyager_list_only');
fileInfo = VoyagerFilesDownload_NAS(char(parser.Results.Date), scratchRoot, ...
    'Spacecraft', parser.Results.Spacecraft, ...
    'Products', parser.Results.Products, ...
    'ListOnly', true, ...
    'PythonScript', char(parser.Results.PythonScript), ...
    'PythonExe', char(parser.Results.PythonExe));

if isempty(fileInfo)
    filenames = {};
elseif isstruct(fileInfo) && isfield(fileInfo, 'filename')
    filenames = {fileInfo.filename}.';
elseif isstruct(fileInfo) && isfield(fileInfo, 'files')
    files = fileInfo.files;
    if isstruct(files) && isfield(files, 'filename')
        filenames = {files.filename}.';
    else
        filenames = {};
    end
else
    filenames = {};
end
end
