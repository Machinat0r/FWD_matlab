function report = Voyager_Download(varargin)
% Voyager_Download Download Voyager 1/2 position, magnetic-field and plasma data.
%
% The default products are NASA/SPDF COHO one-hour merged magnetic-field,
% plasma and ephemeris CDF files plus the independent one-day position CDF.
% Raw CDF files are retained and arranged in an MMS-style directory tree.
%
% Example: full available mission archive for both spacecraft
%   report = Voyager_Download;
%
% Example: a shorter interval for Voyager 2
%   report = Voyager_Download( ...
%       'Date', '2018-01-01/2020-12-31', ...
%       'Spacecraft', 2, ...
%       'DataRoot', 'Z:\SPART-WORK\Data\Voyager');
%
% Example: recommended high-resolution products
%   report = Voyager_Download('Products', ...
%       'mag2s,mag48s_vim,plasma_fine,position1hr');
%
% Example: both post-1991 unreviewed 1.92-second MAG sensors
%   report = Voyager_Download('Products', 'mag2s_unreviewed', ...
%       'Date', '1991-01-01/2030-12-31');
%
% Name-value inputs
%   Date        Inclusive date interval, START/END (default: mission to today)
%   Spacecraft  1, 2, or [1 2] (default: [1 2])
%   Products    Product names or both/highres/mag2s_unreviewed/all
%               (default: both)
%   DataRoot    Archive root (default: Z:\SPART-WORK\Data\Voyager)
%   StageDir    Local staging directory (default: tempdir/Voyager_staging)
%   ManifestName JSON manifest filename under DataRoot
%   Threads     Number of simultaneous file downloads (default: 5)
%   CheckSize   Compare local size with HTTP Content-Length (default: true)
%   Force       Redownload already complete files (default: false)
%   ListOnly    Discover files without downloading (default: false)
%   PythonExe   Explicit Python 3 executable (default: auto-detect)
%
% Important data caveat
%   Voyager 1 PLS stopped returning usable plasma data in 1980. Later
%   COHO CDF records therefore retain official fill values for density,
%   speed, and temperature where measurements are unavailable.
%   NASA marks hires1991_2030 primary/secondary MAG CDFs as generally not
%   science quality. Their point-quality flags must be retained and used.

parser = inputParser;
parser.CaseSensitive = false;

todayUtc = datetime('today', 'TimeZone', 'UTC');
defaultDate = ['1977-08-20/', char(string(todayUtc, 'yyyy-MM-dd'))];

addParameter(parser, 'Date', defaultDate, @(x) ischar(x) || isstring(x));
addParameter(parser, 'Spacecraft', [1 2], @isnumeric);
addParameter(parser, 'Products', 'coho1hr,position1day', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(parser, 'DataRoot', 'Z:\SPART-WORK\Data\Voyager', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'StageDir', fullfile(tempdir, 'Voyager_staging'), ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'ManifestName', 'voyager_download_manifest.json', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'Threads', 5, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 1 && x <= 5 && fix(x) == x);
addParameter(parser, 'CheckSize', true, @isBooleanScalar);
addParameter(parser, 'Force', false, @isBooleanScalar);
addParameter(parser, 'ListOnly', false, @isBooleanScalar);
addParameter(parser, 'PythonExe', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

opts = parser.Results;
spacecraft = unique(double(opts.Spacecraft(:).'), 'stable');
if isempty(spacecraft) || any(~ismember(spacecraft, [1 2]))
    error('Voyager:InvalidSpacecraft', ...
        'Spacecraft must contain only 1 and/or 2.');
end

products = normalizeProducts(opts.Products);
dataRoot = char(opts.DataRoot);
stageDir = char(opts.StageDir);
manifestName = validateManifestName(opts.ManifestName);

if ~logical(opts.ListOnly)
    if ~isfolder(dataRoot)
        [ok, message] = mkdir(dataRoot);
        if ~ok
            error('Voyager:CreateDataRootFailed', ...
                'Cannot create data root %s: %s', dataRoot, message);
        end
    end
    if ~isfolder(stageDir)
        [ok, message] = mkdir(stageDir);
        if ~ok
            error('Voyager:CreateStageDirFailed', ...
                'Cannot create staging directory %s: %s', stageDir, message);
        end
    end
end

report = VoyagerFilesDownload_NAS(char(opts.Date), dataRoot, ...
    'Spacecraft', spacecraft, ...
    'Products', products, ...
    'StageDir', stageDir, ...
    'ManifestName', manifestName, ...
    'Threads', opts.Threads, ...
    'CheckSize', opts.CheckSize, ...
    'Force', opts.Force, ...
    'ListOnly', opts.ListOnly, ...
    'PythonExe', char(opts.PythonExe));
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

function value = normalizeProducts(products)
if iscell(products)
    value = strjoin(cellfun(@char, products, 'UniformOutput', false), ',');
else
    value = char(products);
end
value = regexprep(strtrim(value), '[;\s]+', ',');
value = regexprep(value, ',+', ',');
value = regexprep(value, '^,|,$', '');
if isempty(value)
    error('Voyager:InvalidProducts', 'Products cannot be empty.');
end
end
