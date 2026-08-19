function product = Voyager_Read_CDF_Product(sourceFile, profile, varargin)
%Voyager_Read_CDF_Product Read selected Voyager CDF variables safely.
%   PRODUCT = Voyager_Read_CDF_Product(FILE, PROFILE) uses the bundled
%   Python bridge to convert CDF values to a little-endian binary cache,
%   then returns the values as a MATLAB structure. Time variables are UTC
%   datetimes. No interpolation, resampling, or gap filling is performed.

parser = inputParser;
parser.CaseSensitive = false;
programRoot = fileparts(mfilename('fullpath'));
addParameter(parser, 'CacheRoot', fullfile(tempdir, ...
    'Voyager_monthly_cdf_cache'), @isTextScalar);
addParameter(parser, 'PythonExe', '', @isTextScalar);
parse(parser, varargin{:});
opts = parser.Results;

sourceFile = char(sourceFile);
profile = char(profile);
cacheRoot = char(opts.CacheRoot);
if ~isfile(sourceFile)
    error('VoyagerCDF:SourceMissing', 'CDF file is unavailable: %s', sourceFile);
end
if ~isfolder(cacheRoot)
    [ok, message] = mkdir(cacheRoot);
    if ~ok, error('VoyagerCDF:CacheCreateFailed', '%s', message); end
end

% Reuse a previously validated binary conversion directly. This avoids
% launching Python once per monthly CDF during short batch plotting jobs.
metadataFile = findExistingCache(cacheRoot, sourceFile, profile);
if ~isempty(metadataFile)
    product = readBinaryCache(metadataFile);
    return
end

bridgeFile = fullfile(programRoot, 'voyager_cdf_bridge.py');
if ~isfile(bridgeFile)
    error('VoyagerCDF:BridgeMissing', 'Missing CDF bridge: %s', bridgeFile);
end
pythonCommand = resolvePython(opts.PythonExe, bridgeFile);
command = sprintf('%s %s --source %s --cache-root %s --profile %s', ...
    pythonCommand, quoteArgument(bridgeFile), quoteArgument(sourceFile), ...
    quoteArgument(cacheRoot), quoteArgument(profile));
[status, output] = system(command);
if status ~= 0
    error('VoyagerCDF:BridgeFailed', ...
        'CDF bridge failed for %s.\n%s', sourceFile, output);
end
match = regexp(output, 'CACHE_JSON\|([^\r\n]+)', 'tokens', 'once');
if isempty(match)
    error('VoyagerCDF:BridgeOutput', ...
        'CDF bridge did not return a cache path. Output:\n%s', output);
end
product = readBinaryCache(strtrim(match{1}));
end

function metadataFile = findExistingCache(cacheRoot, sourceFile, profile)
persistent indexedRoot cacheIndex
convertedRoot = fullfile(cacheRoot, 'converted');
metadataFile = '';
if ~isfolder(convertedRoot)
    return
end
if isempty(indexedRoot) || ~strcmp(indexedRoot, convertedRoot)
    cacheIndex = containers.Map('KeyType', 'char', 'ValueType', 'char');
    entries = dir(fullfile(convertedRoot, '*.json'));
    [~, order] = sort([entries.datenum], 'descend');
    entries = entries(order);
    for ii = 1:numel(entries)
        candidate = fullfile(entries(ii).folder, entries(ii).name);
        try
            item = jsondecode(fileread(candidate));
            if ~isfield(item, 'source_file') || ~isfield(item, 'profile') || ...
                    ~isfield(item, 'source_size') || ...
                    ~isfield(item, 'binary_file') || ...
                    ~isfield(item, 'binary_size')
                continue
            end
            binaryFile = fullfile(entries(ii).folder, item.binary_file);
            binaryInfo = dir(binaryFile);
            if isempty(binaryInfo) || ...
                    binaryInfo(1).bytes ~= double(item.binary_size)
                continue
            end
            key = [char(item.profile), '|', char(item.source_file)];
            if ~isKey(cacheIndex, key)
                cacheIndex(key) = candidate;
            end
        catch
            % Ignore incomplete or unrelated cache metadata.
        end
    end
    indexedRoot = convertedRoot;
end
key = [profile, '|', sourceFile];
if isKey(cacheIndex, key)
    candidate = cacheIndex(key);
    try
        item = jsondecode(fileread(candidate));
        sourceInfo = dir(sourceFile);
        binaryInfo = dir(fullfile(fileparts(candidate), item.binary_file));
        if ~isempty(sourceInfo) && ~isempty(binaryInfo) && ...
                sourceInfo(1).bytes == double(item.source_size) && ...
                binaryInfo(1).bytes == double(item.binary_size)
            metadataFile = candidate;
        end
    catch
        metadataFile = '';
    end
end
end

function product = readBinaryCache(metadataFile)
metadata = jsondecode(fileread(metadataFile));
binaryFile = fullfile(fileparts(metadataFile), metadata.binary_file);
fid = fopen(binaryFile, 'rb', 'ieee-le');
if fid < 0
    error('VoyagerCDF:BinaryOpen', 'Cannot open cache: %s', binaryFile);
end
cleanup = onCleanup(@() fclose(fid));
product = struct('available', true, 'source_file', metadata.source_file, ...
    'profile', metadata.profile, 'variable_meta', struct);
for ii = 1:numel(metadata.variables)
    item = metadata.variables(ii);
    field = matlab.lang.makeValidName(item.name);
    fseek(fid, double(item.offset_bytes), 'bof');
    value = fread(fid, double(item.count), 'double=>double');
    shape = double(item.shape(:).');
    if isempty(shape)
        shape = [1 1];
    elseif isscalar(shape)
        shape = [shape 1]; %#ok<AGROW>
    end
    value = reshape(value, shape);
    if logical(item.is_time)
        value = datetime(value, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    end
    product.(field) = value;
    product.variable_meta.(field) = item;
end
clear cleanup
end

function command = resolvePython(explicit, bridgeFile)
candidates = {};
explicit = strtrim(char(explicit));
if ~isempty(explicit), candidates{end + 1} = quoteArgument(explicit); end %#ok<AGROW>
environment = strtrim(getenv('VOYAGER_PYTHON'));
if ~isempty(environment), candidates{end + 1} = quoteArgument(environment); end %#ok<AGROW>
if ispc
    bundled = fullfile(getenv('USERPROFILE'), '.cache', 'codex-runtimes', ...
        'codex-primary-runtime', 'dependencies', 'python', 'python.exe');
    if isfile(bundled), candidates{end + 1} = quoteArgument(bundled); end %#ok<AGROW>
    candidates = [candidates, {'py -3', 'python', 'python3'}]; %#ok<AGROW>
else
    candidates = [candidates, {'python3', 'python'}]; %#ok<AGROW>
end
for ii = 1:numel(candidates)
    test = sprintf('%s %s --help', candidates{ii}, quoteArgument(bridgeFile));
    [status, ~] = system(test);
    if status == 0
        command = candidates{ii};
        return
    end
end
error('VoyagerCDF:PythonUnavailable', ...
    'A Python interpreter with numpy and cdflib is required.');
end

function output = quoteArgument(input)
dq = char(34);
output = [dq, strrep(char(input), dq, [dq, dq]), dq];
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end
