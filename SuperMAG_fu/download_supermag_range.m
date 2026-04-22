function savedFiles = download_supermag_range(userid, station, tStart, tEnd, outDir, varargin)
%DOWNLOAD_SUPERMAG_RANGE Download SuperMAG magnetometer data using fetchSuperMAG.
%
% savedFiles = download_supermag_range(userid, station, tStart, tEnd, outDir, ...)
%
% Required:
%   userid  : your SuperMAG login name (string/char)
%   station : IAGA code (e.g., 'VIC') or list (cellstr/string array)
%   tStart  : datetime OR char/string in 'yyyy-MM-dd''T''HH:mm' (UTC recommended)
%   tEnd    : datetime OR char/string in 'yyyy-MM-dd''T''HH:mm'
%   outDir  : folder to save .mat files
%
% Name-Value options:
%   'Flags'        : flag string passed to fetchSuperMAG (default 'ALL')
%   'ChunkSeconds' : seconds per request (default 86400 = 1 day)
%   'PauseSeconds' : pause between requests (default 1)
%   'MaxRetries'   : retry count per chunk (default 3)
%   'CheckInventory' : true/false, call inventory first (default false)

% ---- Parse inputs
p = inputParser;
p.addParameter('Flags', 'ALL', @(x)ischar(x) || isstring(x));
p.addParameter('ChunkSeconds', 86400, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('PauseSeconds', 1, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('MaxRetries', 3, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CheckInventory', false, @(x)islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

% ---- Normalize times
tStart = local_to_datetime(tStart);
tEnd   = local_to_datetime(tEnd);
if tEnd <= tStart
    error('tEnd must be later than tStart.');
end

% ---- Normalize station
stationArg = station_to_arg(station);

% ---- Ensure output dir
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ---- Optional inventory check (helps validate station/time coverage)
if opt.CheckInventory
    startStr = dt_to_supermag_str(tStart);
    extentSec = seconds(min(tEnd, tStart + seconds(opt.ChunkSeconds)) - tStart);
    inv = call_fetchSuperMAG('inventory', userid, startStr, extentSec); %#ok<NASGU>
    % 说明：inventory 返回结构/表格的具体形态取决于客户端版本；
    % 这里只是触发一次调用，确保账号/网络/路径都通。
end

% ---- Download loop
chunk = seconds(opt.ChunkSeconds);
t0 = tStart;
savedFiles = {};

while t0 < tEnd
    t1 = min(t0 + chunk, tEnd);
    extentSec = seconds(t1 - t0);
    startStr = dt_to_supermag_str(t0);

    % Filename: include station + start time
    stationTag = regexprep(stationArg, '[^A-Za-z0-9,]+', '_');
    stationTag = strrep(stationTag, ',', '-');
    fname = sprintf('supermag_%s_%s_%ds.mat', stationTag, datestr(t0,'yyyymmddTHHMM'), round(extentSec));
    fpath = fullfile(outDir, fname);

    % Retry wrapper
    lastErr = [];
    for k = 1:opt.MaxRetries
        try
            % MATLAB client doc indicates:
            % sm_data = fetchSuperMAG('data', userid, start, extent, flags, station)
            % Some variants may also return a status as first output.
            [status, sm_data] = call_fetchSuperMAG('data', userid, startStr, extentSec, opt.Flags, stationArg); %#ok<ASGLU>

            meta.userid  = userid;
            meta.station = stationArg;
            meta.flags   = opt.Flags;
            meta.t0      = t0;
            meta.t1      = t1;

            save(fpath, 'sm_data', 'meta', '-v7.3');
            savedFiles{end+1,1} = fpath; %#ok<AGROW>
            lastErr = [];
            break; % success
        catch ME
            lastErr = ME;
            pause( max(0.5, opt.PauseSeconds) ); % small backoff
        end
    end

    if ~isempty(lastErr)
        warning('Failed chunk %s -> %s. Last error: %s', startStr, dt_to_supermag_str(t1), lastErr.message);
    end

    pause(opt.PauseSeconds);
    t0 = t1;
end

end

% ===== Helper functions =====

function dt = local_to_datetime(x)
    if isdatetime(x)
        dt = x;
    elseif ischar(x) || isstring(x)
        dt = datetime(x, 'InputFormat', "yyyy-MM-dd'T'HH:mm", 'TimeZone', 'UTC');
    else
        error('tStart/tEnd must be datetime or a string like yyyy-MM-ddTHH:mm');
    end

    % If timezone not set, assume UTC (recommended for SuperMAG requests)
    if isempty(dt.TimeZone)
        dt.TimeZone = 'UTC';
    end
end

function s = dt_to_supermag_str(dt)
    % SuperMAG examples use ISO-like 'YYYY-MM-DDTHH:MM'
    dt.TimeZone = 'UTC';
    s = char(datetime(dt, 'Format', "yyyy-MM-dd'T'HH:mm"));
end

function stationArg = station_to_arg(station)
    if ischar(station) || isstring(station)
        stationArg = char(station);
    elseif iscell(station)
        stationArg = strjoin(station, ',');
    elseif isstring(station)
        stationArg = strjoin(cellstr(station), ',');
    else
        error('station must be char/string or a cell array of station codes.');
    end
end

function varargout = call_fetchSuperMAG(varargin)
    % Handles both 1-output and 2-output variants.
    %
    % Possible doc forms:
    %   stations = fetchSuperMAG('inventory', userid, start, extent)
    %   sm_data  = fetchSuperMAG('data', userid, start, extent, flags, station)
    %
    % Some clients may return [status, data] as two outputs.
    try
        [varargout{1:nargout}] = fetchSuperMAG(varargin{:});
    catch ME
        if contains(ME.message, 'Too many output arguments')
            % Fall back to single output
            out = fetchSuperMAG(varargin{:});
            varargout = cell(1, max(1,nargout));
            varargout{1} = [];
            if nargout >= 2
                varargout{2} = out;
            else
                varargout{1} = out;
            end
        else
            rethrow(ME);
        end
    end
end
