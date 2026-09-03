function product = Voyager_Read_CDF_Product(sourceFile, profile)
%Voyager_Read_CDF_Product Read a Voyager CDF directly with IRFU-MATLAB.
%   PRODUCT = Voyager_Read_CDF_Product(FILE, PROFILE) reads the original
%   CDF through dataobj. Exact FILLVAL values and values outside VALIDMIN /
%   VALIDMAX are changed to NaN. No interpolation, resampling, smoothing,
%   averaging, extrapolation, or replacement by zero is performed.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-02

%% input check
sourceFile = char(sourceFile);
profile = lower(char(profile));
if ~isfile(sourceFile)
    error('VoyagerCDF:SourceMissing', ...
        'CDF file is unavailable: %s', sourceFile);
end
if exist('dataobj', 'file') ~= 2
    error('VoyagerCDF:IRFUUnavailable', ...
        'IRFU-MATLAB dataobj must be on the MATLAB path.');
end

%% select variables
cdfObj = dataobj(sourceFile);
available = string(cdfObj.Variables(:,1));
requested = profileVariables(profile, available);

product = struct('available', true, 'source_file', sourceFile, ...
    'profile', profile, 'reader', 'IRFU dataobj direct CDF read', ...
    'variable_meta', struct);
product.global_attributes = cdfObj.GlobalAttributes;

%% read and quality-screen each variable
for ii = 1:numel(requested)
    sourceName = char(requested(ii));
    raw = getv(cdfObj, sourceName);
    if isempty(raw), continue, end

    field = matlab.lang.makeValidName(sourceName);
    if strcmpi(sourceName, 'Epoch') || isEpochType(raw)
        value = readEpoch(raw);
        meta = buildMeta(raw, sourceName, 0, 0, 0);
        meta.is_time = true;
    else
        value = raw.data;
        if isnumeric(value) || islogical(value)
            [value, nFill, nBelow, nAbove] = screenNumeric(value, raw);
        else
            nFill = 0; nBelow = 0; nAbove = 0;
        end
        meta = buildMeta(raw, sourceName, nFill, nBelow, nAbove);
        meta.is_time = false;
    end

    product.(field) = value;
    product.variable_meta.(field) = meta;
end
end

function requested = profileVariables(profile, available)
switch profile
    case 'mag48s'
        requested = ["Epoch", "Epoch_ephem", "spacecraftID", ...
            "F1", "BR", "BT", "BN", "dF", "dBR", "dBT", ...
            "dBN", "Radius", "hg_lat", "hg_lon", "hgi_lon"];
    case 'coho'
        requested = ["Epoch", "heliocentricDistance", ...
            "heliographicLatitude", "heliographicLongitude", ...
            "ABS_B", "F", "BR", "BT", "BN", "V", ...
            "elevAngle", "azimuthAngle", "protonDensity", ...
            "protonTemp"];
        patterns = ["^protonFlux\d+_LECP$", "^protonFlux\d+_CRS$"];
        for ii = 1:numel(patterns)
            requested = [requested, ... %#ok<AGROW>
                available(~cellfun('isempty', ...
                regexp(cellstr(available), patterns(ii), 'once')))'];
        end
    case 'pls'
        requested = ["Epoch", "chi2", "V_rtn", "xV_rtn", ...
            "V", "xV", "dens", "xdens", "w", "xw", "T", "xT"];
    case 'mag2s'
        requested = ["Epoch", "Epoch2", "F1", "F2", "B1", ...
            "B2", "B3", "RMSB1", "RMSB2", "RMSB3", ...
            "quality_flag", "dataConfidence", "numDetailPoints", ...
            "scDistance"];
    case 'lecp_sector_daily'
        requested = ["Epoch", "DeltaT", "FHDU_SectoredFluxes", ...
            "FHDU_SectoredFluxUncertainties", "FHDU_SectoredRates", ...
            "FHDU_SectoredRateUncertainties", "FHDU_Energy", ...
            "FHDU_EnergyRange", "SectorIterator", "Hydrogen_Channels", ...
            "Hydrogen_Channels_Label"];
    otherwise
        error('VoyagerCDF:Profile', 'Unknown CDF profile: %s', profile);
end
requested = unique(requested(ismember(requested, available)), 'stable');
end

function tf = isEpochType(raw)
tf = false;
if isfield(raw, 'type')
    tf = contains(lower(string(raw.type)), ["epoch", "tt2000"]);
end
end

function value = readEpoch(raw)
if isdatetime(raw.data)
    value = raw.data;
    value.TimeZone = 'UTC';
    return
end

type = "";
if isfield(raw, 'type'), type = lower(string(raw.type)); end
if contains(type, "tt2000")
    epochUnix = EpochTT(int64(raw.data(:))).epochUnix;
elseif contains(type, "epoch16")
    error('VoyagerCDF:Epoch16', ...
        'CDF_EPOCH16 is not expected for the selected Voyager products.');
else
    % IRFU dataobj stores CDF_EPOCH records internally as POSIX seconds.
    epochUnix = double(raw.data(:));
end
value = datetime(epochUnix, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
end

function [value, nFill, nBelow, nAbove] = screenNumeric(value, raw)
value = double(value);
nFill = 0; nBelow = 0; nAbove = 0;

if isfield(raw, 'FILLVAL') && isnumeric(raw.FILLVAL) && isscalar(raw.FILLVAL)
    mask = value == double(raw.FILLVAL);
    nFill = nnz(mask);
    value(mask) = NaN;
end
if isfield(raw, 'VALIDMIN') && isnumeric(raw.VALIDMIN) && ...
        isscalar(raw.VALIDMIN)
    mask = isfinite(value) & value < double(raw.VALIDMIN);
    nBelow = nnz(mask);
    value(mask) = NaN;
end
if isfield(raw, 'VALIDMAX') && isnumeric(raw.VALIDMAX) && ...
        isscalar(raw.VALIDMAX)
    mask = isfinite(value) & value > double(raw.VALIDMAX);
    nAbove = nnz(mask);
    value(mask) = NaN;
end

% Very large CDF sentinel values are kept out even when metadata is absent.
mask = isfinite(value) & abs(value) >= 1e30;
nFill = nFill + nnz(mask);
value(mask) = NaN;
value(~isfinite(value)) = NaN;
end

function meta = buildMeta(raw, sourceName, nFill, nBelow, nAbove)
meta = struct;
meta.name = sourceName;
meta.shape = size(raw.data);
meta.count = numel(raw.data);
meta.attributes = struct;
attributeNames = {'FIELDNAM', 'LABLAXIS', 'CATDESC', 'UNITS', ...
    'VALIDMIN', 'VALIDMAX', 'FILLVAL', 'DEPEND_0', 'DEPEND_1', ...
    'DEPEND_2', 'DEPEND_3'};
for ii = 1:numel(attributeNames)
    name = attributeNames{ii};
    if isfield(raw, name)
        meta.attributes.(name) = raw.(name);
    end
end
meta.fill_value_rejected = nFill;
meta.below_valid_min_rejected = nBelow;
meta.above_valid_max_rejected = nAbove;
end
