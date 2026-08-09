function data = Voyager_Read_CDF(filename)
% Voyager_Read_CDF Read a downloaded Voyager science CDF.
%
% Supported products:
%   * COHO merged magnetic field, plasma, and position at 1-hour cadence
%   * Legacy calibrated MAG data at 1.92-second cadence
%   * Post-1991 primary/secondary unreviewed MAG data at 1.92-second cadence
%   * Reviewed Voyager Interstellar Mission MAG data at 48-second cadence
%   * HelioWeb position data at either 1-hour or 1-day cadence
%
% CDF fill values and values outside VALIDMIN/VALIDMAX are converted to
% NaN. Products with separate science and ephemeris/header time axes retain
% both axes explicitly.
%
% Example
%   data = Voyager_Read_CDF('voyager2_2s_mag_19770908_v01.cdf');
%   plot(data.time, data.B_rtn_nT)

filename = char(filename);
if ~isfile(filename)
    error('Voyager:CDFNotFound', 'CDF file does not exist: %s', filename);
end

info = cdfinfo(filename);
variableNames = info.Variables(:, 1);
data = struct;
data.source_file = filename;

if all(ismember({'heliocentricDistance', 'BR', 'BT', 'BN', ...
        'V', 'protonDensity', 'protonTemp'}, variableNames))
    data = readCoho(data, filename, info);

elseif all(ismember({'Epoch', 'Epoch2', 'F1', 'F2', ...
        'B1', 'B2', 'B3'}, variableNames))
    data = readLegacyMag2s(data, filename, info);

elseif all(ismember({'Epoch', 'Epoch_ephem', 'F1', 'BR', 'BT', 'BN', ...
        'Radius', 'hg_lat', 'hg_lon', 'hgi_lon'}, variableNames))
    data = readVimMag48s(data, filename, info);

elseif all(ismember({'Epoch', 'RAD_AU', 'HGI_LAT', 'HGI_LON'}, ...
        variableNames))
    data = readHelioPosition(data, filename, info);

else
    error('Voyager:UnsupportedCDF', ...
        ['The file is not a supported Voyager COHO, MAG 1.92 s, ', ...
        'VIM MAG 48 s, or HelioWeb position CDF: %s'], filename);
end

data.variable_units = unitsStruct(info, variableNames);
end

function data = readCoho(data, filename, info)
data.product = 'coho1hr_merged_mag_plasma';
data.time = readVariable(filename, info, 'Epoch', true);
data.position_time = data.time;
data.magnetic_cadence_seconds = 3600;
data.position_cadence_seconds = 3600;
data.position_r_au = readVariable(filename, info, ...
    'heliocentricDistance', false);
data.position_hgi_lat_deg = readVariable(filename, info, ...
    'heliographicLatitude', false);
data.position_hgi_lon_deg = readVariable(filename, info, ...
    'heliographicLongitude', false);

br = readVariable(filename, info, 'BR', false);
bt = readVariable(filename, info, 'BT', false);
bn = readVariable(filename, info, 'BN', false);
data.B_rtn_nT = [br, bt, bn];
data.B_abs_nT = readVariable(filename, info, 'ABS_B', false);
data.B_vector_magnitude_nT = readVariable(filename, info, 'F', false);

data.flow_speed_kms = readVariable(filename, info, 'V', false);
data.flow_elevation_deg = readVariable(filename, info, ...
    'elevAngle', false);
data.flow_azimuth_deg = readVariable(filename, info, ...
    'azimuthAngle', false);
theta = data.flow_elevation_deg .* (pi / 180);
phi = data.flow_azimuth_deg .* (pi / 180);
vr = data.flow_speed_kms .* cos(theta) .* cos(phi);
vt = data.flow_speed_kms .* cos(theta) .* sin(phi);
vn = data.flow_speed_kms .* sin(theta);
data.flow_velocity_rtn_kms = [vr, vt, vn];

data.proton_density_cm3 = readVariable(filename, info, ...
    'protonDensity', false);
data.proton_temperature_K = readVariable(filename, info, ...
    'protonTemp', false);
data.proton_temperature_eV = data.proton_temperature_K ./ 11604.51812;

[x, y, z] = sphericalPosition(data.position_r_au, ...
    data.position_hgi_lat_deg, data.position_hgi_lon_deg);
data.position_hgi_au = [x, y, z];
end

function data = readLegacyMag2s(data, filename, info)
lowerFilename = lower(filename);
if contains(lowerFilename, '_2s_mag-pri_')
    data.product = 'mag2s_unreviewed_post1991';
    data.sensor = 'primary';
    data.sensor_location = 'out-board';
    data.science_quality = 'generally_not_science_quality';
elseif contains(lowerFilename, '_2s_mag-sec_')
    data.product = 'mag2s_unreviewed_post1991';
    data.sensor = 'secondary';
    data.sensor_location = 'in-board';
    data.science_quality = 'generally_not_science_quality';
else
    data.product = 'mag2s_calibrated_legacy';
    data.sensor = '';
    data.sensor_location = '';
    data.science_quality = 'calibrated_unreviewed';
end
data.magnetic_cadence_seconds = 1.92;

% Epoch2 is the 1.92-second science time axis. Epoch indexes each 48-second
% source record and is also the dependency for status and ephemeris fields.
data.time = readVariable(filename, info, 'Epoch2', true);
data.time_block = readVariable(filename, info, 'Epoch', true);
data.position_time = data.time_block;

b1 = readVariable(filename, info, 'B1', false);
b2 = readVariable(filename, info, 'B2', false);
b3 = readVariable(filename, info, 'B3', false);
data.B_rtn_nT = [b1, b2, b3];
data.B_abs_nT = readVariable(filename, info, 'F1', false);
data.B_vector_magnitude_nT = readVariable(filename, info, 'F2', false);

rms1 = readOptionalVariable(filename, info, 'RMSB1', false);
rms2 = readOptionalVariable(filename, info, 'RMSB2', false);
rms3 = readOptionalVariable(filename, info, 'RMSB3', false);
if ~isempty(rms1) && ~isempty(rms2) && ~isempty(rms3)
    data.B_rms_rtn_nT = [rms1, rms2, rms3];
end

data.position_r_au = readVariable(filename, info, 'scDistance', false);
data.position_hgi_lon_deg = readVariable(filename, info, ...
    'scLon', false) .* (180 / pi);
data.position_hgi_lat_deg = readVariable(filename, info, ...
    'scLat', false) .* (180 / pi);
x = readVariable(filename, info, 'x', false);
y = readVariable(filename, info, 'y', false);
z = readVariable(filename, info, 'z', false);
data.position_hgi_au = [x, y, z];

blockCount = numel(data.time_block);
quality = struct;
quality.time = data.time_block;
quality.data_confidence_raw = recordsFirst( ...
    readOptionalVariable(filename, info, 'dataConfidence', false), ...
    blockCount);
if ~isempty(quality.data_confidence_raw)
    quality.all_confidence_flags_zero = ...
        all(isfinite(quality.data_confidence_raw), 2) & ...
        all(quality.data_confidence_raw == 0, 2);
else
    quality.all_confidence_flags_zero = false(0, 1);
end
quality.mag_status_raw = recordsFirst( ...
    readOptionalVariable(filename, info, 'magStatus', false), blockCount);
quality.mag_status_bits_raw = recordsFirst( ...
    readOptionalVariable(filename, info, 'magStatusBYTE', false), ...
    blockCount);
quality.num_primary_data = recordsFirst( ...
    readOptionalVariable(filename, info, 'numPrimaryData', false), ...
    blockCount);
quality.num_secondary_data = recordsFirst( ...
    readOptionalVariable(filename, info, 'numSecondaryData', false), ...
    blockCount);
quality.num_detail_points = recordsFirst( ...
    readOptionalVariable(filename, info, 'numDetailPoints', false), ...
    numel(data.time));
quality.data_confidence_convention = ...
    'CDF dataConfidence flags: 0=OK; retain raw flags when filtering.';
data.quality = quality;
end

function data = readVimMag48s(data, filename, info)
data.product = 'mag48s_vim_reviewed';
data.magnetic_cadence_seconds = 48;
data.position_cadence_seconds = 86400;

data.time = readVariable(filename, info, 'Epoch', true);
data.spacecraft_id = readOptionalVariable(filename, info, ...
    'spacecraftID', false);

br = readVariable(filename, info, 'BR', false);
bt = readVariable(filename, info, 'BT', false);
bn = readVariable(filename, info, 'BN', false);
data.B_rtn_nT = [br, bt, bn];
data.B_abs_nT = readVariable(filename, info, 'F1', false);

data.B_abs_uncertainty_nT = readOptionalVariable(filename, info, ...
    'dF', false);
dbr = readOptionalVariable(filename, info, 'dBR', false);
dbt = readOptionalVariable(filename, info, 'dBT', false);
dbn = readOptionalVariable(filename, info, 'dBN', false);
if ~isempty(dbr) && ~isempty(dbt) && ~isempty(dbn)
    data.B_rtn_uncertainty_nT = [dbr, dbt, dbn];
end

data.time_ephemeris = readVariable(filename, info, 'Epoch_ephem', true);
data.position_time = data.time_ephemeris;
data.position_r_au = readVariable(filename, info, 'Radius', false);
data.position_hg_lat_deg = readVariable(filename, info, 'hg_lat', false);
data.position_hg_lon_deg = readVariable(filename, info, 'hg_lon', false);
data.position_hgi_lat_deg = data.position_hg_lat_deg;
data.position_hgi_lon_deg = readVariable(filename, info, ...
    'hgi_lon', false);
[x, y, z] = sphericalPosition(data.position_r_au, ...
    data.position_hgi_lat_deg, data.position_hgi_lon_deg);
data.position_hgi_au = [x, y, z];
end

function data = readHelioPosition(data, filename, info)
data.time = readVariable(filename, info, 'Epoch', true);
[data.product, cadenceSeconds] = positionProduct(info, filename, data.time);
data.position_cadence_seconds = cadenceSeconds;
data.position_time = data.time;
data.position_r_au = readVariable(filename, info, 'RAD_AU', false);
data.position_hgi_lat_deg = readVariable(filename, info, ...
    'HGI_LAT', false);
data.position_hgi_lon_deg = readVariable(filename, info, ...
    'HGI_LON', false);
data.position_hg_lat_deg = readVariable(filename, info, 'HG_LAT', false);
data.position_hg_lon_deg = readVariable(filename, info, 'HG_LON', false);
data.position_se_lat_deg = readVariable(filename, info, 'SE_LAT', false);
data.position_se_lon_deg = readVariable(filename, info, 'SE_LON', false);
[x, y, z] = sphericalPosition(data.position_r_au, ...
    data.position_hgi_lat_deg, data.position_hgi_lon_deg);
data.position_hgi_au = [x, y, z];
end

function value = readVariable(filename, info, variableName, isEpoch)
variableIndex = find(strcmp(info.Variables(:, 1), variableName), 1);
if isempty(variableIndex)
    error('Voyager:VariableMissing', ...
        'Variable %s is absent from %s.', variableName, filename);
end

recordCount = info.Variables{variableIndex, 3};
if isempty(recordCount) || all(recordCount == 0)
    if isEpoch
        value = NaT(0, 1, 'TimeZone', 'UTC');
    else
        value = zeros(0, 1);
    end
    return
end

if isEpoch
    try
        value = cdfread(filename, 'Variables', {variableName}, ...
            'CombineRecords', true, 'DatetimeType', 'datetime');
    catch
        value = cdfread(filename, 'Variables', {variableName}, ...
            'CombineRecords', true);
    end
else
    value = cdfread(filename, 'Variables', {variableName}, ...
        'CombineRecords', true);
end
value = normalizeRaw(value);

if isEpoch
    if isnumeric(value)
        value = datetime(value, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    elseif isdatetime(value)
        value.TimeZone = 'UTC';
    end
    value = value(:);
    return
end

value = double(value);
if isvector(value)
    value = value(:);
end

fillValue = getVariableAttribute(info, 'FILLVAL', variableName);
validMin = getVariableAttribute(info, 'VALIDMIN', variableName);
validMax = getVariableAttribute(info, 'VALIDMAX', variableName);

if ~isempty(fillValue) && isnumeric(fillValue)
    for ii = 1:numel(fillValue)
        value(value == double(fillValue(ii))) = NaN;
    end
end
if ~isempty(validMin) && isnumeric(validMin)
    value(value < double(validMin(1))) = NaN;
end
if ~isempty(validMax) && isnumeric(validMax)
    value(value > double(validMax(1))) = NaN;
end
value(~isfinite(value) | abs(value) >= 1e30) = NaN;
end

function value = readOptionalVariable(filename, info, variableName, isEpoch)
if any(strcmp(info.Variables(:, 1), variableName))
    value = readVariable(filename, info, variableName, isEpoch);
else
    value = [];
end
end

function value = recordsFirst(value, recordCount)
% CDF vector-valued records can be returned as components-by-records.
if isempty(value) || recordCount == 0
    return
end
if size(value, 1) ~= recordCount && size(value, 2) == recordCount
    value = value.';
end
end

function [product, cadenceSeconds] = positionProduct(info, filename, time)
metadata = lower(strjoin({filename, ...
    globalAttributeText(info, 'Logical_source'), ...
    globalAttributeText(info, 'Logical_file_id'), ...
    globalAttributeText(info, 'Data_type'), ...
    globalAttributeText(info, 'Time_resolution')}, ' '));

if contains(metadata, 'helio1hr') || contains(metadata, '1 hour')
    product = 'helio1hr_position';
    cadenceSeconds = 3600;
elseif contains(metadata, 'helio1day') || contains(metadata, '1 day')
    product = 'helio1day_position';
    cadenceSeconds = 86400;
else
    cadenceSeconds = inferCadenceSeconds(time);
    if isfinite(cadenceSeconds) && cadenceSeconds < 12 * 3600
        product = 'helio1hr_position';
    else
        product = 'helio1day_position';
    end
end
end

function cadenceSeconds = inferCadenceSeconds(time)
cadenceSeconds = NaN;
if numel(time) < 2
    return
end
sampleEnd = min(numel(time), 1001);
intervals = seconds(diff(time(1:sampleEnd)));
intervals = intervals(isfinite(intervals) & intervals > 0);
if ~isempty(intervals)
    cadenceSeconds = median(intervals);
end
end

function text = globalAttributeText(info, attributeName)
text = '';
if ~isfield(info, 'GlobalAttributes')
    return
end
fields = fieldnames(info.GlobalAttributes);
fieldIndex = find(strcmpi(fields, attributeName), 1);
if isempty(fieldIndex)
    return
end
raw = info.GlobalAttributes.(fields{fieldIndex});
if isempty(raw)
    return
end
try
    parts = string(raw(:));
    parts = parts(strlength(parts) > 0);
    text = char(strjoin(parts, ' '));
catch
    if ischar(raw)
        text = raw;
    end
end
end

function value = normalizeRaw(value)
if ~iscell(value)
    return
end
if isempty(value)
    value = [];
    return
end
if isscalar(value)
    value = value{1};
    return
end
try
    value = vertcat(value{:});
catch
    try
        value = cell2mat(value);
    catch
        error('Voyager:CDFShapeUnsupported', ...
            'Cannot combine CDF records into a numeric array.');
    end
end
end

function value = getVariableAttribute(info, attributeName, variableName)
value = [];
fields = fieldnames(info.VariableAttributes);
fieldIndex = find(strcmpi(fields, attributeName), 1);
if isempty(fieldIndex)
    return
end
entries = info.VariableAttributes.(fields{fieldIndex});
for ii = 1:size(entries, 1)
    if strcmp(entries{ii, 1}, variableName)
        value = entries{ii, 2};
        return
    end
end
end

function output = unitsStruct(info, variableNames)
output = struct;
for ii = 1:numel(variableNames)
    name = variableNames{ii};
    unit = getVariableAttribute(info, 'UNITS', name);
    if isempty(unit)
        unit = '';
    end
    if isstring(unit)
        unit = char(unit);
    elseif isnumeric(unit)
        unit = num2str(unit);
    end
    output.(matlab.lang.makeValidName(name)) = unit;
end
end

function [x, y, z] = sphericalPosition(radius, latitudeDeg, longitudeDeg)
latitude = latitudeDeg .* (pi / 180);
longitude = longitudeDeg .* (pi / 180);
x = radius .* cos(latitude) .* cos(longitude);
y = radius .* cos(latitude) .* sin(longitude);
z = radius .* sin(latitude);
end
