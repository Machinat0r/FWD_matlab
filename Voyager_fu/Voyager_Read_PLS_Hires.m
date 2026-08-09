function data = Voyager_Read_PLS_Hires(filename)
% Voyager_Read_PLS_Hires Read Voyager PLS high-resolution CDF data.
%
% DATA = Voyager_Read_PLS_Hires(FILENAME) reads either of these NASA/SPDF
% CDAWeb products:
%
%   VOYAGER1_PLS_HIRES_PLASMA_DATA
%   VOYAGER2_PLS_HIRES_PLASMA_DATA
%   VOYAGER2_PLS_HIRES_PLASMA_DATA_HSH
%
% The first two products contain solar-wind proton moments. The HSH
% product contains Voyager 2 heliosheath fits and their uncertainties.
% Sampling intervals vary from 12 to 192 seconds.
%
% The returned structure has a common layout for all three products:
%
%   source_file       Input CDF path
%   spacecraft        'voyager1', 'voyager2', or 'unknown'
%   product           'pls_hires_solar_wind' or
%                     'pls_hires_heliosheath'
%   Epoch             Nx1 UTC datetime
%   V                 Nx1 proton bulk speed, km/s
%   dens              Nx1 proton number density, cm^-3
%   V_thermal         Nx1 proton thermal speed, km/s
%   V_rtn             Nx3 proton velocity [Vr Vt Vn], km/s
%   T_eV              Nx1 proton temperature, eV
%   T_K               Nx1 proton temperature, K
%   chi2              Nx1 fit chi-square, when available
%   xV                Nx1 1-sigma bulk-speed uncertainty, km/s
%   xdens             Nx1 1-sigma density uncertainty, cm^-3
%   xV_thermal        Nx1 1-sigma thermal-speed uncertainty, km/s
%   xV_rtn            Nx3 1-sigma RTN velocity uncertainty, km/s
%   xT_eV             Nx1 1-sigma temperature uncertainty, eV
%   xT_K              Nx1 1-sigma temperature uncertainty, K
%   timetable         The same common variables in one timetable
%   source_variables  Original CDF variable used for each common field
%   variable_units    Units copied from the original CDF metadata
%   missing_variables Common fields that could not be read or derived
%
% Missing CDF variables are represented by NaN arrays with the appropriate
% number of rows. CDF FILLVAL values, values outside VALIDMIN/VALIDMAX, and
% values whose absolute magnitude is at least 1e30 are converted to NaN.
%
% The regular solar-wind files store thermal speed rather than temperature.
% For those files the official Voyager conversion is used:
%
%   T_eV = 0.0052 * V_thermal.^2
%   T_K  = 60.5   * V_thermal.^2
%
% Example
%   file = 'voyager2_pls_hires_plasma_data_20060101_v01.cdf';
%   d = Voyager_Read_PLS_Hires(file);
%   plot(d.Epoch, d.dens)
%   grid on
%   ylabel('Proton density (cm^{-3})')

filename = char(filename);
if ~isfile(filename)
    error('Voyager:PLSHiresFileNotFound', ...
        'CDF file does not exist: %s', filename);
end

info = cdfinfo(filename);
variableNames = info.Variables(:, 1);

epochName = findVariableName(variableNames, {'Epoch', 'epoch', 'Time'});
if isempty(epochName)
    error('Voyager:PLSHiresEpochMissing', ...
        'No supported epoch variable was found in %s.', filename);
end
epoch = readCdfVariable(filename, info, epochName, true);
numberOfRecords = numel(epoch);
if numberOfRecords == 0
    error('Voyager:PLSHiresEmptyFile', ...
        'The epoch variable contains no records: %s', filename);
end

[bulkSpeed, bulkSpeedName] = readOptionalScalar(filename, info, ...
    variableNames, {'V', 'speed', 'V_bulk'}, numberOfRecords);
[density, densityName] = readOptionalScalar(filename, info, ...
    variableNames, {'dens', 'density', 'protonDensity'}, numberOfRecords);
[thermalSpeed, thermalSpeedName] = readOptionalScalar(filename, info, ...
    variableNames, {'V_thermal', 'w', 'thermalSpeed'}, numberOfRecords);
[velocityRtn, velocityRtnName] = readOptionalVector(filename, info, ...
    variableNames, {'V_rtn', 'V_RTN'}, numberOfRecords, 3);

if isempty(velocityRtnName)
    [vr, vrName] = readOptionalScalar(filename, info, variableNames, ...
        {'Vr', 'VR', 'V_R'}, numberOfRecords);
    [vt, vtName] = readOptionalScalar(filename, info, variableNames, ...
        {'Vt', 'VT', 'V_T'}, numberOfRecords);
    [vn, vnName] = readOptionalScalar(filename, info, variableNames, ...
        {'Vn', 'VN', 'V_N'}, numberOfRecords);
    if ~isempty(vrName) && ~isempty(vtName) && ~isempty(vnName)
        velocityRtn = [vr, vt, vn];
        velocityRtnName = sprintf('%s,%s,%s', vrName, vtName, vnName);
    end
end

if isempty(bulkSpeedName) && any(isfinite(velocityRtn(:)))
    bulkSpeed = sqrt(sum(velocityRtn .^ 2, 2));
    bulkSpeedName = 'derived_from_V_rtn';
end

[chi2, chi2Name] = readOptionalScalar(filename, info, variableNames, ...
    {'chi2', 'chisq', 'chi_square'}, numberOfRecords);
[xDensity, xDensityName] = readOptionalScalar(filename, info, ...
    variableNames, {'xdens', 'x_density', 'dens_sigma'}, numberOfRecords);
[xThermalSpeed, xThermalSpeedName] = readOptionalScalar(filename, info, ...
    variableNames, {'xw', 'xV_thermal', 'thermalSpeed_sigma'}, ...
    numberOfRecords);
[xVelocityRtn, xVelocityRtnName] = readOptionalVector(filename, info, ...
    variableNames, {'xV_rtn', 'xV_RTN', 'V_rtn_sigma'}, ...
    numberOfRecords, 3);
[xBulkSpeed, xBulkSpeedName] = readOptionalScalar(filename, info, ...
    variableNames, {'xV', 'x_speed', 'V_sigma'}, numberOfRecords);

if isempty(xBulkSpeedName) && ~isempty(xVelocityRtnName)
    xBulkSpeed = vectorMagnitudeUncertainty(velocityRtn, xVelocityRtn);
    xBulkSpeedName = 'derived_from_xV_rtn';
end

[temperatureRaw, temperatureName] = readOptionalScalar(filename, info, ...
    variableNames, {'T', 'T_eV', 'T_K', 'protonTemp'}, numberOfRecords);
[xTemperatureRaw, xTemperatureName] = readOptionalScalar(filename, info, ...
    variableNames, {'xT', 'xT_eV', 'xT_K', 'T_sigma'}, numberOfRecords);

if isempty(temperatureName)
    temperatureEv = 0.0052 .* thermalSpeed .^ 2;
    temperatureK = 60.5 .* thermalSpeed .^ 2;
    if any(isfinite(temperatureEv))
        temperatureName = 'derived_from_V_thermal';
    end
else
    temperatureUnit = getVariableUnit(info, temperatureName);
    [temperatureEv, temperatureK] = convertTemperature( ...
        temperatureRaw, temperatureUnit, temperatureName);
end

if isempty(xTemperatureName)
    xTemperatureEv = 2 .* 0.0052 .* abs(thermalSpeed) .* xThermalSpeed;
    xTemperatureK = 2 .* 60.5 .* abs(thermalSpeed) .* xThermalSpeed;
    if any(isfinite(xTemperatureEv))
        xTemperatureName = 'propagated_from_xV_thermal';
    end
else
    xTemperatureUnit = getVariableUnit(info, xTemperatureName);
    if isempty(xTemperatureUnit) && ~isempty(temperatureName) && ...
            ~startsWith(temperatureName, 'derived_')
        xTemperatureUnit = getVariableUnit(info, temperatureName);
    end
    [xTemperatureEv, xTemperatureK] = convertTemperature( ...
        xTemperatureRaw, xTemperatureUnit, xTemperatureName);
end

data = struct;
data.source_file = filename;
data.spacecraft = identifySpacecraft(filename, info);
if ~isempty(chi2Name) || ~isempty(xTemperatureName) || ...
        contains(lower(filename), '_hsh_')
    data.product = 'pls_hires_heliosheath';
else
    data.product = 'pls_hires_solar_wind';
end

data.Epoch = epoch;
data.V = bulkSpeed;
data.dens = density;
data.V_thermal = thermalSpeed;
data.V_rtn = velocityRtn;
data.T_eV = temperatureEv;
data.T_K = temperatureK;
data.chi2 = chi2;
data.xV = xBulkSpeed;
data.xdens = xDensity;
data.xV_thermal = xThermalSpeed;
data.xV_rtn = xVelocityRtn;
data.xT_eV = xTemperatureEv;
data.xT_K = xTemperatureK;

data.source_variables = struct( ...
    'Epoch', epochName, ...
    'V', bulkSpeedName, ...
    'dens', densityName, ...
    'V_thermal', thermalSpeedName, ...
    'V_rtn', velocityRtnName, ...
    'T', temperatureName, ...
    'chi2', chi2Name, ...
    'xV', xBulkSpeedName, ...
    'xdens', xDensityName, ...
    'xV_thermal', xThermalSpeedName, ...
    'xV_rtn', xVelocityRtnName, ...
    'xT', xTemperatureName);
data.variable_units = unitsStruct(info, variableNames);
data.missing_variables = missingCommonVariables(data.source_variables);

data.timetable = timetable(data.Epoch, data.V, data.dens, ...
    data.V_thermal, data.V_rtn, data.T_eV, data.T_K, data.chi2, ...
    data.xV, data.xdens, data.xV_thermal, data.xV_rtn, ...
    data.xT_eV, data.xT_K, ...
    'VariableNames', {'V', 'dens', 'V_thermal', 'V_rtn', ...
    'T_eV', 'T_K', 'chi2', 'xV', 'xdens', 'xV_thermal', ...
    'xV_rtn', 'xT_eV', 'xT_K'});
data.timetable.Properties.DimensionNames{1} = 'Epoch';
end

function [value, actualName] = readOptionalScalar(filename, info, ...
        variableNames, aliases, numberOfRecords)
actualName = findVariableName(variableNames, aliases);
if isempty(actualName)
    value = nan(numberOfRecords, 1);
    return
end
raw = readCdfVariable(filename, info, actualName, false);
value = shapeRecords(raw, numberOfRecords, 1, actualName);
end

function [value, actualName] = readOptionalVector(filename, info, ...
        variableNames, aliases, numberOfRecords, width)
actualName = findVariableName(variableNames, aliases);
if isempty(actualName)
    value = nan(numberOfRecords, width);
    return
end
raw = readCdfVariable(filename, info, actualName, false);
value = shapeRecords(raw, numberOfRecords, width, actualName);
end

function actualName = findVariableName(variableNames, aliases)
actualName = '';
for ii = 1:numel(aliases)
    index = find(strcmpi(variableNames, aliases{ii}), 1);
    if ~isempty(index)
        actualName = variableNames{index};
        return
    end
end
end

function value = readCdfVariable(filename, info, variableName, isEpoch)
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
        value = datetime(value, 'ConvertFrom', 'datenum', ...
            'TimeZone', 'UTC');
    elseif isdatetime(value)
        value.TimeZone = 'UTC';
    end
    value = value(:);
    return
end

value = double(value);
fillValue = getVariableAttribute(info, 'FILLVAL', variableName);
validMin = getVariableAttribute(info, 'VALIDMIN', variableName);
validMax = getVariableAttribute(info, 'VALIDMAX', variableName);

if isnumeric(fillValue)
    for ii = 1:numel(fillValue)
        value(value == double(fillValue(ii))) = NaN;
    end
end
if isnumeric(validMin) && ~isempty(validMin)
    value(value < double(validMin(1))) = NaN;
end
if isnumeric(validMax) && ~isempty(validMax)
    value(value > double(validMax(1))) = NaN;
end
value(~isfinite(value) | abs(value) >= 1e30) = NaN;
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
        error('Voyager:PLSHiresShapeUnsupported', ...
            'Cannot combine CDF records into a numeric array.');
    end
end
end

function value = shapeRecords(value, numberOfRecords, width, variableName)
if isempty(value)
    value = nan(numberOfRecords, width);
    return
end

if width == 1
    if numel(value) ~= numberOfRecords
        error('Voyager:PLSHiresRecordMismatch', ...
            ['Variable %s contains %d values, while Epoch contains ' ...
            '%d records.'], variableName, numel(value), numberOfRecords);
    end
    value = reshape(value, numberOfRecords, 1);
    return
end

if size(value, 1) == numberOfRecords && size(value, 2) == width
    return
end
if size(value, 2) == numberOfRecords && size(value, 1) == width
    value = value.';
    return
end
if numel(value) == numberOfRecords * width
    value = reshape(value, width, numberOfRecords).';
    return
end
error('Voyager:PLSHiresRecordMismatch', ...
    ['Variable %s has size %s; expected %d-by-%d records or its ' ...
    'transpose.'], variableName, mat2str(size(value)), ...
    numberOfRecords, width);
end

function uncertainty = vectorMagnitudeUncertainty(vector, vectorSigma)
speed = sqrt(sum(vector .^ 2, 2));
uncertainty = nan(size(speed));
valid = isfinite(speed) & speed > 0 & ...
    all(isfinite(vector), 2) & all(isfinite(vectorSigma), 2);
weights = vector(valid, :) ./ speed(valid);
uncertainty(valid) = sqrt(sum((weights .* vectorSigma(valid, :)) .^ 2, 2));
end

function [temperatureEv, temperatureK] = convertTemperature( ...
        temperature, unit, variableName)
unitText = lower(strtrim(char(unit)));
nameText = lower(variableName);
isKelvin = strcmp(unitText, 'k') || contains(unitText, 'kelvin') || ...
    endsWith(nameText, '_k');
isEv = contains(unitText, 'ev') || endsWith(nameText, '_ev');

if isKelvin && ~isEv
    temperatureK = temperature;
    temperatureEv = temperature ./ 11604.51812;
else
    % The HSH T and xT variables are in eV. Treat an unlabelled T-like
    % variable as eV to remain compatible with those official files.
    temperatureEv = temperature;
    temperatureK = temperature .* 11604.51812;
end
end

function value = getVariableAttribute(info, attributeName, variableName)
value = [];
if ~isfield(info, 'VariableAttributes') || ...
        isempty(info.VariableAttributes)
    return
end
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

function unit = getVariableUnit(info, variableName)
unit = getVariableAttribute(info, 'UNITS', variableName);
if isempty(unit)
    unit = '';
elseif isstring(unit)
    unit = char(unit);
elseif isnumeric(unit)
    unit = num2str(unit);
end
end

function output = unitsStruct(info, variableNames)
output = struct;
for ii = 1:numel(variableNames)
    name = variableNames{ii};
    output.(matlab.lang.makeValidName(name)) = ...
        getVariableUnit(info, name);
end
end

function spacecraft = identifySpacecraft(filename, info)
lowerName = lower(filename);
if contains(lowerName, 'voyager1') || contains(lowerName, 'voyager_1')
    spacecraft = 'voyager1';
    return
end
if contains(lowerName, 'voyager2') || contains(lowerName, 'voyager_2')
    spacecraft = 'voyager2';
    return
end

spacecraft = 'unknown';
if ~isfield(info, 'GlobalAttributes')
    return
end
fields = fieldnames(info.GlobalAttributes);
for ii = 1:numel(fields)
    values = info.GlobalAttributes.(fields{ii});
    if iscell(values)
        try
            textValue = lower(strjoin(cellfun(@char, values(:), ...
                'UniformOutput', false), ' '));
        catch
            continue
        end
    else
        textValue = lower(char(values));
    end
    if contains(textValue, 'voyager 1') || contains(textValue, 'voyager1')
        spacecraft = 'voyager1';
        return
    end
    if contains(textValue, 'voyager 2') || contains(textValue, 'voyager2')
        spacecraft = 'voyager2';
        return
    end
end
end

function missing = missingCommonVariables(sourceVariables)
fields = fieldnames(sourceVariables);
missing = {};
for ii = 1:numel(fields)
    if isempty(sourceVariables.(fields{ii}))
        missing{end + 1, 1} = fields{ii}; %#ok<AGROW>
    end
end
end
