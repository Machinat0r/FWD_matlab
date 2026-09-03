function product = Case1_Read_LECP_CDFs(sourceFiles)
%Case1_Read_LECP_CDFs Read original annual LECP CDFs directly with IRFU.
%   Concatenation changes no measurement. Conflicting duplicate epochs,
%   channel ordering, or units raise an error. No intermediate flux CSV.

%% original files
sourceFiles = string(sourceFiles);
sourceFiles = unique(sourceFiles(:), 'stable');
if isempty(sourceFiles) || any(~isfile(sourceFiles))
    error('VoyagerCase1:NativeCDFMissing', 'Original LECP CDF files are missing.');
end
recordFields = {'Epoch', 'DeltaT', 'FHDU_SectoredFluxes', ...
    'FHDU_SectoredFluxUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
parts = cell(numel(sourceFiles), 1);
sha256 = strings(numel(sourceFiles), 1);
records = zeros(numel(sourceFiles), 1);
cadence = strings(numel(sourceFiles), 1);
for ii = 1:numel(sourceFiles)
    parts{ii} = Voyager_Read_CDF_Product(sourceFiles(ii), 'lecp_sector_daily');
    p = parts{ii};
    resolution = lower(join(string(p.global_attributes.Time_resolution), ' '));
    if contains(resolution, 'hourly')
        cadence(ii) = "hour";
    elseif contains(resolution, 'daily')
        cadence(ii) = "day";
    else
        error('VoyagerCase1:LECPIdentity', 'Unrecognized original LECP product identity.');
    end
    if ~all(isfield(p, recordFields)) || ...
            size(p.FHDU_SectoredFluxes, 2) ~= 16 || ...
            size(p.FHDU_SectoredFluxes, 3) ~= 8 || ...
            size(p.FHDU_EnergyRange, 2) ~= 2 || ...
            size(p.FHDU_EnergyRange, 3) ~= 16 || ...
            ~isequal(p.SectorIterator(:).', 1:8)
        error('VoyagerCase1:LECPNativeSchema', 'Unexpected schema: %s', sourceFiles(ii));
    end
    % P1 is hydrogen channel 10 in the audited final-CDF product.
    labels = strtrim(string(p.Hydrogen_Channels_Label));
    if ~isequal(find(labels == "P1"), 10)
        error('VoyagerCase1:P1Label', 'Original P1 channel label/order changed.');
    end
    good = isfinite(p.FHDU_Energy(:, 10));
    if any(abs(p.FHDU_Energy(good, 10) - 0.73) > 1e-5)
        error('VoyagerCase1:P1Energy', 'P1 energy changed: %s', sourceFiles(ii));
    end
    if ii > 1
        a = parts{1}.variable_meta.FHDU_SectoredFluxes.attributes;
        b = p.variable_meta.FHDU_SectoredFluxes.attributes;
        if ~isequal(a.UNITS, b.UNITS) || ...
                ~isequal(parts{1}.Hydrogen_Channels, p.Hydrogen_Channels) || ...
                cadence(ii) ~= cadence(1)
            error('VoyagerCase1:LECPUnits', 'LECP units or channel order differ.');
        end
    end
    sha256(ii) = Case1_File_SHA256(char(sourceFiles(ii)));
    records(ii) = numel(p.Epoch);
end

%% concatenate record-varying arrays only
% Keep an explicit field list so no unmerged first-file record array leaks.
product = struct('available', true, 'profile', parts{1}.profile, ...
    'reader', parts{1}.reader, 'variable_meta', parts{1}.variable_meta, ...
    'SectorIterator', parts{1}.SectorIterator, ...
    'Hydrogen_Channels', parts{1}.Hydrogen_Channels, ...
    'Hydrogen_Channels_Label', parts{1}.Hydrogen_Channels_Label);
for jj = 1:numel(recordFields)
    field = recordFields{jj};
    values = cellfun(@(p) p.(field), parts, 'UniformOutput', false);
    product.(field) = cat(1, values{:});
end
[epoch, order] = sort(product.Epoch);
sourceFileIndex = repelem((1:numel(sourceFiles)).', records);
sourceRecord = cell(numel(sourceFiles), 1);
for ii = 1:numel(sourceFiles), sourceRecord{ii} = (1:records(ii)).'; end
sourceRecord = vertcat(sourceRecord{:});
sourceFileIndex = sourceFileIndex(order);
sourceRecord = sourceRecord(order);
for jj = 1:numel(recordFields)
    field = recordFields{jj};
    value = product.(field);
    product.(field) = value(order, :, :);
end
duplicate = find(diff(epoch) == seconds(0));
for ii = duplicate(:).'
    for jj = 2:numel(recordFields)
        field = recordFields{jj};
        value = product.(field);
        if ~isequaln(value(ii, :, :), value(ii+1, :, :))
            error('VoyagerCase1:ConflictingEpoch', ...
                'Conflicting LECP records at %s; no automatic priority.', char(epoch(ii)));
        end
    end
end
[~, keep] = unique(epoch, 'stable');
for jj = 1:numel(recordFields)
    field = recordFields{jj};
    value = product.(field);
    product.(field) = value(keep, :, :);
end
product.source_file = cellstr(sourceFiles);
product.AccumulationStartUTC = product.Epoch;
product.AccumulationEndUTC = product.Epoch + seconds(product.DeltaT(:));
product.SourceManifest = table(sourceFiles, sha256, records, cadence, ...
    'VariableNames', {'SourceFile', 'SHA256', 'Records', 'Cadence'});
product.SourceMetadata = cellfun(@(p) struct('Global', p.global_attributes, ...
    'Variables', p.variable_meta), parts, 'UniformOutput', false);
product.ProductCadence = cadence(1);
product.SourceFileIndex = sourceFileIndex(keep);
product.SourceRecordNumber = sourceRecord(keep);
product.DuplicateIdenticalRecordsRemoved = numel(epoch) - numel(keep);
end
