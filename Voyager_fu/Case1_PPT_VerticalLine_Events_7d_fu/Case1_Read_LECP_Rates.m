function product = Case1_Read_LECP_Rates(sourceFiles)
%Case1_Read_LECP_Rates Read original native LECP rate CDFs with IRFU.
%   No conversion or averaging occurs here. Original Epoch, DeltaT, rates,
%   uncertainties and energy metadata are retained after the common CDF
%   FILLVAL/VALIDMIN/VALIDMAX screening. P1 is channel 10; S8 is diagnostic.

%% original files and audited rate schema
sourceFiles = string(sourceFiles);
sourceFiles = unique(sourceFiles(:), 'stable');
if isempty(sourceFiles) || any(~isfile(sourceFiles))
    error('VoyagerCase1:L1CDFMissing', 'Original LECP rate CDFs are missing.');
end
fields = {'Epoch', 'DeltaT', 'FHDU_SectoredRates', ...
    'FHDU_SectoredRateUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
parts = cell(numel(sourceFiles), 1);
sha256 = strings(numel(sourceFiles), 1);
records = zeros(numel(sourceFiles), 1);
for ii = 1:numel(sourceFiles)
    p = Voyager_Read_CDF_Product(sourceFiles(ii), 'lecp_sector_daily');
    if ~all(isfield(p, [fields, {'SectorIterator', ...
            'Hydrogen_Channels', 'Hydrogen_Channels_Label'}]))
        error('VoyagerCase1:L1Schema', 'Missing rate CDF variables: %s', sourceFiles(ii));
    end
    n = numel(p.Epoch);
    if ~isdatetime(p.Epoch) || any(isnat(p.Epoch)) || ...
            numel(p.DeltaT) ~= n || ...
            ~isequal(size(p.FHDU_SectoredRates), [n, 16, 8]) || ...
            ~isequal(size(p.FHDU_SectoredRateUncertainties), [n, 16, 8]) || ...
            ~isequal(size(p.FHDU_Energy), [n, 16]) || ...
            ~isequal(size(p.FHDU_EnergyRange), [n, 2, 16]) || ...
            ~isequal(p.SectorIterator(:).', 1:8)
        error('VoyagerCase1:L1Schema', 'Unexpected rate CDF dimensions: %s', sourceFiles(ii));
    end
    labels = strtrim(string(p.Hydrogen_Channels_Label));
    if ~isequal(find(labels == "P1"), 10)
        error('VoyagerCase1:L1P1Label', 'P1 channel label/order changed.');
    end
    energy = p.FHDU_Energy(:, 10);
    if ~any(isfinite(energy)) || any(abs(energy(isfinite(energy))-0.73) > 1e-5)
        error('VoyagerCase1:L1P1Energy', 'P1 energy changed: %s', sourceFiles(ii));
    end
    for name = ["FHDU_SectoredRates", "FHDU_SectoredRateUncertainties"]
        a = p.variable_meta.(name).attributes;
        if ~isfield(a, 'UNITS') || ...
                ~strcmpi(strtrim(string(a.UNITS)), "Counts/second")
            error('VoyagerCase1:L1Units', 'Expected Counts/second: %s', sourceFiles(ii));
        end
    end
    if ii > 1 && (~isequal(parts{1}.Hydrogen_Channels, p.Hydrogen_Channels) || ...
            ~isequal(parts{1}.Hydrogen_Channels_Label, p.Hydrogen_Channels_Label))
        error('VoyagerCase1:L1ChannelOrder', 'Rate CDF channel order differs.');
    end
    p.Epoch = p.Epoch(:); p.DeltaT = p.DeltaT(:);
    parts{ii} = p;
    records(ii) = n;
    sha256(ii) = Case1_File_SHA256(char(sourceFiles(ii)));
end

%% concatenate exact record arrays; duplicate payloads must agree
first = parts{1};
product = struct('available', true, 'profile', first.profile, ...
    'reader', first.reader, 'variable_meta', first.variable_meta, ...
    'SectorIterator', first.SectorIterator, ...
    'Hydrogen_Channels', first.Hydrogen_Channels, ...
    'Hydrogen_Channels_Label', first.Hydrogen_Channels_Label);
for jj = 1:numel(fields)
    values = cellfun(@(p) p.(fields{jj}), parts, 'UniformOutput', false);
    product.(fields{jj}) = cat(1, values{:});
end
fileIndex = repelem((1:numel(sourceFiles)).', records);
recordNumber = cell(numel(sourceFiles), 1);
for ii = 1:numel(sourceFiles), recordNumber{ii} = (1:records(ii)).'; end
recordNumber = vertcat(recordNumber{:});
[epoch, order] = sort(product.Epoch);
fileIndex = fileIndex(order); recordNumber = recordNumber(order);
for jj = 1:numel(fields)
    value = product.(fields{jj});
    product.(fields{jj}) = value(order, :, :);
end
duplicate = find(diff(epoch) == seconds(0));
for ii = duplicate(:).'
    for jj = 2:numel(fields)
        value = product.(fields{jj});
        if ~isequaln(value(ii, :, :), value(ii+1, :, :))
            error('VoyagerCase1:ConflictingL1Epoch', ...
                'Conflicting rate records at %s; no automatic priority.', char(epoch(ii)));
        end
    end
end
[~, keep, group] = unique(epoch, 'stable');
% Retain all source identities even when identical duplicate measurements
% are represented once. A duplicate is never counted twice in an average.
recordAliases = cell(numel(keep), 1);
% A singleton's original identity already lives in SourceFileIndex and
% SourceRecordNumber. Only duplicate groups need an additional alias table.
groupCount = accumarray(group, 1, [numel(keep), 1]);
for ii = find(groupCount > 1).'
    rows = find(group == ii);
    recordAliases{ii} = table(sourceFiles(fileIndex(rows)), recordNumber(rows), ...
        'VariableNames', {'SourceCDF', 'SourceCDFRecord'});
end
for jj = 1:numel(fields)
    value = product.(fields{jj});
    product.(fields{jj}) = value(keep, :, :);
end
product.source_file = cellstr(sourceFiles);
product.SourceManifest = table(sourceFiles, sha256, records, ...
    repmat("native", numel(sourceFiles), 1), 'VariableNames', ...
    {'SourceFile', 'SHA256', 'Records', 'Cadence'});
product.SourceMetadata = cellfun(@(p) struct('Global', p.global_attributes, ...
    'Variables', p.variable_meta), parts, 'UniformOutput', false);
product.SourceFileIndex = fileIndex(keep);
product.SourceRecordNumber = recordNumber(keep);
product.SourceRecordAliases = recordAliases;
product.ProductCadence = "native";
product.AccumulationStartUTC = product.Epoch;
product.AccumulationEndUTC = product.Epoch+seconds(product.DeltaT);
product.DuplicateIdenticalRecordsRemoved = numel(epoch)-numel(keep);
end
