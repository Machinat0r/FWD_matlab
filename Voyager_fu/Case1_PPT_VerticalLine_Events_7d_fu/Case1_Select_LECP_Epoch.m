function [keep, audit] = Case1_Select_LECP_Epoch(data, startTime, endTime)
%Case1_Select_LECP_Epoch Apply the user-approved original-Epoch policy.
%   Reject DeltaT < 0. Keep all other records whose original Epoch is in
%   [startTime,endTime). No duration limit, UTC rebinning, averaging or split.
%   CDF fill values have already been mapped to NaN by the existing reader.

%% select original rows without changing the data
if ~all(isfield(data, {'Epoch', 'DeltaT'}))
    error('VoyagerCase1:EpochMetadataMissing', 'Original Epoch and DeltaT are required.');
end
epoch = data.Epoch(:);
deltaT = data.DeltaT(:);
if numel(epoch) ~= numel(deltaT)
    error('VoyagerCase1:EpochMetadataSize', 'Epoch and DeltaT lengths differ.');
end
inWindow = epoch >= startTime & epoch < endTime;
negative = deltaT < 0;
keep = inWindow & ~negative;
rejected = find(inWindow & negative);

%% record every excluded source row
excluded = table(rejected, epoch(rejected), deltaT(rejected), ...
    'VariableNames', {'SourceRow', 'EpochUTC', 'DeltaT_s'});
if all(isfield(data, {'SourceFileIndex', 'SourceRecordNumber', 'SourceManifest'}))
    sourceFile = data.SourceManifest.SourceFile(data.SourceFileIndex(rejected));
    sourceRecord = data.SourceRecordNumber(rejected);
    % Scalar source tables can return a 0-by-0 result for an empty index.
    % The exclusion table requires 0-by-1 columns; source values stay exact.
    excluded.SourceCDF = sourceFile(:);
    excluded.SourceCDFRecord = sourceRecord(:);
end
audit = struct;
audit.Policy = 'epoch_drop_negative_deltat';
audit.InputRecords = nnz(inWindow);
audit.NegativeDeltaTRejected = numel(rejected);
audit.RetainedRecords = nnz(keep);
audit.ExcludedRecords = excluded;
audit.KeptSourceRows = find(keep);
audit.MissingDeltaTRetained = nnz(keep & isnan(deltaT));
audit.Description = ['Official sector flux and uncertainty unchanged; ', ...
    'anchor at original Epoch; discard DeltaT<0 only; no upper duration ', ...
    'limit, no cross-boundary rejection, no flux interpolation or rebinning.'];
end
