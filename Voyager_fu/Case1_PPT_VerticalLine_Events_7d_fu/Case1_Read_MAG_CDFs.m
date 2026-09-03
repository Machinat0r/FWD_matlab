function [samples, sourceManifest] = Case1_Read_MAG_CDFs(sourceFiles)
%Case1_Read_MAG_CDFs Read genuine hourly magnetic vectors with provenance.
%   Original COHO CDF values are read using the existing IRFU reader. No
%   averaging, interpolation or gap filling is performed by this function.
%   All records remain in the candidate table, including source fill values
%   decoded as NaN. The nearest-match helper checks complete valid vectors.

%% original files and their exact record identities
sourceFiles = unique(string(sourceFiles(:)), 'stable');
samples = table;
sha256 = strings(numel(sourceFiles), 1);
records = zeros(numel(sourceFiles), 1);
validVectors = zeros(numel(sourceFiles), 1);
firstEpoch = NaT(numel(sourceFiles), 1, 'TimeZone', 'UTC');
lastEpoch = firstEpoch;
for ii = 1:numel(sourceFiles)
    p = Voyager_Read_CDF_Product(sourceFiles(ii), 'coho');
    if ~all(isfield(p, {'Epoch','BR','BT','BN'}))
        error('VoyagerCase1:MAGVariables', ...
            'Original COHO vector variables unavailable: %s', sourceFiles(ii));
    end
    n = numel(p.Epoch);
    part = table(p.Epoch(:), p.BR(:), p.BT(:), p.BN(:), ...
        repmat(sourceFiles(ii), n, 1), (1:n).', ...
        'VariableNames', {'EpochUTC','BR','BT','BN','SourceCDF','SourceCDFRecord'});
    samples = [samples; part]; %#ok<AGROW>
    sha256(ii) = string(Case1_File_SHA256(char(sourceFiles(ii))));
    records(ii) = n;
    B = part{:, {'BR','BT','BN'}};
    validVectors(ii) = nnz(all(isfinite(B), 2) & vecnorm(B, 2, 2) > 0);
    if n > 0
        firstEpoch(ii) = min(part.EpochUTC);
        lastEpoch(ii) = max(part.EpochUTC);
    end
end
if isempty(samples)
    samples = table(NaT(0,1,'TimeZone','UTC'), zeros(0,1), zeros(0,1), ...
        zeros(0,1), strings(0,1), zeros(0,1), ...
        'VariableNames', {'EpochUTC','BR','BT','BN','SourceCDF','SourceCDFRecord'});
else
    samples = sortrows(samples, 'EpochUTC');
end
sourceManifest = table(sourceFiles, sha256, records, validVectors, firstEpoch, ...
    lastEpoch, 'VariableNames', {'SourceFile','SHA256','Records','ValidVectors', ...
    'FirstEpochUTC','LastEpochUTC'});
end
