function output = Case1_Apply_L1_Fallback(l2, rates, startTime, endTime, cadence, sourcePriority)
%Case1_Apply_L1_Fallback Combine whole L1/L2 bins with explicit priority.
%   Default sourcePriority='l2_first' preserves every L2 row in a bin with
%   any complete positive S1-S7 L2 record. Otherwise a complete L1 mean
%   may replace that bin's incomplete L2 records or populate an empty bin.
%   sourcePriority='l1_first' prefers a complete positive S1-S7 L1 bin
%   mean, replacing all L2 rows in that bin, including complete L2 rows.
%   If L1 is incomplete or absent, the selected original L2 rows remain.
%   Sectors from different levels are never mixed. No magnetic field or
%   attitude is calculated, interpolated or filled by this function.
%
%   L1 uses finite nonnegative rates independently per sector, with zero
%   included. Sigma is sqrt(sum(sigma.^2))/N only if every contributing rate
%   has a finite nonnegative sigma. J = mean(R)/(0.44*(1.78-0.57)). No S8
%   subtraction. L1 is explicitly derived, timed at bin center; its DeltaT
%   is NaN, and every actual original DeltaT is preserved in the audit.

%% input and original-Epoch selections
if nargin < 6, sourcePriority = 'l2_first'; end
sourcePriority = string(sourcePriority);
if ~isscalar(sourcePriority) || ~any(sourcePriority == ["l2_first", "l1_first"])
    error('VoyagerCase1:L1SourcePriority', 'Source priority must be l2_first or l1_first.');
end
cadence = string(cadence);
if ~isscalar(cadence) || ~any(cadence == ["day", "hour"])
    error('VoyagerCase1:L1Cadence', 'Cadence must be day or hour.');
end
if endTime <= startTime
    error('VoyagerCase1:L1Window', 'A positive half-open UTC time window is required.');
end
if cadence == "day", width = days(1); else, width = hours(1); end
factor = 1/(0.44*(1.78-0.57));
[l2Keep, l2Selection] = Case1_Select_LECP_Epoch(l2, startTime, endTime);
[l1Keep, l1Selection] = Case1_Select_LECP_Epoch(rates, startTime, endTime);
l2Rows = find(l2Keep); l1Rows = find(l1Keep);
l2Bin = dateshift(l2.Epoch, 'start', char(cadence));
l1Bin = dateshift(rates.Epoch, 'start', char(cadence));
bins = unique([l2Bin(l2Rows); l1Bin(l1Rows)]);
l2P1 = reshape(l2.FHDU_SectoredFluxes(:, 10, :), numel(l2.Epoch), 8);
l1P1 = reshape(rates.FHDU_SectoredRates(:, 10, :), numel(rates.Epoch), 8);
l1Sigma = reshape(rates.FHDU_SectoredRateUncertainties(:, 10, :), numel(rates.Epoch), 8);
completeL2 = all(isfinite(l2P1(:, 1:7)) & l2P1(:, 1:7) > 0, 2);

%% compute candidates without modifying either source product
nBin = numel(bins);
meanRate = nan(nBin, 8); meanSigma = nan(nBin, 8);
sampleCount = zeros(nBin, 8); applied = false(nBin, 1);
decision = strings(nBin, 1); sourceRecords = cell(nBin, 1);
candidateL1Rows = cell(nBin, 1); candidateL2Rows = cell(nBin, 1);
binEnergy = nan(nBin, 16); binEnergyRange = nan(nBin, 2, 16);
for ii = 1:nBin
    rr = l1Rows(l1Bin(l1Rows) == bins(ii));
    ll = l2Rows(l2Bin(l2Rows) == bins(ii));
    candidateL1Rows{ii} = rr; candidateL2Rows{ii} = ll;
    sourceRecords{ii} = rateRecordTable(rates, rr, l1P1, l1Sigma);
    for sector = 1:8
        valid = isfinite(l1P1(rr, sector)) & l1P1(rr, sector) >= 0;
        useRows = rr(valid); n = numel(useRows);
        sampleCount(ii, sector) = n;
        if n == 0, continue, end
        meanRate(ii, sector) = mean(l1P1(useRows, sector));
        sigma = l1Sigma(useRows, sector);
        if all(isfinite(sigma) & sigma >= 0)
            meanSigma(ii, sector) = sqrt(sum(sigma.^2))/n;
        end
    end
    if sourcePriority == "l2_first" && any(completeL2(ll))
        decision(ii) = "preserve_all_L2_rows_complete_L2_in_bin";
    elseif isempty(rr)
        decision(ii) = "no_L1_Epoch_in_bin";
        if sourcePriority == "l1_first"
            decision(ii) = "l1_first_no_L1_Epoch_keep_original_L2";
        end
    elseif ~all(isfinite(meanRate(ii, 1:7)) & meanRate(ii, 1:7) > 0)
        decision(ii) = "L1_mean_missing_or_nonpositive_sector";
        if sourcePriority == "l1_first"
            decision(ii) = "l1_first_incomplete_L1_mean_keep_original_L2";
        end
    elseif bins(ii)+width/2 < startTime || bins(ii)+width/2 >= endTime
        decision(ii) = "derived_bin_center_outside_requested_window";
        if sourcePriority == "l1_first"
            decision(ii) = "l1_first_L1_bin_center_outside_window_keep_original_L2";
        end
    else
        applied(ii) = true;
        if isempty(ll), decision(ii) = "add_complete_L1_bin_mean";
        else, decision(ii) = "replace_incomplete_L2_with_complete_L1_bin_mean"; end
        if sourcePriority == "l1_first"
            if isempty(ll), decision(ii) = "l1_first_add_complete_L1_bin_mean";
            else, decision(ii) = "l1_first_replace_all_L2_with_complete_L1_bin_mean"; end
        end
        % Preserve the source P1 energy metadata; the conversion bandwidth
        % is separately documented and does not rewrite the original CDF.
        finiteRows = rr(isfinite(rates.FHDU_Energy(rr, 10)));
        if isempty(finiteRows)
            error('VoyagerCase1:L1BinEnergy', 'No finite P1 energy for a used rate bin.');
        end
        binEnergy(ii, 10) = rates.FHDU_Energy(finiteRows(1), 10);
        ranges = reshape(rates.FHDU_EnergyRange(rr, :, 10), numel(rr), 2);
        ranges = unique(ranges(all(isfinite(ranges), 2), :), 'rows');
        if size(ranges, 1) ~= 1
            error('VoyagerCase1:L1BinEnergy', 'P1 energy range missing or conflicting in a used rate bin.');
        end
        binEnergyRange(ii, :, 10) = ranges;
        l2Keep(ll) = false;
    end
end

%% copy unchanged L2 records; append explicitly derived L1 rows
fields = {'Epoch', 'DeltaT', 'FHDU_SectoredFluxes', ...
    'FHDU_SectoredFluxUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
keepRows = find(l2Keep); useBins = find(applied);
nKeep = numel(keepRows); nNew = numel(useBins); nOut = nKeep+nNew;
output = l2;
for ii = 1:numel(fields)
    value = l2.(fields{ii});
    output.(fields{ii}) = value(keepRows, :, :);
end
newFlux = nan(nNew, 16, 8); newSigma = nan(nNew, 16, 8);
newFlux(:, 10, :) = reshape(meanRate(useBins, :)*factor, nNew, 1, 8);
newSigma(:, 10, :) = reshape(meanSigma(useBins, :)*factor, nNew, 1, 8);
output.Epoch = [output.Epoch; bins(useBins)+width/2];
output.DeltaT = [output.DeltaT; nan(nNew, 1)];
output.FHDU_SectoredFluxes = cat(1, output.FHDU_SectoredFluxes, newFlux);
output.FHDU_SectoredFluxUncertainties = cat(1, output.FHDU_SectoredFluxUncertainties, newSigma);
output.FHDU_Energy = [output.FHDU_Energy; binEnergy(useBins, :)];
output.FHDU_EnergyRange = cat(1, output.FHDU_EnergyRange, binEnergyRange(useBins, :, :));
manifestOffset = height(l2.SourceManifest);
output.SourceManifest = [l2.SourceManifest; rates.SourceManifest];
output.SourceMetadata = [l2.SourceMetadata(:); rates.SourceMetadata(:)];
output.source_file = cellstr(output.SourceManifest.SourceFile);
newFileIndex = zeros(nNew, 1);
for ii = 1:nNew
    % The scalar compatibility index names the first contributing CDF.
    % L1SourceRecords carries every file/record, including duplicate aliases.
    rr = candidateL1Rows{useBins(ii)};
    newFileIndex(ii) = manifestOffset+rates.SourceFileIndex(rr(1));
end
output.SourceFileIndex = [l2.SourceFileIndex(keepRows); newFileIndex];
output.SourceRecordNumber = [l2.SourceRecordNumber(keepRows); nan(nNew, 1)];
output.OriginalL2SourceRow = [keepRows; nan(nNew, 1)];
output.SourceProduct = [repmat("L2", nKeep, 1); repmat("L1_UTC_mean", nNew, 1)];
output.SourceToDifferentialFluxFactor = [ones(nKeep, 1); repmat(factor, nNew, 1)];
output.SectorSampleCount = [double(isfinite(l2P1(keepRows, :))); sampleCount(useBins, :)];
output.L1RawRate = [nan(nKeep, 8); meanRate(useBins, :)];
output.L1SourceRecords = [cell(nKeep, 1); sourceRecords(useBins)];
output.L1CandidateIndex = [nan(nKeep, 1); useBins];
output.BinStartUTC = dateshift(output.Epoch, 'start', char(cadence));
output.BinEndUTC = output.BinStartUTC+width;
[~, order] = sort(output.Epoch);
rowFields = [fields, {'SourceFileIndex', 'SourceRecordNumber', ...
    'OriginalL2SourceRow', 'SourceProduct', 'SourceToDifferentialFluxFactor', ...
    'SectorSampleCount', 'L1RawRate', 'L1SourceRecords', 'L1CandidateIndex', ...
    'BinStartUTC', 'BinEndUTC'}];
for ii = 1:numel(rowFields)
    value = output.(rowFields{ii});
    output.(rowFields{ii}) = value(order, :, :);
end
assert(numel(output.Epoch) == nOut);
output.ProductCadence = cadence;
output.AccumulationStartUTC = output.Epoch;
output.AccumulationEndUTC = output.Epoch+seconds(output.DeltaT);
output.meta = output.variable_meta;

%% full provenance, candidates and displaced original L2 payload
removedRows = setdiff(l2Rows, keepRows, 'stable');
replaced = struct('OriginalRows', removedRows);
for ii = 1:numel(fields)
    value = l2.(fields{ii}); replaced.(fields{ii}) = value(removedRows, :, :);
end
replaced.SourceFileIndex = l2.SourceFileIndex(removedRows);
replaced.SourceRecordNumber = l2.SourceRecordNumber(removedRows);
replaced.SourceCDF = l2.SourceManifest.SourceFile(replaced.SourceFileIndex);
candidates = table(bins, bins+width, bins+width/2, applied, decision, ...
    meanRate, meanSigma, sampleCount, candidateL1Rows, candidateL2Rows, sourceRecords, ...
    'VariableNames', {'BinStartUTC', 'BinEndUTC', 'EpochUTC', 'Applied', ...
    'Decision', 'MeanRate', 'MeanRateSigma', 'SectorSampleCount', ...
    'L1Rows', 'L2Rows', 'SourceRecords'});
audit = struct;
audit.SourcePriority = sourcePriority;
audit.Method = 'L2-first complete-bin L1 P1 rate-mean fallback';
priorityDescription = ['Any complete positive L2 S1-S7 record protects all L2 ', ...
    'records in its bin. Whole-bin L1 substitute requires seven positive mean rates. '];
if sourcePriority == "l1_first"
    audit.Method = 'L1-first complete-bin P1 rate means with original L2 fallback';
    priorityDescription = ['A complete positive S1-S7 L1 bin mean takes priority ', ...
        'and replaces all selected L2 records in that bin, including complete L2. ', ...
        'If L1 is absent or incomplete, all selected original L2 records remain. ', ...
        'The source-priority decision uses sector flux only and does not depend on B. '];
end
audit.MethodDescription = ['Native Epoch selects the half-open UTC bin; ', ...
    'negative DeltaT rejected, unknown/zero retained; no duration weighting ', ...
    'or splitting. Finite nonnegative rate mean independently per sector; ', ...
    'zero is included and missing remains missing. All contributing sigmas ', ...
    'required for sqrt(sum(sigma^2))/N. ', priorityDescription, ...
    'No cross-level sector mixing, ', ...
    'background subtraction, magnetic filling, geometry or threshold change.'];
audit.StartUTC = startTime; audit.EndUTC = endTime; audit.Cadence = cadence;
audit.L2RecordSelection = l2Selection; audit.L1RecordSelection = l1Selection;
audit.SourceManifest = output.SourceManifest;
audit.L1SourceMetadata = rates.SourceMetadata;
audit.L2SourceMetadata = l2.SourceMetadata;
audit.Candidates = candidates; audit.AppliedBins = candidates(applied, :);
audit.ReplacedL2 = replaced; audit.L2PreservedSourceRows = keepRows;
audit.L2PreservedRecords = nKeep; audit.L2ReplacedRecords = numel(removedRows);
audit.L1AddedRecords = nNew;
audit.RateToFluxFactor = factor;
audit.GeometryFactor_cm2_sr = 0.44;
audit.ConversionEnergyRange_MeV = [0.57, 1.78];
audit.ConversionIsUserApprovedLegacyL1Path = true;
audit.SourceCDFEnergyMetadataUnchanged = true;
audit.L1RepresentativeEpoch = 'UTC bin center: day 12:00 or hour :30';
audit.L1DeltaT = 'NaN for derived row; original per-record DeltaT in SourceRecords';
audit.S8Policy = 'Diagnostic only; never subtracted or included in seven-sector completeness';
output.L1FallbackAudit = audit;
end

function result = rateRecordTable(rates, rows, rawRate, rawSigma)
% Every selected native record is listed, including individually empty
% sectors and its exact original duration. Per-sector masks identify use.
sourceFile = rates.SourceManifest.SourceFile(rates.SourceFileIndex(rows));
sourceRecord = rates.SourceRecordNumber(rows);
result = table(rows, sourceFile(:), sourceRecord(:), rates.Epoch(rows), ...
    rates.DeltaT(rows), rawRate(rows, :), rawSigma(rows, :), ...
    isfinite(rawRate(rows, :)) & rawRate(rows, :) >= 0, ...
    'VariableNames', {'SourceRow', 'SourceCDF', 'SourceCDFRecord', 'EpochUTC', ...
    'DeltaT_s', 'Rate_S1_to_S8', 'Sigma_S1_to_S8', 'Contributes_S1_to_S8'});
if isfield(rates, 'SourceRecordAliases')
    result.IdenticalRecordAliases = rates.SourceRecordAliases(rows);
end
end
