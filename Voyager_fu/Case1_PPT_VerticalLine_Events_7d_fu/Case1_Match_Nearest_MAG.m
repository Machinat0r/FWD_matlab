function result = Case1_Match_Nearest_MAG(epoch, candidates, maxGapHours)
%Case1_Match_Nearest_MAG Match genuine source magnetic vectors by nearest Epoch.
%   result = Case1_Match_Nearest_MAG(epoch, candidates, maxGapHours)
%   candidates is a table with EpochUTC, BR, BT, BN, SourceCDF and
%   SourceCDFRecord. EpochUTC and epoch must be zoned datetime arrays.
%   maxGapHours is explicit, inclusive and may be Inf for an unlimited audit.
%   A source vector is eligible when its Epoch and all three components are
%   finite and the vector is nonzero. Values are copied from a single source
%   row without averaging, component-wise matching or magnetic interpolation.
%   At an equal time distance, choose the earlier source Epoch. Identical
%   duplicate vectors keep their earliest input row. Conflicting duplicate
%   vectors at the same Epoch raise an explicit error for provenance review.
%   The caller controls which missing records are eligible for this fallback.
%   This helper never replaces an existing mean unless the caller requests it.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% input check
narginchk(3, 3);
validateattributes(maxGapHours, {'numeric'}, ...
    {'real', 'scalar', 'nonnegative', 'nonnan'}, mfilename, 'maxGapHours');
if ~isdatetime(epoch) || (~isvector(epoch) && ~isempty(epoch)) || isempty(epoch.TimeZone)
    error('VoyagerCase1:MAGQueryEpoch', 'epoch must be a zoned datetime vector.');
end
required = {'EpochUTC', 'BR', 'BT', 'BN', 'SourceCDF', 'SourceCDFRecord'};
if ~istable(candidates) || ~all(ismember(required, candidates.Properties.VariableNames))
    error('VoyagerCase1:MAGCandidateFields', 'The magnetic source table is missing required fields.');
end
if ~isdatetime(candidates.EpochUTC) || size(candidates.EpochUTC, 2) ~= 1 || ...
        isempty(candidates.EpochUTC.TimeZone)
    error('VoyagerCase1:MAGCandidateEpoch', 'candidates.EpochUTC must be a zoned datetime column.');
end
for ii = 2:4
    validateattributes(candidates.(required{ii}), {'numeric'}, ...
        {'real', 'size', [height(candidates), 1]}, mfilename, required{ii});
end
validateattributes(candidates.SourceCDFRecord, {'numeric'}, ...
    {'real', 'size', [height(candidates), 1]}, mfilename, 'SourceCDFRecord');
if size(candidates.SourceCDF, 2) ~= 1
    error('VoyagerCase1:MAGSourceCDF', 'SourceCDF must have one path per source row.');
end
epoch = epoch(:);
epoch.TimeZone = 'UTC';
sourceTime = candidates.EpochUTC;
sourceTime.TimeZone = 'UTC';
sourceB = [double(candidates.BR), double(candidates.BT), double(candidates.BN)];

%% output and audit fields
n = numel(epoch);
result = struct;
result.QueryEpochUTC = epoch;
result.SourceRow = NaN(n, 1);
result.SourceEpochUTC = NaT(n, 1, 'TimeZone', 'UTC');
result.DeltaSeconds = NaN(n, 1);
result.AbsDeltaHours = NaN(n, 1);
result.Accepted = false(n, 1);
result.B_RTN = NaN(n, 3);
result.CandidateB_RTN = NaN(n, 3);
result.SourceCDF = strings(n, 1);
result.SourceCDFRecord = NaN(n, 1);
result.Status = repmat("invalid_query_epoch", n, 1);
result.Status(~isnat(epoch)) = "no_valid_mag_candidate";
result.MaxGapHours = double(maxGapHours);
result.ValidSourceRows = zeros(0, 1);
result.ConflictingDuplicateSourceRows = zeros(0, 1);
result.DuplicateIdenticalSourceRows = zeros(0, 1);
result.Policy = ['Nearest original magnetic product Epoch; finite complete ', ...
    'nonzero source vector copied unchanged; explicit inclusive maximum ', ...
    'absolute gap; equal distance selects earlier Epoch; no averaging or ', ...
    'interpolation of magnetic values. Identical duplicate source vectors ', ...
    'keep the earliest input row; conflicting duplicate Epochs raise an error. ', ...
    'Source candidate metadata remain available even when the gap is rejected.'];

%% eligible original source rows and duplicate-time audit
valid = ~isnat(sourceTime) & all(isfinite(sourceB), 2) & any(sourceB ~= 0, 2);
sourceRows = find(valid);
if isempty(sourceRows), return; end
[~, order] = sort(sourceTime(sourceRows));
sourceRows = sourceRows(order);
first = [1; find(diff(sourceTime(sourceRows)) ~= seconds(0)) + 1];
last = [first(2:end) - 1; numel(sourceRows)];
keep = false(numel(sourceRows), 1);
for ii = 1:numel(first)
    group = first(ii):last(ii);
    rows = sourceRows(group);
    same = all(all(sourceB(rows, :) == sourceB(rows(1), :), 2));
    if same
        % Explicit input-row ordering also makes duplicate handling stable.
        [~, earliest] = min(rows);
        keep(group(earliest)) = true;
        if numel(rows) > 1
            rows(earliest) = [];
            result.DuplicateIdenticalSourceRows = ...
                [result.DuplicateIdenticalSourceRows; rows(:)]; %#ok<AGROW>
        end
    else
        error('VoyagerCase1:MAGConflictingDuplicate', ...
            'Conflicting magnetic vectors at %s UTC in candidate rows %s.', ...
            char(sourceTime(rows(1))), mat2str(rows(:).'));
    end
end
sourceRows = sourceRows(keep);
result.ValidSourceRows = sourceRows;
queryRows = find(~isnat(epoch));
if isempty(sourceRows) || isempty(queryRows), return; end

%% reuse IRFU nearest matching on row indices, then copy raw magnetic values
referenceTime = sourceTime(sourceRows(1));
sourceSeconds = seconds(sourceTime(sourceRows) - referenceTime);
querySeconds = seconds(epoch(queryRows) - referenceTime);
if numel(sourceRows) == 1
    selected = ones(numel(queryRows), 1);
else
    if exist('irf_resamp', 'file') ~= 2
        error('VoyagerCase1:IRFUNearestMissing', ...
            'Initialize IRFU with Case1_Add_IRFU_Path before magnetic matching.');
    end
    indexSeries = [sourceSeconds, (1:numel(sourceRows)).'];
    match = irf_resamp(indexSeries, querySeconds, 'nearest');
    selected = round(match(:, 2));
    % INTERP1 nearest selects the later point at a midpoint. Apply our
    % documented earlier-Epoch tie rule using the original time distances.
    previous = max(selected - 1, 1);
    selectedDistance = abs(querySeconds - sourceSeconds(selected));
    previousDistance = abs(querySeconds - sourceSeconds(previous));
    usePrevious = previous < selected & previousDistance <= selectedDistance;
    selected(usePrevious) = previous(usePrevious);
end
rows = sourceRows(selected);
delta = seconds(sourceTime(rows) - epoch(queryRows));
accepted = abs(delta)/3600 <= maxGapHours;
result.SourceRow(queryRows) = rows;
result.SourceEpochUTC(queryRows) = sourceTime(rows);
result.DeltaSeconds(queryRows) = delta;
result.AbsDeltaHours(queryRows) = abs(delta)/3600;
result.CandidateB_RTN(queryRows, :) = sourceB(rows, :);
result.SourceCDF(queryRows) = string(candidates.SourceCDF(rows));
result.SourceCDFRecord(queryRows) = double(candidates.SourceCDFRecord(rows));
result.Accepted(queryRows) = accepted;
result.Status(queryRows) = "gap_exceeds_limit";
result.Status(queryRows(accepted)) = "nearest_valid_source";
result.B_RTN(queryRows(accepted), :) = sourceB(rows(accepted), :);
end
