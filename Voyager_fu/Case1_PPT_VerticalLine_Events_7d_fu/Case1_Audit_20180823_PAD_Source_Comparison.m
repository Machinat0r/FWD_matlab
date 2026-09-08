function Report = Case1_Audit_20180823_PAD_Source_Comparison
% Read the existing 2018-08-23 hourly PAD audit only.
% Preserve each visible S1--S7 value independently and compare stored L1
% S4 means with replaced L2 values. No new scientific selection, averaging,
% data correction, plotting or source modification is applied.
% Output rounding comparisons describe numerical agreement only.

%% existing source audit
SourceMAT = 'Z:/SPART-WORK/Data/Voyager/voyager1/lecp/1h/derived/pitch_angle/2013-2021/predicted_ck/V1_Case1-S02-L06_20180823_20180823_COHO1h_raw_LECP_P1_pitch_angle_predictedCK_1h_nativeCDF_Epoch.mat';
Source = load(SourceMAT);
PAD = Source.pitchAngleTable;
Audit = Source.l1FallbackAudit;
Flux = PAD{:,cellstr(compose('Flux_S%d_1h',1:7))};
Pitch = PAD{:,cellstr(compose('PA_S%d_deg',1:7))};
Visible = find(PAD.PADUsable);

%% current production display records, one row per original sector
Groups = table;
for ii = 1:numel(Visible)
    Row = Visible(ii);
    Valid = isfinite(Pitch(Row,:)) & isfinite(Flux(Row,:)) & Flux(Row,:)>0;
    Sectors = find(Valid);
    [Angles,Order] = sort(Pitch(Row,Valid));
    Values = Flux(Row,Valid);
    Values = Values(Order);
    Sectors = Sectors(Order);
    for jj = 1:numel(Angles)
        Entry = table(Row,PAD.EpochUTC(Row),jj,1, ...
            string(sprintf('S%d',Sectors(jj))), ...
            Angles(jj),Values(jj),Values(jj),{Angles(jj)},{Values(jj)}, ...
            'VariableNames',{'PADRow','EpochUTC','Group','NumberOfSectors', ...
            'Sectors','DisplayPA_deg','DisplayFlux','LargestOriginalFlux', ...
            'OriginalPA_deg','OriginalFlux'});
        Groups = [Groups;Entry]; %#ok<AGROW>
    end
end
VisibleFlux = Flux(Visible,:);
[OriginalMax,Linear] = max(VisibleFlux,[],'all','omitnan');
[LocalRow,Sector] = ind2sub(size(VisibleFlux),Linear);
OriginalRow = Visible(LocalRow);
[DisplayMax,GroupRow] = max(Groups.DisplayFlux);

%% stored L1 means against the exact original L2 payload in replaced bins
Replaced = Audit.ReplacedL2;
Bins = dateshift(Replaced.Epoch,'start','hour');
[Found,CandidateRows] = ismember(Bins,Audit.Candidates.BinStartUTC);
assert(all(Found),'Every replaced original L2 row must have a stored candidate.');
assert(all(Audit.Candidates.Applied(CandidateRows)), ...
    'Comparison is restricted to already-applied L1 replacement bins.');
L1Mean = Audit.Candidates.MeanRate(CandidateRows,4);
L2Flux = reshape(Replaced.FHDU_SectoredFluxes(:,10,4),[],1);
Difference = L2Flux-L1Mean;
ExactEquality = L2Flux==L1Mean;
% Two decimals are used only to report agreement with displayed source
% precision. This comparison never controls source or scientific selection.
RoundedToTwoDecimalsEqual = round(L2Flux,2)==round(L1Mean,2);
DerivedL1 = L1Mean*Audit.RateToFluxFactor;
Comparison = table(Bins,Replaced.Epoch,Replaced.SourceRecordNumber,L1Mean, ...
    L2Flux,Difference,ExactEquality,RoundedToTwoDecimalsEqual,DerivedL1, ...
    DerivedL1./L2Flux,'VariableNames',{'BinStartUTC','OriginalL2EpochUTC', ...
    'OriginalL2CDFRecord','StoredL1S4MeanRate','OriginalL2S4Flux', ...
    'L2MinusL1Rate','ExactlyEqual','EqualAfterTwoDecimalRounding', ...
    'DerivedL1S4Flux','DerivedL1OverOriginalL2'});

%% return complete read-only results
Report = struct;
Report.SourceMAT = SourceMAT;
Report.WindowStartUTC = Audit.StartUTC;
Report.WindowEndUTC = Audit.EndUTC;
Report.SourcePriority = Audit.SourcePriority;
Report.SectorMergeApplied = false;
if isfield(Source, 'opts') && isfield(Source.opts, 'PitchMergeToleranceDeg')
    Report.DeprecatedStoredPitchMergeTolerance_deg = ...
        Source.opts.PitchMergeToleranceDeg;
end
Report.VisiblePADRecords = numel(Visible);
Report.OriginalVisibleSectorValues = numel(VisibleFlux);
Report.DisplayGroupCount = height(Groups);
Report.MergedGroupCount = nnz(Groups.NumberOfSectors>1);
Report.RecordsContainingMerge = numel(unique(Groups.PADRow(Groups.NumberOfSectors>1)));
Report.SectorValuesRemovedByGrouping = numel(VisibleFlux)-height(Groups);
Report.OriginalVisibleMaximum = OriginalMax;
Report.OriginalVisibleMaximumEpochUTC = PAD.EpochUTC(OriginalRow);
Report.OriginalVisibleMaximumSector = Sector;
Report.DisplayMaximum = DisplayMax;
Report.DisplayMaximumGroup = Groups(GroupRow,:);
Report.DisplayGroups = Groups;
Report.MergedDisplayGroups = Groups(Groups.NumberOfSectors>1,:);
Report.ReplacedS4Comparison = Comparison;
Report.ReplacedS4Rows = height(Comparison);
Report.ExactS4EqualCount = nnz(ExactEquality);
Report.S4EqualAfterTwoDecimalRoundingCount = nnz(RoundedToTwoDecimalsEqual);
Report.S4MaximumAbsoluteDifference = max(abs(Difference));
Report.S4MedianAbsoluteDifference = median(abs(Difference));
Report.S4MinimumDerivedToL2Ratio = min(DerivedL1./L2Flux);
Report.S4MaximumDerivedToL2Ratio = max(DerivedL1./L2Flux);
Report.RateToFluxFactor = Audit.RateToFluxFactor;
Report.InterpretationLimit = ['Numerical agreement between different product fields ', ...
    'does not establish physical equivalence or validate the legacy conversion.'];
disp(rmfield(Report,{'DisplayGroups','MergedDisplayGroups','ReplacedS4Comparison'}));
disp('Merged display groups (expected empty under current policy):');
disp(Report.MergedDisplayGroups);
disp('Stored S4 comparison:');
disp(Comparison);
end
