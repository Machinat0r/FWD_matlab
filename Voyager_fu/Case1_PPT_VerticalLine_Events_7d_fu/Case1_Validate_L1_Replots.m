function result = Case1_Validate_L1_Replots(baselineFile,cadenceFilter)
%Case1_Validate_L1_Replots Independently validate requested overwritten audits.
%   Compare only the event/cadence files explicitly saved in BASELINEFILE.
%   Existing usable L2 rows are identified by CDF, CDF record and Epoch;
%   numeric science fields must remain exact (PA permits 1e-10 degree).
%   New L1 UTC means are rebuilt from their original CDF records for audit.
%   This function creates no figures and modifies no production artifact.

%% baseline, source reader and result containers
if nargin<1
    baselineFile=['Z:/SPART-WORK/Data/Voyager/voyager1/lecp/validation/', ...
        'l1_fallback/baseline_20260903_094142_200.mat'];
end
B=load(baselineFile,'baseline'); baseline=B.baseline;
originalItemCount=numel(baseline.Items);
if nargin<2, cadenceFilter=["1d","1h"]; end
cadenceFilter=string(cadenceFilter);
assert(all(ismember(cadenceFilter,["1d","1h"])), ...
    'Cadence filter must contain only 1d and/or 1h.');
baseline.Items=baseline.Items(ismember(string({baseline.Items.Cadence}),cadenceFilter));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
result=struct('BaselineFile',string(baselineFile), ...
    'BaselineSHA256',string(Case1_File_SHA256(baselineFile)), ...
    'StartedUTC',datetime('now','TimeZone','UTC'), ...
    'Scope','Only the original event/cadence audits listed in baseline.Items', ...
    'Events',table,'Summary',table,'Peaks',{{}},'FatalError',"");
result.CadenceFilter=cadenceFilter;
result.OriginalBaselineItemCount=originalItemCount;
result.ExpectedAuditCount=numel(baseline.Items);
result.Checks=table('Size',[0 5], ...
    'VariableTypes',{'string','string','string','logical','string'}, ...
    'VariableNames',{'EventID','Cadence','Check','Passed','Detail'});
result.Changes=table('Size',[0 8], ...
    'VariableTypes',{'string','string','double','double','string','double','logical','string'}, ...
    'VariableNames',{'EventID','Cadence','OldRow','NewRow','Field', ...
    'MaxAbsDifference','WithinAllowedTolerance','Detail'});
sourceCache=containers.Map('KeyType','char','ValueType','any');
hashCache=containers.Map('KeyType','char','ValueType','char');
binInventoryCache=containers.Map('KeyType','char','ValueType','any');
currentID=""; currentCadence="";
factor=1/(0.44*(1.78-0.57));
checkRows=cell(32768,5); checkPassed=false(32768,1); checkCount=0;

%% inspect only the existing counterparts of the saved baseline
for ik=1:numel(baseline.Items)
    item=baseline.Items(ik);
    currentID=string(item.EventID); currentCadence=string(item.Cadence);
    firstCheck=checkCount+1;
    old=item.Saved.pitchAngleTable;
    oldUsable=usable(old); new=table; newUsable=false(0,1);
    added=false(0,1); changed=false; errorText=""; newHash="";
    preserved=0; remainingB=0; oldPeak="not_applicable"; newPeak=oldPeak;
    oldCount=NaN; newCount=NaN; oldSelected=NaN; newSelected=NaN;
    peakChanged=false; peakRecord=struct;
    try
        addCheck('counterpart_exists',isfile(item.AuditFile),string(item.AuditFile));
        assert(isfile(item.AuditFile),'Saved audit counterpart is missing.');
        S=load(char(item.AuditFile)); new=S.pitchAngleTable;
        newUsable=usable(new); newHash=string(Case1_File_SHA256(item.AuditFile));
        changed=newHash~=string(item.AuditSHA256);
        if ismember('SourceProduct',new.Properties.VariableNames)
            added=string(new.SourceProduct)=="L1_UTC_mean";
        else
            added=false(height(new),1);
        end

        %% unchanged original usable L2 records, independent of merged row index
        originalRows=find(oldUsable); matchingRows=nan(numel(originalRows),1);
        for im=1:numel(originalRows)
            io=originalRows(im);
            needed={'SourceCDF','SourceCDFRecord','EpochUTC'};
            assert(all(ismember(needed,new.Properties.VariableNames)), ...
                'New table lacks original-source identity fields.');
            match=find(string(new.SourceCDF)==string(old.SourceCDF(io)) & ...
                new.SourceCDFRecord==old.SourceCDFRecord(io) & ...
                new.EpochUTC==old.EpochUTC(io));
            addCheck('original_usable_record_identity',numel(match)==1, ...
                sprintf('old row %d has %d source-key matches',io,numel(match)));
            if numel(match)~=1, continue, end
            j=match(1); matchingRows(im)=j; preserved=preserved+1;
            addCheck('original_usable_still_usable',newUsable(j) && ~added(j), ...
                sprintf('old row %d -> new row %d',io,j));
        end
        found=isfinite(matchingRows); originalRows=originalRows(found);
        matchingRows=matchingRows(found);
        fields=old.Properties.VariableNames;
        for iv=1:numel(fields)
            name=fields{iv}; value=old.(name);
            % SourceRow is deliberately re-indexed in the merged product.
            if strcmp(name,'SourceRow')||isempty(originalRows), continue, end
            if ~(isnumeric(value)||islogical(value)||isdatetime(value)), continue, end
            exists=ismember(name,new.Properties.VariableNames);
            addCheck('original_numeric_field_present',exists,name);
            if ~exists, continue, end
            a=value(originalRows,:); b=new.(name)(matchingRows,:);
            tolerance=0;
            if startsWith(name,'PA_S'), tolerance=1e-10; end
            [same,delta,exact]=sameValue(a,b,tolerance);
            addCheck('original_usable_numeric_preserved',same, ...
                sprintf('%d matched original rows: %s; max delta %.12g',numel(originalRows),name,delta));
            if ~exact
                for im=1:numel(originalRows)
                    [rowSame,rowDelta,rowExact]=sameValue(a(im,:),b(im,:),tolerance);
                    if rowExact, continue, end
                    result.Changes(end+1,:)={currentID,currentCadence,originalRows(im), ...
                        matchingRows(im),string(name),rowDelta,rowSame,"Original usable L2 record"};
                end
            end
        end

        %% explicit new L1 provenance and independent bin calculations
        cadence='day'; suffix='1d'; width=days(1);
        if currentCadence=="1h", cadence='hour'; suffix='1h'; width=hours(1); end
        if any(added)
            A=findFallbackAudit(S,new);
            assert(isstruct(A)&&isfield(A,'AppliedBins'), ...
                'L1 means require an AppliedBins audit.');
            addCheck('derived_row_count_matches_audit', ...
                nnz(added)==height(A.AppliedBins),sprintf('%d table / %d audit', ...
                nnz(added),height(A.AppliedBins)));
            addCheck('legacy_conversion_declared', ...
                isfield(A,'RateToFluxFactor') && A.RateToFluxFactor==factor, ...
                'Explicit approved legacy factor, not an official L2 calibration.');
            for j=find(added).'
                begin=dateshift(new.EpochUTC(j),'start',cadence);
                addCheck('derived_center_not_native_epoch', ...
                    new.EpochUTC(j)==begin+width/2 && isnan(new.SourceCDFRecord(j)) && ...
                    isnan(new.DeltaT_s(j)),sprintf('new row %d',j));
                addCheck('derived_factor_exact', ...
                    new.SourceToDifferentialFluxFactor(j)==factor,sprintf('new row %d',j));
                candidate=find(A.AppliedBins.EpochUTC==new.EpochUTC(j));
                assert(numel(candidate)==1,'Derived row does not identify one applied bin.');
                records=A.AppliedBins.SourceRecords{candidate};
                assert(istable(records)&&height(records)>0,'Missing original L1 records.');
                if ismember('L1SourceRecords',new.Properties.VariableNames)
                    addCheck('row_and_audit_source_records_identical', ...
                        isequaln(new.L1SourceRecords{j},records),sprintf('new row %d',j));
                end
                addCheck('original_epochs_within_bin_and_deltat', ...
                    all(records.EpochUTC>=begin & records.EpochUTC<begin+width) && ...
                    ~any(records.DeltaT_s<0),sprintf('new row %d',j));
                expectedKeys=expectedBinKeys(A.SourceManifest,begin,cadence);
                actualKeys=string(records.SourceCDF)+"#"+string(records.SourceCDFRecord);
                if ismember('IdenticalRecordAliases',records.Properties.VariableNames)
                    for ia=1:height(records)
                        aliases=records.IdenticalRecordAliases{ia};
                        if istable(aliases)&&~isempty(aliases)
                            actualKeys=[actualKeys;string(aliases.SourceCDF)+"#"+ ...
                                string(aliases.SourceCDFRecord)]; %#ok<AGROW>
                        end
                    end
                end
                addCheck('all_native_bin_records_accounted_for', ...
                    isequal(sort(unique(actualKeys)),sort(unique(expectedKeys))) && ...
                    numel(unique(records.EpochUTC))==height(records),sprintf('new row %d',j));
                raw=nan(height(records),8); sigma=raw;
                for ir=1:height(records)
                    sourceName=char(records.SourceCDF(ir));
                    if ~isKey(sourceCache,sourceName)
                        sourceCache(sourceName)=Voyager_Read_CDF_Product(sourceName,'lecp_sector_daily');
                        hashCache(sourceName)=Case1_File_SHA256(sourceName);
                    end
                    P=sourceCache(sourceName); sr=records.SourceCDFRecord(ir);
                    assert(isfinite(sr)&&sr>=1&&sr==fix(sr)&&sr<=numel(P.Epoch), ...
                        'Invalid native CDF record index.');
                    raw(ir,:)=reshape(P.FHDU_SectoredRates(sr,10,:),1,8);
                    sigma(ir,:)=reshape(P.FHDU_SectoredRateUncertainties(sr,10,:),1,8);
                    addCheck('native_L1_payload_exact', ...
                        P.Epoch(sr)==records.EpochUTC(ir) && ...
                        isequaln(P.DeltaT(sr),records.DeltaT_s(ir)) && ...
                        isequaln(raw(ir,:),records.Rate_S1_to_S8(ir,:)) && ...
                        isequaln(sigma(ir,:),records.Sigma_S1_to_S8(ir,:)), ...
                        sprintf('new row %d source record %d',j,sr));
                    mi=find(string(A.SourceManifest.SourceFile)==string(sourceName));
                    addCheck('native_L1_hash_matches_manifest',numel(mi)==1 && ...
                        string(A.SourceManifest.SHA256(mi))==string(hashCache(sourceName)),sourceName);
                end
                expectedRate=nan(1,8); expectedSigma=nan(1,8); count=zeros(1,8);
                for is=1:8
                    good=isfinite(raw(:,is))&raw(:,is)>=0;
                    count(is)=nnz(good);
                    if count(is)==0, continue, end
                    expectedRate(is)=mean(raw(good,is));
                    q=sigma(good,is);
                    if all(isfinite(q)&q>=0)
                        expectedSigma(is)=sqrt(sum(q.^2))/count(is);
                    end
                end
                flux=sectorArray(new(j,:),'RawFlux',suffix);
                actualRate=sectorArray(new(j,:),'RawRate',suffix);
                actualSigma=sectorArray(new(j,:),'FluxUncertainty',suffix);
                actualCount=sectorArray(new(j,:),'Samples',suffix);
                addCheck('all_seven_derived_flux_positive', ...
                    all(isfinite(flux(1:7))&flux(1:7)>0),sprintf('new row %d',j));
                addCheck('derived_rate_mean_exact',isequaln(actualRate,expectedRate),sprintf('new row %d',j));
                addCheck('derived_flux_from_rate_exact', ...
                    isequaln(flux,actualRate*factor) && isequaln(flux,expectedRate*factor),sprintf('new row %d',j));
                addCheck('derived_sigma_exact', ...
                    isequaln(actualSigma,expectedSigma*factor),sprintf('new row %d',j));
                addCheck('derived_sample_counts_exact', ...
                    isequaln(actualCount,count) && ...
                    isequaln(A.AppliedBins.SectorSampleCount(candidate,:),count),sprintf('new row %d',j));
                if ismember('Contributes_S1_to_S8',records.Properties.VariableNames)
                    addCheck('contributing_sector_masks_exact', ...
                        isequal(records.Contributes_S1_to_S8,isfinite(raw)&raw>=0),sprintf('new row %d',j));
                end
                if ~isempty(old)
                    oldFlux=sectorArray(old,'RawFlux',suffix);
                    protected=all(isfinite(oldFlux(:,1:7))&oldFlux(:,1:7)>0,2) & ...
                        dateshift(old.EpochUTC,'start',cadence)==begin;
                    addCheck('no_derived_mean_in_complete_L2_bin', ...
                        ~any(protected),sprintf('new row %d',j));
                end
            end
        end
        if ~isempty(new)
            rawFlux=sectorArray(new,'RawFlux',suffix);
            angles=angleArray(new);
            expectedUsable=all(isfinite(rawFlux(:,1:7))&rawFlux(:,1:7)>0,2) & ...
                all(isfinite(angles(:,1:7)),2);
            addCheck('PAD_gate_unchanged',isequal(newUsable,expectedUsable), ...
                'Only seven finite positive fluxes and finite PAs.');
            bnames={sprintf('BR_%s_nT',cadenceName(cadence)), ...
                sprintf('BT_%s_nT',cadenceName(cadence)),sprintf('BN_%s_nT',cadenceName(cadence))};
            V=new{:,bnames};
            remainingB=nnz(all(isfinite(rawFlux(:,1:7))&rawFlux(:,1:7)>0,2) & ...
                (~all(isfinite(V),2)|vecnorm(V,2,2)==0));
        end

        %% actual five-panel selection, rendering disabled
        if currentCadence=="1h"
            er=baseline.Catalog(string(baseline.Catalog.EventID)==currentID,:);
            assert(height(er)==1,'Baseline event identity is ambiguous.');
            er.PlotStartUTC=er.StartUTC-days(item.Saved.opts.ContextDays);
            er.PlotEndUTCExclusive=er.EndUTCExclusive+days(item.Saved.opts.ContextDays);
            options=item.Saved.opts; options.RenderFigure=false; options.PADCadence='hour';
            pold=Case1_Plot_Peak_PAD(old,er,'',options);
            pnew=Case1_Plot_Peak_PAD(new,er,'',options);
            assert(~pold.FigureCreated&&~pnew.FigureCreated,'Unexpected figure creation.');
            stamp=char(datetime(er.StartUTC,'Format','yyyyMMdd'));
            peakName=['V1_',char(currentID),'_',stamp,'_P1_peak5_hourly.mat'];
            peakFile=fullfile(item.Saved.opts.PeakPADDataFolder,peakName);
            addCheck('actual_peak_audit_exists',isfile(peakFile),string(peakFile));
            assert(isfile(peakFile),'Actual production peak MAT is missing.');
            savedPeak=load(peakFile,'peakAudit'); actual=savedPeak.peakAudit;
            peakFields={'SelectedEpochUTC','SelectedCount','PlottedCount', ...
                'SelectedTableRows','SelectedPADUsable','PeakEpochUTC','PeakFlux', ...
                'PeakPADUsable','RawFlux','RawSigma','PA_deg','NormalizationFlux', ...
                'NormalizedFlux','NormalizedSigma','Status'};
            for ip=1:numel(peakFields)
                name=peakFields{ip};
                addCheck('actual_peak_payload_matches_recalculation', ...
                    isfield(actual,name)&&isequaln(actual.(name),pnew.(name)),name);
            end
            peakFigure=string(actual.FigureFile);
            addCheck('actual_peak_PNG_exists',isfile(peakFigure)&&actual.FigureCreated,peakFigure);
            assert(isfile(peakFigure),'Actual production peak PNG is missing.');
            picture=imfinfo(peakFigure);
            addCheck('actual_peak_PNG_decodes',strcmpi(picture.Format,'png') && ...
                picture.Width>0&&picture.Height>0, ...
                sprintf('%d x %d',picture.Width,picture.Height));
            peakMATHash=string(Case1_File_SHA256(peakFile));
            peakPNGHash=string(Case1_File_SHA256(peakFigure));
            oldPeak=string(pold.Status); newPeak=string(pnew.Status);
            oldCount=pold.PlottedCount; newCount=pnew.PlottedCount;
            oldSelected=pold.SelectedCount; newSelected=pnew.SelectedCount;
            peakChanged=~isequaln(pold.PeakEpochUTC,pnew.PeakEpochUTC) || ...
                ~isequaln(pold.PeakFlux,pnew.PeakFlux);
            peakRecord=struct('EventID',currentID,'OldStatus',oldPeak,'NewStatus',newPeak, ...
                'OldPeakEpochUTC',pold.PeakEpochUTC,'NewPeakEpochUTC',pnew.PeakEpochUTC, ...
                'OldPeakFlux',pold.PeakFlux,'NewPeakFlux',pnew.PeakFlux, ...
                'OldSelectedEpochUTC',pold.SelectedEpochUTC,'NewSelectedEpochUTC',pnew.SelectedEpochUTC, ...
                'OldPlottedCount',oldCount,'NewPlottedCount',newCount, ...
                'ActualPeakMAT',string(peakFile),'ActualPeakMATSHA256',peakMATHash, ...
                'ActualPeakPNG',peakFigure,'ActualPeakPNGSHA256',peakPNGHash, ...
                'ActualPNGDimensions',[picture.Width,picture.Height]);
            result.Peaks{end+1,1}=peakRecord;
        end
    catch ME
        errorText=string(getReport(ME,'extended','hyperlinks','off'));
        addCheck('event_validation_completed',false,string(ME.message));
    end
    passed=all(checkPassed(firstCheck:checkCount))&&strlength(errorText)==0;
    row=table(currentID,currentCadence,string(item.AuditFile),changed, ...
        string(item.AuditSHA256),newHash,height(old),nnz(oldUsable), ...
        height(new),nnz(newUsable),nnz(added),nnz(added&newUsable), ...
        preserved,remainingB,height(old)==0,height(old)==0&&any(newUsable), ...
        oldPeak,newPeak,oldCount,newCount,oldSelected,newSelected,peakChanged,passed,errorText, ...
        'VariableNames',{'EventID','Cadence','AuditFile','AuditFileChanged', ...
        'BaselineSHA256','NewSHA256','OldRows','OldUsable','NewRows','NewUsable', ...
        'AddedL1Rows','AddedL1Usable','OriginalUsableRecordsFound', ...
        'RemainingCompleteFluxMissingB','OriginallyNoL2Rows','NoL2EventRecovered', ...
        'OldPeakStatus','NewPeakStatus','OldPeakPlottedCount','NewPeakPlottedCount', ...
        'OldPeakSelectedCount','NewPeakSelectedCount','PeakIdentityChanged','Passed','Error'});
    result.Events=[result.Events;row]; %#ok<AGROW>
    fprintf('L1_VALIDATE %s %s old=%d/%d new=%d/%d L1=%d pass=%d\n', ...
        currentID,currentCadence,nnz(oldUsable),height(old),nnz(newUsable),height(new),nnz(added),passed);
end

%% compact summary and classified audit output
result.Checks=cell2table(checkRows(1:checkCount,:), ...
    'VariableNames',{'EventID','Cadence','Check','Passed','Detail'});
result.RowCountMeaning='Counts are occurrences in each requested event window; overlapping windows are not deduplicated.';
for cadence=["1d","1h"]
    T=result.Events(result.Events.Cadence==cadence,:);
    if isempty(T), continue, end
    row=table(cadence,height(T),nnz(T.AuditFileChanged),sum(T.OldRows), ...
        sum(T.OldUsable),sum(T.NewRows),sum(T.NewUsable), ...
        sum(T.AddedL1Rows),sum(T.AddedL1Usable),nnz(T.OriginallyNoL2Rows), ...
        nnz(T.NoL2EventRecovered),sum(T.RemainingCompleteFluxMissingB), ...
        nnz(T.Passed),'VariableNames',{'Cadence','Events','ChangedAudits', ...
        'OldRows','OldUsable','NewRows','NewUsable','AddedL1Rows','AddedL1Usable', ...
        'OriginallyNoL2Events','RecoveredNoL2Events','RemainingMissingB','PassedEvents'});
    result.Summary=[result.Summary;row]; %#ok<AGROW>
end
result.CompletedUTC=datetime('now','TimeZone','UTC');
result.Passed=all(result.Checks.Passed)&&all(result.Events.Passed)&& ...
    height(result.Events)==numel(baseline.Items);
result.CodeFile=[mfilename('fullpath'),'.m'];
result.CodeSHA256=Case1_File_SHA256(result.CodeFile);
folder='Z:/SPART-WORK/Data/Voyager/voyager1/lecp/validation/l1_fallback';
if ~isfolder(folder), mkdir(folder); end
result.AuditFile=fullfile(folder,['l1_replots_validation_', ...
    char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS')),'.mat']);
save(result.AuditFile,'result','-v7.3');
disp(result.Summary);
if any(~result.Checks.Passed), disp(result.Checks(~result.Checks.Passed,:)); end
fprintf('L1_REPLOT_VALIDATION passed=%d checks=%d/%d audit=%s\n', ...
    result.Passed,nnz(result.Checks.Passed),height(result.Checks),result.AuditFile);

    function addCheck(name,passed,detail)
        checkCount=checkCount+1;
        if checkCount>size(checkRows,1)
            capacity=2*size(checkRows,1);
            checkRows{capacity,5}=[];
            checkPassed(capacity)=false;
        end
        checkRows(checkCount,:)={currentID,currentCadence,string(name),logical(passed),string(detail)};
        checkPassed(checkCount)=logical(passed);
    end

    function keys=expectedBinKeys(manifest,begin,cadence)
        % Independent inventory includes every retained native record and
        % aliases, so an omitted rate cannot falsely pass a sample-count test.
        nativeFiles=string(manifest.SourceFile(string(manifest.Cadence)=="native"));
        inventoryKey=char(string(cadence)+"|"+join(nativeFiles,"|"));
        if ~isKey(binInventoryCache,inventoryKey)
            allBin=zeros(0,1); allKeys=strings(0,1);
            for fi=1:numel(nativeFiles)
                filename=char(nativeFiles(fi));
                if ~isKey(sourceCache,filename)
                    sourceCache(filename)=Voyager_Read_CDF_Product(filename,'lecp_sector_daily');
                    hashCache(filename)=Case1_File_SHA256(filename);
                end
                P=sourceCache(filename); retained=find(~(P.DeltaT<0));
                bins=posixtime(dateshift(P.Epoch(retained),'start',cadence));
                allBin=[allBin;bins]; %#ok<AGROW>
                allKeys=[allKeys;repmat(nativeFiles(fi),numel(retained),1)+"#"+string(retained)]; %#ok<AGROW>
            end
            [u,~,g]=unique(allBin);
            grouped=accumarray(g,(1:numel(g))',[],@(r){allKeys(r)});
            inventory=containers.Map('KeyType','double','ValueType','any');
            for gi=1:numel(u), inventory(u(gi))=grouped{gi}; end
            binInventoryCache(inventoryKey)=inventory;
        end
        inventory=binInventoryCache(inventoryKey); bin=posixtime(begin);
        keys=strings(0,1);
        if isKey(inventory,bin), keys=inventory(bin); end
    end
end

function tf=usable(T)
tf=false(height(T),1);
if ismember('PADUsable',T.Properties.VariableNames), tf=logical(T.PADUsable); end
end

function X=sectorArray(T,prefix,suffix)
X=nan(height(T),8);
for ii=1:8, X(:,ii)=T.(sprintf('%s_S%d_%s',prefix,ii,suffix)); end
end

function X=angleArray(T)
X=nan(height(T),8);
for ii=1:8, X(:,ii)=T.(sprintf('PA_S%d_deg',ii)); end
end

function out=cadenceName(cadence)
if strcmp(cadence,'day'), out='daily'; else, out='hourly'; end
end

function A=findFallbackAudit(S,T)
A=struct;
if isfield(S,'l1FallbackAudit'), A=S.l1FallbackAudit; end
if ~isfield(A,'AppliedBins') && isstruct(T.Properties.UserData) && ...
        isfield(T.Properties.UserData,'L1FallbackAudit')
    A=T.Properties.UserData.L1FallbackAudit;
end
end

function [same,delta,exact]=sameValue(a,b,tolerance)
exact=isequaln(a,b); same=exact; delta=0;
if exact, return, end
if ~isequal(size(a),size(b)), delta=Inf; return, end
if isdatetime(a)
    if any(isnat(a)~=isnat(b),'all'), delta=Inf; return, end
    d=abs(seconds(a-b));
elseif isnumeric(a)||islogical(a)
    if any(isnan(a)~=isnan(b),'all')||any(isinf(a)~=isinf(b),'all')
        delta=Inf; return
    end
    d=abs(double(a)-double(b));
else
    delta=Inf; return
end
v=d(isfinite(d));
if isempty(v), delta=Inf; else, delta=max(v); end
same=tolerance>0&&delta<=tolerance;
end
