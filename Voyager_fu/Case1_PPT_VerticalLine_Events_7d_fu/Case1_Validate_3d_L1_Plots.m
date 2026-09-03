function result = Case1_Validate_3d_L1_Plots(baselineFile,cadenceFilter)
%Case1_Validate_3d_L1_Plots Independent requested-event L1-first regression.
%   Validate [D-3,D+4) UTC, both PAD cadences, original-CDF calculations,
%   three-dimensional pointing, and actual overview/peak PNG/MAT artifacts.
%   BASELINEFILE identifies the 45 requested events and existing audit paths;
%   no old L2-preservation condition is imposed under the new L1 priority.
%   Optional CADENCEFILTER is '1d' or '1h'; default checks all 90 audits.
%   No figures, production data, source CDFs or old audits are modified.

%% fixed scope and existing source reader
assert(nargin>=1 && isfile(baselineFile),'Supply the saved 90-audit baseline.');
if nargin<2, cadenceFilter=["1d","1h"]; end
cadenceFilter=unique(string(cadenceFilter(:).'),'stable');
assert(all(ismember(cadenceFilter,["1d","1h"])),'Use 1d and/or 1h.');
loaded=load(baselineFile,'baseline'); baseline=loaded.baseline;
items=baseline.Items(ismember(string({baseline.Items.Cadence}),cadenceFilter));
catalog=baseline.Catalog;
assert(height(catalog)==45 && all(catalog.Spacecraft==1), ...
    'This test is restricted to the 45 requested V1 events.');
assert(numel(items)==45*numel(unique(cadenceFilter)),'Incomplete baseline scope.');
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
result=struct('BaselineFile',string(baselineFile), ...
    'BaselineSHA256',string(Case1_File_SHA256(baselineFile)), ...
    'CadenceFilter',cadenceFilter,'ExpectedAudits',numel(items), ...
    'StartedUTC',datetime('now','TimeZone','UTC'),'Events',table, ...
    'Summary',table,'Artifacts',table,'PeakReferences',{{}}, ...
    'Scope','Only baseline-listed V1 events; [D-3,D+4) UTC; L1-first', ...
    'CountMeaning','Record counts are appearances per event window, not deduplicated across overlapping windows.');
result.Artifacts=table('Size',[0 7], ...
    'VariableTypes',{'string','string','string','string','string','double','double'}, ...
    'VariableNames',{'EventID','Cadence','Kind','File','SHA256','Width','Height'});
sourceCache=containers.Map('KeyType','char','ValueType','any');
poolCache=containers.Map('KeyType','char','ValueType','any');
magCache=containers.Map('KeyType','char','ValueType','any');
hashCache=containers.Map('KeyType','char','ValueType','char');
roleCache=containers.Map('KeyType','char','ValueType','char');
checkRows=cell(32768,5); checkCount=0; checkPassed=false(32768,1);
currentID=""; currentCadence="";
factor=1/(0.44*(1.78-0.57));

%% exact counterpart files; never scan unrelated event folders
for ik=1:numel(items)
    item=items(ik); currentID=string(item.EventID); currentCadence=string(item.Cadence);
    firstCheck=checkCount+1; eventError="";
    nRows=0; nL1=0; nL2=0; nUsable=0; fluxOnly=0; bOnly=0; both=0; other=0;
    attitudeMissing=0; maxUN=NaN; peakStatus="not_applicable"; peakCount=NaN;
    peakGap="not_applicable"; expectedL1Count=0; expectedL2Count=0;
    cadence='day'; suffix='1d'; width=days(1);
    if currentCadence=="1h", cadence='hour'; suffix='1h'; width=hours(1); end
    er=catalog(string(catalog.EventID)==currentID,:);
    assert(height(er)==1,'Ambiguous requested event identity.');
    D=dateshift(er.StartUTC,'start','day'); begin=D-days(3); finish=D+days(4);
    try
        check('audit_exists',isfile(item.AuditFile),string(item.AuditFile));
        saved=load(char(item.AuditFile)); T=saved.pitchAngleTable; opts=saved.opts;
        A=saved.l1FallbackAudit; nRows=height(T);
        check('context_and_priority',opts.ContextDays==3 && ...
            strcmp(opts.LECPSourcePriority,'l1_first') && ...
            strcmp(A.SourcePriority,'l1_first') && opts.LECPLevel1Fallback, ...
            'ContextDays=3 and explicit l1_first in production options and source audit.');
        check('exact_seven_day_window',A.StartUTC==begin && A.EndUTC==finish && ...
            seconds(finish-begin)==7*86400 && strcmp(string(A.Cadence),cadence), ...
            string(begin)+" / "+string(finish));
        check('approved_unchanged_scientific_options', ...
            strcmp(opts.AccumulationPolicy,'epoch_drop_negative_deltat') && ...
            strcmp(opts.LECPBackgroundMode,'none') && strcmp(opts.PitchAngleMethod,'predicted_ck') && ...
            opts.PredictedAttitudeApproved && opts.NominalLECPGeometryApproved && ...
            opts.HourlyAttitudeApproved,'No B/flux interpolation, background subtraction or new validity gate.');
        check('declared_legacy_L1_conversion',A.RateToFluxFactor==factor && ...
            isequal(A.ConversionEnergyRange_MeV,[0.57,1.78]) && ...
            A.ConversionIsUserApprovedLegacyL1Path, ...
            'Legacy conversion is explicit; official L2-equivalent calibration is not assumed.');
        verifyCode(saved.codeManifest);
        manifest=A.SourceManifest;
        check('source_manifest_propagated',isequaln(saved.sourceLECP,manifest),'Overview sources match method audit.');
        for im=1:height(manifest)
            check('CDF_source_hash',string(manifest.SHA256(im))==fileHash(manifest.SourceFile(im),'CDF'), ...
                string(manifest.SourceFile(im)));
        end
        rates=getPool(manifest,"native","rate");
        l2=getPool(manifest,string(cadence),"flux");
        rKeep=rates.EpochUTC>=begin & rates.EpochUTC<finish & ~(rates.DeltaT_s<0);
        lKeep=l2.EpochUTC>=begin & l2.EpochUTC<finish & ~(l2.DeltaT_s<0);
        check('L1_negative_DeltaT_selection',A.L1RecordSelection.RetainedRecords==nnz(rKeep) && ...
            A.L1RecordSelection.NegativeDeltaTRejected==nnz(rates.EpochUTC>=begin & ...
            rates.EpochUTC<finish & rates.DeltaT_s<0),'Original Epoch anchoring; negative DeltaT alone excluded.');
        check('L2_negative_DeltaT_selection',A.L2RecordSelection.RetainedRecords==nnz(lKeep) && ...
            A.L2RecordSelection.NegativeDeltaTRejected==nnz(l2.EpochUTC>=begin & ...
            l2.EpochUTC<finish & l2.DeltaT_s<0),'L2 remains at its original Epoch.');
        check('candidate_windows_in_scope',all(A.Candidates.BinStartUTC>=begin & ...
            A.Candidates.BinEndUTC<=finish),'No source mean outside the requested seven days.');
        if nRows>0
            check('all_table_epochs_in_scope',all(T.EpochUTC>=begin & T.EpochUTC<finish), ...
                sprintf('%d plotted input records',nRows));
            check('known_source_products',all(ismember(string(T.SourceProduct),["L1_UTC_mean","L2"])), ...
                'Only explicitly derived L1 means or original L2.');
            if ismember('SourcePriority',T.Properties.VariableNames)
                check('row_source_priority',all(string(T.SourcePriority)=="l1_first"),'Per-row priority.');
            end
            isL1=string(T.SourceProduct)=="L1_UTC_mean";
            nL1=nnz(isL1); nL2=nRows-nL1;
            rawFlux=sectorArray(T,'RawFlux',suffix);
            sigma=sectorArray(T,'FluxUncertainty',suffix);
            pa=sectorArray(T,'PA','deg');
            displayFlux=rawFlux; displayFlux(~isfinite(displayFlux)|displayFlux<=0)=NaN;
            displayFlux(:,8)=rawFlux(:,8);
            check('plotted_flux_matches_source_payload', ...
                isequaln(sectorArray(T,'Flux',suffix),displayFlux), ...
                'No extra scaling or background correction in the displayed sectors.');
            check('table_bin_boundaries',isequal(T.BinStartUTC,dateshift(T.EpochUTC,'start',cadence)) && ...
                isequal(T.BinEndUTC,T.BinStartUTC+width),'Half-open UTC bin boundaries.');
            check('source_display_energy_separated',all(T.P1DisplayEnergyLower_MeV==0.57 & ...
                T.P1DisplayEnergyUpper_MeV==1.78),'Display band only; source bounds verified per row.');
            check('no_background_correction',~any(T.BackgroundCorrectionApplied), 'S8 is diagnostic only.');
        else
            isL1=false(0,1); rawFlux=nan(0,8); sigma=rawFlux; pa=rawFlux;
        end

        %% independently rebuild every UTC bin, including bins with no output
        expectedReplaced=zeros(0,1);
        bins=(begin:width:finish-width).';
        for ib=1:numel(bins)
            b=bins(ib); stop=b+width;
            rr=find(rKeep & rates.EpochUTC>=b & rates.EpochUTC<stop);
            ll=find(lKeep & l2.EpochUTC>=b & l2.EpochUTC<stop);
            [meanRate,meanSigma,counts]=independentMean(rates.Value(rr,:),rates.Sigma(rr,:));
            complete=all(isfinite(meanRate(1:7))&meanRate(1:7)>0);
            out=zeros(0,1);
            if nRows>0, out=find(T.EpochUTC>=b & T.EpochUTC<stop); end
            ci=find(A.Candidates.BinStartUTC==b);
            if ~isempty(rr)||~isempty(ll)
                check('one_candidate_for_input_bin',isscalar(ci),string(b));
                assert(isscalar(ci),'Missing or duplicate source candidate bin.');
                check('candidate_mean_sigma_count_exact', ...
                    isequaln(A.Candidates.MeanRate(ci,:),meanRate) && ...
                    isequaln(A.Candidates.MeanRateSigma(ci,:),meanSigma) && ...
                    isequal(A.Candidates.SectorSampleCount(ci,:),counts),string(b));
                check('L1_first_decision_independent',logical(A.Candidates.Applied(ci))==complete,string(b));
                verifyRecords(A.Candidates.SourceRecords{ci},rates,rr,begin,finish);
            else
                check('no_invented_empty_candidate',isempty(ci),string(b));
            end
            if complete
                expectedL1Count=expectedL1Count+1;
                expectedReplaced=[expectedReplaced;ll]; %#ok<AGROW>
                check('complete_L1_wins_entire_bin',numel(out)==1 && all(isL1(out)),string(b));
                if numel(out)~=1||~all(isL1(out)), continue, end
                j=out(1);
                check('L1_center_and_record_semantics',T.EpochUTC(j)==b+width/2 && ...
                    isnan(T.SourceCDFRecord(j)) && isnan(T.DeltaT_s(j)) && ...
                    string(T.SourceCDF(j))==rates.SourceCDF(rr(1)),sprintf('row %d',j));
                check('L1_native_rate_mean_exact', ...
                    isequaln(sectorArray(T(j,:),'RawRate',suffix),meanRate),sprintf('row %d',j));
                check('L1_flux_sigma_count_exact', ...
                    isequaln(rawFlux(j,:),meanRate*factor) && ...
                    isequaln(sigma(j,:),meanSigma*factor) && ...
                    isequaln(sectorArray(T(j,:),'Samples',suffix),counts) && ...
                    T.SourceToDifferentialFluxFactor(j)==factor,sprintf('row %d',j));
                check('L1_source_trace_matches_candidate', ...
                    isequaln(T.L1SourceRecords{j},A.Candidates.SourceRecords{ci}),sprintf('row %d',j));
                band=unique(rates.EnergyRange(rr,:),'rows');
                check('L1_original_energy_bounds',size(band,1)==1 && ...
                    isequaln(T{j,{'P1EnergyLower_MeV','P1EnergyUpper_MeV'}},band),sprintf('row %d',j));
            else
                expectedL2Count=expectedL2Count+numel(ll);
                check('incomplete_L1_keeps_only_original_L2',numel(out)==numel(ll) && ...
                    ~any(isL1(out)),string(b));
                for il=ll(:).'
                    matches=zeros(0,1);
                    if nRows>0
                        matches=find(T.EpochUTC==l2.EpochUTC(il) & ...
                            string(T.SourceCDF)==l2.SourceCDF(il) & ...
                            T.SourceCDFRecord==l2.SourceCDFRecord(il));
                    end
                    check('L2_source_identity',isscalar(matches),string(l2.EpochUTC(il)));
                    if ~isscalar(matches), continue, end
                    j=matches(1);
                    check('L2_original_flux_sigma_DeltaT_exact', ...
                        isequaln(rawFlux(j,:),l2.Value(il,:)) && ...
                        isequaln(sigma(j,:),l2.Sigma(il,:)) && ...
                        isequaln(T.DeltaT_s(j),l2.DeltaT_s(il)) && ...
                        T.SourceToDifferentialFluxFactor(j)==1 && ...
                        isequaln(T{j,{'P1EnergyLower_MeV','P1EnergyUpper_MeV'}},l2.EnergyRange(il,:)), ...
                        sprintf('row %d',j));
                end
            end
        end
        check('complete_output_accounting',nL1==expectedL1Count && nL2==expectedL2Count && ...
            A.L1AddedRecords==nL1,sprintf('L1=%d L2=%d',nL1,nL2));
        verifyReplaced(A.ReplacedL2,l2,expectedReplaced);

        %% independent magnetic means and full three-component PA
        if nRows>0
            M=getMag(saved.sourceMAG);
            expectedB=nan(nRows,3); expectedSamples=zeros(nRows,1);
            for j=1:nRows
                b=dateshift(T.EpochUTC(j),'start',cadence);
                v=M.B(M.EpochUTC>=b & M.EpochUTC<b+width,:);
                v=v(all(isfinite(v),2),:);
                expectedSamples(j)=size(v,1);
                if ~isempty(v), expectedB(j,:)=mean(v,1); end
            end
            names={['BR_',longCadence(cadence),'_nT'],['BT_',longCadence(cadence),'_nT'], ...
                ['BN_',longCadence(cadence),'_nT']};
            actualB=T{:,names};
            check('same_UTC_bin_magnetic_mean_exact',isequaln(actualB,expectedB) && ...
                isequal(T.MAGVectorSampleCount,expectedSamples),'No nearest-B replacement or interpolation.');
            [expectedPA,expectedU,found]=verifyPointing(T,actualB);
            check('three_component_PA_matches',closeValue(pa,expectedPA,1e-10), ...
                'acos(dot(B,u_RTN)/|B|), all three components retained.');
            maxUN=max(abs(expectedU(:,:,3)),[],'all','omitnan');
            attitudeMissing=nnz(~found);
            completeFlux=all(isfinite(rawFlux(:,1:7))&rawFlux(:,1:7)>0,2);
            bGood=all(isfinite(actualB),2)&vecnorm(actualB,2,2)>0;
            completePA=all(isfinite(pa(:,1:7)),2);
            expectedUsable=completeFlux&completePA;
            check('only_seven_sector_PAD_gate',isequal(logical(T.PADUsable),expectedUsable), ...
                'No count, coverage, span, RMS or uncertainty thresholds.');
            nUsable=nnz(expectedUsable);
            bOnly=nnz(completeFlux&~bGood);
            fluxOnly=nnz(~completeFlux&bGood);
            both=nnz(~completeFlux&~bGood);
            other=nnz(completeFlux&bGood&~completePA);
        end

        %% actual artifacts; no renderer is called by this validator
        [~,stem]=fileparts(item.AuditFile);
        overviewPNG=fullfile(opts.OutputFolder,string(stem)+".png");
        verifyPNG(overviewPNG,"overview");
        if currentCadence=="1h"
            stamp=char(datetime(D,'Format','yyyyMMdd'));
            peakMAT=fullfile(opts.PeakPADDataFolder, ...
                ['V1_',char(currentID),'_',stamp,'_P1_peak5_hourly.mat']);
            check('peak_MAT_exists',isfile(peakMAT),string(peakMAT));
            p=load(peakMAT,'peakAudit'); actual=p.peakAudit;
            reference=independentPeak(T,begin,finish,suffix);
            check('peak_window_and_priority',actual.WindowStartUTC==begin && ...
                actual.WindowEndUTCExclusive==finish && ...
                strcmp(actual.L1FallbackAudit.SourcePriority,'l1_first'),'Peak search uses the same seven-day L1-first table.');
            fields=fieldnames(reference);
            for ip=1:numel(fields)
                field=fields{ip};
                check('independent_peak_payload',isfield(actual,field) && ...
                    isequaln(actual.(field),reference.(field)),field);
            end
            check('peak_full_source_context',isequaln(actual.SourceLECP,saved.sourceLECP) && ...
                isequaln(actual.PointingAudit,T.Properties.UserData),'Actual saved overview context.');
            code=table(string(actual.CodeFiles(:)),string(actual.CodeSHA256(:)), ...
                'VariableNames',{'SourceFile','SHA256'});
            check('peak_code_manifest_matches_overview',isequaln(code,saved.codeManifest), ...
                'Actual production code hashes shared by overview and peak.');
            verifyCode(code); verifyPNG(actual.FigureFile,"peak5");
            reference.PeakMAT=string(peakMAT);
            reference.PeakMATSHA256=fileHash(peakMAT,'audit');
            reference.EventID=currentID;
            result.PeakReferences{end+1,1}=reference;
            peakStatus=string(actual.Status); peakCount=actual.PlottedCount;
            peakGap="none";
            if strcmp(actual.Status,'peak_pad_unavailable')
                row=actual.PeakTableRow;
                f=all(isfinite(actual.RawFlux(3,:))&actual.RawFlux(3,:)>0);
                V=T{row,{['BR_',longCadence(cadence),'_nT'],['BT_',longCadence(cadence),'_nT'], ...
                    ['BN_',longCadence(cadence),'_nT']}};
                b=all(isfinite(V))&&norm(V)>0;
                if f&&~b, peakGap="complete_flux_missing_B";
                elseif ~f&&~b, peakGap="incomplete_flux_and_missing_B";
                elseif ~f&&b, peakGap="incomplete_flux_B_valid";
                else, peakGap="complete_flux_B_valid_PA_missing"; end
            end
        end
    catch ME
        eventError=string(getReport(ME,'extended','hyperlinks','off'));
        check('event_completed',false,string(ME.message));
    end
    passed=all(checkPassed(firstCheck:checkCount))&&strlength(eventError)==0;
    row=table(currentID,currentCadence,begin,finish,nRows,nL1,nL2,nUsable, ...
        bOnly,fluxOnly,both,other,attitudeMissing,maxUN,peakStatus,peakCount,peakGap,passed,eventError, ...
        'VariableNames',{'EventID','Cadence','WindowStartUTC','WindowEndUTCExclusive', ...
        'Rows','L1Rows','L2Rows','Usable','CompleteFluxMissingB','IncompleteFluxBValid', ...
        'IncompleteFluxMissingB','CompleteFluxBValidMissingPA','AttitudeMissing', ...
        'MaxAbsParticleUN','PeakStatus','PeakPlottedCount','PeakGap','Passed','Error'});
    result.Events=[result.Events;row]; %#ok<AGROW>
    fprintf('3D_L1_FIRST_VALIDATE %s %s rows=%d L1=%d L2=%d usable=%d passed=%d\n', ...
        currentID,currentCadence,nRows,nL1,nL2,nUsable,passed);
    if ~passed
        failed=find(~checkPassed(firstCheck:checkCount))+firstCheck-1;
        for jf=failed(:).'
            fprintf('  FAILED %s: %s\n',string(checkRows{jf,3}),string(checkRows{jf,5}));
        end
    end
end

%% classified audit and complete expected PNG inventory
result.Checks=cell2table(checkRows(1:checkCount,:), ...
    'VariableNames',{'EventID','Cadence','Check','Passed','Detail'});
result.Checks.Passed=checkPassed(1:checkCount);
result.SourceHashes=table(string(keys(hashCache))',string(values(hashCache))', ...
    string(values(roleCache,keys(hashCache)))','VariableNames',{'SourceFile','SHA256','Role'});
for c=unique(cadenceFilter)
    E=result.Events(result.Events.Cadence==c,:);
    row=table(c,height(E),nnz(E.Usable>0),sum(E.Rows),sum(E.L1Rows),sum(E.L2Rows), ...
        sum(E.Usable),sum(E.CompleteFluxMissingB),sum(E.IncompleteFluxBValid), ...
        sum(E.IncompleteFluxMissingB),sum(E.AttitudeMissing),nnz(E.Passed), ...
        'VariableNames',{'Cadence','Events','EventsWithPAD','Rows','L1Rows','L2Rows','Usable', ...
        'CompleteFluxMissingB','IncompleteFluxBValid','IncompleteFluxMissingB','AttitudeMissing','PassedEvents'});
    result.Summary=[result.Summary;row]; %#ok<AGROW>
end
result.ExpectedPNGCount=45*nnz(cadenceFilter=="1d")+90*nnz(cadenceFilter=="1h");
result.PNGInventoryComplete=height(result.Artifacts)==result.ExpectedPNGCount && ...
    numel(unique(result.Artifacts.File))==result.ExpectedPNGCount;
result.CompletedUTC=datetime('now','TimeZone','UTC');
result.Passed=all(result.Checks.Passed)&&all(result.Events.Passed)&& ...
    height(result.Events)==result.ExpectedAudits&&result.PNGInventoryComplete;
result.CodeFile=[mfilename('fullpath'),'.m'];
result.CodeSHA256=Case1_File_SHA256(result.CodeFile);
folder='Z:/SPART-WORK/Data/Voyager/voyager1/lecp/validation/context3d_l1_first';
if ~isfolder(folder), mkdir(folder); end
result.AuditFile=fullfile(folder,['validate_3d_l1_first_', ...
    char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS')),'.mat']);
save(result.AuditFile,'result','-v7.3');
disp(result.Summary);
if any(~result.Checks.Passed), disp(result.Checks(~result.Checks.Passed,:)); end
fprintf('3D_L1_FIRST_VALIDATION passed=%d checks=%d/%d PNG=%d/%d audit=%s\n', ...
    result.Passed,nnz(result.Checks.Passed),height(result.Checks),height(result.Artifacts), ...
    result.ExpectedPNGCount,result.AuditFile);

    function check(name,passed,detail)
        assert(isscalar(passed),'Independent check must produce one logical scalar.');
        checkCount=checkCount+1;
        if checkCount>size(checkRows,1)
            capacity=2*size(checkRows,1); checkRows{capacity,5}=[]; checkPassed(capacity)=false;
        end
        checkRows(checkCount,:)={currentID,currentCadence,string(name),logical(passed),string(detail)};
        checkPassed(checkCount)=logical(passed);
    end

    function hash=fileHash(filename,role)
        filename=char(filename);
        if ~isKey(hashCache,filename)
            hashCache(filename)=Case1_File_SHA256(filename); roleCache(filename)=char(role);
        end
        hash=string(hashCache(filename));
    end

    function P=source(filename,profile)
        filename=char(filename);
        if ~isKey(sourceCache,filename)
            sourceCache(filename)=Voyager_Read_CDF_Product(filename,profile);
            fileHash(filename,'CDF');
        end
        P=sourceCache(filename);
    end

    function pool=getPool(manifest,cadence,kind)
        files=string(manifest.SourceFile(string(manifest.Cadence)==cadence));
        assert(~isempty(files),'Expected source product is absent from the manifest.');
        key=char(kind+"|"+join(files,"|"));
        if isKey(poolCache,key), pool=poolCache(key); return, end
        parts=cell(numel(files),1);
        if kind=="rate", valueField='FHDU_SectoredRates'; sigmaField='FHDU_SectoredRateUncertainties';
        else, valueField='FHDU_SectoredFluxes'; sigmaField='FHDU_SectoredFluxUncertainties'; end
        for fi=1:numel(files)
            P=source(files(fi),'lecp_sector_daily'); n=numel(P.Epoch);
            labels=strtrim(string(P.Hydrogen_Channels_Label));
            check('P1_channel_and_sector_numbering',isequal(find(labels=="P1"),10) && ...
                isequal(P.SectorIterator(:)',1:8),files(fi));
            parts{fi}=table(repmat(files(fi),n,1),(1:n)',P.Epoch(:),P.DeltaT(:), ...
                reshape(P.(valueField)(:,10,:),n,8),reshape(P.(sigmaField)(:,10,:),n,8), ...
                reshape(P.FHDU_EnergyRange(:,:,10),n,2),'VariableNames', ...
                {'SourceCDF','SourceCDFRecord','EpochUTC','DeltaT_s','Value','Sigma','EnergyRange'});
        end
        allRows=vertcat(parts{:}); allRows=sortrows(allRows,'EpochUTC');
        duplicate=find(diff(allRows.EpochUTC)==seconds(0));
        for di=duplicate(:)'
            check('duplicate_native_payload_agrees', ...
                isequaln(allRows{di,{'DeltaT_s','Value','Sigma','EnergyRange'}}, ...
                allRows{di+1,{'DeltaT_s','Value','Sigma','EnergyRange'}}),string(allRows.EpochUTC(di)));
        end
        [~,keep,group]=unique(allRows.EpochUTC,'stable');
        alias=cell(numel(keep),1);
        for gi=unique(group(duplicate))'
            at=find(group==gi);
            alias{gi}=allRows.SourceCDF(at)+"#"+string(allRows.SourceCDFRecord(at));
        end
        pool=allRows(keep,:); pool.AliasKeys=alias;
        poolCache(key)=pool;
    end

    function verifyRecords(records,pool,rows,windowStart,windowStop)
        check('source_record_count',height(records)==numel(rows),'One contribution per distinct original Epoch.');
        if height(records)~=numel(rows), return, end
        check('native_record_payload_and_scope', ...
            isequal(string(records.SourceCDF),pool.SourceCDF(rows)) && ...
            isequal(records.SourceCDFRecord,pool.SourceCDFRecord(rows)) && ...
            isequal(records.EpochUTC,pool.EpochUTC(rows)) && ...
            isequaln(records.DeltaT_s,pool.DeltaT_s(rows)) && ...
            isequaln(records.Rate_S1_to_S8,pool.Value(rows,:)) && ...
            isequaln(records.Sigma_S1_to_S8,pool.Sigma(rows,:)) && ...
            all(records.EpochUTC>=windowStart & records.EpochUTC<windowStop) && ...
            ~any(records.DeltaT_s<0),'Every original L1 rate, sigma, Epoch and DeltaT checked.');
        check('per_sector_contribution_mask', ...
            isequal(records.Contributes_S1_to_S8,isfinite(pool.Value(rows,:))&pool.Value(rows,:)>=0), ...
            'Zero contributes; NaN and negative values do not.');
        for ai=1:numel(rows)
            expected=pool.AliasKeys{rows(ai)};
            if isempty(expected), continue, end
            actual=records.IdenticalRecordAliases{ai};
            check('duplicate_aliases_retained',istable(actual) && ...
                isequal(sort(string(actual.SourceCDF)+"#"+string(actual.SourceCDFRecord)),sort(expected)), ...
                'Aliases retained but counted once.');
        end
    end

    function verifyReplaced(replaced,pool,rows)
        check('replaced_L2_count',numel(replaced.SourceRecordNumber)==numel(rows), ...
            'Complete L1 bins replace all L2 records, including complete L2.');
        if numel(replaced.SourceRecordNumber)~=numel(rows), return, end
        expected=pool(rows,:);
        check('replaced_L2_original_payload_retained', ...
            isequal(string(replaced.SourceCDF),expected.SourceCDF) && ...
            isequal(replaced.SourceRecordNumber,expected.SourceCDFRecord) && ...
            isequal(replaced.Epoch,expected.EpochUTC) && ...
            isequaln(replaced.DeltaT,expected.DeltaT_s) && ...
            isequaln(reshape(replaced.FHDU_SectoredFluxes(:,10,:),numel(rows),8),expected.Value) && ...
            isequaln(reshape(replaced.FHDU_SectoredFluxUncertainties(:,10,:),numel(rows),8),expected.Sigma), ...
            'Displaced L2 remains auditable and unchanged.');
    end

    function M=getMag(files)
        files=string(files(:));
        if isempty(files)
            M=table(NaT(0,1,'TimeZone','UTC'),nan(0,3),'VariableNames',{'EpochUTC','B'});
            return
        end
        key=char(join(files,"|"));
        if isKey(magCache,key), M=magCache(key); return, end
        parts=cell(numel(files),1);
        for mi=1:numel(files)
            P=source(files(mi),'coho');
            parts{mi}=table(P.Epoch(:),[P.BR(:),P.BT(:),P.BN(:)], ...
                'VariableNames',{'EpochUTC','B'});
        end
        M=sortrows(vertcat(parts{:}),'EpochUTC'); [~,keep]=unique(M.EpochUTC,'stable');
        M=M(keep,:); magCache(key)=M;
    end

    function [expectedPA,U,found]=verifyPointing(T,B)
        P=T.Properties.UserData; G=P.Geometry; A=P.Attitude; n=height(T);
        axis=[-sind(15);cosd(15);0]; hga=[0;0;-1]; theta=22.5+(0:7)*45;
        look=hga*cosd(theta)+cross(axis,hga)*sind(theta); particle=-look;
        check('official_nominal_geometry',isequal(G.Sector,1:8) && isequal(G.ActiveSectors,1:7) && ...
            G.BlockedSector==8 && closeValue(G.AxisSC,axis,1e-14) && ...
            closeValue(G.CenterFromHGA_deg,theta,1e-14) && ...
            closeValue(G.LookSC,look,1e-14) && closeValue(G.ParticleSC,particle,1e-14), ...
            'Nominal installation axis, numbered centers and v=-look.');
        check('geometry_official_sources_recorded',numel(G.SourceURL)>=3 && ...
            any(contains(string(G.SourceURL),'voyager.ftecs.com')) && ...
            any(contains(string(G.SourceURL),'naif.jpl.nasa.gov')) && ...
            any(contains(string(G.SourceURL),'pds-ppi.igpp.ucla.edu')),string(G.Model));
        check('attitude_evaluated_at_product_Epoch',isequal(P.TimeUTC,T.EpochUTC) && ...
            isequal(A.TimeUTC,T.EpochUTC) && size(A.C_SC_to_RTN,3)==n, ...
            'L1 centers and original L2 Epochs, never a fixed RT geometry.');
        for si=1:height(A.Sources)
            check('attitude_kernel_hash',string(A.Sources.SHA256(si))== ...
                fileHash(A.Sources.SourceFile(si),'attitude'),string(A.Sources.SourceFile(si)));
        end
        found=logical(A.Found(:)); U=nan(n,8,3); expectedPA=nan(n,8); expectedMu=nan(n,8);
        for jt=1:n
            if ~found(jt), continue, end
            C=A.C_SC_to_RTN(:,:,jt);
            check('full_right_handed_rotation',norm(C*C'-eye(3),'fro')<1e-10 && ...
                abs(det(C)-1)<1e-10 && ...
                closeValue(C,A.C_J2000_to_RTN(:,:,jt)*A.C_SC_to_J2000(:,:,jt),1e-12),sprintf('row %d',jt));
            u=C*particle; U(jt,:,:)=reshape(u',1,8,3);
            if all(isfinite(B(jt,:)))&&norm(B(jt,:))>0
                expectedMu(jt,:)=max(-1,min(1,(B(jt,:)/norm(B(jt,:)))*u));
                expectedPA(jt,:)=acosd(expectedMu(jt,:));
            end
        end
        actualU=nan(n,8,3);
        for si=1:8
            actualU(:,si,1)=T.(sprintf('ParticleUR_S%d',si));
            actualU(:,si,2)=T.(sprintf('ParticleUT_S%d',si));
            actualU(:,si,3)=T.(sprintf('ParticleUN_S%d',si));
            check('pitch_cosine_full_3D',closeValue(T.(sprintf('Mu_S%d',si)),expectedMu(:,si),1e-12),sprintf('S%d',si));
        end
        check('particle_RTN_from_full_attitude',closeValue(actualU,U,1e-12) && ...
            closeValue(P.ParticleRTN,U,1e-12) && closeValue(P.LookRTN,-U,1e-12), ...
            'All R/T/N components from C times the documented body geometry.');
        if any(found)
            check('N_component_not_forced_zero',any(abs(U(found,1:7,3))>1e-8,'all'), ...
                'Computed central directions contain nonzero N components.');
        end
    end

    function verifyCode(code)
        required=["Voyager_Case1_Plot_Events.m","Case1_Config.m","Case1_Read_LECP_Rates.m", ...
            "Case1_Apply_L1_Fallback.m","Case1_Select_LECP_Epoch.m", ...
            "Case1_LECP_Geometry.m","Case1_Read_Predicted_Attitude.m", ...
            "Case1_Predicted_LECP_Pointing.m","Case1_Plot_Peak_PAD.m"];
        files=string(code.SourceFile);
        check('core_code_manifest_complete',all(arrayfun(@(s)any(endsWith(files,s)),required)), ...
            'Core selection, CDF, geometry, attitude and rendering files recorded.');
        for ic=1:height(code)
            check('production_code_hash',string(code.SHA256(ic))==fileHash(files(ic),'code'),files(ic));
        end
    end

    function verifyPNG(filename,kind)
        filename=string(filename); check('PNG_exists',isfile(filename),filename);
        info=imfinfo(filename);
        check('PNG_decodes',strcmpi(info.Format,'png')&&info.Width>0&&info.Height>0,filename);
        row=table(currentID,currentCadence,string(kind),filename,fileHash(filename,'PNG'), ...
            info.Width,info.Height,'VariableNames',{'EventID','Cadence','Kind','File','SHA256','Width','Height'});
        result.Artifacts=[result.Artifacts;row];
    end
end

function [rate,sigma,count]=independentMean(raw,uncertainty)
rate=nan(1,8); sigma=rate; count=zeros(1,8);
for si=1:8
    valid=isfinite(raw(:,si))&raw(:,si)>=0; count(si)=nnz(valid);
    if count(si)==0, continue, end
    rate(si)=mean(raw(valid,si)); q=uncertainty(valid,si);
    if all(isfinite(q)&q>=0), sigma(si)=sqrt(sum(q.^2))/count(si); end
end
end

function value=sectorArray(T,prefix,suffix)
value=nan(height(T),8);
for si=1:8, value(:,si)=T.(sprintf('%s_S%d_%s',prefix,si,suffix)); end
end

function name=longCadence(cadence)
if strcmp(cadence,'day'), name='daily'; else, name='hourly'; end
end

function same=closeValue(a,b,tolerance)
if ~isequal(size(a),size(b))||any(isnan(a)~=isnan(b),'all')||any(isinf(a)~=isinf(b),'all')
    same=false; return
end
d=abs(a-b); same=all(d(isfinite(d))<=tolerance);
end

function R=independentPeak(T,begin,finish,suffix)
n=height(T); selected=nan(5,1); peak=NaN; peakFlux=NaN;
peakTime=NaT(1,1,'TimeZone','UTC'); candidate=zeros(0,1); rowMaximum=nan(n,1);
flux=nan(n,7); sigma=flux; pa=flux; complete=false(n,1);
if n>0
    flux=sectorArray(T,'RawFlux',suffix); flux=flux(:,1:7);
    sigma=sectorArray(T,'FluxUncertainty',suffix); sigma=sigma(:,1:7);
    pa=sectorArray(T,'PA','deg'); pa=pa(:,1:7);
    for j=1:n
        values=flux(j,isfinite(flux(j,:))&flux(j,:)>0);
        if ~isempty(values), rowMaximum(j)=max(values); end
    end
    candidate=find(isfinite(rowMaximum)&T.EpochUTC>=begin&T.EpochUTC<finish);
    complete=logical(T.PADUsable)&all(isfinite(flux)&flux>0,2)&all(isfinite(pa),2);
    if ~isempty(candidate)
        peakFlux=max(rowMaximum(candidate)); ties=candidate(rowMaximum(candidate)==peakFlux);
        [~,first]=min(T.EpochUTC(ties)); peak=ties(first); peakTime=T.EpochUTC(peak); selected(3)=peak;
        eligible=find(complete&T.EpochUTC>=begin&T.EpochUTC<finish);
        [~,order]=sort(T.EpochUTC(eligible)); eligible=eligible(order);
        before=eligible(T.EpochUTC(eligible)<peakTime); after=eligible(T.EpochUTC(eligible)>peakTime);
        before=before(max(1,end-1):end); after=after(1:min(2,numel(after)));
        selected(3-numel(before):2)=before; selected(4:3+numel(after))=after;
    end
end
times=NaT(5,1,'TimeZone','UTC'); valid=false(5,1); f=nan(5,7); u=f; a=f; scale=nan(5,1);
normalized=f; normalizedSigma=f; rows=cell(5,1); normalizing=cell(5,1); band=nan(5,2);
for ip=1:5
    j=selected(ip); if ~isfinite(j), continue, end
    times(ip)=T.EpochUTC(j); valid(ip)=complete(j); f(ip,:)=flux(j,:); u(ip,:)=sigma(j,:); a(ip,:)=pa(j,:);
    scale(ip)=rowMaximum(j); normalized(ip,:)=f(ip,:)/scale(ip); normalizedSigma(ip,:)=u(ip,:)/scale(ip);
    normalizing{ip}=find(f(ip,:)==scale(ip)); rows{ip}=T(j,:); rows{ip}.Properties.UserData=[];
    band(ip,:)=T{j,{'P1EnergyLower_MeV','P1EnergyUpper_MeV'}};
end
status="no_retained_records";
if n>0, status="no_positive_sector_flux"; end
peakUsable=isfinite(peak)&&complete(peak);
if isfinite(peak)
    if nnz(valid)==5, status="complete";
    elseif ~peakUsable, status="peak_pad_unavailable";
    else, status="partial_neighbors"; end
end
R=struct('Status',status,'PeakTableRow',peak,'PeakEpochUTC',peakTime,'PeakFlux',peakFlux, ...
    'PeakPADUsable',peakUsable,'SelectedTableRows',selected,'SelectedCount',nnz(isfinite(selected)), ...
    'PlottedCount',nnz(valid),'SelectedEpochUTC',times,'SelectedPADUsable',valid, ...
    'SelectedRows',{rows},'RawFlux',f,'RawSigma',u,'PA_deg',a,'NormalizationFlux',scale, ...
    'NormalizingSectors',{normalizing},'NormalizedFlux',normalized,'NormalizedSigma',normalizedSigma, ...
    'SourceEnergyMeV',band,'CandidateTableRows',candidate,'RowMaximumSectorFlux',rowMaximum, ...
    'Sectors',1:7,'PanelOffsets',(-2:2)');
end
