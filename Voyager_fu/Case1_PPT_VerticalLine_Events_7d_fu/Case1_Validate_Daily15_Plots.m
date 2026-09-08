function report = Case1_Validate_Daily15_Plots(result,baseline)
%Case1_Validate_Daily15_Plots Check daily-only +/-15-day replacement.
%   The caller captures BASELINE before rendering. Each Items.Saved contains
%   the complete former daily MAT. ProtectedFiles lists all non-daily PNGs
%   and hourly/peak/selected-time MATs that must remain byte-for-byte intact.
%   No production figure or science data are modified by this function.

%% explicit 47-event scope and check ledger
R = result.Run.ReportV1;
catalog = result.Catalog;
report = struct;
report.CreatedUTC = datetime('now','TimeZone','UTC');
report.Scope = "47 V1 daily overviews only; [D-15,D+16) UTC; existing hourly and PAD artifacts protected.";
report.Comparison = "Old daily payload matched by Epoch, product, source file and CDF record; paired NaN record numbers identify L1 means. SourceRow is intentionally excluded.";
report.AbsoluteTolerance = 1e-10;
report.Checks = table('Size',[0 4], ...
    'VariableTypes',{'string','string','logical','string'}, ...
    'VariableNames',{'EventID','Check','Passed','Detail'});
report.Events = table;
report.Artifacts = table;
report.ProtectedFiles = baseline.ProtectedFiles;
report.ProtectedFiles.CurrentSHA256 = strings(height(baseline.ProtectedFiles),1);
report.ProtectedFiles.Unchanged = false(height(baseline.ProtectedFiles),1);
report.Errors = strings(0,1);
currentID = "all";
check('exact_47_V1_catalog',height(catalog)==47 && ...
    numel(unique(string(catalog.EventID)))==47 && all(catalog.Spacecraft==1), ...
    sprintf('%d catalog rows',height(catalog)));
check('exact_47_reports',height(R)==47 && ...
    numel(unique(string(R.EventID)))==47 && ...
    isequal(sort(string(R.EventID)),sort(string(catalog.EventID))), ...
    sprintf('%d distinct report events',numel(unique(string(R.EventID)))));
check('all_reports_ok',all(string(R.Status)=="ok"),join(unique(string(R.Status)),', '));
oldIDs = reshape(string({baseline.Items.EventID}),[],1);
check('complete_old_daily_baseline',numel(baseline.Items)==47 && ...
    isequal(sort(oldIDs),sort(string(catalog.EventID))), ...
    sprintf('%d old daily audits',numel(baseline.Items)));

%% exact event windows, regenerated PNGs and unchanged overlapping values
for ie = 1:height(R)
    currentID = string(R.EventID(ie));
    firstCheck = height(report.Checks)+1;
    E = struct('EventID',currentID,'Rows',0,'L1Rows',0,'L2Rows',0, ...
        'PADUsableRows',0,'MissingFluxRows',0,'MissingPARows',0, ...
        'CompleteFluxMissingB',0,'OldRowsCompared',0,'Added2012Rows',0,'Passed',false,'Error',"");
    try
        id = find(string(catalog.EventID)==currentID);
        assert(isscalar(id),'Expected one catalog row per event.');
        D = dateshift(catalog.StartUTC(id),'start','day');
        begin = D-days(15); finish = D+days(16);
        check('exact_31_day_report_window',R.PlotStartUTC(ie)==begin && ...
            R.PlotEndUTCExclusive(ie)==finish && seconds(finish-begin)==31*86400, ...
            string(begin)+" / "+string(finish));
        matFile = string(R.PitchAngleAuditFile(ie));
        saved = load(char(matFile),'pitchAngleTable','opts','l1FallbackAudit');
        T = saved.pitchAngleTable; opts = saved.opts;
        check('daily_15day_L1_options',opts.ContextDays==15 && ...
            strcmp(opts.PADCadence,'day') && strcmp(opts.LECPSourcePriority,'l1_first') && ...
            ~opts.ExportPeakPAD,'Daily context=15, L1-first, peak renderer disabled.');
        check('method_audit_window',saved.l1FallbackAudit.StartUTC==begin && ...
            saved.l1FallbackAudit.EndUTC==finish && ...
            strcmp(saved.l1FallbackAudit.Cadence,'day') && ...
            strcmp(saved.l1FallbackAudit.SourcePriority,'l1_first'), ...
            'Saved source-selection audit uses the same complete UTC days.');
        check('all_epochs_inside_window',isempty(T) || ...
            all(T.EpochUTC>=begin & T.EpochUTC<finish),sprintf('%d rows',height(T)));

        pngFile = string(R.FigureFile(ie));
        info = imfinfo(char(pngFile)); pixels = imread(char(pngFile));
        check('PNG_decodes',strcmpi(info.Format,'png') && ...
            size(pixels,1)==info.Height && size(pixels,2)==info.Width && ...
            info.Width>0 && info.Height>0,pngFile);
        newHash = string(Case1_File_SHA256(char(pngFile)));
        oldPng = find(strcmpi(string(baseline.DailyFiles.SourceFile),pngFile));
        check('PNG_is_replacement',isscalar(oldPng) && ...
            all(newHash~=string(baseline.DailyFiles.SHA256(oldPng))), ...
            'Current PNG SHA256 differs from the captured former daily PNG.');
        artifact = table(currentID,pngFile,newHash,info.Width,info.Height,matFile, ...
            string(Case1_File_SHA256(char(matFile))), ...
            'VariableNames',{'EventID','PNGFile','PNGSHA256','Width','Height','MATFile','MATSHA256'});
        report.Artifacts = [report.Artifacts;artifact];

        oldItem = find(string({baseline.Items.EventID})==currentID);
        assert(isscalar(oldItem),'Expected one captured old daily audit.');
        old = baseline.Items(oldItem).Saved.pitchAngleTable;
        E.OldRowsCompared = height(old);
        compareOldRows(old,T);
        oldBegin = D-days(3); oldFinish = D+days(4);
        overlapCount = 0;
        if ~isempty(T), overlapCount = nnz(T.EpochUTC>=oldBegin & T.EpochUTC<oldFinish); end
        check('old_seven_days_same_row_count',overlapCount==height(old), ...
            sprintf('old=%d new-overlap=%d',height(old),overlapCount));

        E.Rows = height(T);
        if ~isempty(T)
            F = sectorValues(T,'RawFlux','1d');
            PA = sectorValues(T,'PA','deg');
            fluxOK = all(isfinite(F)&F>0,2);
            paOK = all(isfinite(PA),2);
            B = T{:,{'BR_daily_nT','BT_daily_nT','BN_daily_nT'}};
            bOK = all(isfinite(B),2)&vecnorm(B,2,2)>0;
            check('only_seven_sector_PAD_gate',isequal(logical(T.PADUsable),fluxOK&paOK), ...
                'Only finite positive S1-S7 flux and finite PA; uncertainty is not a gate.');
            check('source_priority_per_row',all(string(T.SourcePriority)=="l1_first"), ...
                'All saved records retain the approved priority.');
            check('known_source_products',all(ismember(string(T.SourceProduct),["L1_UTC_mean","L2"])), ...
                'No new source-product type.');
            E.L1Rows = nnz(string(T.SourceProduct)=="L1_UTC_mean");
            E.L2Rows = nnz(string(T.SourceProduct)=="L2");
            E.PADUsableRows = nnz(fluxOK&paOK);
            E.MissingFluxRows = nnz(~fluxOK); E.MissingPARows = nnz(~paOK);
            E.CompleteFluxMissingB = nnz(fluxOK&~bOK);
            if currentID == "Case1-S09-L01"
                added = T.EpochUTC < datetime(2013,1,1,'TimeZone','UTC');
                E.Added2012Rows = nnz(added);
                check('new_2012_window_uses_original_2012_CDF', ...
                    any(added) && all(contains(string(T.SourceCDF(added)),'20120101')), ...
                    sprintf('%d added-year rows trace to the downloaded 2012 source CDFs.', nnz(added)));
            end
        else
            check('only_seven_sector_PAD_gate',true,'No particle input rows in this window.');
        end
    catch ME
        E.Error = string(getReport(ME,'extended','hyperlinks','off'));
        report.Errors(end+1,1) = currentID+": "+string(ME.message);
        check('event_validation_completed',false,ME.message);
    end
    E.Passed = all(report.Checks.Passed(firstCheck:end));
    report.Events = [report.Events;struct2table(E)];
    fprintf('DAILY15_QA %d/%d %s rows=%d L1=%d L2=%d usable=%d old=%d passed=%d\n', ...
        ie,height(R),currentID,E.Rows,E.L1Rows,E.L2Rows,E.PADUsableRows,E.OldRowsCompared,E.Passed);
end

%% protected non-daily artifacts remain byte-for-byte unchanged
currentID = "protected";
for ip = 1:height(report.ProtectedFiles)
    filename = string(report.ProtectedFiles.SourceFile(ip));
    try
        report.ProtectedFiles.CurrentSHA256(ip) = string(Case1_File_SHA256(char(filename)));
        report.ProtectedFiles.Unchanged(ip) = ...
            report.ProtectedFiles.CurrentSHA256(ip)==string(report.ProtectedFiles.SHA256(ip));
        check('non_daily_artifact_unchanged',report.ProtectedFiles.Unchanged(ip),filename);
    catch ME
        check('non_daily_artifact_unchanged',false,filename+": "+string(ME.message));
    end
end
currentID = "all";
check('complete_new_artifact_inventory',height(report.Artifacts)==47 && ...
    numel(unique(report.Artifacts.PNGFile))==47 && ...
    numel(unique(report.Artifacts.MATFile))==47,'47 distinct daily PNG/MAT pairs.');
report.PassedCount = nnz(report.Checks.Passed);
report.CheckCount = height(report.Checks);
report.Passed = all(report.Checks.Passed) && height(report.Events)==47;
report.CompletedUTC = datetime('now','TimeZone','UTC');
report.CodeFile = string([mfilename('fullpath'),'.m']);
report.CodeSHA256 = string(Case1_File_SHA256(char(report.CodeFile)));
report.CountMeaning = "Rows are event-window appearances; overlapping windows are not deduplicated.";
folder = char(result.AuditFolder);
assert(startsWith(lower(strrep(folder,'\','/')),'z:/spart-work/data/voyager/'), ...
    'Save validation data beneath the Voyager data root.');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(folder,['validate_daily15_',stamp,'.mat']));
save(report.AuditFile,'report','-v7.3');
fprintf('DAILY15_VALIDATION passed=%d checks=%d/%d\n%s\n', ...
    report.Passed,report.PassedCount,report.CheckCount,report.AuditFile);
if ~report.Passed, disp(report.Checks(~report.Checks.Passed,:)); end

    function check(name,passed,detail)
        assert(isscalar(passed),'A check must return one logical value.');
        report.Checks(end+1,:) = {currentID,string(name),logical(passed),string(detail)};
    end

    function compareOldRows(old,new)
        names = {'BR_daily_nT','BT_daily_nT','BN_daily_nT','MAGVectorSampleCount','PADUsable'};
        for sector = 1:7
            for prefix = {'RawFlux','FluxUncertainty','RawRate','Samples'}
                names{end+1} = sprintf('%s_S%d_1d',prefix{1},sector); %#ok<AGROW>
            end
            names{end+1} = sprintf('PA_S%d_deg',sector); %#ok<AGROW>
            for component = 'RTN'
                names{end+1} = sprintf('ParticleU%c_S%d',component,sector); %#ok<AGROW>
            end
        end
        if isempty(old)
            check('old_daily_payload_unchanged',true,'Old table has no records.'); return
        end
        check('comparison_columns_present',all(ismember(names,old.Properties.VariableNames)) && ...
            all(ismember(names,new.Properties.VariableNames)),'All flux, uncertainty, rate, count, PA, direction and B columns.');
        map = nan(height(old),1);
        for row = 1:height(old)
            record = new.SourceCDFRecord==old.SourceCDFRecord(row) | ...
                (isnan(new.SourceCDFRecord)&isnan(old.SourceCDFRecord(row)));
            match = find(new.EpochUTC==old.EpochUTC(row) & ...
                string(new.SourceProduct)==string(old.SourceProduct(row)) & ...
                strcmpi(string(new.SourceCDF),string(old.SourceCDF(row))) & record);
            if isscalar(match), map(row)=match; end
        end
        check('every_old_record_matched_once',all(isfinite(map)) && ...
            numel(unique(map(isfinite(map))))==height(old), ...
            sprintf('%d/%d unique old records matched',nnz(isfinite(map)),height(old)));
        if any(~isfinite(map)), return, end
        for in = 1:numel(names)
            name = names{in};
            check('old_daily_payload_unchanged',closeValue(old.(name),new.(name)(map,:),1e-10),name);
        end
    end
end

function V = sectorValues(T,prefix,suffix)
V = nan(height(T),7);
for sector = 1:7
    V(:,sector) = T.(sprintf('%s_S%d_%s',prefix,sector,suffix));
end
end

function same = closeValue(a,b,tolerance)
same = isequal(size(a),size(b));
if ~same, return, end
same = isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b));
if ~same, return, end
infinite = isinf(a);
same = isequal(a(infinite),b(infinite));
finite = isfinite(a)&isfinite(b);
same = same && all(abs(double(a(finite))-double(b(finite)))<=tolerance);
end
