function finalize_kh_results()
%FINALIZE_KH_RESULTS Build event inventory and consolidated SDC audit tables.

outputRoot = 'C:\Users\Administrator\Documents\KH';
events = kh_event_catalog();
status = readtable(fullfile(outputRoot,'run_status_all.csv'), ...
    'TextType','string','VariableNamingRule','preserve');

n = height(events);
figuresOK = zeros(n,1);
figuresPartial = zeros(n,1);
figuresNoData = zeros(n,1);
figuresFailed = zeros(n,1);
figuresUsable = zeros(n,1);
resultClass = strings(n,1);
figureDirectory = strings(n,1);

for k = 1:n
    rows = status(status.EventID==events.EventID(k),:);
    figuresOK(k) = sum(rows.Status=="ok");
    figuresPartial(k) = sum(rows.Status=="partial");
    figuresNoData(k) = sum(rows.Status=="no_data");
    figuresFailed(k) = sum(rows.Status=="failed");
    figuresUsable(k) = figuresOK(k)+figuresPartial(k);
    if figuresFailed(k) > 0
        resultClass(k) = "failed";
    elseif figuresNoData(k) == 4
        resultClass(k) = "source_absent_4sc";
    elseif figuresPartial(k) > 0
        resultClass(k) = "4sc_available_with_partial_panels";
    else
        resultClass(k) = "4sc_complete_9panel";
    end
    figureDirectory(k) = string(fullfile(outputRoot,'figures', ...
        char(events.EventID(k)+"_"+string(events.StartUTC(k),'yyyyMMdd_HHmmss'))));
end

inventory = events;
inventory.FiguresGenerated = repmat(4,n,1);
inventory.FiguresOK = figuresOK;
inventory.FiguresPartial = figuresPartial;
inventory.FiguresNoData = figuresNoData;
inventory.FiguresFailed = figuresFailed;
inventory.FiguresUsable = figuresUsable;
inventory.FourSpacecraftOverviewReady = figuresUsable==4;
inventory.ResultClass = resultClass;
inventory.FigureDirectory = figureDirectory;
writetable(inventory,fullfile(outputRoot,'MMS_KH_event_result_inventory.csv'));

audit = status(status.Status=="partial" | status.Status=="no_data",:);
audit.SDCPublicConclusion = repmat("target-window public source absent",height(audit),1);
audit.DownloadAction = repmat("not downloaded: no public target-window file available",height(audit),1);
audit.AuditMethod = repmat("MMS SDC public file_names API; priority audit and rounds 1-3",height(audit),1);
audit.AuditDateUTC = repmat("2026-07-18",height(audit),1);
writetable(audit,fullfile(outputRoot,'MMS_SDC_public_missing_audit_all.csv'));

fprintf('KH_FINALIZE inventory=%d audit_rows=%d usable=%d ok=%d partial=%d no_data=%d failed=%d\n', ...
    height(inventory),height(audit),sum(figuresUsable),sum(figuresOK), ...
    sum(figuresPartial),sum(figuresNoData),sum(figuresFailed));
end
