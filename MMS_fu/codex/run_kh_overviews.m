function summary = run_kh_overviews(mode,maxEvents)
%RUN_KH_OVERVIEWS Batch four-spacecraft MMS KH overview production.
%
%   run_kh_overviews('smoke',1)    first priority event only
%   run_kh_overviews('priority')   all published rolled-up cases
%   run_kh_overviews('all')        complete de-duplicated catalog

if nargin < 1 || isempty(mode), mode = 'priority'; end
if nargin < 2 || isempty(maxEvents), maxEvents = inf; end

dataRoot = 'Z:\SPART-WORK\Data\MMS';
outputRoot = 'C:\Users\Administrator\Documents\KH';
if ~exist(outputRoot,'dir'), mkdir(outputRoot); end
if ~exist(fullfile(outputRoot,'figures'),'dir'), mkdir(fullfile(outputRoot,'figures')); end

events = kh_event_catalog();
writetable(events,fullfile(outputRoot,'MMS_KH_published_event_catalog.csv'));

switch lower(string(mode))
    case {"smoke","priority"}
        selected = events(startsWith(events.Tier,"A_"),:);
    case "all"
        selected = [events(startsWith(events.Tier,"A_"),:); events(~startsWith(events.Tier,"A_"),:)];
    otherwise
        error('Unknown mode: %s',mode);
end
if isfinite(maxEvents)
    selected = selected(1:min(height(selected),maxEvents),:);
end

mms.db_init('local_file_db',[dataRoot filesep]);
plotCode = dir(which('plot_mms_kh_overview'));

logPath = fullfile(outputRoot,sprintf('MATLAB_run_%s.log',lower(char(string(mode)))));
diary(logPath);
diary on;
cleanupDiary = onCleanup(@() diary('off'));

status = emptyStatusTable();
windows = table;
statusPath = fullfile(outputRoot,sprintf('run_status_%s.csv',lower(char(string(mode)))));
windowPath = fullfile(outputRoot,sprintf('selected_burst_windows_%s.csv',lower(char(string(mode)))));
progressPath = fullfile(outputRoot,sprintf('progress_%s.txt',lower(char(string(mode)))));

fprintf('KH_BATCH_START mode=%s events=%d figures_requested=%d\n',char(string(mode)),height(selected),4*height(selected));
for k = 1:height(selected)
    ev = selected(k,:);
    fprintf('KH_EVENT_START %d/%d %s %s to %s\n',k,height(selected),ev.EventID, ...
        char(string(ev.StartUTC)),char(string(ev.EndUTC)));
    try
        win = kh_find_burst_window(ev,dataRoot,20);
        windows = appendTable(windows,win);
        writetable(windows,windowPath);
    catch ME
        fprintf(2,'KH_WINDOW_ERROR %s %s\n',ev.EventID,ME.message);
        win = fallbackWindow(ev);
        windows = appendTable(windows,win);
        writetable(windows,windowPath);
    end

    for sc = 1:4
        eventDir = fullfile(outputRoot,'figures',char(ev.EventID + "_" + string(ev.StartUTC,'yyyyMMdd_HHmmss')));
        expectedName = sprintf('%s_%s_MMS%d.png',ev.EventID, ...
            char(string(win.PlotStartUTC,'yyyyMMdd_HHmmss')),sc);
        expectedPath = fullfile(eventDir,expectedName);
        existing = dir(expectedPath);
        isCurrent = ~isempty(existing) && existing(1).bytes > 50000 && ...
            ~isempty(plotCode) && existing(1).datenum >= plotCode(1).datenum;
        % Smoke mode deliberately redraws its four figures so code changes
        % are visually verified instead of accepting stale test images.
        if isCurrent && lower(string(mode)) ~= "smoke"
            row = makeStatus(ev.EventID,sc,win.PlotStartUTC,win.PlotEndUTC, ...
                string(expectedPath),"existing","current code/window output reused");
            fprintf('KH_FIGURE_EXISTING %s MMS%d\n',ev.EventID,sc);
        else
            try
                row = plot_mms_kh_overview(ev,win,sc,outputRoot);
                fprintf('KH_FIGURE_OK %s MMS%d %s\n',ev.EventID,sc,row.OutputFile);
            catch ME
                close all force;
                msg = string(getReport(ME,'extended','hyperlinks','off'));
                row = makeStatus(ev.EventID,sc,win.PlotStartUTC,win.PlotEndUTC,"","failed",msg);
                fprintf(2,'KH_FIGURE_FAILED %s MMS%d %s\n',ev.EventID,sc,ME.message);
            end
        end
        status = [status; row]; %#ok<AGROW>
        writetable(status,statusPath);
        writeProgress(progressPath,mode,k,height(selected),ev.EventID,sc,row.Status);
        drawnow;
    end
    fprintf('KH_EVENT_DONE %s\n',ev.EventID);
end

complete = sum(status.Status=="ok");
partial = sum(status.Status=="partial");
noData = sum(status.Status=="no_data");
reused = sum(status.Status=="existing");
ok = complete + partial + reused;
failed = sum(status.Status=="failed");
summary = table(string(mode),height(selected),4*height(selected),ok,complete, ...
    partial,noData,reused,failed, ...
    'VariableNames',{'Mode','Events','FiguresRequested','FiguresAvailable', ...
    'FiguresComplete','FiguresPartial','FiguresNoData','FiguresReused','FiguresFailed'});
writetable(summary,fullfile(outputRoot,sprintf('summary_%s.csv',lower(char(string(mode))))));
fprintf('KH_BATCH_DONE mode=%s available=%d failed=%d\n',char(string(mode)),ok,failed);
clear cleanupDiary
diary off;
end

function t = emptyStatusTable()
t = table(strings(0,1),zeros(0,1),NaT(0,1,'TimeZone','UTC'),NaT(0,1,'TimeZone','UTC'), ...
    strings(0,1),strings(0,1),strings(0,1), ...
    'VariableNames',{'EventID','Spacecraft','PlotStartUTC','PlotEndUTC','OutputFile','Status','Warnings'});
end

function row = makeStatus(id,sc,t0,t1,file,status,warnings)
row = table(string(id),sc,t0,t1,string(file),string(status),string(warnings), ...
    'VariableNames',{'EventID','Spacecraft','PlotStartUTC','PlotEndUTC','OutputFile','Status','Warnings'});
end

function out = appendTable(out,row)
if isempty(out), out = row; else, out = [out; row]; end
end

function win = fallbackWindow(ev)
t0 = ev.PreferredStartUTC;
t1 = ev.PreferredEndUTC;
if isnat(t0) || isnat(t1)
    t0 = ev.StartUTC;
    t1 = min(ev.EndUTC,t0+minutes(20));
end
win = table(string(ev.EventID),ev.StartUTC,ev.EndUTC,t0,t1,"fallback",0,"", ...
    'VariableNames',{'EventID','EventStartUTC','EventEndUTC','PlotStartUTC','PlotEndUTC','Mode','CoverageOf16','EmptyProducts'});
end

function writeProgress(path,mode,k,n,id,sc,status)
fid = fopen(path,'w');
if fid < 0, return; end
c = onCleanup(@() fclose(fid));
fprintf(fid,'mode=%s\nevent=%d/%d\nid=%s\nspacecraft=%d\nstatus=%s\nupdated=%s\n', ...
    char(string(mode)),k,n,char(string(id)),sc,char(string(status)),char(datetime('now','TimeZone','UTC')));
end
