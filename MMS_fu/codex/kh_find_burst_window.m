function info = kh_find_burst_window(eventRow,dataRoot,maxMinutes)
%KH_FIND_BURST_WINDOW Select a representative local four-spacecraft window.
%
% The selector indexes only the four products needed by the overview:
% FGM brst L2, EDP brst L2 dce, FPI brst des-moms and dis-moms.  It first
% seeks file start times common to all 16 spacecraft/product combinations.
% If strict intersection is empty, it chooses the best-covered timestamps
% and reports that fact in INFO.Mode.  No CDF data are read here.

if nargin < 3 || isempty(maxMinutes), maxMinutes = 20; end
if nargin < 2 || isempty(dataRoot), dataRoot = 'Z:\SPART-WORK\Data\MMS'; end

eventStart = eventRow.StartUTC;
eventEnd = eventRow.EndUTC;
preferredStart = eventRow.PreferredStartUTC;
preferredEnd = eventRow.PreferredEndUTC;
if isnat(preferredStart) || isnat(preferredEnd)
    targetStart = eventStart;
    targetEnd = eventEnd;
else
    targetStart = max(eventStart,preferredStart);
    targetEnd = min(eventEnd,preferredEnd);
end
targetCenter = targetStart + (targetEnd-targetStart)/2;

products = {
    fullfile('fgm','brst','l2')
    fullfile('edp','brst','l2','dce')
    fullfile('fpi','brst','l2','des-moms')
    fullfile('fpi','brst','l2','dis-moms')};

sets = cell(16,1);
labels = strings(16,1);
q = 0;
for sc = 1:4
    for p = 1:numel(products)
        q = q + 1;
        labels(q) = sprintf('mms%d/%s',sc,products{p});
        sets{q} = localFileStarts(dataRoot,sc,products{p},eventStart-minutes(3),eventEnd);
    end
end

strict = sets{1};
for q = 2:numel(sets)
    strict = intersect(strict,sets{q});
end

if ~isempty(strict)
    candidate = strict;
    candidateCount = 16*ones(size(candidate));
    mode = "strict_4sc_all_products";
else
    unionTimes = NaT(0,1,'TimeZone','UTC');
    for q = 1:numel(sets), unionTimes = union(unionTimes,sets{q}); end
    candidate = unionTimes;
    candidateCount = zeros(size(candidate));
    for j = 1:numel(candidate)
        for q = 1:numel(sets)
            candidateCount(j) = candidateCount(j) + any(sets{q}==candidate(j));
        end
    end
    if ~isempty(candidateCount)
        bestCount = max(candidateCount);
        candidate = candidate(candidateCount==bestCount);
        candidateCount = bestCount*ones(size(candidate));
        mode = "best_effort_" + bestCount + "of16";
    else
        mode = "no_local_file_index";
    end
end

candidate = sort(candidate);
if isempty(candidate)
    plotStart = targetStart;
    plotEnd = min(targetEnd,plotStart+minutes(maxMinutes));
    coverage = 0;
else
    % Split timestamps into runs; MMS burst files usually start 110-120 s
    % apart.  A 180 s threshold tolerates short irregular segments.
    breaks = [true; seconds(diff(candidate)) > 180];
    runId = cumsum(breaks);
    ids = unique(runId);
    score = -inf(size(ids));
    runStart = NaT(size(ids),'TimeZone','UTC');
    runEnd = NaT(size(ids),'TimeZone','UTC');
    for j = 1:numel(ids)
        tt = candidate(runId==ids(j));
        runStart(j) = max(eventStart,tt(1));
        runEnd(j) = min(eventEnd,tt(end)+seconds(120));
        overlap = max(0,seconds(min(runEnd(j),targetEnd)-max(runStart(j),targetStart)));
        centerDistance = abs(seconds((runStart(j)+(runEnd(j)-runStart(j))/2)-targetCenter));
        durationScore = seconds(runEnd(j)-runStart(j));
        score(j) = 1e6*double(overlap>0) + 100*overlap + durationScore - centerDistance;
    end
    [~,best] = max(score);
    availableStart = runStart(best);
    availableEnd = runEnd(best);

    % If a publication supplied a preferred interval, preserve as much of
    % it as the local run permits.  Otherwise center on the event midpoint.
    desiredStart = max(availableStart,targetStart);
    desiredEnd = min(availableEnd,targetEnd);
    if desiredEnd <= desiredStart
        desiredCenter = min(max(targetCenter,availableStart),availableEnd);
        desiredStart = max(availableStart,desiredCenter-minutes(maxMinutes/2));
        desiredEnd = min(availableEnd,desiredStart+minutes(maxMinutes));
    end
    if minutes(desiredEnd-desiredStart) > maxMinutes
        desiredCenter = min(max(targetCenter,desiredStart),desiredEnd);
        plotStart = max(desiredStart,desiredCenter-minutes(maxMinutes/2));
        plotEnd = plotStart+minutes(maxMinutes);
        if plotEnd > desiredEnd
            plotEnd = desiredEnd;
            plotStart = plotEnd-minutes(maxMinutes);
        end
    else
        plotStart = desiredStart;
        plotEnd = desiredEnd;
    end
    coverage = max(candidateCount);
end

if plotEnd <= plotStart
    plotStart = eventStart;
    plotEnd = min(eventEnd,eventStart+minutes(maxMinutes));
end
if plotEnd <= plotStart
    error('KH:InvalidWindow','No positive-duration plot window for %s.',eventRow.EventID);
end

info = table(string(eventRow.EventID),eventStart,eventEnd,plotStart,plotEnd, ...
    string(mode),coverage,strjoin(labels(cellfun(@isempty,sets)),','), ...
    'VariableNames',{'EventID','EventStartUTC','EventEndUTC','PlotStartUTC','PlotEndUTC','Mode','CoverageOf16','EmptyProducts'});
end

function times = localFileStarts(dataRoot,sc,product,eventStart,eventEnd)
persistent cache
if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end

days = dateshift(eventStart,'start','day'):caldays(1):dateshift(eventEnd,'start','day');
times = NaT(0,1,'TimeZone','UTC');
for day = days
    dayDir = fullfile(dataRoot,sprintf('mms%d',sc),product, ...
        char(string(day,'yyyy')),char(string(day,'MM')),char(string(day,'dd')));
    key = lower(dayDir);
    if isKey(cache,key)
        dayTimes = cache(key);
    else
        files = dir(fullfile(dayDir,'*.cdf'));
        dayTimes = NaT(0,1,'TimeZone','UTC');
        for k = 1:numel(files)
            tok = regexp(files(k).name,'_(\d{14})_v','tokens','once');
            if ~isempty(tok)
                dayTimes(end+1,1) = datetime(tok{1},'InputFormat','yyyyMMddHHmmss','TimeZone','UTC'); %#ok<AGROW>
            end
        end
        dayTimes = unique(dayTimes);
        cache(key) = dayTimes;
    end
    times = union(times,dayTimes);
end
times = times(times>=eventStart & times<=eventEnd);
end
