%% Juno Orbit 7 interval shown in the supplied Schok et al. figure crop
% Displayed values are read directly from archived source fields:
%   MAG:   official Sun-State/JSO 60-s BX, BY, BZ records
%   JADE:  official Level-5 electron moment N_CC records
%   Waves: calibrated electric spectral-density records
% No coordinate rotation, vector magnitude, plasma moment, interpolation,
% channel average, or spectral-peak proxy is calculated in this script.

clear; clc;
if ~strcmp(getenv('JUNO_SKIP_MATLAB_RENDER'),'1')
    close all;
end

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
databaseRoot = 'Z:\SPART-WORK\Data\Juno';
assert(isfolder(databaseRoot), 'Juno database is unavailable: %s', databaseRoot);

% Visible interval of the supplied crop.  The author-supplied uncropped
% Figure 2 extends farther back, to 2017-06-14 15:59:05 UTC.
timeRange = [datenum(2017,6,27,8,0,0), datenum(2017,7,9,12,0,0)];
dateTicks = datenum(2017,6,29):2:datenum(2017,7,9);
waveDisplaySeconds = 300;
waveFrequencyRangeHz = [73.244 14929];
waveColorLimits = [1e-14 1e-10];

% Region intervals copied from the author-released Figure 2 object.
% Region codes: 1=intermediate, 2=outer flux pileup, 3=plasmadisc.
paperRegions = [
    datenum(2017,6,19,10,0,0),  datenum(2017,6,27,10,0,0),  1
    datenum(2017,6,29,15,10,0), datenum(2017,7,3,15,3,0),   2
    datenum(2017,7,4,0,1,0),    datenum(2017,7,5,19,16,0),  1
    datenum(2017,7,5,23,46,0),  datenum(2017,7,9,0,47,0),   3
    ];

waveFiles = discoverWaveFiles(databaseRoot, timeRange);
magFiles = discoverMagFiles(databaseRoot, timeRange);
momentFiles = discoverElectronMomentFiles(databaseRoot, timeRange);

assert(numel(waveFiles) == 13, 'Expected 13 Waves files; found %d.', numel(waveFiles));
assert(numel(magFiles) == 13, 'Expected 13 MAG files; found %d.', numel(magFiles));
assert(numel(momentFiles) == 13, ...
    'Expected 13 electron-moment files; found %d.', numel(momentFiles));

fprintf('Reading %d official JSO/SS MAG r60s files ...\n', numel(magFiles));
[magTime, magB] = processMag(magFiles, timeRange);

fprintf('Reading %d official JADE electron-moment files ...\n', numel(momentFiles));
[densityTime, electronDensity, densitySigma, densityRadiusRj] = ...
    processElectronMoments(momentFiles, timeRange);

fprintf('Selecting direct Waves source records from %d files ...\n', numel(waveFiles));
[waveEdges, waveSourceTime, wavePsd, frequencyHz] = ...
    processWaves(waveFiles, timeRange, waveDisplaySeconds);

[radiusTickSourceTime, radiusTickRj] = nearestDirectRadius( ...
    densityTime, densityRadiusRj, dateTicks, 8*3600);

% Headless fallback for systems where MATLAB R2025a+ WebGL cannot create a
% figure with the active display adapter.  MATLAB still performs every
% source-file read and direct-record selection; the MAT file contains only
% the plotted source values and fixed display metadata.
cacheFile = fullfile(outputDir, ...
    'juno_orbit7_schok_crop_direct_source_cache.mat');
if strcmp(getenv('JUNO_SKIP_MATLAB_RENDER'),'1')
    save(cacheFile,'timeRange','dateTicks','paperRegions', ...
        'magTime','magB','densityTime','electronDensity','densitySigma', ...
        'densityRadiusRj','radiusTickSourceTime','radiusTickRj', ...
        'waveEdges','waveSourceTime','wavePsd','frequencyHz', ...
        'waveFrequencyRangeHz','waveColorLimits','-v7');
    fprintf('Saved direct-source MATLAB cache: %s\n',cacheFile);
    fprintf('Set JUNO_SKIP_MATLAB_RENDER=0 to render in MATLAB.\n');
    return
end

fig = createFigure(timeRange, dateTicks, radiusTickRj, paperRegions, ...
    magTime, magB, densityTime, electronDensity, ...
    waveEdges, wavePsd, frequencyHz, waveFrequencyRangeHz, waveColorLimits);

baseName = 'juno_orbit7_schok_crop_mag_ne_waves_20170627_20170709';
pngFile = fullfile(outputDir, [baseName '.png']);
pdfFile = fullfile(outputDir, [baseName '.pdf']);
exportgraphics(fig, pngFile, 'Resolution',260);
exportgraphics(fig, pdfFile, 'ContentType','image', 'Resolution',260);

fprintf('Saved %s\n', pngFile);
fprintf('Saved %s\n', pdfFile);
fprintf('Displayed interval: %s to %s UTC\n', ...
    datestr(timeRange(1),'yyyy-mm-dd HH:MM:SS'), ...
    datestr(timeRange(2),'yyyy-mm-dd HH:MM:SS'));
fprintf('Files used: MAG=%d, electron moments=%d, Waves=%d\n', ...
    numel(magFiles), numel(momentFiles), numel(waveFiles));
fprintf('Direct records shown: MAG=%d, electron density=%d, Waves cells=%d\n', ...
    size(magB,1), nnz(isfinite(electronDensity)), ...
    nnz(any(isfinite(wavePsd),2)));
fprintf('JADE direct fields: N_CC and N_SIGMA_CC; no density is derived.\n');
fprintf('Electron-density finite range: %.4g to %.4g cm^-3\n', ...
    min(electronDensity,[],'omitnan'), max(electronDensity,[],'omitnan'));
fprintf('Electron-density uncertainty finite range: %.4g to %.4g cm^-3\n', ...
    min(densitySigma,[],'omitnan'), max(densitySigma,[],'omitnan'));
fprintf('Waves selected-source offset: maximum %.2f s\n', ...
    max(abs(waveSourceTime - 0.5*(waveEdges(1:end-1)+waveEdges(2:end))') ...
    *86400,[],'omitnan'));
for iTick = 1:numel(dateTicks)
    if isfinite(radiusTickRj(iTick))
        fprintf('  %s UTC: R=%.3f R_J from source %s UTC\n', ...
            datestr(dateTicks(iTick),'yyyy-mm-dd'), radiusTickRj(iTick), ...
            datestr(radiusTickSourceTime(iTick),'yyyy-mm-dd HH:MM:SS.FFF'));
    else
        fprintf('  %s UTC: no direct SC_POS_R record within 8 h\n', ...
            datestr(dateTicks(iTick),'yyyy-mm-dd'));
    end
end

%% Local functions
function files = discoverWaveFiles(databaseRoot, timeRange)
    folder = fullfile(databaseRoot, 'Waves_data', '2017165_ORBIT_07');
    entries = dir(fullfile(folder, 'WAV_*_E_V*.CSV'));
    files = selectDailyFiles(entries, ...
        '^WAV_(\d{4})(\d{3})T(\d{2})(\d{2})(\d{2})_E_V(\d+)\.CSV$', ...
        timeRange, true);
end

function files = discoverMagFiles(databaseRoot, timeRange)
    folder = fullfile(databaseRoot, 'MAG', '2017', 'JUPITER', 'SS', 'PERI-07');
    entries = dir(fullfile(folder, 'fgm_jno_l3_2017*ss_r60s_v*.sts'));
    files = selectDailyFiles(entries, ...
        '^fgm_jno_l3_(\d{4})(\d{3})ss_r60s_v(\d+)\.sts$', ...
        timeRange, false);
end

function files = discoverElectronMomentFiles(databaseRoot, timeRange)
    firstDay = floor(timeRange(1));
    lastDay = floor(timeRange(2)-eps);
    files = cell(lastDay-firstDay+1,1);
    for iDay = 1:numel(files)
        dayNumber = firstDay+iDay-1;
        parts = datevec(dayNumber);
        yearValue = parts(1);
        doyValue = dayNumber-datenum(yearValue,1,1)+1;
        dayCode = sprintf('%04d%03d',yearValue,doyValue);
        folder = fullfile(databaseRoot,'JADE_data',sprintf('%04d',yearValue), ...
            dayCode,'ELECTRONS');
        entries = dir(fullfile(folder, ...
            sprintf('JAD_L50_HLS_ELC_MOM_ISO_2D_ELECTRONS_%s_V*.CSV',dayCode)));
        assert(~isempty(entries), 'Missing official electron moments for %s.', dayCode);
        versions = NaN(numel(entries),1);
        for iEntry = 1:numel(entries)
            token = regexp(entries(iEntry).name,'_V(\d+)\.CSV$','tokens','once');
            assert(~isempty(token), 'Unexpected moment filename: %s',entries(iEntry).name);
            versions(iEntry) = str2double(token{1});
        end
        [~,best] = max(versions);
        files{iDay} = fullfile(entries(best).folder,entries(best).name);
    end
end

function files = selectDailyFiles(entries, expression, timeRange, hasClock)
    paths = {};
    starts = [];
    versions = [];
    keys = {};
    for iEntry = 1:numel(entries)
        token = regexp(entries(iEntry).name, expression, 'tokens','once');
        if isempty(token)
            continue
        end
        if hasClock
            fileStart = datenum(str2double(token{1}),1,1, ...
                str2double(token{3}),str2double(token{4}),str2double(token{5})) + ...
                str2double(token{2})-1;
            version = str2double(token{6});
        else
            fileStart = datenum(str2double(token{1}),1,1) + str2double(token{2})-1;
            version = str2double(token{3});
        end
        if fileStart < timeRange(2) && fileStart+1 > timeRange(1)
            paths{end+1,1} = fullfile(entries(iEntry).folder,entries(iEntry).name); %#ok<AGROW>
            starts(end+1,1) = fileStart; %#ok<AGROW>
            versions(end+1,1) = version; %#ok<AGROW>
            keys{end+1,1} = sprintf('%04d%03d',str2double(token{1}),str2double(token{2})); %#ok<AGROW>
        end
    end
    [uniqueKeys,~,group] = unique(keys,'stable');
    selected = NaN(numel(uniqueKeys),1);
    for iKey = 1:numel(uniqueKeys)
        candidates = find(group==iKey);
        [~,localBest] = max(versions(candidates));
        selected(iKey) = candidates(localBest);
    end
    paths = paths(selected);
    starts = starts(selected);
    [~,order] = sort(starts);
    files = paths(order);
end

function [timeNum, magneticField] = processMag(files, timeRange)
    timeNum = [];
    magneticField = [];
    for iFile = 1:numel(files)
        fprintf('  MAG %2d/%2d: %s\n',iFile,numel(files),files{iFile});
        [thisTime,thisField] = readJunoMag(files{iFile});
        use = thisTime>=timeRange(1) & thisTime<timeRange(2);
        timeNum = [timeNum;thisTime(use)]; %#ok<AGROW>
        magneticField = [magneticField;thisField(use,:)]; %#ok<AGROW>
    end
    [timeNum,order] = sort(timeNum);
    magneticField = magneticField(order,:);
    [timeNum,uniqueIndex] = unique(timeNum,'last');
    magneticField = magneticField(uniqueIndex,:);
end

function [timeNum,magneticField] = readJunoMag(fileName)
    assert(~isempty(regexp(char(fileName),'ss_r60s_v\d+\.sts$','once')), ...
        'Expected official SS/JSO r60s MAG: %s',fileName);
    fid = fopen(fileName,'r');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    dataStart = -1;
    while ~feof(fid)
        lineStart = ftell(fid);
        thisLine = fgetl(fid);
        if ischar(thisLine) && ~isempty(regexp(thisLine, ...
                '^\s*20\d{2}\s+\d{1,3}\s+\d+','once'))
            dataStart = lineStart;
            break
        end
    end
    assert(dataStart>=0,'No numeric MAG records found in %s',fileName);
    fseek(fid,dataStart,'bof');
    values = textscan(fid,repmat('%f',1,14),'CollectOutput',true, ...
        'ReturnOnError',false);
    values = values{1};
    assert(size(values,2)==14,'Unexpected MAG row format in %s',fileName);
    timeNum = datenum(values(1,1),1,1)+values(:,7)-1;
    magneticField = values(:,8:10);
    magneticField(abs(magneticField)>1e5) = NaN;
end

function [timeNum,density,densitySigma,radiusRj] = ...
        processElectronMoments(files,timeRange)
    timeNum = [];
    density = [];
    densitySigma = [];
    radiusRj = [];
    for iFile = 1:numel(files)
        fprintf('  JADE %2d/%2d: %s\n',iFile,numel(files),files{iFile});
        [thisTime,thisDensity,thisSigma,thisRadius] = ...
            readElectronMoments(files{iFile});
        use = thisTime>=timeRange(1) & thisTime<timeRange(2);
        timeNum = [timeNum;thisTime(use)]; %#ok<AGROW>
        density = [density;thisDensity(use)]; %#ok<AGROW>
        densitySigma = [densitySigma;thisSigma(use)]; %#ok<AGROW>
        radiusRj = [radiusRj;thisRadius(use)]; %#ok<AGROW>
    end
    [timeNum,order] = sort(timeNum);
    density = density(order);
    densitySigma = densitySigma(order);
    radiusRj = radiusRj(order);
    [timeNum,uniqueIndex] = unique(timeNum,'last');
    density = density(uniqueIndex);
    densitySigma = densitySigma(uniqueIndex);
    radiusRj = radiusRj(uniqueIndex);
end

function [timeNum,density,densitySigma,radiusRj] = readElectronMoments(fileName)
    fid = fopen(fileName,'r');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    values = textscan(fid,['%q' repmat('%f',1,23)],'Delimiter',',', ...
        'CollectOutput',true,'ReturnOnError',false,'EmptyValue',NaN);
    utc = values{1};
    numericValues = values{2};
    assert(size(numericValues,2)==23,'Unexpected moments row format in %s',fileName);
    timeNum = parseDoyUtc(utc);
    radiusRj = numericValues(:,9);      % Direct SC_POS_R field, R_J
    density = numericValues(:,17);      % Direct N_CC field, cm^-3
    densitySigma = numericValues(:,18); % Direct N_SIGMA_CC field, cm^-3
    radiusRj(radiusRj<0 | radiusRj>1000) = NaN;
    density(density<=0 | density>1e8) = NaN;
    densitySigma(densitySigma<0 | densitySigma>1e8) = NaN;
end

function timeNum = parseDoyUtc(utc)
    utcChars = char(utc);
    years = str2double(cellstr(utcChars(:,1:4)));
    doys = str2double(cellstr(utcChars(:,6:8)));
    hours = str2double(cellstr(utcChars(:,10:11)));
    minutes = str2double(cellstr(utcChars(:,13:14)));
    seconds = str2double(cellstr(utcChars(:,16:end)));
    timeNum = datenum(years,ones(size(years)),ones(size(years)), ...
        hours,minutes,seconds)+doys-1;
end

function [timeEdges,sourceTime,gridPsd,frequencyHz] = ...
        processWaves(files,timeRange,displaySeconds)
    nBin = round(diff(timeRange)*86400/displaySeconds);
    assert(abs(nBin*displaySeconds-diff(timeRange)*86400)<1e-3, ...
        'Time range must contain an integer number of display intervals.');
    timeEdges = timeRange(1)+(0:nBin)*displaySeconds/86400;
    binCenters = 0.5*(timeEdges(1:end-1)+timeEdges(2:end));
    sourceTime = NaN(nBin,1);
    gridPsd = NaN(nBin,126);
    sourceDistance = inf(nBin,1);
    rawFrequency = [];
    for iFile = 1:numel(files)
        fprintf('  Waves %2d/%2d: %s\n',iFile,numel(files),files{iFile});
        [thisBin,thisTime,thisPsd,thisFrequency] = ...
            readJunoWavesSurveyDirect(files{iFile},timeRange, ...
            displaySeconds,nBin,binCenters);
        if isempty(rawFrequency)
            rawFrequency = thisFrequency(:);
        else
            assert(all(abs(thisFrequency(:)-rawFrequency)<= ...
                max(1e-9,1e-10*abs(rawFrequency))), ...
                'Waves frequency grids differ: %s',files{iFile});
        end
        for iRecord = 1:numel(thisBin)
            targetBin = thisBin(iRecord);
            distance = abs(thisTime(iRecord)-binCenters(targetBin));
            if distance<sourceDistance(targetBin)
                sourceTime(targetBin) = thisTime(iRecord);
                gridPsd(targetBin,:) = thisPsd(iRecord,:);
                sourceDistance(targetBin) = distance;
            end
        end
    end
    [frequencyHz,frequencyOrder] = sort(rawFrequency);
    gridPsd = gridPsd(:,frequencyOrder);
    validFrequency = isfinite(frequencyHz) & frequencyHz>0;
    frequencyHz = frequencyHz(validFrequency);
    gridPsd = gridPsd(:,validFrequency);
    gridPsd(gridPsd<=0) = NaN;
end

function [binIndex,timeNum,psd,frequencyHz] = ...
        readJunoWavesSurveyDirect(fileName,timeRange,binSeconds,nBin,binCenters)
    % Stream the large survey CSV.  Only the single source record nearest
    % each display-bin center is parsed in full; no spectral/time average
    % or interpolation is performed.
    fid = fopen(fileName,'r');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    header = cell(5,1);
    for iLine = 1:5
        header{iLine} = fgetl(fid);
    end
    headerFields = split(string(header{4}),',');
    headerFields = erase(headerFields,'"');
    assert(numel(headerFields)>=154,'Incomplete Waves header in %s',fileName);
    frequencyHz = str2double(headerFields(29:154)).';
    assert(numel(frequencyHz)==126 && all(isfinite(frequencyHz)), ...
        'Invalid Waves frequency grid in %s',fileName);

    binIndex = zeros(nBin,1);
    timeNum = NaN(nBin,1);
    psd = NaN(nBin,126);
    outputCount = 0;
    activeBin = NaN;
    bestDistance = Inf;
    bestTime = NaN;
    bestLine = '';
    while ~feof(fid)
        thisLine = fgetl(fid);
        if ~ischar(thisLine) || isempty(thisLine)
            continue
        end
        firstComma = find(thisLine==',',1,'first');
        if isempty(firstComma)
            continue
        end
        secondRelative = find(thisLine(firstComma+1:end)==',',1,'first');
        if isempty(secondRelative)
            continue
        end
        secondComma = firstComma+secondRelative;
        scet = thisLine(firstComma+1:secondComma-1);
        if numel(scet)<21
            continue
        end
        thisTime = parseSingleDoyUtc(scet);
        if thisTime<timeRange(1) || thisTime>=timeRange(2)
            continue
        end
        thisBin = floor((thisTime-timeRange(1))*86400/binSeconds)+1;
        if thisBin<1 || thisBin>nBin
            continue
        end
        if ~isfinite(activeBin) || thisBin~=activeBin
            if isfinite(activeBin) && ~isempty(bestLine)
                [rowPsd,rowOkay] = parseDirectWaveLine(bestLine);
                if rowOkay
                    outputCount = outputCount+1;
                    binIndex(outputCount) = activeBin;
                    timeNum(outputCount) = bestTime;
                    psd(outputCount,:) = rowPsd;
                end
            end
            activeBin = thisBin;
            bestDistance = Inf;
            bestTime = NaN;
            bestLine = '';
        end
        thisDistance = abs(thisTime-binCenters(thisBin));
        if thisDistance<bestDistance
            bestDistance = thisDistance;
            bestTime = thisTime;
            bestLine = thisLine;
        end
    end
    if isfinite(activeBin) && ~isempty(bestLine)
        [rowPsd,rowOkay] = parseDirectWaveLine(bestLine);
        if rowOkay
            outputCount = outputCount+1;
            binIndex(outputCount) = activeBin;
            timeNum(outputCount) = bestTime;
            psd(outputCount,:) = rowPsd;
        end
    end
    binIndex = binIndex(1:outputCount);
    timeNum = timeNum(1:outputCount);
    psd = psd(1:outputCount,:);
end

function [rowPsd,rowOkay] = parseDirectWaveLine(thisLine)
    fields = strsplit(thisLine,',');
    rowOkay = numel(fields)==154;
    rowPsd = NaN(1,126);
    if ~rowOkay
        return
    end
    rowPsd = str2double(fields(29:154));
    qualityFlag = erase(fields{4},'"');
    if flagIsSet(qualityFlag,2)
        rowPsd(1:43) = NaN;
    end
    if flagIsSet(qualityFlag,3)
        rowPsd(44:61) = NaN;
    end
end

function timeNum = parseSingleDoyUtc(scet)
    digits = double(scet)-double('0');
    yearValue = 1000*digits(1)+100*digits(2)+10*digits(3)+digits(4);
    doyValue = 100*digits(6)+10*digits(7)+digits(8);
    hourValue = 10*digits(10)+digits(11);
    minuteValue = 10*digits(13)+digits(14);
    secondValue = 10*digits(16)+digits(17)+str2double(scet(18:end));
    timeNum = datenum(yearValue,1,1)+doyValue-1+ ...
        (hourValue*3600+minuteValue*60+secondValue)/86400;
end

function tf = flagIsSet(flagString,flagNumber)
    parts = strsplit(flagString,':');
    tf = numel(parts)>=flagNumber && strcmp(parts{flagNumber},'1');
end

function [sourceTime,radiusRj] = nearestDirectRadius( ...
        timeNum,sourceRadiusRj,tickValues,maxOffsetSeconds)
    valid = isfinite(timeNum) & isfinite(sourceRadiusRj);
    timeNum = timeNum(valid);
    sourceRadiusRj = sourceRadiusRj(valid);
    sourceTime = NaN(size(tickValues));
    radiusRj = NaN(size(tickValues));
    for iTick = 1:numel(tickValues)
        [distance,index] = min(abs(timeNum-tickValues(iTick)));
        if distance*86400<=maxOffsetSeconds
            sourceTime(iTick) = timeNum(index);
            radiusRj(iTick) = sourceRadiusRj(index);
        end
    end
end

function fig = createFigure(timeRange,dateTicks,radiusTickRj,paperRegions, ...
        magTime,magB,densityTime,electronDensity, ...
        waveEdges,wavePsd,frequencyHz,waveFrequencyRangeHz,waveColorLimits)
    fig = figure('Color','w','Position',[25 25 1800 1250], ...
        'Renderer','opengl');
    left = 0.085;
    width = 0.805;
    axMag = axes(fig,'Position',[left 0.705 width 0.205]);
    axDensity = axes(fig,'Position',[left 0.455 width 0.195]);
    axWave = axes(fig,'Position',[left 0.105 width 0.295]);

    styleAxis(axMag);
    xlim(axMag,timeRange);
    ylim(axMag,[-35 35]);
    shadePaperRegions(axMag,paperRegions);
    hold(axMag,'on');
    hBx = plot(axMag,magTime,magB(:,1),'Color',[0 0.35 0.85],'LineWidth',0.7);
    hBy = plot(axMag,magTime,magB(:,2),'Color',[0.05 0.58 0.18],'LineWidth',0.7);
    hBz = plot(axMag,magTime,magB(:,3),'Color',[0.85 0.12 0.12],'LineWidth',0.7);
    ylabel(axMag,'B_{JSO} (nT)');
    magLegend = legend(axMag,[hBx hBy hBz], ...
        {'B_{x,JSO}','B_{y,JSO}','B_{z,JSO}'}, ...
        'Location','northwest','NumColumns',3,'Box','off','FontSize',9);
    magLegend.AutoUpdate = 'off';

    styleAxis(axDensity);
    xlim(axDensity,timeRange);
    set(axDensity,'YScale','log');
    ylim(axDensity,[1e-4 1]);
    shadePaperRegions(axDensity,paperRegions);
    hold(axDensity,'on');
    plot(axDensity,densityTime,electronDensity,'.','Color',[0.05 0.35 0.78], ...
        'MarkerSize',4,'LineStyle','none');
    ylabel(axDensity,'n_e (cm^{-3})');

    plotWaveSpectrum(axWave,waveEdges,wavePsd,frequencyHz, ...
        waveFrequencyRangeHz,waveColorLimits);
    xlim(axWave,timeRange);

    set([axMag axDensity axWave],'XTick',dateTicks);
    axMag.XTickLabel = [];
    axDensity.XTickLabel = [];
    axWave.XTickLabel = cellstr(datestr(dateTicks,'yyyy-mm-dd'));
    axWave.XTickLabelRotation = 18;

    linkaxes([axMag axDensity axWave],'x');
    drawnow;

    radiusLabels = repmat({'--'},size(dateTicks));
    for iTick = 1:numel(dateTicks)
        if isfinite(radiusTickRj(iTick))
            radiusLabels{iTick} = sprintf('%.1f',radiusTickRj(iTick));
        end
    end
    axRadius = axes(fig,'Position',axMag.Position,'Color','none', ...
        'XAxisLocation','top','YAxisLocation','right','YColor','none', ...
        'XLim',timeRange,'XTick',dateTicks,'XTickLabel',radiusLabels, ...
        'YTick',[],'TickDir','out','Box','off','FontSize',9, ...
        'HandleVisibility','off');
    xlabel(axRadius,'Radial distance (R_J)');
    linkaxes([axMag axDensity axWave axRadius],'x');

    colormap(fig,parula(256));
    wavePosition = axWave.Position;
    cb = colorbar(axWave,'eastoutside');
    axWave.Position = wavePosition;
    cb.Position = [0.915 wavePosition(2) 0.015 wavePosition(4)];
    cb.Ticks = 10.^(-14:-10);
    cb.Label.String = 'E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 9;
end

function shadePaperRegions(ax,paperRegions)
    colors = [0.88 0.82 0.91; 0.72 0.72 0.72; 0.98 0.72 0.68];
    alphaValues = [0.48 0.48 0.48];
    hold(ax,'on');
    yLimits = ylim(ax);
    for iRegion = 1:size(paperRegions,1)
        t1 = paperRegions(iRegion,1);
        t2 = paperRegions(iRegion,2);
        if t2<=ax.XLim(1) || t1>=ax.XLim(2)
            continue
        end
        code = paperRegions(iRegion,3);
        h = patch(ax,[t1 t2 t2 t1], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)],colors(code,:), ...
            'FaceAlpha',alphaValues(code),'EdgeColor','none', ...
            'HandleVisibility','off');
        uistack(h,'bottom');
    end
end

function plotWaveSpectrum(ax,timeEdges,psd,frequencyHz,frequencyRange,colorLimits)
    useFrequency = frequencyHz>=frequencyRange(1)/1.2 & ...
        frequencyHz<=frequencyRange(2)*1.2;
    frequencyHz = frequencyHz(useFrequency);
    directPsd = psd(:,useFrequency).';
    frequencyMidpoints = sqrt(frequencyHz(1:end-1).*frequencyHz(2:end));
    frequencyEdges = [frequencyHz(1)^2/frequencyMidpoints(1); ...
        frequencyMidpoints; frequencyHz(end)^2/frequencyMidpoints(end)];
    [xGrid,yGrid] = meshgrid(timeEdges,frequencyEdges);
    colorData = NaN(size(xGrid));
    colorData(1:end-1,1:end-1) = directPsd;
    alphaData = zeros(size(xGrid));
    alphaData(1:end-1,1:end-1) = isfinite(directPsd) & directPsd>0;
    h = surface(ax,xGrid,yGrid,zeros(size(xGrid)),colorData, ...
        'FaceColor','flat','EdgeColor','none');
    set(h,'AlphaData',alphaData,'FaceAlpha','flat','AlphaDataMapping','none');
    view(ax,2);
    set(ax,'YScale','log','YDir','normal','ColorScale','log', ...
        'Layer','top','FontSize',10,'TickDir','out','Box','on', ...
        'Color',[0.88 0.88 0.88],'XGrid','off','YGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off');
    ylim(ax,frequencyRange);
    yticks(ax,[1e2 1e3 1e4]);
    yticklabels(ax,{'100 Hz','1 kHz','10 kHz'});
    ylabel(ax,'Frequency');
    caxis(ax,colorLimits);
end

function styleAxis(ax)
    set(ax,'FontSize',10,'TickDir','out','Box','on','Layer','top', ...
        'Color','w','XGrid','off','YGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off');
end
