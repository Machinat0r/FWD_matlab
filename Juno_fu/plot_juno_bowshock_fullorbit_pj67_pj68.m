%% Juno PJ67-to-PJ68 full revolution: direct JSO MAG and Waves measurements
% This physical perijove-to-perijove revolution contains both bow-shock
% crossings reported on 2024-12-10/11.  The PDS orbit number changes at the
% intervening apojove (2024-12-10 21:12:42 UTC), so Waves ORBIT 67/68 and
% MAG PERI-67/68 must be combined.
%
% Data products:
%   JNO-J-3-FGM-CAL-V1.0, official Sun-State/JSO 60-s MAG
%   JNO-E/J/SS-WAV-3-CDR-SRVFULL-V2.0, calibrated E spectra
%   JNO-J/SW-JAD-3-CALIBRATED-V1.0, official SC_POS_R metadata only
%
% Every displayed physical value is copied from a source record.  Waves uses
% one actual source record per 600-s display interval; MAG reads the official
% r60s source records directly.  Each x tick uses the nearest JADE SC_POS_R
% source record, already expressed in R_J.  No coordinate rotation, plasma
% moment, spectral-peak proxy, bin statistic, or vector magnitude is calculated.

clear; close all; clc;

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
databaseRoot = 'Z:\SPART-WORK\Data\Juno';
assert(isfolder(databaseRoot), 'Juno database is unavailable: %s', databaseRoot);

% Actual PJ times from the final trajectory geometry.  A perijove-to-
% perijove interval is used so both December bow-shock cases remain in one
% complete physical revolution.
tOrbit = [datenum(2024,11,24,13,5,28), ...
          datenum(2024,12,27, 5,22,30)];
shockRampRanges = [
    datenum(2024,12,10,17,24,51.757), datenum(2024,12,10,17,24,54.507)
    datenum(2024,12,11, 3, 6, 1.461), datenum(2024,12,11, 3, 6, 4.211)
    ];
shockTimes = mean(shockRampRanges,2);
assert(all(shockRampRanges(:) >= tOrbit(1) & shockRampRanges(:) < tOrbit(2)), ...
    'A verified shock ramp lies outside the full-revolution interval.');

waveBinSeconds = 600;

waveFiles = discoverWaveFiles(databaseRoot, tOrbit);
magFiles = discoverMagFiles(databaseRoot, tOrbit);
dateTicks = ceil(tOrbit(1)):2:floor(tOrbit(2));
jadePositionFiles = discoverJadePositionFiles(databaseRoot, dateTicks);

assert(numel(waveFiles) == 35, ...
    'Expected 35 Waves E files for PJ67-to-PJ68; found %d.', numel(waveFiles));
assert(numel(magFiles) == 34, ...
    'Expected 34 MAG daily files for PJ67-to-PJ68; found %d.', numel(magFiles));
assert(numel(jadePositionFiles) == numel(dateTicks), ...
    'Expected one JADE SC_POS_R file per x tick; found %d for %d ticks.', ...
    numel(jadePositionFiles), numel(dateTicks));

fprintf('Selecting direct source records from %d Waves files at %d-s spacing ...\n', ...
    numel(waveFiles), waveBinSeconds);
[waveTime, wavePsd, frequencyHz, receiverBoundaries] = ...
    processWaves(waveFiles, tOrbit, waveBinSeconds);

fprintf('Reading direct official JSO/SS records from %d MAG r60s files ...\n', ...
    numel(magFiles));
[magTime, magB] = processMag(magFiles, tOrbit);

fprintf('Reading direct JADE SC_POS_R records for %d x ticks ...\n', ...
    numel(dateTicks));
[positionSourceTime, tickRadiusRj] = ...
    readJadeRadiusAtTicks(jadePositionFiles, dateTicks);

finiteWave = wavePsd(isfinite(wavePsd) & wavePsd > 0);
logColorLimits = robustLimits(log10(finiteWave), 0.005, 0.997);
logColorLimits = [floor(2*logColorLimits(1))/2, ...
    ceil(2*logColorLimits(2))/2];
if diff(logColorLimits) < 5
    centerValue = mean(logColorLimits);
    logColorLimits = centerValue + [-2.5 2.5];
end
colorLimits = 10.^logColorLimits;

fig = createFullOrbitFigure(tOrbit, magTime, magB, ...
    waveTime, wavePsd, frequencyHz, ...
    receiverBoundaries, colorLimits, shockTimes, ...
    dateTicks, tickRadiusRj);

baseName = 'juno_pj67_to_pj68_mag_waves_fullband_20241124_20241227';
pngFile = fullfile(outputDir, [baseName '.png']);
pdfFile = fullfile(outputDir, [baseName '.pdf']);
exportgraphics(fig, pngFile, 'Resolution',260);
exportgraphics(fig, pdfFile, 'ContentType','image', 'Resolution',260);

fprintf('Saved %s\n', pngFile);
fprintf('Saved %s\n', pdfFile);
fprintf('Physical revolution: %s to %s UTC (%.3f days)\n', ...
    datestr(tOrbit(1),'yyyy-mm-dd HH:MM:SS'), ...
    datestr(tOrbit(2),'yyyy-mm-dd HH:MM:SS'), diff(tOrbit));
fprintf('Files used: Waves=%d, MAG=%d, JADE position=%d\n', ...
    numel(waveFiles), numel(magFiles), numel(jadePositionFiles));
fprintf('Direct source samples: Waves=%d/%d, MAG=%d/%d\n', ...
    nnz(any(isfinite(wavePsd),2)), size(wavePsd,1), ...
    nnz(any(isfinite(magB),2)), size(magB,1));
fprintf('Frequency centers: %.4g Hz to %.4g Hz (%d channels)\n', ...
    min(frequencyHz), max(frequencyHz), numel(frequencyHz));
fprintf('Receiver separators: %.7g, %.7g, %.7g Hz\n', receiverBoundaries);
fprintf('Color limits: %.4g to %.4g (V/m)^2/Hz (log color scale)\n', ...
    colorLimits);
fprintf('Shock markers: Case I %s UTC; Case II %s UTC\n', ...
    datestr(shockTimes(1),'yyyy-mm-dd HH:MM:SS.FFF'), ...
    datestr(shockTimes(2),'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('Tick positions use direct JADE SC_POS_R records:\n');
for iTick = 1:numel(dateTicks)
    fprintf('  %s UTC: R=%.3f R_J (%s source)\n', ...
        datestr(dateTicks(iTick),'yyyy-mm-dd'), tickRadiusRj(iTick), ...
        datestr(positionSourceTime(iTick),'HH:MM:SS.FFF'));
end

%% Local functions
function files = discoverWaveFiles(databaseRoot, timeRange)
    orbitFolders = {
        fullfile(databaseRoot, 'Waves_data', '2024313_ORBIT_67')
        fullfile(databaseRoot, 'Waves_data', '2024345_ORBIT_68')
        };
    paths = {};
    starts = [];
    versions = [];
    keys = {};
    for iFolder = 1:numel(orbitFolders)
        entries = dir(fullfile(orbitFolders{iFolder}, 'WAV_*_E_V*.CSV'));
        for iEntry = 1:numel(entries)
            token = regexp(entries(iEntry).name, ...
                '^WAV_(\d{4})(\d{3})T(\d{2})(\d{2})(\d{2})_E_V(\d+)\.CSV$', ...
                'tokens', 'once');
            if isempty(token)
                continue
            end
            fileStart = datenum(str2double(token{1}),1,1, ...
                str2double(token{3}),str2double(token{4}),str2double(token{5})) + ...
                str2double(token{2}) - 1;
            if fileStart < timeRange(2) && fileStart + 1 > timeRange(1)
                paths{end+1,1} = fullfile(entries(iEntry).folder, entries(iEntry).name); %#ok<AGROW>
                starts(end+1,1) = fileStart; %#ok<AGROW>
                versions(end+1,1) = str2double(token{6}); %#ok<AGROW>
                keys{end+1,1} = regexprep(entries(iEntry).name, ...
                    '_V\d+\.CSV$', ''); %#ok<AGROW>
            end
        end
    end
    selected = selectHighestVersions(keys, versions);
    paths = paths(selected);
    starts = starts(selected);
    [~, order] = sort(starts);
    files = paths(order);
end

function files = discoverMagFiles(databaseRoot, timeRange)
    periFolders = {
        fullfile(databaseRoot, 'MAG', '2024', 'JUPITER', 'SS', 'PERI-67')
        fullfile(databaseRoot, 'MAG', '2024', 'JUPITER', 'SS', 'PERI-68')
        };
    paths = {};
    starts = [];
    versions = [];
    keys = {};
    for iFolder = 1:numel(periFolders)
        entries = dir(fullfile(periFolders{iFolder}, ...
            'fgm_jno_l3_2024*ss_r60s_v*.sts'));
        for iEntry = 1:numel(entries)
            token = regexp(entries(iEntry).name, ...
                '^fgm_jno_l3_(\d{4})(\d{3})ss_r60s_v(\d+)\.sts$', ...
                'tokens', 'once');
            if isempty(token)
                continue
            end
            fileStart = datenum(str2double(token{1}),1,1) + ...
                str2double(token{2}) - 1;
            if fileStart < timeRange(2) && fileStart + 1 > timeRange(1)
                paths{end+1,1} = fullfile(entries(iEntry).folder, entries(iEntry).name); %#ok<AGROW>
                starts(end+1,1) = fileStart; %#ok<AGROW>
                versions(end+1,1) = str2double(token{3}); %#ok<AGROW>
                keys{end+1,1} = sprintf('%s%s', token{1}, token{2}); %#ok<AGROW>
            end
        end
    end
    selected = selectHighestVersions(keys, versions);
    paths = paths(selected);
    starts = starts(selected);
    [~, order] = sort(starts);
    files = paths(order);
end

function files = discoverJadePositionFiles(databaseRoot, tickValues)
    files = cell(numel(tickValues),1);
    for iTick = 1:numel(tickValues)
        dateParts = datevec(tickValues(iTick));
        yearValue = dateParts(1);
        doyValue = floor(tickValues(iTick)-datenum(yearValue,1,1))+1;
        dateCode = sprintf('%04d%03d', yearValue, doyValue);
        folder = fullfile(databaseRoot, 'JADE_data', ...
            sprintf('%04d',yearValue), dateCode, 'ION_LOGICALS');
        entries = dir(fullfile(folder, ...
            sprintf('JAD_L30_HLS_ION_LOG_CNT_%s_V*.DAT',dateCode)));
        assert(~isempty(entries), ...
            'Missing JADE SC_POS_R source file for %s.', dateCode);
        versions = NaN(numel(entries),1);
        for iEntry = 1:numel(entries)
            token = regexp(entries(iEntry).name, '_V(\d+)\.DAT$', ...
                'tokens','once');
            assert(~isempty(token), 'Unexpected JADE filename: %s', ...
                entries(iEntry).name);
            versions(iEntry) = str2double(token{1});
        end
        [~,bestIndex] = max(versions);
        files{iTick} = fullfile(entries(bestIndex).folder, ...
            entries(bestIndex).name);
    end
end

function selected = selectHighestVersions(keys, versions)
    if isempty(keys)
        selected = [];
        return
    end
    [uniqueKeys,~,group] = unique(keys, 'stable');
    selected = NaN(numel(uniqueKeys),1);
    for iKey = 1:numel(uniqueKeys)
        candidates = find(group == iKey);
        [~, localBest] = max(versions(candidates));
        selected(iKey) = candidates(localBest);
    end
end

function [gridTime, gridPsd, frequencyHz, receiverBoundaries] = ...
        processWaves(files, timeRange, binSeconds)
    nBin = ceil(diff(timeRange)*86400/binSeconds);
    gridTime = NaN(nBin,1);
    gridPsd = NaN(nBin,126);
    gridDistance = inf(nBin,1);
    rawFrequency = [];

    for iFile = 1:numel(files)
        fprintf('  Waves %2d/%2d: %s\n', iFile, numel(files), files{iFile});
        [timeNum, psd, thisFrequency] = readJunoWavesSurvey(files{iFile});
        if isempty(rawFrequency)
            rawFrequency = thisFrequency(:);
        else
            tolerance = max(1e-9, 1e-10*abs(rawFrequency));
            assert(all(abs(thisFrequency(:)-rawFrequency) <= tolerance), ...
                'Waves frequency grids differ: %s', files{iFile});
        end
        use = timeNum >= timeRange(1) & timeNum < timeRange(2);
        if any(use)
            [gridTime,gridPsd,gridDistance] = selectDirectSamples( ...
                gridTime,gridPsd,gridDistance,timeNum(use),psd(use,:), ...
                timeRange(1),binSeconds);
        end
        clear timeNum psd thisFrequency
    end

    receiverChannels = {1:43, 44:61, 62:88, 89:126};
    assert(numel(rawFrequency) == 126, ...
        'Unexpected Waves survey channel count: %d', numel(rawFrequency));
    receiverBoundaries = NaN(1,3);
    for iReceiver = 1:3
        lowerMaximum = max(rawFrequency(receiverChannels{iReceiver}));
        upperCenters = rawFrequency(receiverChannels{iReceiver+1});
        nextUpperCenter = min(upperCenters(upperCenters > lowerMaximum));
        assert(~isempty(nextUpperCenter), ...
            'Unable to derive receiver boundary %d.', iReceiver);
        receiverBoundaries(iReceiver) = sqrt(lowerMaximum*nextUpperCenter);
    end

    [frequencyHz, frequencyOrder] = sort(rawFrequency);
    gridPsd = gridPsd(:,frequencyOrder);
    validFrequency = isfinite(frequencyHz) & frequencyHz > 0;
    frequencyHz = frequencyHz(validFrequency);
    gridPsd = gridPsd(:,validFrequency);
    gridPsd(gridPsd <= 0) = NaN;
end

function [timeNum, magneticField] = processMag(files, timeRange)
    timeNum = [];
    magneticField = [];
    for iFile = 1:numel(files)
        fprintf('  MAG   %2d/%2d: %s\n', iFile, numel(files), files{iFile});
        [thisTime, thisField] = readJunoMag(files{iFile});
        use = thisTime >= timeRange(1) & thisTime < timeRange(2);
        timeNum = [timeNum; thisTime(use)]; %#ok<AGROW>
        magneticField = [magneticField; thisField(use,:)]; %#ok<AGROW>
    end
    [timeNum, order] = sort(timeNum);
    magneticField = magneticField(order,:);
    [timeNum, uniqueIndex] = unique(timeNum, 'last');
    magneticField = magneticField(uniqueIndex,:);
end

function [sourceTime, radiusRj] = readJadeRadiusAtTicks(files, tickValues)
    sourceTime = NaN(numel(tickValues),1);
    radiusRj = NaN(numel(tickValues),1);
    for iTick = 1:numel(tickValues)
        fprintf('  JADE position %2d/%2d: %s\n', ...
            iTick, numel(tickValues), files{iTick});
        [thisTime, thisRadiusRj] = readJadeScPosR(files{iTick});
        valid = isfinite(thisTime) & isfinite(thisRadiusRj) & ...
            thisRadiusRj >= 0 & thisRadiusRj < 65535;
        thisTime = thisTime(valid);
        thisRadiusRj = thisRadiusRj(valid);
        assert(~isempty(thisTime), 'No valid SC_POS_R records in %s.', ...
            files{iTick});
        [distance,sourceIndex] = min(abs(thisTime-tickValues(iTick)));
        assert(distance <= 300/86400, ...
            'No JADE SC_POS_R record within 300 s of tick %s.', ...
            datestr(tickValues(iTick),'yyyy-mm-dd HH:MM:SS'));
        sourceTime(iTick) = thisTime(sourceIndex);
        radiusRj(iTick) = thisRadiusRj(sourceIndex);
    end
end

function [timeNum, radiusRj] = readJadeScPosR(fileName)
    recordBytes = 45120;
    utcBytesPerRecord = 21;
    scPosRZeroBasedOffset = 88;
    fileInfo = dir(fileName);
    assert(~isempty(fileInfo) && mod(fileInfo.bytes,recordBytes) == 0, ...
        'Unexpected JADE record structure in %s.', fileName);
    recordCount = fileInfo.bytes/recordBytes;

    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fseek(fid,0,'bof');
    [utcBytes,utcCount] = fread(fid,[utcBytesPerRecord recordCount], ...
        '21*uint8=>uint8',recordBytes-utcBytesPerRecord);
    assert(utcCount == utcBytesPerRecord*recordCount, ...
        'Incomplete DIM0_UTC read from %s.', fileName);
    fseek(fid,scPosRZeroBasedOffset,'bof');
    [radiusRj,radiusCount] = fread(fid,recordCount, ...
        'single=>double',recordBytes-4);
    assert(radiusCount == recordCount, ...
        'Incomplete SC_POS_R read from %s.', fileName);

    utcChars = char(utcBytes.');
    years = str2double(cellstr(utcChars(:,1:4)));
    doys = str2double(cellstr(utcChars(:,6:8)));
    hours = str2double(cellstr(utcChars(:,10:11)));
    minutes = str2double(cellstr(utcChars(:,13:14)));
    seconds = str2double(cellstr(utcChars(:,16:21)));
    timeNum = datenum(years,ones(size(years)),ones(size(years)), ...
        hours,minutes,seconds) + doys - 1;
end

function [targetTime,targetValues,targetDistance] = selectDirectSamples( ...
        targetTime,targetValues,targetDistance,timeNum,values,t0,binSeconds)
    if isempty(timeNum)
        return
    end
    binIndex = floor((timeNum-t0)*86400/binSeconds) + 1;
    valid = binIndex >= 1 & binIndex <= size(targetValues,1);
    binIndex = binIndex(valid);
    timeNum = timeNum(valid);
    values = values(valid,:);
    binCenters = t0 + (binIndex-0.5)*binSeconds/86400;
    distance = abs(timeNum-binCenters);
    occupiedBins = unique(binIndex(:)).';
    for iBin = occupiedBins
        candidates = find(binIndex == iBin);
        [candidateDistance,localIndex] = min(distance(candidates));
        if candidateDistance < targetDistance(iBin)
            sourceIndex = candidates(localIndex);
            targetTime(iBin) = timeNum(sourceIndex);
            targetValues(iBin,:) = values(sourceIndex,:);
            targetDistance(iBin) = candidateDistance;
        end
    end
end

function [timeNum, psd, frequencyHz] = readJunoWavesSurvey(fileName)
    fid = fopen(fileName, 'r');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    header = cell(5,1);
    for iLine = 1:5
        header{iLine} = fgetl(fid);
    end
    headerFields = split(string(header{4}), ',');
    headerFields = erase(headerFields, '"');
    frequencyHz = str2double(headerFields(29:154)).';

    frewind(fid);
    formatString = ['%f%q%*q%q' repmat('%*q',1,24) repmat('%f',1,126)];
    values = textscan(fid, formatString, 'Delimiter', ',', ...
        'HeaderLines',5, 'CollectOutput',false, 'EmptyValue',NaN, ...
        'ReturnOnError',false);
    scet = values{2};
    qualityFlags = values{3};
    psd = horzcat(values{4:end});
    if isempty(scet)
        timeNum = [];
        return
    end

    utcChars = char(scet);
    years = str2double(cellstr(utcChars(:,1:4)));
    doys = str2double(cellstr(utcChars(:,6:8)));
    hours = str2double(cellstr(utcChars(:,10:11)));
    minutes = str2double(cellstr(utcChars(:,13:14)));
    seconds = str2double(cellstr(utcChars(:,16:21)));
    timeNum = datenum(years,ones(size(years)),ones(size(years)), ...
        hours,minutes,seconds) + doys - 1;
    badLfrLo = cellfun(@(x) flagIsSet(x,2), qualityFlags);
    badLfrHi = cellfun(@(x) flagIsSet(x,3), qualityFlags);
    psd(badLfrLo,1:43) = NaN;
    psd(badLfrHi,44:61) = NaN;
end

function [timeNum, magneticField] = readJunoMag(fileName)
    assert(~isempty(regexp(char(fileName), 'ss_r60s_v\d+\.sts$', 'once')), ...
        'Expected an official Sun-State/JSO 60-s MAG file: %s', fileName);
    fid = fopen(fileName, 'r');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    dataStart = -1;
    while ~feof(fid)
        lineStart = ftell(fid);
        thisLine = fgetl(fid);
        if ischar(thisLine) && ~isempty(regexp(thisLine, ...
                '^\s*20\d{2}\s+\d{1,3}\s+\d+', 'once'))
            dataStart = lineStart;
            break
        end
    end
    if dataStart < 0
        error('No numeric MAG records found in %s', fileName);
    end
    fseek(fid, dataStart, 'bof');
    values = textscan(fid, repmat('%f',1,14), ...
        'CollectOutput',true, 'ReturnOnError',false);
    values = values{1};
    if size(values,2) ~= 14
        error('Unexpected MAG row format in %s', fileName);
    end
    timeNum = datenum(values(1,1),1,1) + values(:,7) - 1;
    magneticField = values(:,8:10);
    magneticField(abs(magneticField) > 1e5) = NaN;
end

function fig = createFullOrbitFigure(timeRange, magTime, magB, ...
        waveTime, wavePsd, frequencyHz, ...
        receiverBoundaries, colorLimits, shockTimes, ...
        dateTicks, tickRadiusRj)
    fig = figure('Color','w', 'Position',[30 30 1800 1260], ...
        'Renderer','opengl');
    layout = tiledlayout(fig, 10, 1, 'TileSpacing','compact', ...
        'Padding','compact');

    % Symmetric-log is an axis mapping only; plotted values originate from
    % the BX/BY/BZ SUN-STATE fields in the official JSO source records.
    linearThreshold = 1;
    transformedB = symmetricLog(magB, linearThreshold);
    finiteMagnitude = abs(magB(:));
    finiteMagnitude = finiteMagnitude(isfinite(finiteMagnitude));
    maxMag = max(finiteMagnitude);
    plotLimit = 1.04*symmetricLog(maxMag, linearThreshold);

    axMag = nexttile(layout, 1, [3 1]);
    hold(axMag, 'on');
    hBx = plot(axMag, magTime, transformedB(:,1), ...
        'Color',[0.19 0.39 0.78], 'LineWidth',0.72);
    hBy = plot(axMag, magTime, transformedB(:,2), ...
        'Color',[0.14 0.55 0.25], 'LineWidth',0.72);
    hBz = plot(axMag, magTime, transformedB(:,3), ...
        'Color',[0.78 0.19 0.19], 'LineWidth',0.72);
    styleTimeSeriesAxis(axMag);
    xlim(axMag, timeRange);
    ylim(axMag, [-plotLimit plotLimit]);
    physicalTicks = 10.^(0:max(0,floor(log10(maxMag))));
    tickValues = [-fliplr(physicalTicks) 0 physicalTicks];
    yticks(axMag, symmetricLog(tickValues, linearThreshold));
    yticklabels(axMag, cellstr(compose('%g',tickValues)));
    ylabel(axMag, {'B_{JSO} (nT)','symmetric log'});
    magLegend = legend(axMag, [hBx hBy hBz], ...
        {'B_{x,JSO}','B_{y,JSO}','B_{z,JSO}'}, 'Location','northeast', ...
        'NumColumns',3, 'Box','off', 'FontSize',9);
    magLegend.AutoUpdate = 'off';

    axWave = nexttile(layout, 4, [6 1]);
    plotSpectrum(axWave, waveTime, wavePsd, frequencyHz, ...
        receiverBoundaries, colorLimits);
    xlim(axWave, timeRange);

    shockAxes = [axMag axWave];
    shockLabels = {
        'Shock I'
        'Shock II'
        };
    shockVerticalAlign = {'top','bottom'};
    shockHorizontalAlign = {'right','left'};
    shockColor = [0.56 0.05 0.68];
    for iShock = 1:numel(shockTimes)
        for iAxis = 1:numel(shockAxes)
            if iAxis == 1
                thisLabel = shockLabels{iShock};
            else
                thisLabel = '';
            end
            xline(shockAxes(iAxis), shockTimes(iShock), '--', thisLabel, ...
                'Color',shockColor, 'LineWidth',1.15, ...
                'LabelOrientation','horizontal', ...
                'LabelVerticalAlignment',shockVerticalAlign{iShock}, ...
                'LabelHorizontalAlignment',shockHorizontalAlign{iShock}, ...
                'FontSize',9, 'FontWeight','bold', ...
                'HandleVisibility','off');
        end
    end

    linkaxes([axMag axWave], 'x');
    set([axMag axWave], 'XTick',dateTicks);
    axMag.XTickLabel = [];
    axWave.XTickLabel = [];
    axWave.XTickLabelRotation = 0;
    axWave.FontSize = 9;
    addDateAndRadiusLabels(axWave, timeRange, dateTicks, tickRadiusRj);

    axBottomSpace = nexttile(layout, 10);
    axis(axBottomSpace, 'off');

    colormap(fig, turbo(256));
    cb = colorbar(axWave, 'eastoutside');
    cb.Label.String = 'E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 10;
end

function transformed = symmetricLog(values, threshold)
    transformed = sign(values).*log10(1 + abs(values)/threshold);
end

function addDateAndRadiusLabels(ax, timeRange, tickValues, radiusRj)
    assert(numel(tickValues) == numel(radiusRj), ...
        'Tick times and SC_POS_R values have different lengths.');
    xNormalized = (tickValues-timeRange(1))/diff(timeRange);
    for iTick = 1:numel(tickValues)
        horizontalAlignment = 'center';
        if iTick == 1
            horizontalAlignment = 'left';
        elseif iTick == numel(tickValues)
            horizontalAlignment = 'right';
        end
        text(ax,xNormalized(iTick),-0.042, ...
            datestr(tickValues(iTick),'mm/dd'), ...
            'Units','normalized','HorizontalAlignment',horizontalAlignment, ...
            'VerticalAlignment','top','FontSize',8.5, ...
            'Clipping','off','Interpreter','none');
        text(ax,xNormalized(iTick),-0.088, ...
            sprintf('R/R_J %.1f',radiusRj(iTick)), ...
            'Units','normalized','HorizontalAlignment',horizontalAlignment, ...
            'VerticalAlignment','top','FontSize',8.5, ...
            'Clipping','off','Interpreter','tex');
    end
end

function plotSpectrum(ax, timeNum, psd, frequencyHz, ...
        receiverBoundaries, colorLimits)
    validTime = isfinite(timeNum);
    timeNum = timeNum(validTime);
    psd = psd(validTime,:);
    [timeNum, timeOrder] = sort(timeNum(:));
    psd = psd(timeOrder,:);
    [frequencyHz, frequencyOrder] = sort(frequencyHz(:));
    directPsd = psd(:,frequencyOrder).';

    timeEdges = [timeNum(1)-0.5*(timeNum(2)-timeNum(1)); ...
        0.5*(timeNum(1:end-1)+timeNum(2:end)); ...
        timeNum(end)+0.5*(timeNum(end)-timeNum(end-1))];
    logFrequency = log10(frequencyHz);
    frequencyMidpoints = 0.5*(logFrequency(1:end-1)+logFrequency(2:end));
    frequencyEdges = [2*logFrequency(1)-frequencyMidpoints(1); ...
        frequencyMidpoints; ...
        2*logFrequency(end)-frequencyMidpoints(end)];
    frequencyEdges(1) = log10(50);
    frequencyEdges(end) = log10(41e6);

    [xGrid,yGrid] = meshgrid(timeEdges,frequencyEdges);
    colorData = NaN(size(xGrid));
    colorData(1:end-1,1:end-1) = directPsd;
    alphaData = zeros(size(xGrid));
    alphaData(1:end-1,1:end-1) = isfinite(directPsd) & directPsd > 0;
    h = surface(ax,xGrid,yGrid,zeros(size(xGrid)),colorData, ...
        'FaceColor','flat','EdgeColor','none');
    set(h,'AlphaData',alphaData,'FaceAlpha','flat', ...
        'AlphaDataMapping','none');
    view(ax,2);
    set(ax,'YDir','normal','Layer','top','FontSize',10, ...
        'TickDir','out','Box','on','Color',[0.88 0.88 0.88], ...
        'ColorScale','log');
    ylim(ax,log10([50 41e6]));
    yticks(ax,log10([50 1e2 1e3 1e4 1e5 1e6 1e7 4.1e7]));
    yticklabels(ax,{'50 Hz','100 Hz','1 kHz','10 kHz','100 kHz', ...
        '1 MHz','10 MHz','41 MHz'});
    ylabel(ax,'Frequency');
    caxis(ax,colorLimits);
    set(ax,'XGrid','off','YGrid','off','XMinorGrid','off', ...
        'YMinorGrid','off');
    for value = receiverBoundaries
        yline(ax,log10(value),'--k','LineWidth',0.9);
    end
end

function tf = flagIsSet(flagString, flagNumber)
    parts = strsplit(flagString, ':');
    tf = numel(parts) >= flagNumber && strcmp(parts{flagNumber}, '1');
end

function limits = robustLimits(values, lowerFraction, upperFraction)
    values = sort(values(:));
    nValue = numel(values);
    if nValue == 0
        limits = [-16 -9];
        return
    end
    lowerIndex = max(1,min(nValue,round(1+lowerFraction*(nValue-1))));
    upperIndex = max(1,min(nValue,round(1+upperFraction*(nValue-1))));
    limits = [values(lowerIndex),values(upperIndex)];
end

function styleTimeSeriesAxis(ax)
    set(ax,'FontSize',10,'TickDir','out','Box','on', ...
        'Layer','top','Color','w','XGrid','off','YGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off');
end
