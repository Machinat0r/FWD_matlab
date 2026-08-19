%% Juno current-style direct JADE JSS-RTP field and Waves plotter
% The default preset is the PJ67-to-PJ68 physical revolution containing the
% two bow-shock crossings on 2024-12-10/11.  Set environment variable
% JUNO_PLOT_PRESET=RAPTIS_20231001 to draw the event analysed by Raptis et al.
% (Nature, 2026) on 2023-10-01 between 08:00 and 20:00 UTC.
%
% Data products:
%   JNO-J/SW-JAD-5-CALIBRATED-V1.0, direct MAG_VECTOR_JSSRTP
%   JNO-E/J/SS-WAV-3-CDR-SRVFULL-V2.0, calibrated E spectra
%
% Br/Btheta/Bphi are copied directly from the JADE Level-5 binary field
% MAG_VECTOR_JSSRTP at its paired MAG_UTC.  No magnetic-coordinate
% transformation is performed.  At the user's explicit request, total field
% is calculated as |B|=sqrt(Br^2+Btheta^2+Bphi^2).  Waves uses one actual
% source record per display interval.  No plasma moment, spectral-peak
% proxy, or bin statistic is calculated.  Tick positions use the direct
% SC_POS_JUPITER_J2000XYZ and SC_POS_R fields at DIM0_UTC.  XYZ is converted
% from km to R_J with the label's fixed 1 R_J = 71492 km conversion; R is
% read directly in R_J.

clear; close all; clc;

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
databaseRoot = 'Z:\SPART-WORK\Data\Juno';
assert(isfolder(databaseRoot), 'Juno database is unavailable: %s', databaseRoot);

plotPreset = getenv('JUNO_PLOT_PRESET');
if isempty(plotPreset)
    plotPreset = 'PJ67_PJ68_2024';
end

purple = [0.56 0.05 0.68];
blue = [0.10 0.35 0.85];
switch upper(plotPreset)
    case 'PJ67_PJ68_2024'
        % Reconstructed physical revolution from PJ67 to PJ68.
        tOrbit = [datenum(2024,11,24,13,5,28), ...
                  datenum(2024,12,27, 5,22,30)];
        shockRampRanges = [
            datenum(2024,12,10,17,24,51.757), ...
                datenum(2024,12,10,17,24,54.507)
            datenum(2024,12,11, 3, 6, 1.461), ...
                datenum(2024,12,11, 3, 6, 4.211)
            ];
        eventTimes = mean(shockRampRanges,2);
        eventColors = repmat(purple,numel(eventTimes),1);
        waveBinSeconds = 600;
        waveOrbitFolderNames = {'2024313_ORBIT_67','2024345_ORBIT_68'};
        dateTicks = ceil(tOrbit(1)):2:floor(tOrbit(2));
        positionProduct = 'L3_ION_LOGICALS';
        expectedWaveFiles = 35;
        baseName = ['juno_pj67_pj68_jade_direct_jss_rtp_linear_btotal_' ...
            'waves_xyzr_20241124_20241227'];
        intervalDescription = 'Physical revolution';
    case 'RAPTIS_20231001'
        % Nature Fig. 1 overview: foreshock transients and inbound shock.
        tOrbit = [datenum(2023,10,1,8,0,0), ...
                  datenum(2023,10,1,20,0,0)];
        eventTimes = [
            datenum(2023,10,1,12,30,0)
            datenum(2023,10,1,12,50,0)
            datenum(2023,10,1,18, 5,0)
            datenum(2023,10,1,18,15,0)
            ];
        eventColors = [purple; purple; blue; blue];
        waveBinSeconds = 30;
        waveOrbitFolderNames = {'2023269_ORBIT_55'};
        dateTicks = tOrbit(1):1/24:tOrbit(2);
        positionProduct = 'L5_ELECTRONS';
        expectedWaveFiles = 1;
        baseName = ['juno_20231001_raptis_relativistic_electron_event_' ...
            'jade_direct_jss_rtp_linear_btotal_waves_xyzr_0800_2000'];
        intervalDescription = 'Nature event interval';
    case 'RAPTIS_20231001_BOWSHOCK'
        % Bow-shock interval with one hour of context on each side.
        tOrbit = [datenum(2023,10,1,17,5,0), ...
                  datenum(2023,10,1,19,15,0)];
        eventTimes = zeros(0,1);
        eventColors = zeros(0,3);
        waveBinSeconds = 30;
        waveOrbitFolderNames = {'2023269_ORBIT_55'};
        dateTicks = tOrbit(1):10/1440:tOrbit(2);
        positionProduct = 'L5_ELECTRONS';
        expectedWaveFiles = 1;
        baseName = ['juno_20231001_raptis_bowshock_crossing_context_' ...
            'jade_direct_jss_rtp_linear_btotal_waves_xyzr_1705_1915'];
        intervalDescription = 'Nature bow-shock crossing context interval';
    otherwise
        error('Unknown JUNO_PLOT_PRESET: %s', plotPreset);
end
assert(all(eventTimes >= tOrbit(1) & eventTimes <= tOrbit(2)), ...
    'A verified event marker lies outside the plotting interval.');

waveFiles = discoverWaveFiles(databaseRoot, tOrbit, waveOrbitFolderNames);
jadeMagFiles = discoverJadeMagFiles(databaseRoot, tOrbit);
if strcmp(positionProduct,'L5_ELECTRONS')
    jadePositionFiles = jadeMagFiles;
    positionProductDescription = 'Level-5 JADE electron DEF';
else
    jadePositionFiles = discoverJadePositionFiles(databaseRoot, dateTicks);
    positionProductDescription = 'Level-3 JADE ION_LOGICALS';
end

expectedDailyFiles = floor(tOrbit(2)-eps)-floor(tOrbit(1))+1;
assert(numel(waveFiles) == expectedWaveFiles, ...
    'Expected %d Waves E files; found %d.', ...
    expectedWaveFiles, numel(waveFiles));
assert(numel(jadeMagFiles) == expectedDailyFiles, ...
    'Expected %d JADE direct-RTP daily files; found %d.', ...
    expectedDailyFiles, numel(jadeMagFiles));
assert(~isempty(jadePositionFiles), 'No direct JADE position source files found.');

fprintf('Selecting direct source records from %d Waves files at %d-s spacing ...\n', ...
    numel(waveFiles), waveBinSeconds);
[waveTime, wavePsd, frequencyHz, receiverBoundaries] = ...
    processWaves(waveFiles, tOrbit, waveBinSeconds);

fprintf(['Reading direct MAG_VECTOR_JSSRTP records from %d JADE Level-5 ' ...
    'DEF files ...\n'], numel(jadeMagFiles));
[magTime, magSpherical] = processJadeMagRtp(jadeMagFiles, tOrbit);
magTotal = sqrt(sum(magSpherical.^2,2));

fprintf(['Reading direct J2000 XYZ and SC_POS_R records from %s ' ...
    'for %d x ticks ...\n'], positionProductDescription, numel(dateTicks));
[positionSourceTime, tickPositionRj] = ...
    readJadePositionAtTicks(jadePositionFiles, dateTicks, tOrbit);

finiteWave = wavePsd(isfinite(wavePsd) & wavePsd > 0);
logColorLimits = robustLimits(log10(finiteWave), 0.005, 0.997);
logColorLimits = [floor(2*logColorLimits(1))/2, ...
    ceil(2*logColorLimits(2))/2];
if diff(logColorLimits) < 5
    centerValue = mean(logColorLimits);
    logColorLimits = centerValue + [-2.5 2.5];
end
colorLimits = 10.^logColorLimits;

fig = createFullOrbitFigure(tOrbit, magTime, magSpherical, magTotal, ...
    waveTime, wavePsd, frequencyHz, ...
    receiverBoundaries, colorLimits, eventTimes, eventColors, ...
    dateTicks, tickPositionRj);

pngFile = fullfile(outputDir, [baseName '.png']);
pdfFile = fullfile(outputDir, [baseName '.pdf']);
exportgraphics(fig, pngFile, 'Resolution',260);
exportgraphics(fig, pdfFile, 'ContentType','image', 'Resolution',260);

fprintf('Saved %s\n', pngFile);
fprintf('Saved %s\n', pdfFile);
fprintf('%s: %s to %s UTC (%.3f days)\n', intervalDescription, ...
    datestr(tOrbit(1),'yyyy-mm-dd HH:MM:SS'), ...
    datestr(tOrbit(2),'yyyy-mm-dd HH:MM:SS'), diff(tOrbit));
fprintf('Files used: Waves=%d, JADE RTP=%d, JADE position=%d\n', ...
    numel(waveFiles), numel(jadeMagFiles), numel(jadePositionFiles));
fprintf('Direct source samples: Waves=%d/%d, JADE RTP=%d/%d\n', ...
    nnz(any(isfinite(wavePsd),2)), size(wavePsd,1), ...
    nnz(any(isfinite(magSpherical),2)), size(magSpherical,1));
fprintf('Total-field samples derived from direct RTP components: %d\n', ...
    nnz(isfinite(magTotal)));
fprintf('Frequency centers: %.4g Hz to %.4g Hz (%d channels)\n', ...
    min(frequencyHz), max(frequencyHz), numel(frequencyHz));
fprintf('Receiver separators: %.7g, %.7g, %.7g Hz\n', receiverBoundaries);
fprintf('Color limits: %.4g to %.4g (V/m)^2/Hz (log color scale)\n', ...
    colorLimits);
if isempty(eventTimes)
    fprintf('Event markers: none\n');
else
    fprintf('Event markers:\n');
    for iEvent = 1:numel(eventTimes)
        fprintf('  %s UTC\n', datestr(eventTimes(iEvent), ...
            'yyyy-mm-dd HH:MM:SS.FFF'));
    end
end
fprintf(['Tick positions use direct %s SC_POS_JUPITER_J2000XYZ and ' ...
    'SC_POS_R records:\n'], positionProductDescription);
for iTick = 1:numel(dateTicks)
    fprintf(['  %s UTC: X=%+.3f, Y=%+.3f, Z=%+.3f, R=%.3f R_J ' ...
        '(%s source)\n'], ...
        datestr(dateTicks(iTick),'yyyy-mm-dd HH:MM'), tickPositionRj(iTick,:), ...
        datestr(positionSourceTime(iTick),'HH:MM:SS.FFF'));
end

%% Local functions
function files = discoverWaveFiles(databaseRoot, timeRange, folderNames)
    orbitFolders = cellfun(@(name) fullfile(databaseRoot, ...
        'Waves_data',name), folderNames(:), 'UniformOutput',false);
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

function files = discoverJadeMagFiles(databaseRoot, timeRange)
    firstDay = floor(timeRange(1));
    lastDay = floor(timeRange(2)-eps);
    files = cell(lastDay-firstDay+1,1);
    for iDay = 1:numel(files)
        dayNumber = firstDay+iDay-1;
        dateParts = datevec(dayNumber);
        yearValue = dateParts(1);
        doyValue = floor(dayNumber-datenum(yearValue,1,1))+1;
        dateCode = sprintf('%04d%03d', yearValue, doyValue);
        folder = fullfile(databaseRoot, 'JADE_data', ...
            sprintf('%04d',yearValue), dateCode, 'ELECTRONS');
        entries = dir(fullfile(folder, ...
            sprintf('JAD_L50_LRS_ELC_ANY_DEF_%s_V*.DAT',dateCode)));
        assert(~isempty(entries), ...
            'Missing JADE MAG_VECTOR_JSSRTP file for %s.', dateCode);
        versions = NaN(numel(entries),1);
        for iEntry = 1:numel(entries)
            token = regexp(entries(iEntry).name, '_V(\d+)\.DAT$', ...
                'tokens','once');
            assert(~isempty(token), 'Unexpected JADE filename: %s', ...
                entries(iEntry).name);
            versions(iEntry) = str2double(token{1});
        end
        [~,bestIndex] = max(versions);
        files{iDay} = fullfile(entries(bestIndex).folder, ...
            entries(bestIndex).name);
    end
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
            'Missing JADE direct XYZ/R source file for %s.', dateCode);
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

function [timeNum, magneticFieldRtp] = processJadeMagRtp(files, timeRange)
    timeNum = [];
    magneticFieldRtp = [];
    for iFile = 1:numel(files)
        fprintf('  JADE RTP %2d/%2d: %s\n', ...
            iFile, numel(files), files{iFile});
        [thisTime, thisField] = readJadeMagRtp(files{iFile});
        use = thisTime >= timeRange(1) & thisTime < timeRange(2);
        timeNum = [timeNum; thisTime(use)]; %#ok<AGROW>
        magneticFieldRtp = [magneticFieldRtp; thisField(use,:)]; %#ok<AGROW>
    end
    [timeNum, order] = sort(timeNum);
    magneticFieldRtp = magneticFieldRtp(order,:);
    [timeNum, uniqueIndex] = unique(timeNum, 'last');
    magneticFieldRtp = magneticFieldRtp(uniqueIndex,:);
    gapIndex = find(diff(timeNum)*86400 > 300);
    if ~isempty(gapIndex)
        gapTime = 0.5*(timeNum(gapIndex)+timeNum(gapIndex+1));
        timeNum = [timeNum; gapTime];
        magneticFieldRtp = [magneticFieldRtp; ...
            NaN(numel(gapTime),size(magneticFieldRtp,2))];
        [timeNum, order] = sort(timeNum);
        magneticFieldRtp = magneticFieldRtp(order,:);
    end
end

function [sourceTime, positionRj] = readJadePositionAtTicks( ...
        files, tickValues, timeRange)
    allTime = [];
    allPositionRj = [];
    for iFile = 1:numel(files)
        fprintf('  JADE position %2d/%2d: %s\n', ...
            iFile, numel(files), files{iFile});
        [thisTime, thisXyzKm, thisRadiusRj] = ...
            readJadeScPosition(files{iFile});
        thisPositionRj = [thisXyzKm/71492, thisRadiusRj];
        valid = isfinite(thisTime) & all(isfinite(thisPositionRj),2) & ...
            thisRadiusRj >= 0 & thisRadiusRj < 65535 & ...
            thisTime >= timeRange(1) & thisTime < timeRange(2);
        allTime = [allTime; thisTime(valid)]; %#ok<AGROW>
        allPositionRj = [allPositionRj; thisPositionRj(valid,:)]; %#ok<AGROW>
    end
    assert(~isempty(allTime), 'No valid direct position records in the interval.');
    [allTime, order] = sort(allTime);
    allPositionRj = allPositionRj(order,:);
    [allTime, uniqueIndex] = unique(allTime, 'last');
    allPositionRj = allPositionRj(uniqueIndex,:);
    sourceTime = NaN(numel(tickValues),1);
    positionRj = NaN(numel(tickValues),4);
    for iTick = 1:numel(tickValues)
        [distance,sourceIndex] = min(abs(allTime-tickValues(iTick)));
        assert(distance <= 10/1440, ...
            'No JADE direct position record within 10 min of tick %s.', ...
            datestr(tickValues(iTick),'yyyy-mm-dd HH:MM:SS'));
        sourceTime(iTick) = allTime(sourceIndex);
        positionRj(iTick,:) = allPositionRj(sourceIndex,:);
    end
end

function [timeNum, xyzKm, radiusRj] = readJadeScPosition(fileName)
    labelFile = regexprep(fileName, '\.DAT$', '.LBL', 'ignorecase');
    assert(isfile(labelFile), 'Missing JADE label: %s', labelFile);
    labelText = fileread(labelFile);
    recordBytes = readPdsLabelInteger(labelText, 'RECORD_BYTES');
    nRecords = readPdsLabelInteger(labelText, 'FILE_RECORDS');
    fileInfo = dir(fileName);
    assert(fileInfo.bytes == recordBytes*nRecords, ...
        'JADE DAT size does not match label dimensions: %s', fileName);
    utcStartByte = readPdsColumnStartByte(labelText, 'DIM0_UTC');
    radiusStartByte = readPdsColumnStartByte(labelText, 'SC_POS_R');
    xyzStartByte = readPdsColumnStartByte(labelText, ...
        'SC_POS_JUPITER_J2000XYZ');

    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    assert(fseek(fid,utcStartByte-1,'bof') == 0, ...
        'Unable to seek to DIM0_UTC in %s', fileName);
    utcBytes = fread(fid,[21 nRecords],'21*uint8=>uint8',recordBytes-21);
    assert(size(utcBytes,2) == nRecords, ...
        'Incomplete DIM0_UTC records in %s', fileName);
    assert(fseek(fid,radiusStartByte-1,'bof') == 0, ...
        'Unable to seek to SC_POS_R in %s', fileName);
    radiusRj = fread(fid,[1 nRecords],'1*single=>double',recordBytes-4, ...
        'ieee-le').';
    assert(numel(radiusRj) == nRecords, ...
        'Incomplete SC_POS_R records in %s', fileName);
    assert(fseek(fid,xyzStartByte-1,'bof') == 0, ...
        'Unable to seek to SC_POS_JUPITER_J2000XYZ in %s', fileName);
    xyzKm = fread(fid,[3 nRecords],'3*single=>double',recordBytes-12, ...
        'ieee-le').';
    assert(size(xyzKm,1) == nRecords, ...
        'Incomplete SC_POS_JUPITER_J2000XYZ records in %s', fileName);

    timeNum = parsePdsDoyUtc(char(utcBytes.'));
    invalid = ~isfinite(radiusRj) | radiusRj == 65535 | ...
        any(~isfinite(xyzKm),2) | any(xyzKm == 65535,2) | ...
        any(abs(xyzKm) > 10008880,2);
    radiusRj(invalid) = NaN;
    xyzKm(invalid,:) = NaN;
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
    badHfrHi = cellfun(@(x) flagIsSet(x,1), qualityFlags);
    psd(badLfrLo,1:43) = NaN;
    psd(badLfrHi,44:61) = NaN;
    psd(badHfrHi,89:126) = NaN;
end

function [timeNum, magneticFieldRtp] = readJadeMagRtp(fileName)
    assert(~isempty(regexp(char(fileName), ...
        'JAD_L50_LRS_ELC_ANY_DEF_\d{7}_V\d+\.DAT$', 'once')), ...
        'Expected a JADE Level-5 LRS electron DEF file: %s', fileName);
    labelFile = regexprep(fileName, '\.DAT$', '.LBL', 'ignorecase');
    assert(isfile(labelFile), 'Missing JADE label: %s', labelFile);
    labelText = fileread(labelFile);
    recordBytes = readPdsLabelInteger(labelText, 'RECORD_BYTES');
    nRecords = readPdsLabelInteger(labelText, 'FILE_RECORDS');
    fileInfo = dir(fileName);
    assert(fileInfo.bytes == recordBytes*nRecords, ...
        'JADE DAT size does not match label dimensions: %s', fileName);
    utcStartByte = readPdsColumnStartByte(labelText, 'MAG_UTC');
    rtpStartByte = readPdsColumnStartByte(labelText, 'MAG_VECTOR_JSSRTP');
    assert(utcStartByte + 21 - 1 <= recordBytes, ...
        'MAG_UTC lies outside a JADE record in %s', labelFile);
    assert(rtpStartByte + 12 - 1 <= recordBytes, ...
        'MAG_VECTOR_JSSRTP lies outside a JADE record in %s', labelFile);

    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    status = fseek(fid, utcStartByte-1, 'bof');
    assert(status == 0, 'Unable to seek to MAG_UTC in %s', fileName);
    utcBytes = fread(fid, [21 nRecords], '21*uint8=>uint8', ...
        recordBytes-21);
    assert(size(utcBytes,2) == nRecords, ...
        'Incomplete MAG_UTC records in %s', fileName);

    status = fseek(fid, rtpStartByte-1, 'bof');
    assert(status == 0, ...
        'Unable to seek to MAG_VECTOR_JSSRTP in %s', fileName);
    rtpValues = fread(fid, [3 nRecords], '3*single=>double', ...
        recordBytes-12, 'ieee-le');
    assert(size(rtpValues,2) == nRecords, ...
        'Incomplete MAG_VECTOR_JSSRTP records in %s', fileName);

    utcChars = char(utcBytes.');
    timeNum = parsePdsDoyUtc(utcChars);
    magneticFieldRtp = rtpValues.';
    invalidRows = any(~isfinite(magneticFieldRtp),2) | ...
        any(magneticFieldRtp == 9990000,2) | ...
        any(abs(magneticFieldRtp) > 1.6e6,2);
    magneticFieldRtp(invalidRows,:) = NaN;
end

function value = readPdsLabelInteger(labelText, keyword)
    token = regexp(labelText, ...
        ['(?m)^\s*' regexptranslate('escape',keyword) ...
         '\s*=\s*(\d+)'], 'tokens', 'once');
    assert(~isempty(token), 'Missing %s in PDS label.', keyword);
    value = str2double(token{1});
end

function startByte = readPdsColumnStartByte(labelText, columnName)
    expression = ['(?s)OBJECT\s*=\s*COLUMN\s+' ...
        'NAME\s*=\s*"?' regexptranslate('escape',columnName) ...
        '"?.*?START_BYTE\s*=\s*(\d+).*?' ...
        'END_OBJECT\s*=\s*COLUMN'];
    token = regexp(labelText, expression, 'tokens', 'once');
    assert(~isempty(token), ...
        'Missing PDS column %s or START_BYTE.', columnName);
    startByte = str2double(token{1});
end

function timeNum = parsePdsDoyUtc(utcChars)
    years = str2double(cellstr(utcChars(:,1:4)));
    doys = str2double(cellstr(utcChars(:,6:8)));
    hours = str2double(cellstr(utcChars(:,10:11)));
    minutes = str2double(cellstr(utcChars(:,13:14)));
    seconds = str2double(cellstr(utcChars(:,16:21)));
    timeNum = datenum(years,ones(size(years)),ones(size(years)), ...
        hours,minutes,seconds) + doys - 1;
end

function fig = createFullOrbitFigure(timeRange, magTime, magB, magTotal, ...
        waveTime, wavePsd, frequencyHz, ...
        receiverBoundaries, colorLimits, eventTimes, eventColors, ...
        dateTicks, tickPositionRj)
    fig = figure('Color','w', 'Position',[30 30 1800 1360], ...
        'Renderer','opengl');
    layout = tiledlayout(fig, 13, 1, 'TileSpacing','compact', ...
        'Padding','compact');

    finiteField = [magB(:); magTotal(:)];
    finiteField = finiteField(isfinite(finiteField));
    fieldLimits = [min(finiteField) max(finiteField)];
    fieldPadding = max(1,0.04*diff(fieldLimits));

    axMag = nexttile(layout, 1, [3 1]);
    hold(axMag, 'on');
    hBr = plot(axMag, magTime, magB(:,1), ...
        'Color',[0.19 0.39 0.78], 'LineWidth',0.72);
    hBtheta = plot(axMag, magTime, magB(:,2), ...
        'Color',[0.95 0.48 0.12], 'LineWidth',0.72);
    hBphi = plot(axMag, magTime, magB(:,3), ...
        'Color',[0.14 0.55 0.25], 'LineWidth',0.72);
    hBtotal = plot(axMag, magTime, magTotal, ...
        'Color',[0 0 0], 'LineWidth',0.9);
    styleTimeSeriesAxis(axMag);
    xlim(axMag, timeRange);
    ylim(axMag,fieldLimits + [-fieldPadding fieldPadding]);
    ylabel(axMag, 'B^{JSS}_{r,\theta,\phi} (nT)');
    magLegend = legend(axMag, [hBr hBtheta hBphi hBtotal], ...
        {'B_r','B_\theta','B_\phi','|B|'}, 'Location','northeast', ...
        'NumColumns',4, 'Box','off', 'FontSize',9);
    magLegend.AutoUpdate = 'off';

    axWave = nexttile(layout, 4, [8 1]);
    plotSpectrum(axWave, waveTime, wavePsd, frequencyHz, ...
        receiverBoundaries, colorLimits);
    xlim(axWave, timeRange);

    eventAxes = [axMag axWave];
    assert(size(eventColors,1) == numel(eventTimes) && ...
        size(eventColors,2) == 3, ...
        'Every event marker must have one RGB color row.');
    for iEvent = 1:numel(eventTimes)
        for iAxis = 1:numel(eventAxes)
            xline(eventAxes(iAxis), eventTimes(iEvent), '--', ...
                'Color',eventColors(iEvent,:), 'LineWidth',1.15, ...
                'HandleVisibility','off');
        end
    end

    linkaxes([axMag axWave], 'x');
    set([axMag axWave], 'XTick',dateTicks);
    axMag.XTickLabel = [];
    axWave.XTickLabel = [];
    axWave.XTickLabelRotation = 0;
    axWave.FontSize = 9;
    addDateAndPositionLabels(axWave, timeRange, dateTicks, tickPositionRj);

    axBottomSpace = nexttile(layout, 12, [2 1]);
    axis(axBottomSpace, 'off');

    colormap(fig, turbo(256));
    cb = colorbar(axWave, 'eastoutside');
    cb.Label.String = 'E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 10;
end

function addDateAndPositionLabels(ax, timeRange, tickValues, positionRj)
    assert(size(positionRj,1) == numel(tickValues) && ...
        size(positionRj,2) == 4, ...
        'Tick times and direct XYZ/R positions have different dimensions.');
    xNormalized = (tickValues-timeRange(1))/diff(timeRange);
    rowY = [-0.076 -0.116 -0.156 -0.196];
    rowLabels = {'X','Y','Z','R'};
    if diff(timeRange) <= 2
        tickLabelFormat = 'HH:MM';
    else
        tickLabelFormat = 'mm/dd';
    end
    text(ax,-0.010,-0.036,'(R_J)', ...
        'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','top','FontSize',8.5, ...
        'Clipping','off','Interpreter','tex');
    for iRow = 1:4
        text(ax,-0.010,rowY(iRow),rowLabels{iRow}, ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','top','FontSize',8.5, ...
            'Clipping','off','Interpreter','tex');
    end
    for iTick = 1:numel(tickValues)
        horizontalAlignment = 'center';
        if iTick == 1
            horizontalAlignment = 'left';
        elseif iTick == numel(tickValues)
            horizontalAlignment = 'right';
        end
        text(ax,xNormalized(iTick),-0.036, ...
            datestr(tickValues(iTick),tickLabelFormat), ...
            'Units','normalized','HorizontalAlignment',horizontalAlignment, ...
            'VerticalAlignment','top','FontSize',8.5, ...
            'Clipping','off','Interpreter','none');
        for iRow = 1:4
            if iRow < 4
                valueText = sprintf('%+.1f',positionRj(iTick,iRow));
            else
                valueText = sprintf('%.1f',positionRj(iTick,iRow));
            end
            text(ax,xNormalized(iTick),rowY(iRow),valueText, ...
                'Units','normalized', ...
                'HorizontalAlignment',horizontalAlignment, ...
                'VerticalAlignment','top','FontSize',8.5, ...
                'Clipping','off','Interpreter','none');
        end
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
