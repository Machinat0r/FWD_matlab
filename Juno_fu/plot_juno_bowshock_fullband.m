%% Juno bow-shock overview: JADE ions, MAG, and full-band Waves spectra
% Data products: JNO-J/SW-JAD-3-CALIBRATED-V1.0,
% JNO-J-3-FGM-CAL-V1.0, and JNO-E/J/SS-WAV-3-CDR-SRVFULL-V2.0.
% Event times: Case I ~2024-12-10 17:25 UTC; Case II ~2024-12-11 03:05 UTC.
% The top survey-channel center is 40.5 MHz; the nominal instrument band
% extends to 41 MHz.  Spectral densities are already calibrated by PDS.
% JADE does not provide a public three-dimensional ion moment for this 2024
% interval.  The speed panel therefore reports the proton-equivalent speed
% derived from the peak energy of the omnidirectional All Stops spectrum.

clear; close all; clc;

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
workspaceDir = fileparts(outputDir);
dataDir = fullfile(workspaceDir, 'work', 'juno_waves_20241210_11');
if ~isfolder(dataDir)
    mkdir(dataDir);
end

baseUrl = ['https://pds-ppi.igpp.ucla.edu/data/' ...
    'JNO-E_J_SS-WAV-3-CDR-SRVFULL-V2.0/'];
relPaths = {
    'DATA/WAVES_SURVEY/2024313_ORBIT_67/WAV_2024345T000000_E_V01.CSV'
    'DATA/WAVES_SURVEY/2024345_ORBIT_68/WAV_2024345T211242_E_V01.CSV'
    'DATA/WAVES_SURVEY/2024345_ORBIT_68/WAV_2024346T000000_E_V01.CSV'
    };

files = strings(numel(relPaths), 1);
for iFile = 1:numel(relPaths)
    [~, name, ext] = fileparts(relPaths{iFile});
    files(iFile) = fullfile(dataDir, [name ext]);
    if ~isfile(files(iFile))
        fprintf('Downloading %s ...\n', relPaths{iFile});
        websave(files(iFile), [baseUrl relPaths{iFile}]);
    end
end

magFiles = strings(2,1);
magFiles(1) = fullfile(dataDir, 'fgm_jno_l3_2024345pc_r1s_v01.sts');
magFiles(2) = fullfile(dataDir, 'fgm_jno_l3_2024346pc_r1s_v01.sts');
jadeFiles = strings(2,1);
jadeFiles(1) = fullfile(dataDir, ...
    'JAD_L30_HLS_ION_LOG_CNT_2024345_V04_subset_1200_2400.DAT');
jadeFiles(2) = fullfile(dataDir, ...
    'JAD_L30_HLS_ION_LOG_CNT_2024346_V04_subset_0000_0702.DAT');
assert(all(isfile(magFiles)), 'The two local 1-s MAG files are required.');
assert(all(isfile(jadeFiles)), ['The two JADE subset files are required. ' ...
    'Run work/download_jade_hls_subset.ps1 first.']);

tOverview = [datenum(2024,12,10,12,0,0), datenum(2024,12,11,7,0,0)];
tShock = [datenum(2024,12,10,17,25,0), datenum(2024,12,11,3,5,0)];

allTime = [];
allPsd = [];
frequencyHz = [];
for iFile = 1:numel(files)
    fprintf('Reading %s ...\n', files(iFile));
    [timeNum, psd, thisFrequency] = readJunoWavesSurvey(files(iFile));
    if isempty(frequencyHz)
        frequencyHz = thisFrequency;
    elseif any(abs(thisFrequency - frequencyHz) > max(1e-9, 1e-10*frequencyHz))
        error('Frequency grids do not match among the input files.');
    end

    use = timeNum >= tOverview(1) & timeNum <= tOverview(2);
    allTime = [allTime; timeNum(use)]; %#ok<AGROW>
    allPsd = [allPsd; psd(use,:)]; %#ok<AGROW>
end

% Sort by time, remove exact duplicated records at the orbit boundary, and
% sort the small receiver-band frequency overlaps into monotonic order.
[allTime, order] = sort(allTime);
allPsd = allPsd(order,:);
[allTime, uniqueIndex] = unique(allTime, 'last');
allPsd = allPsd(uniqueIndex,:);
[frequencyHz, frequencyOrder] = sort(frequencyHz);
allPsd = allPsd(:,frequencyOrder);

validFrequency = isfinite(frequencyHz) & frequencyHz > 0;
frequencyHz = frequencyHz(validFrequency);
allPsd = allPsd(:,validFrequency);
allPsd(allPsd <= 0) = NaN;

% Read 1-s Planetocentric magnetic-field components.
magTime = [];
magB = [];
for iFile = 1:numel(magFiles)
    fprintf('Reading %s ...\n', magFiles(iFile));
    [thisTime, thisB] = readJunoMag1s(magFiles(iFile));
    use = thisTime >= tOverview(1) & thisTime <= tOverview(2);
    magTime = [magTime; thisTime(use)]; %#ok<AGROW>
    magB = [magB; thisB(use,:)]; %#ok<AGROW>
end
[magTime, order] = sort(magTime);
magB = magB(order,:);
[magTime, uniqueIndex] = unique(magTime, 'last');
magB = magB(uniqueIndex,:);
magMagnitude = sqrt(sum(magB.^2, 2));

% Read the 2-min JADE ion-logical spectra and infer the energy of their peak.
jadeTime = [];
jadeSpeed = [];
jadePeakEnergy = [];
jadePeakRate = [];
for iFile = 1:numel(jadeFiles)
    fprintf('Reading %s ...\n', jadeFiles(iFile));
    [thisTime, thisSpeed, thisEnergy, thisRate] = ...
        readJadeIonPeakSpeed(jadeFiles(iFile));
    use = thisTime >= tOverview(1) & thisTime <= tOverview(2);
    jadeTime = [jadeTime; thisTime(use)]; %#ok<AGROW>
    jadeSpeed = [jadeSpeed; thisSpeed(use)]; %#ok<AGROW>
    jadePeakEnergy = [jadePeakEnergy; thisEnergy(use)]; %#ok<AGROW>
    jadePeakRate = [jadePeakRate; thisRate(use)]; %#ok<AGROW>
end
[jadeTime, order] = sort(jadeTime);
jadeSpeed = jadeSpeed(order);
jadePeakEnergy = jadePeakEnergy(order);
jadePeakRate = jadePeakRate(order);

% A common robust color range keeps the three panels quantitatively comparable.
sampleValues = log10(allPsd(1:20:end,:));
sampleValues = sampleValues(isfinite(sampleValues));
colorLimits = robustLimits(sampleValues, 0.005, 0.997);
colorLimits = [floor(2*colorLimits(1))/2, ceil(2*colorLimits(2))/2];
if diff(colorLimits) < 5
    middleValue = mean(colorLimits);
    colorLimits = middleValue + [-2.5 2.5];
end

% Preserve short enhancements in the long overview by taking a 10-s maximum.
[overviewTime, overviewPsd] = binMaximum(allTime, allPsd, 10);

validSpeed = jadeSpeed(isfinite(jadeSpeed));
speedLimits = [floor((min(validSpeed)-20)/50)*50, ...
    ceil((max(validSpeed)+20)/50)*50];
speedLimits(1) = max(0, speedLimits(1));
if diff(speedLimits) < 150
    speedLimits = mean(speedLimits) + [-75 75];
end
magValues = abs([magB(:); magMagnitude]);
magValues = magValues(isfinite(magValues));
magLimit = max(2, ceil(1.05*max(magValues)));

zoomHalfWidth = 30/1440; % 30 minutes in datenum units
use1 = abs(allTime - tShock(1)) <= zoomHalfWidth;
use2 = abs(allTime - tShock(2)) <= zoomHalfWidth;
caseRange1 = tShock(1) + [-1 1]*zoomHalfWidth;
caseRange2 = tShock(2) + [-1 1]*zoomHalfWidth;

figures = gobjects(3,1);
figures(1) = createThreePanelFigure(tOverview, ...
    jadeTime, jadeSpeed, magTime, magB, magMagnitude, ...
    overviewTime, overviewPsd, frequencyHz, colorLimits, ...
    tShock, {'Case I  17:25','Case II  03:05'}, ...
    'Juno at the Jovian bow shock: overview (2024-12-10/11 UTC)', ...
    'Juno/Waves electric-field survey, full interval', ...
    'mm/dd HH:MM', speedLimits, magLimit);

figures(2) = createThreePanelFigure(caseRange1, ...
    jadeTime, jadeSpeed, magTime, magB, magMagnitude, ...
    allTime(use1), allPsd(use1,:), frequencyHz, colorLimits, ...
    tShock(1), {'Bow shock  17:25'}, ...
    'Case I: outbound, magnetosheath \rightarrow solar wind (2024-12-10 UTC)', ...
    'Juno/Waves full-band spectrum, Case I', ...
    'HH:MM', speedLimits, magLimit);

figures(3) = createThreePanelFigure(caseRange2, ...
    jadeTime, jadeSpeed, magTime, magB, magMagnitude, ...
    allTime(use2), allPsd(use2,:), frequencyHz, colorLimits, ...
    tShock(2), {'Bow shock  03:05'}, ...
    'Case II: inbound, solar wind \rightarrow magnetosheath (2024-12-11 UTC)', ...
    'Juno/Waves full-band spectrum, Case II', ...
    'HH:MM', speedLimits, magLimit);

baseNames = {
    'juno_bowshock_overview_mag_speed_waves_20241210_11'
    'juno_bowshock_case1_mag_speed_waves_20241210'
    'juno_bowshock_case2_mag_speed_waves_20241211'
    };
for iFigure = 1:numel(figures)
    pngFile = fullfile(outputDir, [baseNames{iFigure} '.png']);
    pdfFile = fullfile(outputDir, [baseNames{iFigure} '.pdf']);
    exportgraphics(figures(iFigure), pngFile, 'Resolution',260);
    exportgraphics(figures(iFigure), pdfFile, 'ContentType','image', ...
        'Resolution',260);
    fprintf('Saved %s\n', pngFile);
    fprintf('Saved %s\n', pdfFile);
end

fprintf('Rows retained: %d\n', numel(allTime));
fprintf('MAG samples: %d; JADE speed samples: %d (%d finite)\n', ...
    numel(magTime), numel(jadeTime), nnz(isfinite(jadeSpeed)));
fprintf('Frequency centers: %.4g Hz to %.4g Hz (%d channels)\n', ...
    min(frequencyHz), max(frequencyHz), numel(frequencyHz));
fprintf('Color limits: %.2f to %.2f log10[(V/m)^2/Hz]\n', colorLimits);

%% Local functions
function fig = createThreePanelFigure(timeRange, ...
        jadeTime, jadeSpeed, magTime, magB, magMagnitude, ...
        waveTime, wavePsd, frequencyHz, colorLimits, ...
        shockTimes, shockLabels, figureTitle, spectrumTitle, ...
        timeFormat, speedLimits, magLimit)
    fig = figure('Color','w', 'Position',[40 40 1800 1500], ...
        'Renderer','opengl');
    layout = tiledlayout(fig, 4, 1, 'TileSpacing','compact', ...
        'Padding','compact');

    useSpeed = jadeTime >= timeRange(1) & jadeTime <= timeRange(2);
    axSpeed = nexttile(layout, 1);
    plot(axSpeed, jadeTime(useSpeed), jadeSpeed(useSpeed), '-', ...
        'Color',[0.15 0.39 0.72], 'LineWidth',1.15, ...
        'Marker','.', 'MarkerSize',8);
    styleTimeSeriesAxis(axSpeed);
    xlim(axSpeed, timeRange);
    ylim(axSpeed, speedLimits);
    ylabel(axSpeed, {'v_{p,peak}', '(km s^{-1})'});
    title(axSpeed, ['JADE ion spectrum: proton-equivalent peak speed ' ...
        '(80 eV/q to 5 keV/q)'], 'FontWeight','normal');
    datetick(axSpeed, 'x', timeFormat, 'keeplimits', 'keepticks');
    axSpeed.XTickLabel = [];

    useMag = magTime >= timeRange(1) & magTime <= timeRange(2);
    axMag = nexttile(layout, 2);
    hold(axMag, 'on');
    hBx = plot(axMag, magTime(useMag), magB(useMag,1), ...
        'Color',[0.78 0.19 0.19], 'LineWidth',0.75);
    hBy = plot(axMag, magTime(useMag), magB(useMag,2), ...
        'Color',[0.14 0.55 0.25], 'LineWidth',0.75);
    hBz = plot(axMag, magTime(useMag), magB(useMag,3), ...
        'Color',[0.19 0.39 0.78], 'LineWidth',0.75);
    hBt = plot(axMag, magTime(useMag), magMagnitude(useMag), ...
        'k', 'LineWidth',1.05);
    styleTimeSeriesAxis(axMag);
    xlim(axMag, timeRange);
    ylim(axMag, [-magLimit magLimit]);
    ylabel(axMag, 'B (nT)');
    title(axMag, 'Juno/MAG Planetocentric magnetic field (1 s)', ...
        'FontWeight','normal');
    magLegend = legend(axMag, [hBx hBy hBz hBt], ...
        {'B_x','B_y','B_z','|B|'}, 'Location','northeast', ...
        'NumColumns',4, 'Box','off', 'FontSize',9);
    magLegend.AutoUpdate = 'off';
    datetick(axMag, 'x', timeFormat, 'keeplimits', 'keepticks');
    axMag.XTickLabel = [];

    axWave = nexttile(layout, 3, [2 1]);
    plotSpectrum(axWave, waveTime, wavePsd, frequencyHz, colorLimits);
    xlim(axWave, timeRange);
    datetick(axWave, 'x', timeFormat, 'keeplimits', 'keepticks');
    title(axWave, spectrumTitle, 'FontWeight','normal');

    for iShock = 1:numel(shockTimes)
        addShockLineDark(axSpeed, shockTimes(iShock), shockLabels{iShock});
        addShockLineDark(axMag, shockTimes(iShock), '');
        addShockLine(axWave, shockTimes(iShock), shockLabels{iShock});
    end

    linkaxes([axSpeed axMag axWave], 'x');
    colormap(fig, turbo(256));
    cb = colorbar(axWave, 'eastoutside');
    cb.Layout.Tile = 'east';
    cb.Label.String = 'log_{10} E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 10;

    sgtitle(fig, figureTitle, 'FontWeight','bold', 'FontSize',14);
    annotation(fig, 'textbox', [0.055 0.001 0.89 0.024], ...
        'String', sprintf(['JADE L30 HLS ION LOG CNT V04 (2-min); ' ...
        'MAG PC 1-s; Waves SRVFULL V2.0. Highest channel center: %.1f MHz; ' ...
        'nominal upper limit: 41 MHz. Time: UTC.'], max(frequencyHz)/1e6), ...
        'EdgeColor','none', 'HorizontalAlignment','center', ...
        'FontSize',9, 'Color',[0.25 0.25 0.25]);
end

function [timeNum, psd, frequencyHz] = readJunoWavesSurvey(fileName)
    fid = fopen(fileName, 'r');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    if fid < 0
        error('Unable to open %s', fileName);
    end

    header = cell(5,1);
    for iLine = 1:5
        header{iLine} = fgetl(fid);
    end
    headerFields = split(string(header{4}), ',');
    headerFields = erase(headerFields, '"');
    frequencyHz = str2double(headerFields(29:154)).';

    frewind(fid);
    % Fields 1 and 2 provide time, field 4 provides three quality flags,
    % fields 29:154 contain the 126 calibrated electric spectral densities.
    formatString = ['%f%q%*q%q' repmat('%*q',1,24) repmat('%f',1,126)];
    values = textscan(fid, formatString, 'Delimiter', ',', ...
        'HeaderLines', 5, 'CollectOutput', false, 'EmptyValue', NaN, ...
        'ReturnOnError', false);

    sclk = values{1};
    scet = values{2};
    qualityFlags = values{3};
    psd = horzcat(values{4:end});
    if isempty(sclk)
        timeNum = [];
        return
    end

    t0 = parseScetDoy(scet{1});
    timeNum = t0 + (sclk - sclk(1))/86400;

    % PDS recommends excluding power-system-contaminated LFR samples from
    % survey plots.  Flag 2 applies to LFR_LO; flag 3 applies to LFR_HI.
    badLfrLo = cellfun(@(x) flagIsSet(x,2), qualityFlags);
    badLfrHi = cellfun(@(x) flagIsSet(x,3), qualityFlags);
    psd(badLfrLo,1:43) = NaN;
    psd(badLfrHi,44:61) = NaN;
end

function value = parseScetDoy(scet)
    token = regexp(scet, ...
        '^(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$', ...
        'tokens', 'once');
    if isempty(token)
        error('Unexpected SCET format: %s', scet);
    end
    yearValue = str2double(token{1});
    doyValue = str2double(token{2});
    hourValue = str2double(token{3});
    minuteValue = str2double(token{4});
    secondValue = str2double(token{5});
    value = datenum(yearValue,1,1,hourValue,minuteValue,secondValue) + doyValue - 1;
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
    lowerIndex = max(1, min(nValue, round(1 + lowerFraction*(nValue-1))));
    upperIndex = max(1, min(nValue, round(1 + upperFraction*(nValue-1))));
    limits = [values(lowerIndex), values(upperIndex)];
end

function [binnedTime, binnedPsd] = binMaximum(timeNum, psd, binSeconds)
    relativeSeconds = (timeNum - timeNum(1))*86400;
    binNumber = floor(relativeSeconds/binSeconds) + 1;
    firstRow = [1; find(diff(binNumber) > 0) + 1];
    lastRow = [firstRow(2:end) - 1; numel(timeNum)];
    nBin = numel(firstRow);
    binnedTime = NaN(nBin,1);
    binnedPsd = NaN(nBin,size(psd,2));
    for iBin = 1:nBin
        rows = firstRow(iBin):lastRow(iBin);
        binnedTime(iBin) = mean(timeNum(rows));
        binnedPsd(iBin,:) = max(psd(rows,:), [], 1, 'omitnan');
    end
end

function plotSpectrum(ax, timeNum, psd, frequencyHz, colorLimits)
    logPsd = log10(psd.');
    spectrumImage = imagesc(ax, timeNum, log10(frequencyHz), logPsd);
    set(spectrumImage, 'AlphaData', isfinite(logPsd));
    set(ax, 'YDir','normal', 'Layer','top', 'FontSize',10, ...
        'TickDir','out', 'Box','on', 'Color',[0.88 0.88 0.88]);
    ylim(ax, log10([50 41e6]));
    yticks(ax, log10([50 1e2 1e3 1e4 1e5 1e6 1e7 4.1e7]));
    yticklabels(ax, {'50 Hz','100 Hz','1 kHz','10 kHz','100 kHz', ...
        '1 MHz','10 MHz','41 MHz'});
    ylabel(ax, 'Frequency');
    caxis(ax, colorLimits);
    grid(ax, 'off');

    boundaries = [20e3 150e3 3e6];
    for value = boundaries
        yline(ax, log10(value), ':', 'Color',[1 1 1]*0.9, ...
            'LineWidth',0.7, 'Alpha',0.65);
    end
end

function addShockLine(ax, shockTime, labelText)
    xline(ax, shockTime, '--w', labelText, 'LineWidth',1.4, ...
        'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','left', 'FontWeight','bold');
end

function [timeNum, magneticField] = readJunoMag1s(fileName)
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
    formatString = repmat('%f', 1, 14);
    values = textscan(fid, formatString, 'CollectOutput',true, ...
        'ReturnOnError',false);
    values = values{1};
    if size(values,2) ~= 14
        error('Unexpected MAG row format in %s', fileName);
    end

    % Column 7 is decimal day of year, including the subsecond fraction.
    timeNum = datenum(values(1,1),1,1) + values(:,7) - 1;
    magneticField = values(:,8:10);
    magneticField(abs(magneticField) > 1e5) = NaN;
end

function [timeNum, speedKmS, peakEnergyEv, peakRate] = ...
        readJadeIonPeakSpeed(fileName)
    recordBytes = 45120;
    dataOffset = 320;       % PDS START_BYTE 321
    energyOffset = 25920;   % PDS START_BYTE 25921
    nEnergy = 64;
    nLogical = 25;
    allStopsLogical = 16;   % PDS logical counter [15], MATLAB index 16

    info = dir(fileName);
    nRecord = floor(info.bytes/recordBytes);
    if nRecord < 1 || info.bytes ~= nRecord*recordBytes
        error('JADE subset size is not an integer number of records: %s', fileName);
    end

    fid = fopen(fileName, 'r', 'ieee-le');
    if fid < 0
        error('Unable to open %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    timeNum = NaN(nRecord,1);
    speedKmS = NaN(nRecord,1);
    peakEnergyEv = NaN(nRecord,1);
    peakRate = NaN(nRecord,1);
    elementaryCharge = 1.602176634e-19;
    protonMass = 1.67262192369e-27;

    for iRecord = 1:nRecord
        recordStart = (iRecord-1)*recordBytes;
        fseek(fid, recordStart, 'bof');
        utcString = char(fread(fid, 21, '*char').');
        timeNum(iRecord) = parseScetDoy(utcString);

        fseek(fid, recordStart + dataOffset, 'bof');
        rawData = fread(fid, [nLogical nEnergy], ...
            'single=>double', 0, 'ieee-le').';
        fseek(fid, recordStart + energyOffset, 'bof');
        rawEnergy = fread(fid, [nLogical nEnergy], ...
            'single=>double', 0, 'ieee-le').';
        if numel(rawData) ~= nLogical*nEnergy || ...
                numel(rawEnergy) ~= nLogical*nEnergy
            continue
        end

        counts = rawData(:,allStopsLogical);
        energy = rawEnergy(:,allStopsLogical);
        valid = isfinite(counts) & counts > 0 & counts < 2.25e6 & ...
            isfinite(energy) & energy >= 80 & energy <= 5000;
        if nnz(valid) < 4
            continue
        end

        energy = energy(valid);
        counts = counts(valid);
        [energy, order] = sort(energy);
        counts = counts(order);
        smoothLogCounts = movmean(log10(counts), 3, 'Endpoints','shrink');
        [~, peakIndex] = max(smoothLogCounts);

        % Reject spectra with no discernible peak above their positive floor.
        contrast = counts(peakIndex)/max(median(counts), realmin);
        if counts(peakIndex) < 0.02 || contrast < 1.35
            continue
        end

        peakEnergyEv(iRecord) = energy(peakIndex);
        peakRate(iRecord) = counts(peakIndex);
        speedKmS(iRecord) = sqrt(2*elementaryCharge*energy(peakIndex) / ...
            protonMass)/1000;
    end
end

function styleTimeSeriesAxis(ax)
    set(ax, 'FontSize',10, 'TickDir','out', 'Box','on', ...
        'Layer','top', 'Color','w', 'XGrid','on', 'YGrid','on', ...
        'GridAlpha',0.14, 'MinorGridAlpha',0.08);
end

function addShockLineDark(ax, shockTime, labelText)
    xline(ax, shockTime, '--', labelText, 'Color',[0.22 0.22 0.22], ...
        'LineWidth',1.15, 'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','left', 'FontWeight','bold');
end
