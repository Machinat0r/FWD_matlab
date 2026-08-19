%% Juno bow-shock overview: direct JSO MAG and full-band Waves measurements
% Data products: JNO-J-3-FGM-CAL-V1.0 (official Sun-State/JSO MAG) and
% JNO-E/J/SS-WAV-3-CDR-SRVFULL-V2.0.
% Exact shock ramps reported in the paper are 2024-12-10
% 17:24:51.757--17:24:54.507 UTC and 2024-12-11
% 03:06:01.461--03:06:04.211 UTC.  Case windows follow Figures 3 and 4.
% The top survey-channel center is 40.5 MHz; the nominal instrument band
% extends to 41 MHz.  Spectral densities are already calibrated by PDS.
% Every displayed physical value is read directly from a PDS source field.
% No plasma moment, spectral-peak proxy, or magnetic magnitude is calculated.

clear; close all; clc;

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
databaseRoot = 'Z:\SPART-WORK\Data\Juno';
assert(isfolder(databaseRoot), 'Juno database is unavailable: %s', databaseRoot);
fprintf('Database root: %s\n', databaseRoot);

files = strings(3,1);
files(1) = fullfile(databaseRoot, 'Waves_data', '2024313_ORBIT_67', ...
    'WAV_2024345T000000_E_V01.CSV');
files(2) = fullfile(databaseRoot, 'Waves_data', '2024345_ORBIT_68', ...
    'WAV_2024345T211242_E_V01.CSV');
files(3) = fullfile(databaseRoot, 'Waves_data', '2024345_ORBIT_68', ...
    'WAV_2024346T000000_E_V01.CSV');

magFiles = strings(2,1);
magFiles(1) = fullfile(databaseRoot, 'MAG', '2024', 'JUPITER', 'SS', ...
    'PERI-67', 'fgm_jno_l3_2024345ss_r1s_v01.sts');
magFiles(2) = fullfile(databaseRoot, 'MAG', '2024', 'JUPITER', 'SS', ...
    'PERI-68', 'fgm_jno_l3_2024346ss_r1s_v01.sts');

inputFiles = [files; magFiles];
for iInput = 1:numel(inputFiles)
    if ~isfile(inputFiles(iInput))
        error('Required database file is missing: %s', inputFiles(iInput));
    end
end

tOverview = [datenum(2024,12,10,12,0,0), datenum(2024,12,11,7,0,0)];
shockRampRanges = [
    datenum(2024,12,10,17,24,51.757), datenum(2024,12,10,17,24,54.507)
    datenum(2024,12,11, 3, 6, 1.461), datenum(2024,12,11, 3, 6, 4.211)
    ];
caseRange1 = [datenum(2024,12,10,17,0,0), datenum(2024,12,10,18,45,0)];
caseRange2 = [datenum(2024,12,11, 3,0,0), datenum(2024,12,11, 3,10,0)];

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

% Derive one visible separator per receiver transition from the actual
% survey-channel centers.  The first two receiver pairs overlap, so the
% separator is placed after the last channel of the lower receiver.  This
% matches the geometric frequency-bin edges used by plotSpectrum.
receiverChannels = {1:43, 44:61, 62:88, 89:126};
assert(numel(frequencyHz) == 126, ...
    'Unexpected Waves survey channel count: %d', numel(frequencyHz));
receiverBoundaries = NaN(1,3);
for iReceiver = 1:3
    lowerMaximum = max(frequencyHz(receiverChannels{iReceiver}));
    upperCenters = frequencyHz(receiverChannels{iReceiver+1});
    nextUpperCenter = min(upperCenters(upperCenters > lowerMaximum));
    assert(~isempty(nextUpperCenter), ...
        'Unable to derive receiver boundary %d.', iReceiver);
    receiverBoundaries(iReceiver) = sqrt(lowerMaximum*nextUpperCenter);
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

% Read the official 1-s Sun-State components.  The PDS Sun-State axes are
% the JSO axes used here; no coordinate rotation is performed in MATLAB.
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

% A common robust color range keeps the three panels quantitatively comparable.
sampleValues = log10(allPsd(1:20:end,:));
sampleValues = sampleValues(isfinite(sampleValues));
logColorLimits = robustLimits(sampleValues, 0.005, 0.997);
logColorLimits = [floor(2*logColorLimits(1))/2, ...
    ceil(2*logColorLimits(2))/2];
if diff(logColorLimits) < 5
    middleValue = mean(logColorLimits);
    logColorLimits = middleValue + [-2.5 2.5];
end
colorLimits = 10.^logColorLimits;

% The overview uses the source records directly; no temporal bin statistic.
overviewTime = allTime;
overviewPsd = allPsd;

magValues = abs(magB(:));
magValues = magValues(isfinite(magValues));
magLimit = max(2, ceil(1.05*max(magValues)));

use1 = allTime >= caseRange1(1) & allTime <= caseRange1(2);
use2 = allTime >= caseRange2(1) & allTime <= caseRange2(2);
for iCase = 1:size(shockRampRanges,1)
    assert(all(shockRampRanges(iCase,:) >= tOverview(1) & ...
        shockRampRanges(iCase,:) <= tOverview(2)), ...
        'Shock ramp for Case %d lies outside the overview interval.', iCase);
end

figures = gobjects(3,1);
figures(1) = createTwoPanelFigure(tOverview, ...
    magTime, magB, ...
    overviewTime, overviewPsd, frequencyHz, receiverBoundaries, colorLimits, ...
    60, true, magLimit);

figures(2) = createTwoPanelFigure(caseRange1, ...
    magTime, magB, ...
    allTime(use1), allPsd(use1,:), frequencyHz, receiverBoundaries, colorLimits, ...
    15, false, magLimit);

figures(3) = createTwoPanelFigure(caseRange2, ...
    magTime, magB, ...
    allTime(use2), allPsd(use2,:), frequencyHz, receiverBoundaries, colorLimits, ...
    1, false, magLimit);

baseNames = {
    'juno_bowshock_overview_mag_waves_20241210_11'
    'juno_bowshock_case1_mag_waves_20241210'
    'juno_bowshock_case2_mag_waves_20241211'
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
fprintf('MAG source samples: %d\n', numel(magTime));
fprintf('Frequency centers: %.4g Hz to %.4g Hz (%d channels)\n', ...
    min(frequencyHz), max(frequencyHz), numel(frequencyHz));
fprintf('Color limits: %.4g to %.4g (V/m)^2/Hz (log color scale)\n', ...
    colorLimits);
fprintf('Receiver separators: %.6g, %.6g, %.6g Hz\n', ...
    receiverBoundaries);
fprintf(['Verified paper shock ramps: Case I %s--%s UTC; ' ...
    'Case II %s--%s UTC.\n'], ...
    datestr(shockRampRanges(1,1),'HH:MM:SS.FFF'), ...
    datestr(shockRampRanges(1,2),'HH:MM:SS.FFF'), ...
    datestr(shockRampRanges(2,1),'HH:MM:SS.FFF'), ...
    datestr(shockRampRanges(2,2),'HH:MM:SS.FFF'));

%% Local functions
function fig = createTwoPanelFigure(timeRange, ...
        magTime, magB, ...
        waveTime, wavePsd, frequencyHz, receiverBoundaries, colorLimits, ...
        tickMinutes, showDateChanges, magLimit)
    fig = figure('Color','w', 'Position',[40 40 1800 1200], ...
        'Renderer','opengl');
    layout = tiledlayout(fig, 3, 1, 'TileSpacing','compact', ...
        'Padding','compact');

    useMag = magTime >= timeRange(1) & magTime <= timeRange(2);
    axMag = nexttile(layout, 1);
    hold(axMag, 'on');
    hBx = plot(axMag, magTime(useMag), magB(useMag,1), ...
        'Color',[0.19 0.39 0.78], 'LineWidth',0.75);
    hBy = plot(axMag, magTime(useMag), magB(useMag,2), ...
        'Color',[0.14 0.55 0.25], 'LineWidth',0.75);
    hBz = plot(axMag, magTime(useMag), magB(useMag,3), ...
        'Color',[0.78 0.19 0.19], 'LineWidth',0.75);
    styleTimeSeriesAxis(axMag);
    xlim(axMag, timeRange);
    ylim(axMag, [-magLimit magLimit]);
    ylabel(axMag, 'B_{JSO} (nT)');
    magLegend = legend(axMag, [hBx hBy hBz], ...
        {'B_{x,JSO}','B_{y,JSO}','B_{z,JSO}'}, 'Location','northeast', ...
        'NumColumns',3, 'Box','off', 'FontSize',9);
    magLegend.AutoUpdate = 'off';

    axWave = nexttile(layout, 2, [2 1]);
    plotSpectrum(axWave, waveTime, wavePsd, frequencyHz, ...
        receiverBoundaries, colorLimits);
    xlim(axWave, timeRange);

    linkaxes([axMag axWave], 'x');
    [tickValues, tickLabels] = makeTimeTicks(timeRange, tickMinutes, ...
        showDateChanges);
    set([axMag axWave], 'XTick',tickValues);
    axMag.XTickLabel = [];
    axWave.XTickLabel = tickLabels;
    colormap(fig, turbo(256));
    cb = colorbar(axWave, 'eastoutside');
    cb.Label.String = 'E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 10;
end

function [tickValues, tickLabels] = makeTimeTicks(timeRange, ...
        tickMinutes, showDateChanges)
    stepDays = tickMinutes/1440;
    % The tolerance keeps an exact endpoint tick despite datenum roundoff.
    firstIndex = ceil(timeRange(1)/stepDays - 1e-4);
    lastIndex = floor(timeRange(2)/stepDays + 1e-4);
    tickValues = (firstIndex:lastIndex)*stepDays;
    endpointTolerance = 0.1/86400; % 0.1 s in datenum units
    if ~isempty(tickValues) && ...
            abs(tickValues(1)-timeRange(1)) < endpointTolerance
        tickValues(1) = timeRange(1);
    end
    if ~isempty(tickValues) && ...
            abs(tickValues(end)-timeRange(2)) < endpointTolerance
        tickValues(end) = timeRange(2);
    end
    tickLabels = cellstr(datestr(tickValues, 'HH:MM'));

    if showDateChanges && ~isempty(tickValues)
        dateChange = [true, diff(floor(tickValues)) ~= 0];
        dateLabels = cellstr(datestr(tickValues(dateChange), ...
            'mm/dd HH:MM'));
        tickLabels(dateChange) = dateLabels;
    end
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

    % PDS recommends excluding power-system-contaminated LFR samples from
    % survey plots.  Flag 2 applies to LFR_LO; flag 3 applies to LFR_HI.
    badLfrLo = cellfun(@(x) flagIsSet(x,2), qualityFlags);
    badLfrHi = cellfun(@(x) flagIsSet(x,3), qualityFlags);
    psd(badLfrLo,1:43) = NaN;
    psd(badLfrHi,44:61) = NaN;
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

function plotSpectrum(ax, timeNum, psd, frequencyHz, ...
        receiverBoundaries, colorLimits)
    if numel(timeNum) < 2 || numel(frequencyHz) < 2
        error('At least two time and frequency samples are required.');
    end

    [timeNum, timeOrder] = sort(timeNum(:));
    psd = psd(timeOrder,:);
    [frequencyHz, frequencyOrder] = sort(frequencyHz(:));
    directPsd = psd(:,frequencyOrder).';

    timeEdges = [timeNum(1) - 0.5*(timeNum(2)-timeNum(1)); ...
        0.5*(timeNum(1:end-1)+timeNum(2:end)); ...
        timeNum(end) + 0.5*(timeNum(end)-timeNum(end-1))];
    logFrequency = log10(frequencyHz);
    logFrequencyMidpoints = 0.5*(logFrequency(1:end-1) + ...
        logFrequency(2:end));
    logFrequencyEdges = [2*logFrequency(1)-logFrequencyMidpoints(1); ...
        logFrequencyMidpoints; ...
        2*logFrequency(end)-logFrequencyMidpoints(end)];
    logFrequencyEdges(1) = log10(50);
    logFrequencyEdges(end) = log10(41e6);

    [xGrid, yGrid] = meshgrid(timeEdges, logFrequencyEdges);
    colorData = NaN(size(xGrid));
    colorData(1:end-1,1:end-1) = directPsd;
    alphaData = zeros(size(xGrid));
    alphaData(1:end-1,1:end-1) = isfinite(directPsd) & directPsd > 0;
    spectrumSurface = surface(ax, xGrid, yGrid, zeros(size(xGrid)), ...
        colorData, 'FaceColor','flat', 'EdgeColor','none');
    set(spectrumSurface, 'AlphaData',alphaData, 'FaceAlpha','flat', ...
        'AlphaDataMapping','none');
    view(ax, 2);
    set(ax, 'YDir','normal', 'Layer','top', 'FontSize',10, ...
        'TickDir','out', 'Box','on', 'Color',[0.88 0.88 0.88], ...
        'ColorScale','log');
    ylim(ax, log10([50 41e6]));
    yticks(ax, log10([50 1e2 1e3 1e4 1e5 1e6 1e7 4.1e7]));
    yticklabels(ax, {'50 Hz','100 Hz','1 kHz','10 kHz','100 kHz', ...
        '1 MHz','10 MHz','41 MHz'});
    ylabel(ax, 'Frequency');
    caxis(ax, colorLimits);
    grid(ax, 'off');

    for value = receiverBoundaries
        yline(ax, log10(value), '--k', 'LineWidth',0.9);
    end
end

function [timeNum, magneticField] = readJunoMag1s(fileName)
    assert(~isempty(regexp(char(fileName), 'ss_r1s_v\d+\.sts$', 'once')), ...
        'Expected an official Sun-State/JSO 1-s MAG file: %s', fileName);
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

function styleTimeSeriesAxis(ax)
    set(ax, 'FontSize',10, 'TickDir','out', 'Box','on', ...
        'Layer','top', 'Color','w', 'XGrid','on', 'YGrid','on', ...
        'GridAlpha',0.14, 'MinorGridAlpha',0.08);
end
