%% Published Juno bow-shock crossings: direct fields and full-band Waves
% Four crossings that were explicitly discussed or analyzed as event cases
% in the papers identified for this project are plotted with a common style.
%
% Science-data rule:
% - JSS Br/Btheta/Bphi are read directly from JADE L5 MAG_VECTOR_JSSRTP.
% - The 2016-07-17 event has no JADE data (JADE was off); that panel reads
%   the direct PC Bx/By/Bz fields used by Hospodarsky et al. (2017).
% - XYZ and R use direct JADE SC_POS_JUPITER_J2000XYZ and SC_POS_R fields.
%   XYZ is only converted from km to R_J using 71492 km/R_J. The PC fallback
%   has direct XYZ but no direct R, so the R row is intentionally shown as --.
% - |B| is calculated from the three displayed components with the user's
%   explicit permission. Waves pixels are selected actual source records;
%   no spectral averaging or interpolation is performed.

clear; close all; clc;

scriptFile = mfilename('fullpath');
outputDir = fileparts(scriptFile);
databaseRoot = 'Z:\SPART-WORK\Data\Juno';
assert(isfolder(databaseRoot),'Juno database is unavailable: %s',databaseRoot);

events = defineEvents(databaseRoot);
statusRows = cell(numel(events),7);

for iEvent = 1:numel(events)
    event = events(iEvent);
    fprintf('\n[%d/%d] %s\n',iEvent,numel(events),event.id);
    try
        [waveTime,wavePsd,frequencyHz,receiverBoundaries] = ...
            processWaves(event.waveFiles,event.timeRange,event.waveBinSeconds);

        if strcmp(event.fieldKind,'JADE_RTP')
            [magTime,magComponents,positionTime,positionXyzKm,positionRadiusRj] = ...
                readJadeContext(event.fieldFile,event.timeRange);
            componentLabels = {'B_r','B_\theta','B_\phi','|B|'};
            magneticAxisLabel = 'B^{JSS}_{r,\theta,\phi} (nT)';
        elseif strcmp(event.fieldKind,'MAG_PC')
            [magTime,magComponents,positionTime,positionXyzKm] = ...
                readJunoMagPc(event.fieldFile,event.timeRange);
            positionRadiusRj = NaN(size(positionTime));
            componentLabels = {'B_x','B_y','B_z','|B|'};
            magneticAxisLabel = 'B^{PC}_{x,y,z} (nT)';
        else
            error('Unknown field kind: %s',event.fieldKind);
        end

        magTotal = sqrt(sum(magComponents.^2,2));
        tickValues = event.timeRange(1):event.tickStepMinutes/1440:event.timeRange(2);
        remainingMinutes = (event.timeRange(2)-tickValues(end))*1440;
        if remainingMinutes >= 0.45*event.tickStepMinutes
            tickValues(end+1) = event.timeRange(2); %#ok<SAGROW>
        end
        [positionSourceTime,tickPositionRj] = positionAtTicks( ...
            positionTime,positionXyzKm,positionRadiusRj,tickValues, ...
            event.maxPositionOffsetMinutes);

        fig = createEventFigure(event,magTime,magComponents,magTotal, ...
            componentLabels,magneticAxisLabel,waveTime,wavePsd, ...
            frequencyHz,receiverBoundaries,tickValues,tickPositionRj);

        pngFile = fullfile(outputDir,[event.outputTag '.png']);
        pdfFile = fullfile(outputDir,[event.outputTag '.pdf']);
        exportgraphics(fig,pngFile,'Resolution',260);
        exportgraphics(fig,pdfFile,'ContentType','image','Resolution',260);
        close(fig);

        fprintf('  Saved %s\n',pngFile);
        fprintf('  Direct magnetic records: %d\n',nnz(any(isfinite(magComponents),2)));
        fprintf('  Direct Waves display records: %d/%d\n', ...
            nnz(any(isfinite(wavePsd),2)),size(wavePsd,1));
        for iTick = 1:numel(tickValues)
            fprintf(['    %s: X=%+.3f Y=%+.3f Z=%+.3f R=%s R_J ' ...
                '(source offset %.1f s)\n'], ...
                datestr(tickValues(iTick),'yyyy-mm-dd HH:MM:SS'), ...
                tickPositionRj(iTick,1:3),formatRadius(tickPositionRj(iTick,4)), ...
                abs(positionSourceTime(iTick)-tickValues(iTick))*86400);
        end

        statusRows(iEvent,:) = {event.id,event.article,event.doi, ...
            datestr(event.shockTime,'yyyy-mm-dd HH:MM:SS.FFF'), ...
            event.direction,pngFile,'OK'};
    catch ME
        statusRows(iEvent,:) = {event.id,event.article,event.doi, ...
            datestr(event.shockTime,'yyyy-mm-dd HH:MM:SS.FFF'), ...
            event.direction,'',ME.message};
        fprintf(2,'  FAILED: %s\n',getReport(ME,'extended'));
    end
end

manifest = cell2table(statusRows,'VariableNames', ...
    {'Event','Article','DOI','ShockUTC','Direction','PNG','Status'});
manifestFile = fullfile(outputDir,'juno_published_bowshock_event_figures.csv');
writetable(manifest,manifestFile);
disp(manifest(:,{'Event','Article','ShockUTC','Status'}));
fprintf('Saved manifest %s\n',manifestFile);
assert(all(strcmp(manifest.Status,'OK')),'One or more event figures failed.');

%% Local functions
function events = defineEvents(root)
    events(1) = makeEvent( ...
        'first_approach_20160624', ...
        [datenum(2016,6,24,7,0,0) datenum(2016,6,24,9,0,0)], ...
        datenum(2016,6,24,8,16,0),30,30,'in', ...
        'McComas et al. (2017); Hospodarsky et al. (2017)', ...
        '10.1002/2017GL072831; 10.1002/2017GL073177', ...
        {fullfile(root,'Waves_data','2016005_JUPITER_APPROACH','2016_176', ...
        'WAV_2016176T000000_E_V02.CSV')}, ...
        'JADE_RTP',fullfile(root,'JADE_data','2016','2016176','ION_TOF', ...
        'JAD_L50_HLS_ION_TOF_DEF_2016176_V01.DAT'), ...
        'juno_bowshock_20160624_0816_mccomas2017_hospodarsky2017');

    events(2) = makeEvent( ...
        'orbit0_20160717', ...
        [datenum(2016,7,17,11,0,0) datenum(2016,7,17,17,0,0)], ...
        datenum(2016,7,17,15,33,0),60,30,'out', ...
        'Hospodarsky et al. (2017)', ...
        '10.1002/2017GL073177', ...
        {fullfile(root,'Waves_data','2016186_ORBIT_00', ...
        'WAV_2016199T000000_E_V02.CSV')}, ...
        'MAG_PC',fullfile(root,'MAG','2016','JUPITER','PC','PERI-00', ...
        'fgm_jno_l3_2016199pc_r1s_v01.sts'), ...
        'juno_bowshock_20160717_1533_hospodarsky2017');

    events(3) = makeEvent( ...
        'joseph_case1_20241210', ...
        [datenum(2024,12,10,17,0,0) datenum(2024,12,10,18,46,0)], ...
        mean([datenum(2024,12,10,17,24,51.757), ...
              datenum(2024,12,10,17,24,54.507)]),15,30,'out', ...
        'Joseph et al. (2026), Case I', ...
        '10.1038/s41467-026-76223-x', ...
        {fullfile(root,'Waves_data','2024313_ORBIT_67', ...
        'WAV_2024345T000000_E_V01.CSV')}, ...
        'JADE_RTP',fullfile(root,'JADE_data','2024','2024345','ELECTRONS', ...
        'JAD_L50_LRS_ELC_ANY_DEF_2024345_V02.DAT'), ...
        'juno_bowshock_20241210_case1_joseph2026');

    events(4) = makeEvent( ...
        'joseph_case2_20241211', ...
        [datenum(2024,12,11,3,0,0) datenum(2024,12,11,3,10,0)], ...
        mean([datenum(2024,12,11,3,6,1.461), ...
              datenum(2024,12,11,3,6,4.211)]),2,10,'in', ...
        'Joseph et al. (2026), Case II', ...
        '10.1038/s41467-026-76223-x', ...
        {fullfile(root,'Waves_data','2024345_ORBIT_68', ...
        'WAV_2024346T000000_E_V01.CSV')}, ...
        'JADE_RTP',fullfile(root,'JADE_data','2024','2024346','ELECTRONS', ...
        'JAD_L50_LRS_ELC_ANY_DEF_2024346_V02.DAT'), ...
        'juno_bowshock_20241211_case2_joseph2026');

    % Use the requested event-specific log-PSD range for the 2016-07-17 case.
    events(2).colorLimits = [10^(-13.5) 1e-10];
    events(2).colorTickExponents = -13.5:0.5:-10;
end

function event = makeEvent(id,timeRange,shockTime,tickStepMinutes, ...
        waveBinSeconds,direction,article,doi,waveFiles,fieldKind, ...
        fieldFile,outputTag)
    event = struct('id',id,'timeRange',timeRange,'shockTime',shockTime, ...
        'tickStepMinutes',tickStepMinutes,'waveBinSeconds',waveBinSeconds, ...
        'direction',direction,'article',article,'doi',doi, ...
        'waveFiles',{waveFiles},'fieldKind',fieldKind, ...
        'fieldFile',fieldFile,'outputTag',outputTag, ...
        'maxPositionOffsetMinutes',10,'colorLimits',[1e-15 1e-8], ...
        'colorTickExponents',{[]});
    for iFile = 1:numel(waveFiles)
        assert(isfile(waveFiles{iFile}),'Missing Waves file: %s',waveFiles{iFile});
    end
    assert(isfile(fieldFile),'Missing field file: %s',fieldFile);
end

function [gridTime,gridPsd,frequencyHz,receiverBoundaries] = ...
        processWaves(files,timeRange,binSeconds)
    nBin = ceil(diff(timeRange)*86400/binSeconds);
    gridTime = NaN(nBin,1);
    gridPsd = NaN(nBin,126);
    gridDistance = inf(nBin,1);
    rawFrequency = [];
    for iFile = 1:numel(files)
        [timeNum,psd,thisFrequency] = readJunoWavesSurvey(files{iFile});
        if isempty(rawFrequency)
            rawFrequency = thisFrequency(:);
        else
            assert(all(abs(thisFrequency(:)-rawFrequency) <= ...
                max(1e-9,1e-10*abs(rawFrequency))), ...
                'Waves frequency grids differ: %s',files{iFile});
        end
        use = timeNum >= timeRange(1) & timeNum < timeRange(2);
        [gridTime,gridPsd,gridDistance] = selectDirectSamples( ...
            gridTime,gridPsd,gridDistance,timeNum(use),psd(use,:), ...
            timeRange(1),binSeconds);
    end
    assert(any(isfinite(gridTime)),'No Waves records in event interval.');
    receiverChannels = {1:43,44:61,62:88,89:126};
    receiverBoundaries = NaN(1,3);
    for iReceiver = 1:3
        lowerMaximum = max(rawFrequency(receiverChannels{iReceiver}));
        upperCenters = rawFrequency(receiverChannels{iReceiver+1});
        nextUpperCenter = min(upperCenters(upperCenters > lowerMaximum));
        receiverBoundaries(iReceiver) = sqrt(lowerMaximum*nextUpperCenter);
    end
    [frequencyHz,order] = sort(rawFrequency);
    gridPsd = gridPsd(:,order);
    gridPsd(gridPsd <= 0) = NaN;
end

function [magTime,magRtp,posTime,xyzKm,radiusRj] = ...
        readJadeContext(fileName,timeRange)
    labelFile = regexprep(fileName,'\.DAT$','.LBL','ignorecase');
    assert(isfile(labelFile),'Missing JADE label: %s',labelFile);
    labelText = fileread(labelFile);
    recordBytes = readPdsLabelInteger(labelText,'RECORD_BYTES');
    nRecords = readPdsLabelInteger(labelText,'FILE_RECORDS');
    info = dir(fileName);
    assert(info.bytes == recordBytes*nRecords, ...
        'JADE DAT size does not match label dimensions: %s',fileName);

    magUtcByte = readPdsColumnStartByte(labelText,'MAG_UTC');
    magRtpByte = readPdsColumnStartByte(labelText,'MAG_VECTOR_JSSRTP');
    posUtcByte = readPdsColumnStartByte(labelText,'DIM0_UTC');
    radiusByte = readPdsColumnStartByte(labelText,'SC_POS_R');
    xyzByte = readPdsColumnStartByte(labelText,'SC_POS_JUPITER_J2000XYZ');

    magUtc = readStridedBytes(fileName,recordBytes,nRecords,magUtcByte,21);
    posUtc = readStridedBytes(fileName,recordBytes,nRecords,posUtcByte,21);
    magRtp = readStridedSingle(fileName,recordBytes,nRecords,magRtpByte,3);
    radiusRj = readStridedSingle(fileName,recordBytes,nRecords,radiusByte,1);
    xyzKm = readStridedSingle(fileName,recordBytes,nRecords,xyzByte,3);
    magTime = parsePdsDoyUtc(char(magUtc.'));
    posTime = parsePdsDoyUtc(char(posUtc.'));

    badMag = any(~isfinite(magRtp),2) | any(magRtp == 9990000,2) | ...
        any(abs(magRtp)>1.6e6,2);
    magRtp(badMag,:) = NaN;
    badPos = ~isfinite(radiusRj) | radiusRj == 65535 | ...
        any(~isfinite(xyzKm),2) | any(xyzKm == 65535,2);
    radiusRj(badPos) = NaN;
    xyzKm(badPos,:) = NaN;

    useMag = magTime >= timeRange(1) & magTime < timeRange(2);
    usePos = posTime >= timeRange(1) & posTime < timeRange(2);
    magTime = magTime(useMag);
    magRtp = magRtp(useMag,:);
    posTime = posTime(usePos);
    xyzKm = xyzKm(usePos,:);
    radiusRj = radiusRj(usePos);
    assert(any(any(isfinite(magRtp))),'No direct JADE RTP records in interval.');
end

function [timeNum,magneticField,positionTime,positionKm] = ...
        readJunoMagPc(fileName,timeRange)
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
    positionKm = values(:,12:14);
    use = timeNum >= timeRange(1) & timeNum < timeRange(2);
    timeNum = timeNum(use);
    magneticField = magneticField(use,:);
    positionKm = positionKm(use,:);
    magneticField(abs(magneticField)>1e5) = NaN;
    positionKm(abs(positionKm)>1e9) = NaN;
    positionTime = timeNum;
end

function bytes = readStridedBytes(fileName,recordBytes,nRecords,startByte,nByte)
    fid = fopen(fileName,'r','ieee-le');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    assert(fseek(fid,startByte-1,'bof')==0,'Unable to seek in %s',fileName);
    bytes = fread(fid,[nByte nRecords],sprintf('%d*uint8=>uint8',nByte), ...
        recordBytes-nByte);
    assert(size(bytes,2)==nRecords,'Incomplete records in %s',fileName);
end

function values = readStridedSingle(fileName,recordBytes,nRecords,startByte,nItem)
    fid = fopen(fileName,'r','ieee-le');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    assert(fseek(fid,startByte-1,'bof')==0,'Unable to seek in %s',fileName);
    values = fread(fid,[nItem nRecords],sprintf('%d*single=>double',nItem), ...
        recordBytes-4*nItem,'ieee-le').';
    assert(size(values,1)==nRecords,'Incomplete records in %s',fileName);
end

function [sourceTime,tickPositionRj] = positionAtTicks( ...
        positionTime,xyzKm,radiusRj,tickValues,maxOffsetMinutes)
    valid = isfinite(positionTime) & all(isfinite(xyzKm),2);
    positionTime = positionTime(valid);
    xyzKm = xyzKm(valid,:);
    radiusRj = radiusRj(valid);
    assert(~isempty(positionTime),'No valid direct XYZ position records.');
    tickPositionRj = NaN(numel(tickValues),4);
    sourceTime = NaN(numel(tickValues),1);
    for iTick = 1:numel(tickValues)
        [distance,index] = min(abs(positionTime-tickValues(iTick)));
        assert(distance <= maxOffsetMinutes/1440, ...
            'No direct position within %.1f min of %s.', ...
            maxOffsetMinutes,datestr(tickValues(iTick)));
        sourceTime(iTick) = positionTime(index);
        tickPositionRj(iTick,1:3) = xyzKm(index,:)/71492;
        if isfinite(radiusRj(index)) && radiusRj(index) ~= 65535
            tickPositionRj(iTick,4) = radiusRj(index);
        end
    end
end

function fig = createEventFigure(event,magTime,magB,magTotal, ...
        componentLabels,magneticAxisLabel,waveTime,wavePsd, ...
        frequencyHz,receiverBoundaries,tickValues,tickPositionRj)
    fig = figure('Color','w','Position',[30 30 1800 1360], ...
        'Renderer','opengl');
    layout = tiledlayout(fig,13,1,'TileSpacing','compact','Padding','compact');
    finiteField = [magB(:);magTotal(:)];
    finiteField = finiteField(isfinite(finiteField));
    fieldLimits = [min(finiteField) max(finiteField)];
    fieldPadding = max(0.05,0.04*diff(fieldLimits));

    axMag = nexttile(layout,1,[3 1]);
    hold(axMag,'on');
    h1 = plot(axMag,magTime,magB(:,1),'Color',[0.19 0.39 0.78],'LineWidth',0.72);
    h2 = plot(axMag,magTime,magB(:,2),'Color',[0.95 0.48 0.12],'LineWidth',0.72);
    h3 = plot(axMag,magTime,magB(:,3),'Color',[0.14 0.55 0.25],'LineWidth',0.72);
    h4 = plot(axMag,magTime,magTotal,'k','LineWidth',0.9);
    styleTimeSeriesAxis(axMag);
    xlim(axMag,event.timeRange);
    ylim(axMag,fieldLimits+[-fieldPadding fieldPadding]);
    ylabel(axMag,magneticAxisLabel);
    legendLocation = 'northeast';
    if diff(event.timeRange)*1440 <= 15
        legendLocation = 'northwest';
    end
    lgd = legend(axMag,[h1 h2 h3 h4],componentLabels,'Location',legendLocation, ...
        'NumColumns',4,'Box','off','FontSize',9);
    lgd.AutoUpdate = 'off';

    axWave = nexttile(layout,4,[8 1]);
    plotSpectrum(axWave,waveTime,wavePsd,frequencyHz, ...
        receiverBoundaries,event.colorLimits);
    xlim(axWave,event.timeRange);

    for ax = [axMag axWave]
        xline(ax,event.shockTime,'--','Color',[0.56 0.05 0.68], ...
            'LineWidth',1.15,'HandleVisibility','off');
    end
    linkaxes([axMag axWave],'x');
    set([axMag axWave],'XTick',tickValues);
    axMag.XTickLabel = [];
    axWave.XTickLabel = [];
    addTimeAndPositionLabels(axWave,event.timeRange,tickValues,tickPositionRj);

    axBottomSpace = nexttile(layout,12,[2 1]);
    axis(axBottomSpace,'off');
    colormap(fig,turbo(256));
    cb = colorbar(axWave,'eastoutside');
    cb.Label.String = 'E PSD  [(V m^{-1})^2 Hz^{-1}]';
    cb.FontSize = 10;
    if ~isempty(event.colorTickExponents)
        cb.Ticks = 10.^event.colorTickExponents;
        cb.TickLabels = arrayfun(@(value) sprintf('10^{%g}',value), ...
            event.colorTickExponents,'UniformOutput',false);
        cb.TickLabelInterpreter = 'tex';
    end
end

function addTimeAndPositionLabels(ax,timeRange,tickValues,positionRj)
    xNorm = (tickValues-timeRange(1))/diff(timeRange);
    rowY = [-0.076 -0.116 -0.156 -0.196];
    rowLabels = {'X','Y','Z','R'};
    text(ax,-0.010,-0.036,'(R_J)','Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',8.5,'Clipping','off','Interpreter','tex');
    for iRow = 1:4
        text(ax,-0.010,rowY(iRow),rowLabels{iRow},'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','top', ...
            'FontSize',8.5,'Clipping','off','Interpreter','none');
    end
    for iTick = 1:numel(tickValues)
        alignment = 'center';
        if xNorm(iTick)<0.03
            alignment = 'left';
        elseif xNorm(iTick)>0.97
            alignment = 'right';
        end
        text(ax,xNorm(iTick),-0.036,datestr(tickValues(iTick),'HH:MM'), ...
            'Units','normalized','HorizontalAlignment',alignment, ...
            'VerticalAlignment','top','FontSize',8.5, ...
            'Clipping','off','Interpreter','none');
        for iRow = 1:4
            if ~isfinite(positionRj(iTick,iRow))
                valueText = '--';
            elseif iRow<4
                valueText = sprintf('%+.1f',positionRj(iTick,iRow));
            else
                valueText = sprintf('%.1f',positionRj(iTick,iRow));
            end
            text(ax,xNorm(iTick),rowY(iRow),valueText,'Units','normalized', ...
                'HorizontalAlignment',alignment,'VerticalAlignment','top', ...
                'FontSize',8.5,'Clipping','off','Interpreter','none');
        end
    end
end

function plotSpectrum(ax,timeNum,psd,frequencyHz,receiverBoundaries,colorLimits)
    validTime = isfinite(timeNum);
    timeNum = timeNum(validTime);
    psd = psd(validTime,:);
    [timeNum,timeOrder] = sort(timeNum(:));
    psd = psd(timeOrder,:);
    [frequencyHz,frequencyOrder] = sort(frequencyHz(:));
    directPsd = psd(:,frequencyOrder).';
    assert(numel(timeNum)>=2,'Too few Waves samples for spectrogram.');
    timeEdges = [timeNum(1)-0.5*(timeNum(2)-timeNum(1)); ...
        0.5*(timeNum(1:end-1)+timeNum(2:end)); ...
        timeNum(end)+0.5*(timeNum(end)-timeNum(end-1))];
    logF = log10(frequencyHz);
    midF = 0.5*(logF(1:end-1)+logF(2:end));
    fEdges = [2*logF(1)-midF(1);midF;2*logF(end)-midF(end)];
    fEdges(1) = log10(50);
    fEdges(end) = log10(41e6);
    [xGrid,yGrid] = meshgrid(timeEdges,fEdges);
    colorData = NaN(size(xGrid));
    colorData(1:end-1,1:end-1) = directPsd;
    alphaData = zeros(size(xGrid));
    alphaData(1:end-1,1:end-1) = isfinite(directPsd) & directPsd>0;
    h = surface(ax,xGrid,yGrid,zeros(size(xGrid)),colorData, ...
        'FaceColor','flat','EdgeColor','none');
    set(h,'AlphaData',alphaData,'FaceAlpha','flat','AlphaDataMapping','none');
    view(ax,2);
    set(ax,'YDir','normal','Layer','top','FontSize',10,'TickDir','out', ...
        'Box','on','Color',[0.88 0.88 0.88],'ColorScale','log');
    ylim(ax,log10([50 41e6]));
    yticks(ax,log10([50 1e2 1e3 1e4 1e5 1e6 1e7 4.1e7]));
    yticklabels(ax,{'50 Hz','100 Hz','1 kHz','10 kHz','100 kHz', ...
        '1 MHz','10 MHz','41 MHz'});
    ylabel(ax,'Frequency');
    caxis(ax,colorLimits);
    set(ax,'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off');
    for value = receiverBoundaries
        yline(ax,log10(value),'--k','LineWidth',0.9);
    end
end

function [targetTime,targetValues,targetDistance] = selectDirectSamples( ...
        targetTime,targetValues,targetDistance,timeNum,values,t0,binSeconds)
    if isempty(timeNum), return; end
    binIndex = floor((timeNum-t0)*86400/binSeconds)+1;
    valid = binIndex>=1 & binIndex<=size(targetValues,1);
    binIndex = binIndex(valid);
    timeNum = timeNum(valid);
    values = values(valid,:);
    centers = t0+(binIndex-0.5)*binSeconds/86400;
    distance = abs(timeNum-centers);
    for iBin = unique(binIndex(:)).'
        candidates = find(binIndex==iBin);
        [candidateDistance,localIndex] = min(distance(candidates));
        if candidateDistance < targetDistance(iBin)
            source = candidates(localIndex);
            targetTime(iBin) = timeNum(source);
            targetValues(iBin,:) = values(source,:);
            targetDistance(iBin) = candidateDistance;
        end
    end
end

function [timeNum,psd,frequencyHz] = readJunoWavesSurvey(fileName)
    fid = fopen(fileName,'r');
    assert(fid>=0,'Unable to open %s',fileName);
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    header = cell(5,1);
    for iLine = 1:5, header{iLine}=fgetl(fid); end
    fields = erase(split(string(header{4}),','),'"');
    frequencyHz = str2double(fields(29:154)).';
    frewind(fid);
    fmt = ['%f%q%*q%q' repmat('%*q',1,24) repmat('%f',1,126)];
    values = textscan(fid,fmt,'Delimiter',',','HeaderLines',5, ...
        'CollectOutput',false,'EmptyValue',NaN,'ReturnOnError',false);
    scet = values{2};
    qualityFlags = values{3};
    psd = horzcat(values{4:end});
    timeNum = parsePdsDoyUtc(char(scet));
    badLo = cellfun(@(x)flagIsSet(x,2),qualityFlags);
    badHi = cellfun(@(x)flagIsSet(x,3),qualityFlags);
    psd(badLo,1:43) = NaN;
    psd(badHi,44:61) = NaN;
end

function value = readPdsLabelInteger(labelText,keyword)
    token = regexp(labelText,['(?m)^\s*' regexptranslate('escape',keyword) ...
        '\s*=\s*(\d+)'],'tokens','once');
    assert(~isempty(token),'Missing %s in PDS label.',keyword);
    value = str2double(token{1});
end

function startByte = readPdsColumnStartByte(labelText,columnName)
    expression = ['(?s)OBJECT\s*=\s*COLUMN\s+' ...
        'NAME\s*=\s*"?' regexptranslate('escape',columnName) ...
        '"?.*?START_BYTE\s*=\s*(\d+).*?END_OBJECT\s*=\s*COLUMN'];
    token = regexp(labelText,expression,'tokens','once');
    assert(~isempty(token),'Missing PDS column %s.',columnName);
    startByte = str2double(token{1});
end

function timeNum = parsePdsDoyUtc(utcChars)
    years = str2double(cellstr(utcChars(:,1:4)));
    doys = str2double(cellstr(utcChars(:,6:8)));
    hours = str2double(cellstr(utcChars(:,10:11)));
    minutes = str2double(cellstr(utcChars(:,13:14)));
    seconds = str2double(cellstr(utcChars(:,16:21)));
    timeNum = datenum(years,ones(size(years)),ones(size(years)), ...
        hours,minutes,seconds)+doys-1;
end

function tf = flagIsSet(flagString,flagNumber)
    parts = strsplit(flagString,':');
    tf = numel(parts)>=flagNumber && strcmp(parts{flagNumber},'1');
end

function styleTimeSeriesAxis(ax)
    set(ax,'FontSize',10,'TickDir','out','Box','on','Layer','top', ...
        'Color','w','XGrid','off','YGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off');
end

function value = formatRadius(radius)
    if isfinite(radius)
        value = sprintf('%.3f',radius);
    else
        value = '--';
    end
end
