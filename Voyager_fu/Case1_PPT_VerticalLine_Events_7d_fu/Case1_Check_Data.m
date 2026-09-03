function sourceManifest = Case1_Check_Data(cfg, catalog)
%Case1_Check_Data Verify every source needed by the 49 Case 1 windows.
%   The returned table records full source paths and the EventIDs that use
%   each source. This function does not modify data or calculate averages.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-02

%% collect required COHO months
sourceFile = strings(0,1);
spacecraft = zeros(0,1);
instrument = strings(0,1);
productLevel = strings(0,1);
cadence = strings(0,1);
eventID = strings(0,1);
required = false(0,1);
note = strings(0,1);

for ii = 1:height(catalog)
    plotStart = catalog.StartUTC(ii) - days(cfg.ContextDays);
    plotEnd = catalog.EndUTCExclusive(ii) + days(cfg.ContextDays);
    monthStart = dateshift(plotStart, 'start', 'month');
    finalMonth = dateshift(plotEnd - seconds(1), 'start', 'month');
    months = (monthStart:calmonths(1):finalMonth)';

    for jj = 1:numel(months)
        yearText = sprintf('%04d', year(months(jj)));
        monthText = sprintf('%02d', month(months(jj)));
        monthFolder = fullfile(cfg.DataRoot, ...
            sprintf('voyager%d', catalog.Spacecraft(ii)), 'coho', ...
            '1hr', 'l2', 'merged_mag_plasma', yearText, monthText);
        entries = dir(fullfile(monthFolder, '*.cdf'));
        if isempty(entries)
            sourceFile(end+1,1) = string(fullfile(monthFolder, '*.cdf')); %#ok<AGROW>
            spacecraft(end+1,1) = catalog.Spacecraft(ii); %#ok<AGROW>
            instrument(end+1,1) = "COHO merged MAG/plasma/particles"; %#ok<AGROW>
            productLevel(end+1,1) = "L2"; %#ok<AGROW>
            cadence(end+1,1) = "1 hr grid"; %#ok<AGROW>
            eventID(end+1,1) = string(catalog.EventID(ii)); %#ok<AGROW>
            required(end+1,1) = true; %#ok<AGROW>
            note(end+1,1) = "Required monthly CDF is missing"; %#ok<AGROW>
        else
            [~, order] = sort({entries.name});
            entries = entries(order);
            for kk = 1:numel(entries)
                sourceFile(end+1,1) = string(fullfile( ... %#ok<AGROW>
                    entries(kk).folder, entries(kk).name));
                spacecraft(end+1,1) = catalog.Spacecraft(ii); %#ok<AGROW>
                instrument(end+1,1) = "COHO merged MAG/plasma/particles"; %#ok<AGROW>
                productLevel(end+1,1) = "L2"; %#ok<AGROW>
                cadence(end+1,1) = "1 hr grid"; %#ok<AGROW>
                eventID(end+1,1) = string(catalog.EventID(ii)); %#ok<AGROW>
                required(end+1,1) = true; %#ok<AGROW>
                note(end+1,1) = ""; %#ok<AGROW>
            end
        end
    end
end

%% add original annual V1 LECP sources; historical subset is comparison only
if any(catalog.Spacecraft == 1)
    nativeCDFs = string(cfg.LECPNativeDailyCDFs);
    productLabel = "daily";
    if strcmp(cfg.PADCadence, 'hour')
        nativeCDFs = string(cfg.LECPNativeHourlyCDFs);
        productLabel = "hourly";
    end
    if isempty(nativeCDFs)
        error('VoyagerCase1:NativeCDFMissing', 'Original annual LECP CDFs are missing.');
    end
    for ii = 1:numel(nativeCDFs)
        sourceFile(end+1,1) = nativeCDFs(ii); %#ok<AGROW>
        spacecraft(end+1,1) = 1; %#ok<AGROW>
        instrument(end+1,1) = "LECP sectored hydrogen flux"; %#ok<AGROW>
        productLevel(end+1,1) = "L2 original annual CDF"; %#ok<AGROW>
        cadence(end+1,1) = productLabel + " product; original Epoch anchor"; %#ok<AGROW>
        eventID(end+1,1) = strjoin(string(catalog.EventID( ...
            catalog.Spacecraft == 1)), ';'); %#ok<AGROW>
        required(end+1,1) = true; %#ok<AGROW>
        note(end+1,1) = "Discard DeltaT<0; preserve other original Epoch/flux/sigma; no flux rebinning"; %#ok<AGROW>
    end
end

if any(catalog.Spacecraft == 1) && cfg.LECPLevel1Fallback
    for ii = 1:numel(cfg.LECPLevel1CDFs)
        sourceFile(end+1,1) = cfg.LECPLevel1CDFs(ii); %#ok<AGROW>
        spacecraft(end+1,1) = 1; %#ok<AGROW>
        instrument(end+1,1) = "LECP sectored hydrogen rate"; %#ok<AGROW>
        productLevel(end+1,1) = "L1 original annual CDF"; %#ok<AGROW>
        cadence(end+1,1) = "Native Epoch; derived UTC " + string(cfg.PADCadence) + " means; priority=" + string(cfg.LECPSourcePriority); %#ok<AGROW>
        eventID(end+1,1) = strjoin(string(catalog.EventID(catalog.Spacecraft == 1)),';'); %#ok<AGROW>
        required(end+1,1) = true; %#ok<AGROW>
        note(end+1,1) = "Approved L1 route; J=mean(R)/[0.44*(1.78-0.57)]; no background correction"; %#ok<AGROW>
    end
end

if strcmp(cfg.PitchAngleMethod, 'predicted_ck')
    kernels = Case1_Attitude_Files(cfg);
    for ii = 1:height(kernels)
        sourceFile(end+1,1) = kernels.SourceFile(ii); %#ok<AGROW>
        spacecraft(end+1,1) = 1; %#ok<AGROW>
        instrument(end+1,1) = "SPICE " + kernels.Role(ii); %#ok<AGROW>
        productLevel(end+1,1) = "official raw kernel/support file"; %#ok<AGROW>
        cadence(end+1,1) = "CK approx 7 d; evaluate at original LECP Epoch"; %#ok<AGROW>
        eventID(end+1,1) = strjoin(string(catalog.EventID( ...
            catalog.Spacecraft == 1)), ';'); %#ok<AGROW>
        required(end+1,1) = kernels.Role(ii) ~= "COMMENT" && ...
            any(catalog.Spacecraft == 1); %#ok<AGROW>
        note(end+1,1) = kernels.SourceURL(ii); %#ok<AGROW>
    end
else
    sourceFile(end+1,1) = string(cfg.LECPSectorPointingFile);
    spacecraft(end+1,1) = 1;
    instrument(end+1,1) = "LECP external sector pointing";
    productLevel(end+1,1) = "derived pointing table";
    cadence(end+1,1) = "daily or finer";
    eventID(end+1,1) = strjoin(string(catalog.EventID(catalog.Spacecraft == 1)), ';');
    required(end+1,1) = false;
    note(end+1,1) = "External pointing interface; no implicit prediction fallback";
end

%% collapse repeated source rows
[uniqueFile, ~, group] = unique(sourceFile, 'stable');
nSource = numel(uniqueFile);
outSpacecraft = zeros(nSource,1);
outInstrument = strings(nSource,1);
outLevel = strings(nSource,1);
outCadence = strings(nSource,1);
outEvents = strings(nSource,1);
outRequired = false(nSource,1);
outExists = false(nSource,1);
outBytes = nan(nSource,1);
outNote = strings(nSource,1);

for ii = 1:nSource
    index = group == ii;
    outSpacecraft(ii) = spacecraft(find(index, 1));
    outInstrument(ii) = instrument(find(index, 1));
    outLevel(ii) = productLevel(find(index, 1));
    outCadence(ii) = cadence(find(index, 1));
    outEvents(ii) = strjoin(unique(eventID(index), 'stable'), ';');
    outRequired(ii) = any(required(index));
    outExists(ii) = isfile(uniqueFile(ii));
    if outExists(ii)
        info = dir(uniqueFile(ii));
        outBytes(ii) = info(1).bytes;
    end
    theseNotes = unique(note(index & note ~= ""), 'stable');
    outNote(ii) = strjoin(theseNotes, '; ');
end

sourceManifest = table(outSpacecraft, outInstrument, outLevel, ...
    outCadence, uniqueFile, outExists, outRequired, outBytes, ...
    outEvents, outNote, 'VariableNames', {'Spacecraft', 'Instrument', ...
    'ProductLevel', 'Cadence', 'SourceFile', 'Exists', 'Required', ...
    'Bytes', 'EventIDs', 'Note'});

sourceManifest.DatasetID = repmat("Voyager COHO one-hour merged", ...
    height(sourceManifest), 1);
sourceManifest.CDFVariables = repmat(strjoin([ ...
    "Epoch", "ABS_B", "F", "BR", "BT", "BN", "V", ...
    "protonDensity", "protonTemp", "protonFlux*_LECP", ...
    "protonFlux*_CRS"], ';'), height(sourceManifest), 1);
isLECP = sourceManifest.Instrument == "LECP sectored hydrogen flux";
sourceManifest.DatasetID(isLECP) = ...
    "VOYAGER-1_LECP_LEV-2-DAILY-AVG";
if strcmp(cfg.PADCadence, 'hour')
    sourceManifest.DatasetID(isLECP) = "VOYAGER-1_LECP_LEV-2-HOURLY-AVG";
end
sourceManifest.CDFVariables(isLECP) = strjoin([ ...
    "Epoch", "DeltaT", "FHDU_SectoredFluxes", ...
    "FHDU_SectoredFluxUncertainties", "FHDU_Energy", "FHDU_EnergyRange", ...
    "SectorIterator", "Hydrogen_Channels", "Hydrogen_Channels_Label"], ';');
isPointing = contains(sourceManifest.Instrument, "sector pointing");
isRate = sourceManifest.Instrument == "LECP sectored hydrogen rate";
sourceManifest.DatasetID(isRate) = "VOYAGER-1_LECP_LEV-1-RATES";
sourceManifest.CDFVariables(isRate) = ...
    "Epoch;DeltaT;FHDU_SectoredRates;FHDU_SectoredRateUncertainties;FHDU_Energy;FHDU_EnergyRange;SectorIterator;Hydrogen_Channels_Label";
sourceManifest.DatasetID(isPointing) = ...
    "Expected SEDR/AACS/NAVMAG-derived pointing";
sourceManifest.CDFVariables(isPointing) = ...
    "EpochUTC;ParticleUR_S1..S8;ParticleUT_S1..S8;ParticleUN_S1..S8";
isKernel = startsWith(sourceManifest.Instrument, "SPICE ");
sourceManifest.DatasetID(isKernel) = "NASA/JPL NAIF Voyager SPICE";
sourceManifest.CDFVariables(isKernel) = ...
    "Non-CDF original CK/FK/SCLK/LSK/SPK/PCK; documented predicted segment";
sourceManifest.SHA256 = strings(height(sourceManifest), 1);
for ii = find(sourceManifest.Exists).'
    sourceManifest.SHA256(ii) = Case1_File_SHA256(char(sourceManifest.SourceFile(ii)));
end
sourceManifest = movevars(sourceManifest, {'DatasetID', 'CDFVariables'}, ...
    'After', 'Cadence');

missingRequired = sourceManifest.Required & ~sourceManifest.Exists;
if any(missingRequired)
    missingList = strjoin(sourceManifest.SourceFile(missingRequired), newline);
    error('VoyagerCase1:RequiredDataMissing', ...
        'Required Voyager source data are missing:%s%s', newline, missingList);
end
end
