function Case1_Write_Audit(fileName, cfg, catalog, sourceManifest, runPlots)
%Case1_Write_Audit Write a reproducible record of inputs and assumptions.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-02

fid = fopen(fileName, 'wt', 'n', 'UTF-8');
if fid < 0
    error('VoyagerCase1:AuditOpen', ...
        'Cannot create audit file: %s', fileName);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Voyager Case 1 figure audit\n');
fprintf(fid, 'Generated UTC: %s\n', ...
    char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd HH:mm:ss XXX')));
fprintf(fid, 'RunPlots: %d\n', logical(runPlots));
fprintf(fid, 'Event count: %d (V1=%d, V2=%d)\n', height(catalog), ...
    nnz(catalog.Spacecraft == 1), nnz(catalog.Spacecraft == 2));
fprintf(fid, 'Data root: %s\n', cfg.DataRoot);
fprintf(fid, 'Output root: %s\n\n', cfg.OutputRoot);

fprintf(fid, '[Time and display]\n');
fprintf(fid, 'Window: [event day - %d d, event day + %d d + 1 d), UTC\n', ...
    cfg.ContextDays, cfg.ContextDays);
fprintf(fid, 'MAG/line-flux visual gap break: %.3g h\n', cfg.MAGGapHours);
fprintf(fid, 'Daily MAG overlay: finite 1 h samples, arithmetic mean, UTC bins\n');
fprintf(fid, 'LECP line channels: original 1 h values; no averaging=%d\n', ...
    ~cfg.LECPDailyAverage);
fprintf(fid, 'No B/flux interpolation, smoothing, extrapolation, or NaN-to-zero\n');
fprintf(fid, 'Only approved predicted attitude uses native CK type-3 rotation interpolation\n');
fprintf(fid, 'Event boundary lines displayed: %d\n', ...
    cfg.ShowEventBoundaries);
fprintf(fid, 'Color percentiles: %.3g and %.3g; rounded outward to 0.25 dex\n\n', ...
    cfg.ColorPercentiles(1), cfg.ColorPercentiles(2));

fprintf(fid, '[CDF quality handling]\n');
fprintf(fid, 'Reader: IRFU-MATLAB dataobj, direct from original CDF\n');
fprintf(fid, 'Exact FILLVAL and values outside VALIDMIN/VALIDMAX become NaN\n');
fprintf(fid, 'Native LECP duplicates: require equal Epoch/DeltaT/flux/sigma/energy before retaining one\n');
fprintf(fid, 'Missing values remain NaN and do not contribute to means\n\n');

fprintf(fid, '[Voyager 1 LECP sector treatment]\n');
fprintf(fid, 'User-approved record policy: %s\n', cfg.AccumulationPolicy);
fprintf(fid, 'Native Epoch is accumulation START; stop=Epoch+DeltaT, not a fixed 24 h interval\n');
fprintf(fid, 'Daily/hourly filename alone does not establish an individual record duration\n');
fprintf(fid, 'Discard DeltaT<0; retain all other in-window records at original Epoch\n');
fprintf(fid, 'No duration upper limit or UTC-hour/day-boundary rejection; no flux rebinning\n');
fprintf(fid, 'Official flux and uncertainty retained without a second temporal average\n');
fprintf(fid, 'PAD cadence: %s; nominal-width display cells centered on original Epoch\n', cfg.PADCadence);
fprintf(fid, 'P1 channel: nearest energy midpoint to 0.73 MeV\n');
fprintf(fid, 'P1 raw CDF bounds: FHDU_EnergyRange, currently 0.57-0.89 MeV; provenance only, physical discrepancy unresolved\n');
fprintf(fid, 'P1 displayed bounds in panels e/f: %.2f-%.2f MeV, user-approved convention following Mosley upper bound\n', cfg.P1DisplayEnergyMeV);
fprintf(fid, 'Display label change applies no energy integration, flux rescaling or CDF metadata mutation\n');
fprintf(fid, 'Explicit user-approved L1 fallback enabled: %d\n',cfg.LECPLevel1Fallback);
fprintf(fid, 'Explicit LECP source priority: %s\n',cfg.LECPSourcePriority);
if strcmp(cfg.LECPSourcePriority,'l1_first')
    fprintf(fid, 'Complete positive seven-sector L1 UTC means take priority over all L2 records in the bin\n');
    fprintf(fid, 'Without a complete L1 bin, retain original L2 rows; PAD completeness remains required\n');
else
    fprintf(fid, 'Complete seven-positive-sector L2 bins take priority; incomplete/missing L2 bins may use complete L1 UTC means\n');
end
fprintf(fid, 'No sector-level product mixing; all replaced L2 payloads retained in audit\n');
fprintf(fid, 'L1 negative DeltaT rejected; samples assigned by Epoch; no duration weighting/splitting or extra thresholds\n');
fprintf(fid, 'L1 J=mean(R)/[0.44*(1.78-0.57)]; historical adopted conversion, not a claim of official L2 calibration\n');
fprintf(fid, 'L1 sigma=sqrt(sum(sigma_i^2))/N if all contributing sigma exist; bin-center Epoch; source DeltaT retained separately\n');
fprintf(fid, 'L1 source records/means/counts/replaced L2 rows preserved in MAT; no background subtraction\n');
fprintf(fid, 'Active sectors: S1-S7; S8: shielded/background diagnostic\n');
fprintf(fid, 'Background mode: %s\n', cfg.LECPBackgroundMode);
fprintf(fid, 'Pitch-angle source policy: %s\n', cfg.PitchAngleMethod);
fprintf(fid, 'Prediction approved: %d; nominal mounting approved: %d\n', ...
    cfg.PredictedAttitudeApproved, cfg.NominalLECPGeometryApproved);
fprintf(fid, 'Attitude evaluated at L2 original Epoch or L1 UTC-bin center (instantaneous prediction)\n');
fprintf(fid, 'CK input cadence approx 7 d; native rotation interpolation, tolerance=0, no extrapolation\n');
fprintf(fid, 'LECP SC clock/cone axis=200/90 deg; S1 center=22.5 deg from HGA; increment=45 deg\n');
fprintf(fid, 'Sector numbering is right-handed about axisSC=(-sin15,cos15,0); HGA=-Z_SC\n');
fprintf(fid, 'P1 low-energy aperture only; particle velocity = minus outward look vector\n');
fprintf(fid, 'RTN R=Sun->V1; T=unit(SunNorth cross R); N=R cross T\n');
fprintf(fid, 'PA=angle(central ray at original Epoch, B vector mean in its UTC day/hour)\n');
fprintf(fid, 'No finite-FOV response integration; missing B or flux stays missing\n');
fprintf(fid, 'As-built sector offsets and actual maneuver/roll corrections unavailable\n');
fprintf(fid, 'User requested no in-figure prediction/nominal-geometry annotation\n');
fprintf(fid, 'Per-figure native MAT audit includes directions, PA, geometry and actual kernel SHA256\n');
fprintf(fid, 'No additional L2 temporal average; L1 uses approved UTC means; no sample-count/coverage threshold\n');
fprintf(fid, 'MAT audit records original CDF path, CDF row, Epoch, DeltaT and excluded rows\n');
fprintf(fid, ['PAD inclusion: every S1-S7 sector has finite positive ', ...
    'flux and finite pitch angle\n']);
fprintf(fid, ['MAG sample count/direction RMS, sector coverage, pitch ', ...
    'span, and relative uncertainty are diagnostic only\n']);
fprintf(fid, 'Display-only pitch center merge tolerance: %.3g deg\n', ...
    cfg.PitchMergeToleranceDeg);
fprintf(fid, 'Pitch-angle table export: %d\n\n', ...
    cfg.ExportPitchAngleTable);

fprintf(fid, '[Five-time hourly PAD figure]\n');
fprintf(fid, 'Export requested for this run: %d\n', cfg.ExportPeakPAD);
fprintf(fid, 'Peak search: maximum positive finite S1-S7 flux in the existing event context window; earliest tie\n');
fprintf(fid, 'Neighbors: nearest two usable PAD records before/after peak, skipping missing PAD; actual Epoch shown\n');
fprintf(fid, 'Missing peak PA leaves center blank; no weaker substitute, outside-window search or gap filling\n');
fprintf(fid, 'Each time: J/Jmax and sigma/Jmax, fixed observed denominator; no denominator error propagation\n');
fprintf(fid, 'PA axis 0-180 deg; seven unmerged sector points; no curve fit or additional data selection\n');
fprintf(fid, 'Peak figure folder: %s\n', cfg.PeakPADFolder);
fprintf(fid, 'Peak MAT audit folder: %s\n\n', cfg.PeakPADDataFolder);

fprintf(fid, '[Source files]\n');
for ii = 1:height(sourceManifest)
    fprintf(fid, ['%s | dataset=%s | cadence=%s | variables=%s | ', ...
        'exists=%d | required=%d | bytes=%.0f | events=%s | sha256=%s | note=%s\n'], ...
        char(sourceManifest.SourceFile(ii)), ...
        char(sourceManifest.DatasetID(ii)), ...
        char(sourceManifest.Cadence(ii)), ...
        char(sourceManifest.CDFVariables(ii)), ...
        sourceManifest.Exists(ii), sourceManifest.Required(ii), ...
        sourceManifest.Bytes(ii), char(sourceManifest.EventIDs(ii)), ...
        char(sourceManifest.SHA256(ii)), char(sourceManifest.Note(ii)));
end
clear cleanup
end
