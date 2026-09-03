function report = Case1_Test_Predicted_Attitude()
%Case1_Test_Predicted_Attitude Numerical regression of the raw-kernel reader.
%   Run in an isolated MATLAB process with an empty SPICE kernel pool.
%   Tests all 3066 UTC-noon daily epochs from 2013-01-02 to 2021-05-25.
%   No flux data, formal pitch-angle product, or figure is generated.
%   The MAT audit is written only beneath the classified Voyager archive.

%% configuration and protection of the calling MATLAB process
cfg = Case1_Config;
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(cfg.IRFURoot, 'contrib', 'mice'));
if exist('mice', 'file') ~= 3
    error('VoyagerCase1:MICEMissing', 'IRFU MICE binary is unavailable.');
end
assertEmptyPool;
files = Case1_Attitude_Files(cfg);
kernelFiles = files.SourceFile(files.Role ~= "COMMENT");
time = (datetime(2013, 1, 2, 12, 0, 0, 'TimeZone', 'UTC'):days(1): ...
    datetime(2021, 5, 25, 12, 0, 0, 'TimeZone', 'UTC')).';
matrixTolerance = 1e-10;
angleToleranceDeg = 1e-10;
report = struct;
report.StartedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = mfilename('fullpath');
report.TestSHA256 = Case1_File_SHA256([report.TestFile, '.m']);
report.MATLABVersion = version;
report.Toolkit = cspice_tkvrsn('TOOLKIT');
report.TimeUTC = time;
report.MatrixTolerance = matrixTolerance;
report.AngleToleranceDeg = angleToleranceDeg;
report.Scope = ['Numerical and historical geometry regression only; ', ...
    'no flux averaging, formal pitch-angle product, or figure.'];
report.Checks = table('Size', [0 5], ...
    'VariableTypes', {'string', 'logical', 'double', 'double', 'string'}, ...
    'VariableNames', {'Check', 'Passed', 'Measured', 'Limit', 'Detail'});
report.FatalError = "";

%% complete daily interval and rotation algebra
try
    addCheck('daily_epoch_count', numel(time) == 3066, numel(time), ...
        3066, 'Inclusive UTC-noon grid; no averaging of attitude.');
    attitude = Case1_Read_Predicted_Attitude(time, cfg);
    report.Sources = attitude.Sources;
    report.Found = attitude.Found;
    addCheck('all_daily_attitudes_found', all(attitude.Found), ...
        sum(attitude.Found), numel(time), 'Every requested event-era day.');
    addCheck('reader_releases_its_kernels', poolIsEmpty, ...
        cspice_ktotal('ALL'), 0, 'Both KEEPER files and pool variables.');
    addCheck('source_hashes_present', ...
        all(strlength(attitude.Sources.SHA256) == 64), ...
        sum(strlength(attitude.Sources.SHA256) == 64), ...
        height(attitude.Sources), 'SHA256 identifies actual source bytes.');

    nTime = numel(time);
    orthogonality = nan(nTime, 3);
    determinant = nan(nTime, 3);
    radialError = nan(nTime, 1);
    roundTripError = nan(nTime, 1);
    angleError = nan(nTime, 1);
    compositionError = nan(nTime, 1);
    firstVector = [1; 2; 3] / norm([1; 2; 3]);
    secondVector = [-2; 3; 1] / norm([-2; 3; 1]);
    originalAngle = vectorAngle(firstVector, secondVector);
    for ii = 1:nTime
        if ~attitude.Found(ii), continue; end
        A = attitude.C_SC_to_J2000(:, :, ii);
        Q = attitude.C_J2000_to_RTN(:, :, ii);
        M = attitude.C_SC_to_RTN(:, :, ii);
        rotations = cat(3, A, Q, M);
        for jj = 1:3
            current = rotations(:, :, jj);
            orthogonality(ii, jj) = norm(current*current.'-eye(3), 'fro');
            determinant(ii, jj) = det(current);
        end
        radial = attitude.PositionJ2000_km(ii, :).';
        radialError(ii) = norm(Q*(radial/norm(radial))-[1; 0; 0]);
        roundTripError(ii) = max(norm(M.'*(M*firstVector)-firstVector), ...
            norm(M.'*(M*secondVector)-secondVector));
        angleError(ii) = abs(vectorAngle(M*firstVector, M*secondVector) ...
            - originalAngle);
        compositionError(ii) = norm(M-Q*A, 'fro');
    end
    report.DailyMetrics = table(time, attitude.Found, orthogonality, ...
        determinant, radialError, roundTripError, angleError, ...
        compositionError, 'VariableNames', {'EpochUTC', 'Found', ...
        'OrthogonalityError', 'Determinant', 'RadialMappingError', ...
        'RoundTripError', 'AngleInvariantErrorDeg', 'CompositionError'});
    addMaximumCheck('rotations_are_orthogonal', orthogonality, ...
        matrixTolerance, 'SC/J2000, J2000/RTN, and SC/RTN matrices.');
    addMaximumCheck('rotations_are_right_handed', abs(determinant-1), ...
        matrixTolerance, 'Determinant must equal +1, including its sign.');
    addMaximumCheck('RTN_radial_mapping', radialError, matrixTolerance, ...
        'Geometric Sun-to-spacecraft position maps to positive R.');
    addMaximumCheck('three_dimensional_round_trip', roundTripError, ...
        matrixTolerance, 'Two arbitrary vectors with all components nonzero.');
    addMaximumCheck('same_frame_angle_invariance', angleError, ...
        angleToleranceDeg, 'Both vectors rotated into the same frame.');
    addMaximumCheck('rotation_composition', compositionError, ...
        matrixTolerance, 'SC-to-RTN equals J2000-to-RTN times SC-to-J2000.');

    %% independent FK/pxform route, using the original CK and FK
    assertEmptyPool;
    kernelCleanup = onCleanup(@() unloadFiles(kernelFiles));
    loadFiles(kernelFiles);
    fkError = nan(nTime, 1);
    for ii = 1:nTime
        if ~attitude.Found(ii), continue; end
        independent = cspice_pxform('VG1_SC_BUS', 'J2000', ...
            attitude.EphemerisTimeTDB(ii));
        fkError(ii) = norm(independent - ...
            attitude.C_SC_to_J2000(:, :, ii), 'fro');
    end
    report.FKComparisonError = fkError;
    addMaximumCheck('FK_pxform_independent_comparison', fkError, ...
        matrixTolerance, 'Uses FK frame chain instead of explicit C transpose.');

    %% independently specified sector geometry and historical sign check
    scanNormalSC = [-sind(15); cosd(15); 0];
    scanBisectorSC = [0; 0; -1];
    theta = 22.5 + (0:7)*45;
    lookSC = scanBisectorSC*cosd(theta) + ...
        cross(scanNormalSC, scanBisectorSC)*sind(theta);
    particleSC = -lookSC;
    geometry = Case1_LECP_Geometry(cfg);
    addMaximumCheck('geometry_matches_independent_formula', ...
        [abs(geometry.LookSC(:)-lookSC(:)); ...
        abs(geometry.ParticleSC(:)-particleSC(:)); ...
        abs(geometry.AxisSC(:)-scanNormalSC); ...
        abs(geometry.HGASC(:)-scanBisectorSC)], matrixTolerance, ...
        'Regression of production geometry against the independent signed construction.');
    adjacentAngle = nan(8, 1);
    for ii = 1:8
        next = mod(ii, 8)+1;
        adjacentAngle(ii) = vectorAngle(lookSC(:, ii), lookSC(:, next));
    end
    bisector = lookSC(:, 1)+lookSC(:, 8);
    bisector = bisector/norm(bisector);
    addMaximumCheck('sector_adjacent_spacing', abs(adjacentAngle-45), ...
        angleToleranceDeg, 'Eight look-axis centres are separated by 45 deg.');
    addMaximumCheck('S1_S8_bisector', norm(bisector-scanBisectorSC), ...
        matrixTolerance, 'The S1/S8 look-axis bisector is SC negative Z.');
    addMaximumCheck('sector_scan_plane', abs(scanNormalSC.'*lookSC), ...
        matrixTolerance, 'Instrument mounting plane is defined in SC coordinates.');
    addMaximumCheck('sector_unit_norms', abs(vecnorm(lookSC)-1), ...
        matrixTolerance, 'No unit-vector component has been set by RTN geometry.');

    % The production x2100 SPK begins after these historical checks. Load
    % archived trajectory/planetary SPKs only for this historical section.
    positionRoot = fullfile(cfg.DataRoot, 'voyager1', 'ephemeris', ...
        'spice', 'kernels');
    historicalFiles = string({ ...
        fullfile(positionRoot, 'Voyager_1.a54206u_V0.2_merged.bsp'); ...
        fullfile(positionRoot, 'de440s.bsp')});
    historicalSHA256 = strings(numel(historicalFiles), 1);
    for ii = 1:numel(historicalFiles)
        if ~isfile(historicalFiles(ii))
            error('VoyagerCase1:HistoricalSPKMissing', ...
                'Historical regression requires existing SPK: %s', historicalFiles(ii));
        end
        historicalSHA256(ii) = Case1_File_SHA256(char(historicalFiles(ii)));
    end
    report.HistoricalAdditionalSPK = table(historicalFiles, historicalSHA256, ...
        'VariableNames', {'SourceFile', 'SHA256'});
    report.HistoricalSPKScope = ...
        'Additional archived SPKs loaded only for 1980 historical geometry checks.';
    historicalCleanup = onCleanup(@() unloadFiles(historicalFiles));
    loadFiles(historicalFiles);
    history = datetime([1980; 1980], [1; 10], [1; 1], ...
        [12; 12], [0; 0], [0; 0], 'TimeZone', 'UTC');
    historicalET = utcToET(history);
    targetLatitude = [30; 31];
    historicalToleranceDeg = 2;
    latitude = nan(2, 1);
    progradeProjection = nan(2, 1);
    historicalLook = nan(2, 3);
    for ii = 1:numel(history)
        scToEcliptic = cspice_pxform('VG1_SC_BUS', 'ECLIPJ2000', ...
            historicalET(ii));
        direction = scToEcliptic*geometry.LookSC(:, 3);
        historicalLook(ii, :) = direction.';
        latitude(ii) = asind(direction(3)/norm(direction));
        position = cspice_spkpos('-31', historicalET(ii), ...
            'ECLIPJ2000', 'NONE', 'SUN');
        prograde = cross([0; 0; 1], position/norm(position));
        prograde = prograde/norm(prograde);
        progradeProjection(ii) = dot(direction, prograde);
    end
    report.Geometry = geometry;
    report.Historical = table(history, targetLatitude, latitude, ...
        historicalLook, progradeProjection, ...
        'VariableNames', {'EpochUTC', 'ExpectedNorthLatitudeDeg', ...
        'LookS3LatitudeDeg', 'LookS3Ecliptic', 'ProgradeProjection'});
    report.HistoricalToleranceDeg = historicalToleranceDeg;
    report.HistoricalProgradeDefinition = ...
        'ECLIPJ2000 prograde tangent = unit(eclipticNorth cross SunToV1).';
    addMaximumCheck('historical_S3_north_ecliptic', ...
        abs(latitude-targetLatitude), historicalToleranceDeg, ...
        'Historical approximate +30/+31 deg checks; no PA data screening.');
    addCheck('historical_S3_counter_prograde', ...
        all(progradeProjection < 0), max(progradeProjection), 0, ...
        'Negative dot product with the explicitly recorded prograde tangent.');
    addCheck('historical_direction_sign_rejects_particle', ...
        all(abs(-latitude-targetLatitude) > historicalToleranceDeg), ...
        min(abs(-latitude-targetLatitude)), historicalToleranceDeg, ...
        'The historical statement describes the outward look; particle velocity is opposite.');
    clear historicalCleanup
    clear kernelCleanup
    assertEmptyPool;

    %% rejected time coverage retains missing values
    outside = datetime([1990; 2028], [4; 1], [15; 1], ...
        [12; 12], [0; 0], [0; 0], 'TimeZone', 'UTC');
    missing = Case1_Read_Predicted_Attitude(outside, cfg);
    missingMatrices = [missing.C_SC_to_J2000(:); ...
        missing.C_J2000_to_RTN(:); missing.C_SC_to_RTN(:)];
    addCheck('outside_predicted_coverage_is_missing', ...
        ~any(missing.Found) && all(isnan(missingMatrices)), ...
        sum(missing.Found), 0, ...
        'Before predicted-era cutoff and after the official CK endpoint.');
    report.OutsideCoverage = struct('TimeUTC', outside, ...
        'Found', missing.Found, 'Status', missing.Status);
    assertEmptyPool;

    %% reader must reject nonempty caller pools without deleting their data
    lskFile = files.SourceFile(files.Role == "LSK");
    [fileProtected, fileDetail] = testFilePool(char(lskFile), time(1), cfg);
    addCheck('caller_kernel_file_protected', fileProtected, ...
        double(fileProtected), 1, fileDetail);
    assertEmptyPool;
    [variableProtected, variableDetail] = testVariablePool(time(1), cfg);
    addCheck('caller_program_variable_protected', variableProtected, ...
        double(variableProtected), 1, variableDetail);
    assertEmptyPool;
catch ME
    report.FatalError = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    addCheck('execution_completed', false, NaN, NaN, string(ME.message));
    if exist('historicalCleanup', 'var'), clear historicalCleanup; end
    if exist('kernelCleanup', 'var'), clear kernelCleanup; end
end

%% reproducible MAT audit; no CSV or formal PA output
report.FinalPoolEmpty = poolIsEmpty;
addCheck('test_releases_its_pool', report.FinalPoolEmpty, ...
    cspice_ktotal('ALL'), 0, 'Only files and sentinel created by this test are removed.');
report.CompletedUTC = datetime('now', 'TimeZone', 'UTC');
report.Passed = all(report.Checks.Passed) && strlength(report.FatalError) == 0;
folder = fullfile(cfg.DataRoot, 'voyager1', 'attitude', 'spice', 'validation');
if ~isfolder(folder), mkdir(folder); end
report.AuditFile = fullfile(folder, ['predicted_attitude_regression_', ...
    char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyyMMdd_HHmmss_SSS')), '.mat']);
save(report.AuditFile, 'report');
fprintf('Predicted-attitude regression: passed=%d, days=%d\n', ...
    report.Passed, numel(time));
fprintf('Test audit: %s\n', report.AuditFile);
disp(report.Checks);

    function addCheck(name, passed, measured, limit, detail)
        report.Checks(end+1, :) = {string(name), logical(passed), ...
            double(measured), double(limit), string(detail)};
    end

    function addMaximumCheck(name, value, limit, detail)
        value = value(:);
        passed = ~isempty(value) && all(isfinite(value)) && all(value <= limit);
        if isempty(value), measured = NaN; else, measured = max(value); end
        addCheck(name, passed, measured, limit, detail);
    end
end

function empty = poolIsEmpty()
[~, hasVariables] = cspice_gnpool('*', 1, 1);
empty = cspice_ktotal('ALL') == 0 && ~hasVariables;
end

function assertEmptyPool()
if ~poolIsEmpty
    error('VoyagerCase1:TestPoolNotEmpty', ...
        'Run this test in an isolated MATLAB process; the caller pool is preserved.');
end
end

function angle = vectorAngle(first, second)
value = dot(first, second)/(norm(first)*norm(second));
angle = acosd(max(-1, min(1, value)));
end

function et = utcToET(time)
time.TimeZone = 'UTC';
time.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSSSSS';
text = cellstr(char(time));
text = cellfun(@(x) [x, ' UTC'], text, 'UniformOutput', false);
et = cspice_str2et(text);
end

function loadFiles(files)
for ii = 1:numel(files)
    cspice_furnsh(char(files(ii)));
end
end

function unloadFiles(files)
for ii = numel(files):-1:1
    try
        cspice_unload(char(files(ii)));
    catch
        % Cleanup only; retain the original test exception for the audit.
    end
end
end

function [protected, detail] = testFilePool(fileName, time, cfg)
assertEmptyPool;
cleanup = onCleanup(@() cspice_unload(fileName));
cspice_furnsh(fileName);
before = cspice_gdpool('DELTET/DELTA_T_A', 1, 1);
identifier = expectPoolRejection(time, cfg);
[~, ~, ~, stillLoaded] = cspice_kinfo(fileName);
[after, variableFound] = cspice_gdpool('DELTET/DELTA_T_A', 1, 1);
protected = identifier == "VoyagerCase1:KernelPoolNotEmpty" && ...
    cspice_ktotal('ALL') == 1 && stillLoaded && variableFound && isequal(before, after);
detail = "Existing LSK and DELTET/DELTA_T_A preserved; exception=" + identifier;
end

function [protected, detail] = testVariablePool(time, cfg)
assertEmptyPool;
% Keep below the current MICE pool-name storage limit. An overlong name
% can be truncated on insertion yet remain untruncated in cspice_dvpool.
name = 'VGR_CASE1_TEST_SENTINEL';
value = [12.345; -67.89];
cleanup = onCleanup(@() cspice_dvpool(name));
cspice_pdpool(name, value);
identifier = expectPoolRejection(time, cfg);
[after, found] = cspice_gdpool(name, 1, 2);
protected = identifier == "VoyagerCase1:KernelPoolNotEmpty" && ...
    cspice_ktotal('ALL') == 0 && found && isequal(value, after(:));
detail = "Program-injected variable preserved; exception=" + identifier;
end

function identifier = expectPoolRejection(time, cfg)
identifier = "no exception";
try
    Case1_Read_Predicted_Attitude(time, cfg);
catch ME
    identifier = string(ME.identifier);
end
end
