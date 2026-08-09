%% Voyager 1/2 highest-resolution supplemental archive
% Downloads the finest broadly available official MAG/PLS products, the
% hourly HelioWeb position product, and JPL/NAIF SPICE trajectory kernels.

clearvars;
clc;

DataRoot = 'Z:\SPART-WORK\Data\Voyager';
StageDir = fullfile(tempdir, 'Voyager_staging_highres');
StartDate = '1977-08-20';
EndDate = char(string(datetime('today', 'TimeZone', 'UTC'), 'yyyy-MM-dd'));

Spacecraft = [1 2];
Products = 'mag2s,mag48s_vim,plasma_fine,position1hr';
IncludeSPICE = ismember(lower(string(getenv('VOYAGER_INCLUDE_SPICE'))), ...
    ["1", "true", "yes"]);
if IncludeSPICE
    Products = [Products, ',spice_spk'];
end
Threads = 5;
CheckSize = true;
Force = false;
ManifestName = 'voyager_highres_download_manifest.json';

Date = [StartDate, '/', EndDate];

report = Voyager_Download( ...
    'Date', Date, ...
    'Spacecraft', Spacecraft, ...
    'Products', Products, ...
    'DataRoot', DataRoot, ...
    'StageDir', StageDir, ...
    'ManifestName', ManifestName, ...
    'Threads', Threads, ...
    'CheckSize', CheckSize, ...
    'Force', Force);

disp(report)
