%% Voyager 1/2 post-1991 unreviewed 1.92-second MAG archive
% Downloads both magnetometers:
%   primary   = out-board sensor
%   secondary = in-board sensor
%
% NASA states that these calibrated files are generally not science
% quality. Keep the two sensors separate and apply the point-quality flags
% stored in each CDF before quantitative analysis.

clearvars;
clc;

DataRoot = 'Z:\SPART-WORK\Data\Voyager';
StageDir = fullfile(tempdir, ...
    'Voyager_staging_mag2s_unreviewed_post1991');
StartDate = '1991-01-01';
EndDate = '2030-12-31';

Spacecraft = [1 2];
Products = 'mag2s_unreviewed';
Threads = 5;
CheckSize = true;
Force = false;
ManifestName = 'voyager_mag2s_unreviewed_post1991_manifest.json';

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
