%% Voyager 1/2 NASA/SPDF data download
% Edit the settings in this section, then run this script in MATLAB.

clearvars;
clc;

DataRoot = 'Z:\SPART-WORK\Data\Voyager';
StageDir = fullfile(tempdir, 'Voyager_staging');
StartDate = '1977-08-20';
EndDate = char(string(datetime('today', 'TimeZone', 'UTC'), 'yyyy-MM-dd'));

Spacecraft = [1 2];
Products = 'coho1hr,position1day';
Threads = 5;
CheckSize = true;
Force = false;

Date = [StartDate, '/', EndDate];

report = Voyager_Download( ...
    'Date', Date, ...
    'Spacecraft', Spacecraft, ...
    'Products', Products, ...
    'DataRoot', DataRoot, ...
    'StageDir', StageDir, ...
    'Threads', Threads, ...
    'CheckSize', CheckSize, ...
    'Force', Force);

disp(report)
