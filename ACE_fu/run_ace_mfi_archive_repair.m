function status = run_ace_mfi_archive_repair(varargin)
%RUN_ACE_MFI_ARCHIVE_REPAIR Safely repair the full ACE/MFI CDF archive.
%
% Downloads are staged on local NTFS, validated with cdflib, then copied
% serially to Z:.  Native-product priority is H3 > H0 > H1 > H2 and lower
% cadence products are fetched when the higher product has hourly gaps.
%
% Example:
%   run_ace_mfi_archive_repair('StartDate','1997-09-02', ...
%       'EndDate',datestr(datetime('today'),'yyyy-mm-dd'),'Threads',4)

    scriptDir = fileparts(mfilename('fullpath'));

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'StartDate', '1997-09-02', @(x)ischar(x) || isstring(x));
    addParameter(p, 'EndDate', datestr(datetime('today'), 'yyyy-mm-dd'), ...
        @(x)ischar(x) || isstring(x));
    addParameter(p, 'DataRoot', 'Z:\SPART-WORK\Data\ACE', ...
        @(x)ischar(x) || isstring(x));
    addParameter(p, 'Threads', 4, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'PythonExe', ...
        'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe', ...
        @(x)ischar(x) || isstring(x));
    addParameter(p, 'StagingDir', fullfile(scriptDir, 'ace_repair_staging'), ...
        @(x)ischar(x) || isstring(x));
    addParameter(p, 'StateDir', fullfile(scriptDir, 'ace_repair_state'), ...
        @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    o = p.Results;

    repairScript = fullfile(scriptDir, 'repair_ace_mfi_archive.py');
    if ~isfile(repairScript)
        error('ACE:RepairScriptMissing', 'Repair script not found: %s', repairScript);
    end

    command = sprintf(['"%s" "%s" --start "%s" --end "%s" ' ...
        '--out-root "%s" --threads %d --staging-dir "%s" --state-dir "%s"'], ...
        char(o.PythonExe), repairScript, char(o.StartDate), char(o.EndDate), ...
        char(o.DataRoot), round(o.Threads), char(o.StagingDir), char(o.StateDir));

    fprintf('Starting ACE/MFI archive repair. Detailed logs will be written to:\n%s\n', ...
        char(o.StateDir));
    status = system(command, '-echo');
    if status ~= 0
        warning('ACE:RepairIncomplete', ...
            ['Repair exited with status %d. This may mean a transfer failed or ' ...
             'NASA has genuine source-data gaps; inspect the newest manifest.'], status);
    end
end
