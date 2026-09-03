function Case1_Add_IRFU_Path(irfuRoot)
%Case1_Add_IRFU_Path Add the local IRFU-MATLAB package to the MATLAB path.
%   Generated paths named .git, .agents, or .codex are excluded.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-02

%% input check
if ~isfolder(irfuRoot)
    error('VoyagerCase1:IRFUMissing', ...
        'IRFU-MATLAB folder is unavailable: %s', irfuRoot);
end

%% add reviewed package paths
pathList = string(strsplit(genpath(irfuRoot), pathsep));
pathList(pathList == "") = [];
blocked = contains(lower(pathList), ["\.git", "\.agents", "\.codex"]);
addpath(char(join(pathList(~blocked), pathsep)));

if exist('dataobj', 'file') ~= 2
    error('VoyagerCase1:IRFUDataobjMissing', ...
        'IRFU dataobj is unavailable after adding: %s', irfuRoot);
end
end
