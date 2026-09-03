function result = Case1_Audit_LECP_Level1_Archive()
%Case1_Audit_LECP_Level1_Archive Audit all downloaded original L1 rate CDFs.
%   Source rates stay unchanged. Reports contain no calibrated flux, average,
%   interpolated data, rebinning, or newly applied scientific quality filter.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% original files and shared audit directory
ParentDir = 'Z:/SPART-WORK/Data/Voyager';
OutputDir = fullfile(ParentDir,'derived','source_audit', ...
    'lecp_official_cdf_archive','20260902_173114');
Manifest = struct([]); Summary = table; Events = table;
k = 0;
for Spacecraft = 1:2
    if Spacecraft == 1, Years = 2013:2021; else, Years = 2019:2020; end
    for Year = Years
        k = k+1;
        R = Case1_Audit_LECP_Level1(ParentDir,OutputDir,Spacecraft,Year);
        if k == 1, Manifest = R.Manifest; else, Manifest(k) = R.Manifest; end %#ok<AGROW>
        S = addvars(R.Summary,repmat(Spacecraft,height(R.Summary),1), ...
            repmat(Year,height(R.Summary),1),'Before',1, ...
            'NewVariableNames',{'Spacecraft','Year'});
        E = addvars(R.Events,repmat(Spacecraft,height(R.Events),1), ...
            repmat(Year,height(R.Events),1),'Before',1, ...
            'NewVariableNames',{'Spacecraft','Year'});
        Summary = [Summary;S]; Events = [Events;E]; %#ok<AGROW>
    end
end

%% auditable lists only; original CDF files remain unchanged
ManifestFile = fullfile(OutputDir,'LECP_Level1_Original_CDF_Manifest.json');
fid = fopen(ManifestFile,'w','n','UTF-8');
assert(fid >= 0,'Cannot write level-1 archive manifest.');
Clean = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(Manifest,'PrettyPrint',true));
writetable(Summary,fullfile(OutputDir,'LECP_Level1_Duration_Summary.csv'));
writetable(Events,fullfile(OutputDir,'LECP_Level1_49_Event_Audit.csv'));
result = struct('OutputDir',OutputDir,'ManifestFile',ManifestFile, ...
    'Files',numel(Manifest),'Bytes',sum([Manifest.bytes]), ...
    'Records',sum([Manifest.records]),'Events',height(Events));
disp(result);
end
