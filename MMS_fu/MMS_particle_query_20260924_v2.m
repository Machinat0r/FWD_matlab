%% 查询本批事件的FEEPS / HPCA，直接调用用户现有SDCFilenames
clear;clc
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
RecordDir='Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2';
addpath(CodeDir);
Date='2026-07-20/2026-09-17';
for Instrument={'feeps','hpca'}
    for Mode={'brst','srvy'}
        Query=struct('instrument',Instrument{1},'mode',Mode{1},'date',Date,'complete',false,'error','');
        try
            fprintf('QUERY %s %s\n',Instrument{1},Mode{1});
            filenames=SDCFilenames(Date,1:4,'inst',Instrument{1},'drm',Mode{1});
            filenames=filenames(contains(filenames,'_l2_'));
            Query.files=filenames;Query.complete=true;
            Query.queriedUTC=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'));
            fprintf('FOUND %s %s: %d\n',Instrument{1},Mode{1},numel(filenames));
        catch ME
            Query.error=getReport(ME,'extended','hyperlinks','off');
            fprintf(2,'%s\n',Query.error);
        end
        fid=fopen(fullfile(RecordDir,['query_' Instrument{1} '_' Mode{1} '.json']),'w','n','UTF-8');
        fprintf(fid,'%s',jsonencode(Query));fclose(fid);
    end
end
fprintf('PARTICLE_QUERY_FINISHED\n');