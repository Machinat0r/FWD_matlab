%% 检查已归档FEEPS survey：只读接口诊断，不修改原数据
addpath('C:\Users\Administrator\Documents\irfu-matlab-master');irf('check_path');
mms.db_init('local_file_db','Z:\SPART-WORK\Data\MMS');
tint=irf.tint('2026-07-28T05:10:00Z/2026-07-28T05:30:00Z');
for particle={'electron','ion'}
    Flux=mms.get_data(['Omniflux' particle{1} '_epd_feeps_srvy_l2'],tint,1);
    if isempty(Flux),fprintf('%s EMPTY\n',particle{1});continue;end
    fprintf('%s rows=%d finite=%d positive=%d units=%s\n',particle{1},size(Flux.data,1),sum(isfinite(Flux.data),'all'),sum(Flux.data>0,'all'),Flux.units);
    disp(Flux.time.start);disp(Flux.time.stop);
    Prefix=['mms1_epd_feeps_srvy_l2_' particle{1}];
    sensors=[3:5 11:12];if strcmp(particle{1},'ion'),sensors=6:8;end
    for sensor=sensors
        A=mms.db_get_ts(['mms1_feeps_srvy_l2_' particle{1}],[Prefix '_top_intensity_sensorid_' num2str(sensor)],tint);
        if isempty(A),continue;end
        fprintf('sensor=%d rows=%d finite=%d positive=%d min=%g max=%g\n',sensor,size(A.data,1),sum(isfinite(A.data),'all'),sum(A.data>0,'all'),min(A.data,[],'all'),max(A.data,[],'all')); disp(A.userData.VALIDMIN);disp(A.userData.VALIDMAX); V=mms.db_get_variable(['mms1_feeps_srvy_l2_' particle{1}],[Prefix '_top_intensity_sensorid_' num2str(sensor)],tint);disp(V.FILLVAL);disp(unique(V.data(:,1))); 
    end
end
fprintf('PARTICLE_DIAGNOSTIC_FINISHED\n');