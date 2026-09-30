%% FEEPS负值与能量轴的只读核查
addpath('C:\Users\Administrator\Documents\irfu-matlab-master');irf('check_path');mms.db_init('local_file_db','Z:\SPART-WORK\Data\MMS');
tint=irf.tint('2026-07-28T05:10:00Z/2026-07-28T05:30:00Z');
for mode={'brst','srvy'}
    for particle={'electron','ion'}
        F=mms.get_data(['Omniflux' particle{1} '_epd_feeps_' mode{1} '_l2'],tint,1);
        if isempty(F),continue;end
        fprintf('%s %s E_min=%g E_max=%g negative=%d min=%g max=%g\n',mode{1},particle{1},min(F.depend{1},[],'all'),max(F.depend{1},[],'all'),sum(F.data<0,'all'),min(F.data,[],'all'),max(F.data,[],'all'));
    end
end
for sensor=[3 4 11]
 V=mms.db_get_variable('mms1_feeps_srvy_l2_electron',['mms1_epd_feeps_srvy_l2_electron_top_intensity_sensorid_' num2str(sensor)],tint);
 fprintf('sensor=%d dataClass=%s fillClass=%s min=%.17g FILL=%.17g\n',sensor,class(V.data),class(V.FILLVAL),min(V.data,[],'all'),V.FILLVAL);
end
fprintf('NEGATIVE_CHECK_FINISHED\n');