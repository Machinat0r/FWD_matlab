%% Raptis et al. (Nature, 2026) Juno bow-shock crossing
% Uses the same direct-JSS-RTP plus full-band Waves plotting program as the
% PJ67-to-PJ68 overview, showing 17:05-19:15 UTC with no event lines.

setenv('JUNO_PLOT_PRESET','RAPTIS_20231001_BOWSHOCK');
try
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'plot_juno_bowshock_pj67_pj68_direct_rtp_fullorbit.m'));
catch exception
    setenv('JUNO_PLOT_PRESET','');
    rethrow(exception)
end
setenv('JUNO_PLOT_PRESET','');
