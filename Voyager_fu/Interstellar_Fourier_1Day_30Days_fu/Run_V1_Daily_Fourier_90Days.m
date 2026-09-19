function result=Run_V1_Daily_Fourier_90Days(varargin)
% Same 2--90 day band and [0,3] color scale as the latest daily Morlet figure.
% Retain the previously approved 90-day Hann windows, 1-day hop and window means.
% The 90-day period is the first nonzero FFT bin; long-period resolution is coarse.
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_1Day_90Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'ReferenceAuditFile','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_90Days/V1_daily_morlet_analysis.mat',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
result=Run_V1_Daily_Fourier('DataRoot',cfg.DataRoot,'OutputRoot',cfg.OutputRoot, ...
    'IRFURoot',cfg.IRFURoot,'ReferenceAuditFile',cfg.ReferenceAuditFile, ...
    'Visible',cfg.Visible,'WindowDays',90,'MaxPeriodDays',90,'PSDColorLimits',[0 3]);
end
