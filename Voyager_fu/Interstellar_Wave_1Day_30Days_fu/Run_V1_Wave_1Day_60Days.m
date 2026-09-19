function result=Run_V1_Wave_1Day_60Days(varargin)
% User requested extending the daily Morlet lowest frequency to a 60-day period.
% Keep the original 79-day linear interpolation, 2-day Nyquist limit and [0,3] colors.
% The shared entry reads original CDF. Prior MAT products are never science input.
p=inputParser;
addParameter(p,'DataRoot','Z:/SPART-WORK/Data/Voyager',@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_60Days',@(x)ischar(x)||isstring(x));
addParameter(p,'IRFURoot','C:/Users/Administrator/Documents/irfu-matlab-master',@(x)ischar(x)||isstring(x));
addParameter(p,'Visible',true,@islogical);
parse(p,varargin{:}); cfg=p.Results;
result=Run_V1_Wave_1Day_30Days('DataRoot',cfg.DataRoot, ...
    'OutputRoot',cfg.OutputRoot,'IRFURoot',cfg.IRFURoot,'Visible',cfg.Visible,'MaxPeriodDays',60);
end
