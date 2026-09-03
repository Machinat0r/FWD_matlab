function pointing = Case1_Predicted_LECP_Pointing(epoch, B_RTN, cfg)
%Case1_Predicted_LECP_Pointing P1 central-ray angles with predicted attitude.
%   Uses the original LECP Epoch and its time-matched magnetic vector B.
%   No interpolation of B/flux, no uN=0 and no finite-aperture deconvolution.

%% retain the original Epoch; do not round to noon or an hour midpoint
if nargin < 3, cfg = Case1_Config; end
validateattributes(B_RTN, {'numeric'}, {'2d', 'ncols', 3});
if size(B_RTN, 1) ~= numel(epoch)
    error('VoyagerCase1:PointingSize', 'B and epoch lengths differ.');
end
time = epoch(:);
time.TimeZone = 'UTC';
if ~isfield(cfg, 'PADCadence'), cfg.PADCadence = 'day'; end
if strcmp(cfg.PADCadence, 'hour')
    if ~isfield(cfg, 'HourlyAttitudeApproved') || ~cfg.HourlyAttitudeApproved
        error('VoyagerCase1:HourlyApproval', ...
            'Hourly predicted-attitude evaluation requires explicit approval.');
    end
end
timeDescription = sprintf('original LECP Epoch direction dot B in the same UTC %s', ...
    cfg.PADCadence);
geometry = Case1_LECP_Geometry(cfg);
attitude = Case1_Read_Predicted_Attitude(time, cfg);

%% transform all three components; retain missing fields as NaN
nTime = numel(time);
particleRTN = nan(nTime, 8, 3);
lookRTN = nan(nTime, 8, 3);
mu = nan(nTime, 8);
alpha = nan(nTime, 8);
for ii = 1:nTime
    if ~attitude.Found(ii), continue; end
    C = attitude.C_SC_to_RTN(:, :, ii);
    u = C * geometry.ParticleSC;
    b = C * geometry.LookSC;
    particleRTN(ii, :, :) = reshape(u.', 1, 8, 3);
    lookRTN(ii, :, :) = reshape(b.', 1, 8, 3);
    B = B_RTN(ii, :).';
    if ~all(isfinite(B)) || norm(B) == 0, continue; end
    % v points into the aperture, while LookSC points out of it.
    mu(ii, :) = max(-1, min(1, (B.'/norm(B))*u));
    alpha(ii, :) = acosd(mu(ii, :));
end

%% audit accompanies the plotted values, independently of the figure label
pointing = struct;
pointing.TimeUTC = time;
pointing.ParticleRTN = particleRTN;
pointing.LookRTN = lookRTN;
pointing.Mu = mu;
pointing.PitchAngle_deg = alpha;
pointing.Available = all(isfinite(alpha(:, 1:7)), 2);
pointing.Attitude = attitude;
pointing.Geometry = geometry;
pointing.Method = sprintf(['Official predicted CK; %s; ', ...
    'official nominal P1 sector centers; v=-look'], timeDescription);
end
