function out = GSM2SM_BR(B, R)
%------written by Wending Fu, Sept.2025 in Beijing------------
% R: [N x 4] (R_E/km)
% out:
%       .mlt   [N x 1] 磁地方时 (小时, 0-24)
%       .mlat  [N x 1] 磁纬度 (度, [-90,90])
%       .L     [N x 1] L-shell (偶极近似)
%       .pos_sm_re [N x 3] 转到 SM 的坐标 (R_E)
%% input
% check input size
assert(size(B,2)==4, 'B must be [N x 4]');
assert(size(R,2)==4, 'R must be [N x 4]');
assert(size(B,1)==size(R,1), 'B & R row counts must match');
if max(abs(R(1,2:4))) >= 1000
    units = irf_units;
    R(:,2:4) = R(:,2:4) ./ units.RE .* 1e3;
end

time = B(:,1);
%%
% GSM → SM
pos_sm_re = irf_gsm2sm(R);
pos_sm_re = pos_sm_re(:,2:4);
B = irf_gsm2sm(B);
B = B(:,2:4);
%% calculate MLT / MLAT / L
x = pos_sm_re(:,1); y = pos_sm_re(:,2); z = pos_sm_re(:,3);
r = sqrt(x.^2 + y.^2 + z.^2);

phi = atan2(y, x);
theta = acos( max(-1,min(1, z ./ max(r,eps))) );

mlat = 90 - theta*180/pi;
mlt = 12 + (phi*180/pi)/15;
mlt = mod(mlt, 24);

L = r ./ (cosd(mlat).^2);

%% calculate Br Btheta Bphi
% e_r     = [sinθ cosφ, sinθ sinφ, cosθ]
% e_theta = [cosθ cosφ, cosθ sinφ, -sinθ]
% e_phi   = [-sinφ, cosφ, 0]

sinth = sin(theta); costh = cos(theta);
cosph = cos(phi);   sinph = sin(phi);

e_r     = [sinth.*cosph,  sinth.*sinph,  costh];
e_theta = [costh.*cosph,  costh.*sinph, -sinth];
e_phi   = [-sinph,         cosph,        0*phi];

Br     = sum(B .* e_r,     2);
Btheta = sum(B .* e_theta, 2);
Bphi   = sum(B .* e_phi,   2);

% ---------- 轴线上的退化处理（x≈0,y≈0） ----------
on_axis = hypot(x,y) < 1e-9;
Bphi(on_axis) = NaN; 
%% output
out = struct();
out.L      = L;
out.Bphi   = Bphi;
out.Br     = Br;
out.Btheta = Btheta;
out.MLT    = mlt;
out.MLat   = mlat;
out.pos_sm = pos_sm_re;
out.B_sm   = B;
end
