function out = calc_electron_canonical_vorticity_phi_3sc( ...
    R1, R2, R3, ...
    ne1, ne2, ne3, ...
    pe1, pe2, pe3, ...
    Ve1, Ve2, Ve3, ...
    B1, B2, B3, ...
    Lvec, Mvec, Nvec, Omega0)
% calc_electron_canonical_vorticity_rotframe_phi_3sc
%
% Electron canonical vorticity diagnostic in a rotating frame using
% three MMS spacecraft.
%
% 重要定义：
%   W_e = omega_e + (qe/me)*B
%
% 注意：
%   这里不定义 W_e,rot = omega_e + 2Omega0 + (qe/me)*B。
%   旋转参考系的贡献作为方程右端的单独项加入：
%
%   dW_e/dt =
%       curl(Ve x W_e)
%     + curl(Ve x 2Omega0)
%     + (grad ne x grad pe)/(me ne^2)
%
% Local coordinates:
%   er   = N
%   ephi = L
%   ez   = M
%
% Local 2D approximation:
%   d/dphi = 0
%
% Electron vorticity:
%   omega_e = curl(Ve)
%
% In local 2D Cartesian approximation:
%   omega_r   = - dVe_phi/dz
%   omega_phi =   dVe_r/dz - dVe_z/dr
%   omega_z   =   dVe_phi/dr
%
% Rotating-frame term:
%   Omega0 vector = Omega0 * ez
%
%   [curl(Ve x 2Omega0)]_phi
%       = 2*Omega0*dVe_phi/dz
%
% Phi component:
%   dW_phi/dt =
%       [curl(Ve x W_e)]_phi
%     + [curl(Ve x 2Omega0)]_phi
%     + Baro_phi
%
% Baro_phi:
%   Baro_phi = 1/(me ne^2) * (dne/dz*dpe/dr - dne/dr*dpe/dz)
%
% Inputs:
%   R1,R2,R3       : [time, Rx, Ry, Rz], km, GSM
%   ne1,ne2,ne3   : [time, ne], cm^-3
%   pe1,pe2,pe3   : [time, pe] or [time, Pxx, Pyy, Pzz], nPa
%   Ve1,Ve2,Ve3   : [time, Vex, Vey, Vez], km/s, GSM
%   B1,B2,B3      : [time, Bx, By, Bz], nT, GSM
%   Lvec,Mvec,Nvec: local basis vectors in GSM
%   Omega0        : rotation rate of the chosen rotating frame, rad/s
%
% Outputs:
%   out.omega_r, out.omega_phi, out.omega_z
%   out.mag_r,   out.mag_phi,   out.mag_z
%   out.W_r,     out.W_phi,     out.W_z
%   out.dVephi_dz
%   out.dWphi_dt
%   out.conv_phi       = [curl(Ve x W_e)]_phi
%   out.rot_phi        = [curl(Ve x 2Omega0)]_phi
%   out.Vexomegae_phi = [curl(Ve x omega_e)]_phi
%   out.baro_phi
%   out.res_phi
%
% Note:
%   This is a local 2D first-order diagnostic. It is not a full 3D
%   closure of the canonical vorticity equation.

%% constants
units = irf_units;
qe = -units.e;
me = units.me;

%% use R1 time as reference
t = R1(:,1);
Nt = length(t);

%% interpolate all data onto R1 time
R1_xyz  = interp_vec_to_t(R1,  t);
R2_xyz  = interp_vec_to_t(R2,  t);
R3_xyz  = interp_vec_to_t(R3,  t);

ne1_s   = interp_scalar_to_t(ne1, t);
ne2_s   = interp_scalar_to_t(ne2, t);
ne3_s   = interp_scalar_to_t(ne3, t);

pe1_s   = interp_pressure_to_t(pe1, t);
pe2_s   = interp_pressure_to_t(pe2, t);
pe3_s   = interp_pressure_to_t(pe3, t);

Ve1_xyz = interp_vec_to_t(Ve1, t);
Ve2_xyz = interp_vec_to_t(Ve2, t);
Ve3_xyz = interp_vec_to_t(Ve3, t);

B1_xyz  = interp_vec_to_t(B1,  t);
B2_xyz  = interp_vec_to_t(B2,  t);
B3_xyz  = interp_vec_to_t(B3,  t);

%% construct local orthonormal basis
er = Nvec(:) ./ norm(Nvec);

ephi = Lvec(:) - dot(Lvec(:), er).*er;
ephi = ephi ./ norm(ephi);

ez_tmp = cross(er, ephi);
ez_tmp = ez_tmp ./ norm(ez_tmp);

Munit = Mvec(:) ./ norm(Mvec);
if dot(ez_tmp, Munit) < 0
    ez_tmp = -ez_tmp;
end
ez = ez_tmp;

% GSM vector A to local vector:
% A_local = A_GSM * Q
Q = [er, ephi, ez];

%% convert units to SI
R1_m = R1_xyz .* 1e3;
R2_m = R2_xyz .* 1e3;
R3_m = R3_xyz .* 1e3;

ne1_m = ne1_s .* 1e6;
ne2_m = ne2_s .* 1e6;
ne3_m = ne3_s .* 1e6;

pe1_Pa = pe1_s .* 1e-9;
pe2_Pa = pe2_s .* 1e-9;
pe3_Pa = pe3_s .* 1e-9;

Ve1_SI = Ve1_xyz .* 1e3;
Ve2_SI = Ve2_xyz .* 1e3;
Ve3_SI = Ve3_xyz .* 1e3;

B1_SI = B1_xyz .* 1e-9;
B2_SI = B2_xyz .* 1e-9;
B3_SI = B3_xyz .* 1e-9;

%% initialize outputs
omega_r   = nan(Nt,1);
omega_phi = nan(Nt,1);
omega_z   = nan(Nt,1);

mag_r   = nan(Nt,1);
mag_phi = nan(Nt,1);
mag_z   = nan(Nt,1);

W_r   = nan(Nt,1);
W_phi = nan(Nt,1);
W_z   = nan(Nt,1);

dVephi_dz = nan(Nt,1);
dVephi_dr = nan(Nt,1);

conv_phi = nan(Nt,1);
rot_phi  = nan(Nt,1);
baro_phi = nan(Nt,1);

Bmean_local  = nan(Nt,3);
Vemean_local = nan(Nt,3);

grad_Ver_rz   = nan(Nt,2);
grad_Vephi_rz = nan(Nt,2);
grad_Vez_rz   = nan(Nt,2);

grad_ne_rz = nan(Nt,2);
grad_pe_rz = nan(Nt,2);

grad_Ar_rz = nan(Nt,2);
grad_Az_rz = nan(Nt,2);

Vexomegae_phi = nan(Nt,1);
grad_Aomega_r_rz = nan(Nt,2);
grad_Aomega_z_rz = nan(Nt,2);

condG = nan(Nt,1);

%% main loop
for it = 1:Nt

    Rtmp = [
        R1_m(it,:);
        R2_m(it,:);
        R3_m(it,:)
    ];

    ne_sc = [
        ne1_m(it);
        ne2_m(it);
        ne3_m(it)
    ];

    pe_sc = [
        pe1_Pa(it);
        pe2_Pa(it);
        pe3_Pa(it)
    ];

    Vetmp = [
        Ve1_SI(it,:);
        Ve2_SI(it,:);
        Ve3_SI(it,:)
    ];

    Btmp = [
        B1_SI(it,:);
        B2_SI(it,:);
        B3_SI(it,:)
    ];

    if any(isnan(Rtmp(:))) || any(isnan(ne_sc)) || any(isnan(pe_sc)) || ...
       any(isnan(Vetmp(:))) || any(isnan(Btmp(:)))
        continue
    end

    if mean(ne_sc,'omitnan') <= 0
        continue
    end

    %% local positions
    Rloc = Rtmp * Q;

    r_sc = Rloc(:,1);
    z_sc = Rloc(:,3);

    Gtmp = [r_sc - mean(r_sc,'omitnan'), z_sc - mean(z_sc,'omitnan')];

    if rank(Gtmp) < 2
        continue
    end

    condG(it) = cond(Gtmp);

    %% local velocity and magnetic field
    Ve_loc = Vetmp * Q;
    B_loc  = Btmp  * Q;

    Ver   = Ve_loc(:,1);
    Vephi = Ve_loc(:,2);
    Vez   = Ve_loc(:,3);

    B0 = mean(B_loc, 1, 'omitnan');

    Vemean_local(it,:) = mean(Ve_loc, 1, 'omitnan');
    Bmean_local(it,:)  = B0;

    %% gradients of velocity components
    grad_Ver   = grad2D_3sc(r_sc, z_sc, Ver);
    grad_Vephi = grad2D_3sc(r_sc, z_sc, Vephi);
    grad_Vez   = grad2D_3sc(r_sc, z_sc, Vez);

    grad_Ver_rz(it,:)   = grad_Ver;
    grad_Vephi_rz(it,:) = grad_Vephi;
    grad_Vez_rz(it,:)   = grad_Vez;

    dVer_dz = grad_Ver(2);

    dVephi_dr(it) = grad_Vephi(1);
    dVephi_dz(it) = grad_Vephi(2);

    dVez_dr = grad_Vez(1);

    %% electron vorticity under local 2D approximation
    %
    % omega_r   = - dVephi/dz
    % omega_phi =   dVer/dz - dVez/dr
    % omega_z   =   dVephi/dr
    %
    omega_r(it)   = -dVephi_dz(it);
    omega_phi(it) =  dVer_dz - dVez_dr;
    omega_z(it)   =  dVephi_dr(it);

    %% magnetic part of canonical vorticity
    %
    % mag = (qe/me)*B
    %
    mag_r(it)   = (qe/me) .* B0(1);
    mag_phi(it) = (qe/me) .* B0(2);
    mag_z(it)   = (qe/me) .* B0(3);

    %% electron canonical vorticity at centroid
    %
    % Still define:
    %   W_e = omega_e + (qe/me)*B
    %
    % Do not include 2Omega0 in W_e.
    %
    W_r(it)   = omega_r(it)   + mag_r(it);
    W_phi(it) = omega_phi(it) + mag_phi(it);
    W_z(it)   = omega_z(it)   + mag_z(it);

    %% Baroclinic term
    grad_ne = grad2D_3sc(r_sc, z_sc, ne_sc);
    grad_pe = grad2D_3sc(r_sc, z_sc, pe_sc);

    grad_ne_rz(it,:) = grad_ne;
    grad_pe_rz(it,:) = grad_pe;

    dne_dr = grad_ne(1);
    dne_dz = grad_ne(2);

    dpe_dr = grad_pe(1);
    dpe_dz = grad_pe(2);

    ne0 = mean(ne_sc, 'omitnan');

    baro_phi(it) = ...
        (dne_dz .* dpe_dr - dne_dr .* dpe_dz) ./ (me .* ne0.^2);

    %% Convective term of canonical vorticity
    %
    % A = Ve x W_e
    %
    % In the strict 2D first-order estimate, omega_e is evaluated at the
    % centroid and assigned to each spacecraft, while magnetic part uses
    % each spacecraft B_i.
    %
    % W_i = omega_centroid + (qe/me)*B_i
    %
    omega_centroid = [omega_r(it), omega_phi(it), omega_z(it)];

    W_sc = zeros(3,3);
    for isc = 1:3
        W_sc(isc,:) = omega_centroid + (qe/me).*B_loc(isc,:);
    end

    A_sc = cross(Ve_loc, W_sc, 2);

    Ar = A_sc(:,1);
    Az = A_sc(:,3);

    grad_Ar = grad2D_3sc(r_sc, z_sc, Ar);
    grad_Az = grad2D_3sc(r_sc, z_sc, Az);

    grad_Ar_rz(it,:) = grad_Ar;
    grad_Az_rz(it,:) = grad_Az;

    dAr_dz = grad_Ar(2);
    dAz_dr = grad_Az(1);

    conv_phi(it) = dAr_dz - dAz_dr;

    %% Convective term of ordinary electron vorticity
    %
    % Aomega = Ve x omega_e
    %
    % In the same 2D first-order estimate as above, omega_e is evaluated at
    % the centroid and assigned to each spacecraft. Therefore the spatial
    % variation of Aomega comes from the measured electron velocity.
    %
    Aomega_sc = cross(Ve_loc, repmat(omega_centroid,3,1), 2);

    Aomega_r = Aomega_sc(:,1);
    Aomega_z = Aomega_sc(:,3);

    grad_Aomega_r = grad2D_3sc(r_sc, z_sc, Aomega_r);
    grad_Aomega_z = grad2D_3sc(r_sc, z_sc, Aomega_z);

    grad_Aomega_r_rz(it,:) = grad_Aomega_r;
    grad_Aomega_z_rz(it,:) = grad_Aomega_z;

    dAomega_r_dz = grad_Aomega_r(2);
    dAomega_z_dr = grad_Aomega_z(1);

    Vexomegae_phi(it) = dAomega_r_dz - dAomega_z_dr;

    %% Rotating-frame correction term
    %
    % Equation uses:
    %   + curl(Ve x 2Omega0)
    %
    % with Omega0 vector = Omega0 * ez.
    %
    % Under local 2D approximation:
    %   [curl(Ve x 2Omega0)]_phi = 2*Omega0*dVe_phi/dz
    %
    % 注意符号：
    %   curl(2Omega0 x Ve)_phi = -2*Omega0*dVe_phi/dz
    %   but the vorticity equation contains curl(Ve x 2Omega0).
    %
    rot_phi(it) = 2 .* Omega0 .* dVephi_dz(it);

end

%% time derivative of W_phi
dWphi_dt = nan(Nt,1);

valid = ~isnan(W_phi) & ~isnan(t);
if sum(valid) >= 3
    Wtmp = W_phi;
    ttmp = t;

    % light smoothing may reduce FPI-scale noise
    % set smooth_N = 1 if you want raw derivative
    smooth_N = 3;

    Wtmp_s = Wtmp;
    Wtmp_s(valid) = movmean(Wtmp(valid), smooth_N, 'omitnan');

    dWphi_dt(valid) = gradient(Wtmp_s(valid), ttmp(valid));
end

%% residual in rotating frame
%
% dW_phi/dt = conv_phi + rot_phi + baro_phi
%
res_phi = dWphi_dt - conv_phi - rot_phi - baro_phi;

%% diagnostic ratios
eps0 = 1e-30;

mag_over_omega_r   = abs(mag_r)   ./ (abs(omega_r)   + eps0);
mag_over_omega_phi = abs(mag_phi) ./ (abs(omega_phi) + eps0);
mag_over_omega_z   = abs(mag_z)   ./ (abs(omega_z)   + eps0);

comp_r   = abs(W_r)   ./ (abs(omega_r)   + abs(mag_r)   + eps0);
comp_phi = abs(W_phi) ./ (abs(omega_phi) + abs(mag_phi) + eps0);
comp_z   = abs(W_z)   ./ (abs(omega_z)   + abs(mag_z)   + eps0);

sign_comp_r   = omega_r   .* mag_r;
sign_comp_phi = omega_phi .* mag_phi;
sign_comp_z   = omega_z   .* mag_z;

%% output
out.t = t;

out.omega_r   = omega_r;
out.omega_phi = omega_phi;
out.omega_z   = omega_z;

out.mag_r   = mag_r;
out.mag_phi = mag_phi;
out.mag_z   = mag_z;

out.W_r   = W_r;
out.W_phi = W_phi;
out.W_z   = W_z;

out.dVephi_dz = dVephi_dz;
out.dVephi_dr = dVephi_dr;

out.dWphi_dt = dWphi_dt;
out.conv_phi = conv_phi;
out.rot_phi  = rot_phi;
out.Vexomegae_phi = Vexomegae_phi;
out.baro_phi = baro_phi;
out.res_phi  = res_phi;

out.conv_plus_rot_phi = conv_phi + rot_phi;
out.VexBmag_phi = conv_phi - Vexomegae_phi;

out.mag_over_omega_r   = mag_over_omega_r;
out.mag_over_omega_phi = mag_over_omega_phi;
out.mag_over_omega_z   = mag_over_omega_z;

out.comp_r   = comp_r;
out.comp_phi = comp_phi;
out.comp_z   = comp_z;

out.sign_comp_r   = sign_comp_r;
out.sign_comp_phi = sign_comp_phi;
out.sign_comp_z   = sign_comp_z;

out.Bmean_local  = Bmean_local;
out.Vemean_local = Vemean_local;

out.grad_Ver_rz   = grad_Ver_rz;
out.grad_Vephi_rz = grad_Vephi_rz;
out.grad_Vez_rz   = grad_Vez_rz;

out.grad_ne_rz = grad_ne_rz;
out.grad_pe_rz = grad_pe_rz;

out.grad_Ar_rz = grad_Ar_rz;
out.grad_Az_rz = grad_Az_rz;

out.grad_Aomega_r_rz = grad_Aomega_r_rz;
out.grad_Aomega_z_rz = grad_Aomega_z_rz;

out.condG = condG;

out.Omega0 = Omega0;

out.er   = er;
out.ephi = ephi;
out.ez   = ez;
out.Q    = Q;

end


function vec_out = interp_vec_to_t(data_in, t_ref)
% data_in: n x 4, [time, x, y, z]
% vec_out: length(t_ref) x 3

if size(data_in,2) < 4
    error('Vector input should be n x 4: [time, x, y, z].');
end

t_in = data_in(:,1);
x_in = data_in(:,2:4);

vec_out = interp1(t_in, x_in, t_ref, 'linear', NaN);

end


function scalar_out = interp_scalar_to_t(data_in, t_ref)
% data_in: n x 2, [time, scalar]

if size(data_in,2) < 2
    error('Scalar input should be n x 2: [time, scalar].');
end

t_in = data_in(:,1);
s_in = data_in(:,2);

scalar_out = interp1(t_in, s_in, t_ref, 'linear', NaN);

end


function p_out = interp_pressure_to_t(data_in, t_ref)
% data_in:
%   n x 2, [time, scalar pressure]
%   n x 4 or more, [time, Pxx, Pyy, Pzz, ...]
%
% output:
%   scalar pressure

if size(data_in,2) == 2
    t_in = data_in(:,1);
    p_in = data_in(:,2);

elseif size(data_in,2) >= 4
    t_in = data_in(:,1);
    p_in = mean(data_in(:,2:4), 2, 'omitnan');

else
    error('Pressure input should be n x 2 or n x 4.');
end

p_out = interp1(t_in, p_in, t_ref, 'linear', NaN);

end


function grad_f = grad2D_3sc(r_sc, z_sc, f_sc)
% Estimate local 2D gradient of scalar f using three spacecraft.
%
% f(r,z) = f0 + df/dr*(r-r0) + df/dz*(z-z0)
%
% output:
%   grad_f = [df/dr, df/dz]

r_sc = r_sc(:);
z_sc = z_sc(:);
f_sc = f_sc(:);

r0 = mean(r_sc, 'omitnan');
z0 = mean(z_sc, 'omitnan');
f0 = mean(f_sc, 'omitnan');

G  = [r_sc - r0, z_sc - z0];
df = f_sc - f0;

if rank(G) < 2
    grad_f = [NaN, NaN];
    return
end

grad_f_col = G \ df;
grad_f = grad_f_col(:).';

end