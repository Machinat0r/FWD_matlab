function out = calc_electron_rotating_vorticity_phi_3sc( ...
    R1, R2, R3, ...
    ne1, ne2, ne3, ...
    pe1, pe2, pe3, ...
    Ve1, Ve2, Ve3, ...
    B1, B2, B3, ...
    E1, E2, E3, ...
    Lvec, Mvec, Nvec, Omega0)
% calc_electron_rotating_vorticity_phi_3sc
%
% Rotating-frame electron vorticity diagnostic using three MMS spacecraft.
%
% This function is designed to address the concern that a Taylor-Proudman
% type discussion should explicitly include the rotating-frame Coriolis term.
%
% Local coordinates:
%   er   = N
%   ephi = L
%   ez   = M
%
% Main equation tested diagnostically:
%
%   [curl(2 Omega0 x Ve)]_phi
%       =
%   [(grad ne x grad pe)/(me ne^2)]_phi
%       +
%   (qe/me) [curl(E + Ve x B)]_phi
%
% under local 2D approximation:
%   d/dphi = 0
%
% Therefore:
%   [curl(2 Omega0 x Ve)]_phi = -2 Omega0 dVe_phi/dz
%
%   Baro_phi = 1/(me ne^2) * (dne/dz*dpe/dr - dne/dr*dpe/dz)
%
%   EM_phi = (qe/me) * [d(E+Ve x B)_r/dz - d(E+Ve x B)_z/dr]
%
% Inputs:
%   R1,R2,R3       [time, x, y, z], km, GSM
%   ne1,ne2,ne3   [time, ne], cm^-3
%   pe1,pe2,pe3   [time, pe] or [time, Pxx, Pyy, Pzz], nPa
%   Ve1,Ve2,Ve3   [time, Vex, Vey, Vez], km/s, GSM
%   B1,B2,B3      [time, Bx, By, Bz], nT, GSM
%   E1,E2,E3      [time, Ex, Ey, Ez], mV/m, GSM
%   Lvec,Mvec,Nvec local basis vectors in GSM
%   Omega0        assumed reference rotation rate, rad/s
%
% Outputs:
%   out.coriolis_curl_phi : [curl(2 Omega0 x Ve)]_phi, s^-2
%   out.baro_phi          : baroclinic term, s^-2
%   out.em_phi            : electromagnetic term, s^-2
%   out.residual_phi      : coriolis_curl_phi - baro_phi - em_phi
%   out.dVephi_dz         : axial gradient of Ve_phi, s^-1
%   out.condG             : condition number of three-spacecraft 2D fit
%
% Important:
%   Omega0 is not independently known for this MMS event. Therefore the
%   result should be described as a rotating-frame diagnostic or sensitivity
%   test, not as a direct proof of a classical Taylor-Proudman state.

%% constants
units = irf_units;
qe = -units.e;
me = units.me;

%% reference time
t = R1(:,1);
Nt = length(t);

%% interpolate all inputs to R1 time
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

E1_xyz  = interp_vec_to_t(E1,  t);
E2_xyz  = interp_vec_to_t(E2,  t);
E3_xyz  = interp_vec_to_t(E3,  t);

%% local orthonormal basis
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

%% SI units
R1_m  = R1_xyz .* 1e3;
R2_m  = R2_xyz .* 1e3;
R3_m  = R3_xyz .* 1e3;

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

E1_SI = E1_xyz .* 1e-3;
E2_SI = E2_xyz .* 1e-3;
E3_SI = E3_xyz .* 1e-3;

%% outputs
coriolis_curl_phi = nan(Nt,1);
baro_phi          = nan(Nt,1);
em_phi            = nan(Nt,1);
residual_phi      = nan(Nt,1);

dVephi_dz = nan(Nt,1);
dVephi_dr = nan(Nt,1);

Ve_phi_mean = nan(Nt,1);
Ve_r_mean   = nan(Nt,1);
Ve_z_mean   = nan(Nt,1);

grad_ne_rz   = nan(Nt,2);
grad_pe_rz   = nan(Nt,2);
grad_Eeff_rz = nan(Nt,4);   % [dEeff_r/dr, dEeff_r/dz, dEeff_z/dr, dEeff_z/dz]

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

    Etmp = [
        E1_SI(it,:);
        E2_SI(it,:);
        E3_SI(it,:)
    ];

    if any(isnan(Rtmp(:))) || any(isnan(ne_sc)) || any(isnan(pe_sc)) || ...
       any(isnan(Vetmp(:))) || any(isnan(Btmp(:))) || any(isnan(Etmp(:)))
        continue
    end

    if mean(ne_sc,'omitnan') <= 0
        continue
    end

    %% local spacecraft positions
    Rloc = Rtmp * Q;
    r_sc = Rloc(:,1);
    z_sc = Rloc(:,3);

    G = [r_sc - mean(r_sc,'omitnan'), z_sc - mean(z_sc,'omitnan')];
    if rank(G) < 2
        continue
    end
    condG(it) = cond(G);

    %% local vectors
    Ve_loc = Vetmp * Q;
    B_loc  = Btmp  * Q;
    E_loc  = Etmp  * Q;

    Ve_r   = Ve_loc(:,1);
    Ve_phi = Ve_loc(:,2);
    Ve_z   = Ve_loc(:,3);

    Ve_r_mean(it)   = mean(Ve_r,   'omitnan');
    Ve_phi_mean(it) = mean(Ve_phi, 'omitnan');
    Ve_z_mean(it)   = mean(Ve_z,   'omitnan');

    %% effective electric field in electron momentum equation
    % Eeff = E + Ve x B, unit V/m
    Eeff_gsm = Etmp + cross(Vetmp, Btmp, 2);
    Eeff_loc = Eeff_gsm * Q;

    Eeff_r = Eeff_loc(:,1);
    Eeff_z = Eeff_loc(:,3);

    %% gradients in local r-z plane
    grad_ne     = grad2D_3sc(r_sc, z_sc, ne_sc);
    grad_pe     = grad2D_3sc(r_sc, z_sc, pe_sc);
    grad_Ve_phi = grad2D_3sc(r_sc, z_sc, Ve_phi);
    grad_Eeff_r = grad2D_3sc(r_sc, z_sc, Eeff_r);
    grad_Eeff_z = grad2D_3sc(r_sc, z_sc, Eeff_z);

    grad_ne_rz(it,:) = grad_ne;
    grad_pe_rz(it,:) = grad_pe;

    grad_Eeff_rz(it,:) = [ ...
        grad_Eeff_r(1), grad_Eeff_r(2), ...
        grad_Eeff_z(1), grad_Eeff_z(2)];

    dne_dr = grad_ne(1);
    dne_dz = grad_ne(2);

    dpe_dr = grad_pe(1);
    dpe_dz = grad_pe(2);

    dVephi_dr(it) = grad_Ve_phi(1);
    dVephi_dz(it) = grad_Ve_phi(2);

    dEeff_r_dz = grad_Eeff_r(2);
    dEeff_z_dr = grad_Eeff_z(1);

    ne0 = mean(ne_sc, 'omitnan');

    %% Coriolis curl phi component
    %
    % [curl(2 Omega0 x Ve)]_phi = -2 Omega0 dVe_phi/dz
    %
    coriolis_curl_phi(it) = -2 .* Omega0 .* dVephi_dz(it);

    %% Baroclinic phi component
    %
    % [(grad ne x grad pe)/(me ne^2)]_phi
    %
    baro_phi(it) = ...
        (dne_dz .* dpe_dr - dne_dr .* dpe_dz) ./ (me .* ne0.^2);

    %% Electromagnetic phi component
    %
    % (qe/me) [curl(E + Ve x B)]_phi
    %
    % under d/dphi = 0:
    % [curl(Eeff)]_phi = dEeff_r/dz - dEeff_z/dr
    %
    em_phi(it) = ...
        (qe ./ me) .* (dEeff_r_dz - dEeff_z_dr);

    %% residual
    residual_phi(it) = coriolis_curl_phi(it) - baro_phi(it) - em_phi(it);

end

%% output
out.t = t;

out.coriolis_curl_phi = coriolis_curl_phi;
out.baro_phi          = baro_phi;
out.em_phi            = em_phi;
out.residual_phi      = residual_phi;

out.dVephi_dz = dVephi_dz;
out.dVephi_dr = dVephi_dr;

out.Ve_phi_mean = Ve_phi_mean;
out.Ve_r_mean   = Ve_r_mean;
out.Ve_z_mean   = Ve_z_mean;

out.grad_ne_rz   = grad_ne_rz;
out.grad_pe_rz   = grad_pe_rz;
out.grad_Eeff_rz = grad_Eeff_rz;

out.condG = condG;

out.Omega0 = Omega0;

out.er   = er;
out.ephi = ephi;
out.ez   = ez;
out.Q    = Q;

end


function vec_out = interp_vec_to_t(data_in, t_ref)
% data_in: n x 4, [time, x, y, z]

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
%   n x 4, [time, Pxx, Pyy, Pzz]
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
% f(r,z) = f0 + df/dr * (r-r0) + df/dz * (z-z0)

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