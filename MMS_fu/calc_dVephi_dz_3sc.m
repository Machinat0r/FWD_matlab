function out = calc_dVephi_dz_3sc_mc( ...
    R1, R2, R3, ...
    Ve1, Ve2, Ve3, ...
    Lvec, Mvec, Nvec, sigma_Ve_kms, Nmc)
% calc_dVephi_dz_3sc_mc
%
% Calculate dVe_phi/dz and its uncertainty using three MMS spacecraft.
%
% sigma_Ve_kms:
%   assumed 1-sigma uncertainty of electron velocity, km/s.
%   If unknown, try values such as 20, 30, or 50 km/s as sensitivity tests.
%
% Nmc:
%   number of Monte Carlo runs, e.g., 1000 or 5000.
%
% Output:
%   out.dVephi_dz       : best estimate, s^-1
%   out.dVephi_dz_std   : Monte Carlo 1-sigma uncertainty, s^-1
%   out.dVephi_dz_sig   : |value|/sigma
%   out.is_zero_1sigma  : true if |value| < 1 sigma
%   out.is_zero_2sigma  : true if |value| < 2 sigma

if nargin < 10 || isempty(sigma_Ve_kms)
    sigma_Ve_kms = 30; 
end

if nargin < 11 || isempty(Nmc)
    Nmc = 2000;
end

%% reference time
t = R1(:,1);
Nt = length(t);

%% interpolate to R1 time
R1_xyz  = interp_vec_to_t(R1,  t);
R2_xyz  = interp_vec_to_t(R2,  t);
R3_xyz  = interp_vec_to_t(R3,  t);

Ve1_xyz = interp_vec_to_t(Ve1, t);
Ve2_xyz = interp_vec_to_t(Ve2, t);
Ve3_xyz = interp_vec_to_t(Ve3, t);

%% local basis
er = Nvec(:) ./ norm(Nvec);

ephi = Lvec(:) - dot(Lvec(:), er).*er;
ephi = ephi ./ norm(ephi);

ez = cross(er, ephi);
ez = ez ./ norm(ez);

Munit = Mvec(:) ./ norm(Mvec);
if dot(ez, Munit) < 0
    ez = -ez;
end

Q = [er, ephi, ez];

%% units
R1_m = R1_xyz .* 1e3;
R2_m = R2_xyz .* 1e3;
R3_m = R3_xyz .* 1e3;

Ve1_SI = Ve1_xyz .* 1e3;
Ve2_SI = Ve2_xyz .* 1e3;
Ve3_SI = Ve3_xyz .* 1e3;

sigma_Ve = sigma_Ve_kms .* 1e3;

%% output arrays
dVephi_dz     = nan(Nt,1);
dVephi_dr     = nan(Nt,1);
dVephi_dz_std = nan(Nt,1);
dVephi_dz_sig = nan(Nt,1);
condX         = nan(Nt,1);

for it = 1:Nt

    Rtmp = [
        R1_m(it,:);
        R2_m(it,:);
        R3_m(it,:)
    ];

    Vetmp = [
        Ve1_SI(it,:);
        Ve2_SI(it,:);
        Ve3_SI(it,:)
    ];

    if any(isnan(Rtmp(:))) || any(isnan(Vetmp(:)))
        continue
    end

    %% local coordinates
    Rloc = Rtmp * Q;
    Veloc = Vetmp * Q;

    r = Rloc(:,1);
    z = Rloc(:,3);

    Vephi = Veloc(:,2);

    %% direct estimate
    X = [ones(3,1), r, z];

    if rank(X) < 3
        continue
    end

    condX(it) = cond(X);

    beta = X \ Vephi;

    dVephi_dr(it) = beta(2);
    dVephi_dz(it) = beta(3);

    %% Monte Carlo uncertainty
    dtmp = nan(Nmc,1);

    for imc = 1:Nmc

        Vephi_mc = Vephi + sigma_Ve .* randn(3,1);

        beta_mc = X \ Vephi_mc;

        dtmp(imc) = beta_mc(3);

    end

    dVephi_dz_std(it) = std(dtmp, 'omitnan');

    if dVephi_dz_std(it) > 0
        dVephi_dz_sig(it) = abs(dVephi_dz(it)) ./ dVephi_dz_std(it);
    end

end

out.t = t;

out.dVephi_dz     = dVephi_dz;
out.dVephi_dr     = dVephi_dr;
out.dVephi_dz_std = dVephi_dz_std;
out.dVephi_dz_sig = dVephi_dz_sig;

out.is_zero_1sigma = abs(dVephi_dz) < dVephi_dz_std;
out.is_zero_2sigma = abs(dVephi_dz) < 2 .* dVephi_dz_std;

out.omega_er = -dVephi_dz;
out.condX = condX;

out.er = er;
out.ephi = ephi;
out.ez = ez;
out.Q = Q;

end


function vec_out = interp_vec_to_t(data_in, t_ref)

if size(data_in,2) < 4
    error('Vector input should be n x 4: [time, x, y, z].');
end

t_in = data_in(:,1);
x_in = data_in(:,2:4);

vec_out = interp1(t_in, x_in, t_ref, 'linear', NaN);

end