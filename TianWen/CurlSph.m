function [curlAphi, curlAtheta, curlAr] = CurlSph(Aphi, Atheta, Ar, phi, theta, r)
%------written by Wending Fu, Jun.2024 in Beijing------------
% In MAVEN-J statistics, the curl_theta/curl_phi may have fair-sized
% deivation due to the low slides on R direction.
%
% The partial derivative is calculated by finite difference
% see also: curl, gradient, calculateJ_sph
%%
narginchk(6,6)
phiGrid = mean(diff(phi,1,1), 'all');
thetaGrid = mean(diff(theta,1,2), 'all');
rGrid = mean(diff(r,1,3), 'all');
%%
% % % [~, dAr_dtheta, dAr_dphi] = gradient(Ar, rGrid, thetaGrid, phiGrid);
% % % [dAtheta_dr, ~, dAtheta_dphi] = gradient(Atheta, rGrid, thetaGrid, phiGrid);
% % % [dAphi_dr, dAphi_dtheta, ~] = gradient(Aphi, rGrid, thetaGrid, phiGrid);

[dAr_dphi, dAr_dtheta, ~] = gradient(Ar,phiGrid, thetaGrid, rGrid);
[dAtheta_dphi, ~, dAtheta_dr] = gradient(Atheta, phiGrid, thetaGrid, rGrid);
[~, dAphi_dtheta, dAphi_dr] = gradient(Aphi, phiGrid, thetaGrid, rGrid);

curlAr = (1 ./ (r .* cos(theta))) .* ...
         (dAphi_dtheta .* cos(theta) - sin(theta) .* Aphi - dAtheta_dphi);

curlAtheta = (1 ./ r) .* ...
             (1 ./ cos(theta) .* dAr_dphi - Aphi - r .* dAphi_dr);

curlAphi = (1 ./ r) .* ...
           (Atheta + r .* dAtheta_dr + tan(theta) .* Ar - dAr_dtheta);
end