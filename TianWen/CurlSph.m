function [curlAphi, curlAtheta, curlAr] = CurlSph(Aphi, Atheta, Ar, phi, theta, r)
%------written by Wending Fu, Jun.2024 in Beijing------------
% In MAVEN-J statistics, the curl_theta/curl_phi may have fair-sized
% deivation due to the low slides on R direction.
%
% The partial derivative is calculated by finite difference
% see also: curl, gradient, calculateJ_sph
%% grid diff
narginchk(6,6)
phiGrid = mean(diff(phi,1,1), 'all');
thetaGrid = mean(diff(theta,1,2), 'all');
rGrid = mean(diff(r,1,3), 'all');
%% medfilt (eliminate nan)
%这个过滤不会对非nan值带来影响
% % % coor = {'r','theta','phi'};
% % % c_eval('NanIdx_? = isnan(A?);', coor);
% % % c_eval('A?(NanIdx_?) = inf;', coor);
% % % 
% % % filt_box = [9,9,1]; % test-- filt box must cover nan region
% % % c_eval('filtA? = medfilt3(A?, filt_box);', coor);
% % % c_eval('A?(NanIdx_?) = filtA?(NanIdx_?);', coor);
%% calculate curl
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