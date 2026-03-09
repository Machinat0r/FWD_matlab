% ---- user inputs ----
ic   = 1;
tint = irf.tint('2018-07-11T14:00:00.00Z/2018-07-11T16:00:00.00Z');

% ---- load PDist omni flux ----
P = mms.get_data('Omnifluxhplus_hpca_srvy_l2', tint, ic);  % PDist

% ---- extract energy axis ----
E = ion_energy.data;      % usually 63 energies (eV)
if isvector(E)
  energy = E(:);             % [63x1]
else
  % in case energy is time-dependent (Nt x 63), take first record
  energy = E(1,:).';
end

% ---- extract flux ----
J = double(P.data);          % [Nt x 63]

% clean fill/invalid values (HPCA fill often very negative; also avoid log of <=0)
J(J < -1e30 | J <= 0) = NaN;

% ---- make sure dimensions match (Nt x Nenergy) ----
if size(J,2) ~= numel(energy) && size(J,1) == numel(energy)
  J = J.';                   % transpose if needed
end

% ---- build specrec for irf_spectrogram ----
spec = struct();
spec.t = P.time.epochUnix;   % [Nt x 1] unix seconds
spec.f = double(energy);             % [63 x 1] eV
spec.p = J;           % color in log10

% spec.p_label = 'log_{10} J_{H^+} (cm^{-2} s^{-1} sr^{-1} eV^{-1})';

% ---- plot ----
h = irf_plot(1,'newfigure');
irf_spectrogram(h, spec);
set(h,'yscale','log');
set(h,'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel(h,'E (keV)');
irf_timeaxis(h,'utc');
title(h, sprintf('MMS%d HPCA H^+ Omni Flux (elev 0–360)', ic));
colormap('jet')

caxis([0 2])
