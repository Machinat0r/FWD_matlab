%% 2026-02-28 inner planets (Mercury-Mars) + Jupiter/Saturn/Uranus/Neptune + Parker spiral
% Heliocentric ecliptic plane (x-y). Low-precision elements (Paul Schlyter).
% Time interval is used only to fetch OMNI solar wind speed; positions use the start time.

clear; clc;

%% User settings
r0_AU    = 0.05;       % inner boundary of spiral (AU)
rMax_AU  = 35;         % spiral out to (AU), needs to exceed Neptune (~30 AU)
nOrbit   = 1200;       % orbit curve resolution

TT = '2026-02-28T00:00:00Z/2026-03-01T00:00:00Z';
dd = regexp(TT,'\d+','match');
Y = str2double(dd{1});  Mo = str2double(dd{2});  Da = str2double(dd{3});  UT_hours = str2double(dd{4});
tintOmni = irf.tint(TT);

%% Fetch solar wind speed (OMNI) for Parker spiral
try
    omniV = irf_get_data_omni(tintOmni, 'v');

    if isa(omniV, 'TSeries')
        v_omni = double(omniV.data);     % km/s
    else
        v_omni = omniV(:,2);             % km/s
    end

    v_omni(~isfinite(v_omni) | v_omni<=0) = NaN;
    Vsw_km_s = nanmedian(v_omni);

    if ~isfinite(Vsw_km_s)
        error('OMNI Vsw is NaN (possibly no data returned).');
    end

catch ME
    warning('Failed to fetch OMNI solar wind speed via irf_get_data_omni: %s\nFallback to 400 km/s.', ME.message);
    Vsw_km_s = 400;
end

%% Constants
AU_km = 149597870.7;                     % km
Omega = 2*pi/(25.38*86400);              % Sun sidereal rotation rate (rad/s)
Vsw_AU_s = Vsw_km_s / AU_km;             % AU/s

%% Day number (Schlyter convention)
d = day_number_schlyter(Y, Mo, Da, UT_hours);

%% Elements at epoch (with linear rates)
E_Mercury = elements_mercury(d);
E_Venus   = elements_venus(d);
E_Mars    = elements_mars(d);
E_Jupiter = elements_jupiter(d);
E_Saturn  = elements_saturn(d);
E_Uranus  = elements_uranus(d);
E_Neptune = elements_neptune(d);

% Earth: Schlyter "Sun" elements give geocentric Sun vector; Earth heliocentric is negation.
E_SunGeo  = elements_sun_for_geocentric(d); % used only to compute Earth position
% For Earth's orbit curve, use same a,e,w with N=i=0 (ecliptic plane).
E_EarthOrbit.N = 0.0; E_EarthOrbit.i = 0.0;
E_EarthOrbit.w = E_SunGeo.w; E_EarthOrbit.a = E_SunGeo.a; E_EarthOrbit.e = E_SunGeo.e;
E_EarthOrbit.M = E_SunGeo.M; % not needed for curve but kept for consistency

%% Positions (heliocentric)
[rMe, lonMe] = planet_xyz(E_Mercury);
[rVe, lonVe] = planet_xyz(E_Venus);

[rSunGeo, ~] = planet_xyz(E_SunGeo);  % "Sun" geocentric vector (AU)
rEa = -rSunGeo;                       % Earth heliocentric (AU)
lonEa = atan2(rEa(2), rEa(1));

[rMa, lonMa] = planet_xyz(E_Mars);
[rJu, lonJu] = planet_xyz(E_Jupiter);
[rSa, lonSa] = planet_xyz(E_Saturn);
[rUr, lonUr] = planet_xyz(E_Uranus);
[rNe, lonNe] = planet_xyz(E_Neptune);

%% Orbit curves (ellipses, not circles)
[orbMeX, orbMeY] = orbit_curve(E_Mercury, nOrbit);
[orbVeX, orbVeY] = orbit_curve(E_Venus,   nOrbit);
[orbEaX, orbEaY] = orbit_curve(E_EarthOrbit, nOrbit);
[orbMaX, orbMaY] = orbit_curve(E_Mars,    nOrbit);
[orbJuX, orbJuY] = orbit_curve(E_Jupiter, nOrbit);
[orbSaX, orbSaY] = orbit_curve(E_Saturn,  nOrbit);
[orbUrX, orbUrY] = orbit_curve(E_Uranus,  nOrbit);
[orbNeX, orbNeY] = orbit_curve(E_Neptune, nOrbit);

%% Parker spiral through Earth
rE_AU = hypot(rEa(1), rEa(2));
phiE  = lonEa; % rad

% phi(r) = phi0 + Omega*(r - r0)/Vsw ; choose phi0 so phi(rE)=phiE
phi0 = phiE - Omega*(rE_AU - r0_AU)/Vsw_AU_s;

rSp   = linspace(r0_AU, rMax_AU, 2500);
phiSp = phi0 + Omega*(rSp - r0_AU)/Vsw_AU_s;
xSp   = rSp .* cos(phiSp);
ySp   = rSp .* sin(phiSp);

%% Plot
figure('Color','w'); hold on;

% Orbits
plot(orbMeX, orbMeY, '--', 'DisplayName','Mercury orbit');
plot(orbVeX, orbVeY, '--', 'DisplayName','Venus orbit');
plot(orbEaX, orbEaY, '--', 'DisplayName','Earth orbit');
% plot(orbMaX, orbMaY, '--', 'DisplayName','Mars orbit');
plot(orbJuX, orbJuY, '--', 'DisplayName','Jupiter orbit');
plot(orbSaX, orbSaY, '--', 'DisplayName','Saturn orbit');
plot(orbUrX, orbUrY, '--', 'DisplayName','Uranus orbit');
plot(orbNeX, orbNeY, '--', 'DisplayName','Neptune orbit');

% Parker spiral (uncomment to show)
% plot(xSp, ySp, 'LineWidth', 1.3, 'DisplayName','Parker spiral (through Earth)');

% Sun + planets
scatter(0, 0, 90, '*', 'DisplayName','Sun');
scatter(rMe(1), rMe(2), 45, 'filled', 'DisplayName','Mercury');
scatter(rVe(1), rVe(2), 45, 'filled', 'DisplayName','Venus');
scatter(rEa(1), rEa(2), 45, 'filled', 'DisplayName','Earth');
% scatter(rMa(1), rMa(2), 45, 'filled', 'DisplayName','Mars');
scatter(rJu(1), rJu(2), 55, 'filled', 'DisplayName','Jupiter');
scatter(rSa(1), rSa(2), 55, 'filled', 'DisplayName','Saturn');
scatter(rUr(1), rUr(2), 55, 'filled', 'DisplayName','Uranus');
scatter(rNe(1), rNe(2), 55, 'filled', 'DisplayName','Neptune');

% Labels (place slightly outward)
label_point(rMe, lonMe, 'Mercury');
label_point(rVe, lonVe, 'Venus');
% label_point(rEa, lonEa, 'Earth');
% label_point(rMa, lonMa, 'Mars');
label_point(rJu, lonJu, 'Jupiter');
label_point(rSa, lonSa, 'Saturn');
label_point(rUr, lonUr, 'Uranus');
label_point(rNe, lonNe, 'Neptune');
% text(0.15, 0.15, 'Sun');

% Spiral label (uncomment if spiral is plotted)
% idxLab = find(rSp > 2.5, 1, 'first');
% text(xSp(idxLab)*1.03, ySp(idxLab)*1.03, sprintf('Parker spiral (V_{sw}=%g km/s)', Vsw_km_s));

axis equal; grid on;
xlabel('x (AU)'); ylabel('y (AU)');
title(sprintf('Heliocentric ecliptic positions + Parker spiral (schematic)\n%04d-%02d-%02d %02.0f:00 UTC', Y, Mo, Da, UT_hours));

lim = 35;
xlim([-lim, lim]); ylim([-lim, lim]);
legend('Location','northeastoutside');

%% Numeric summary
fprintf('%04d-%02d-%02d %02.0f:00 UTC (schematic)\n', Y, Mo, Da, UT_hours);
print_summary('Mercury', rMe, lonMe);
print_summary('Venus',   rVe, lonVe);
print_summary('Earth',   rEa, lonEa);
% print_summary('Mars',    rMa, lonMa);
print_summary('Jupiter', rJu, lonJu);
print_summary('Saturn',  rSa, lonSa);
print_summary('Uranus',  rUr, lonUr);
print_summary('Neptune', rNe, lonNe);

%% ===================== Local functions =====================

function d = day_number_schlyter(y, m, D, UT_hours)
d = 367*y ...
  - 7 * floor((y + floor((m+9)/12))/4) ...
  + floor(275*m/9) ...
  + D - 730530 ...
  + UT_hours/24;
end

function E = kepler_E(M, e)
E = M + e*sin(M)*(1 + e*cos(M));
for k = 1:60
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f/fp;
    E = E + dE;
    if abs(dE) < 1e-12, break; end
end
end

function [rVec, lon] = planet_xyz(elem)
N = deg2rad(mod(elem.N,360));
i = deg2rad(mod(elem.i,360));
w = deg2rad(mod(elem.w,360));
a = elem.a; e = elem.e;
M = deg2rad(mod(elem.M,360));

E = kepler_E(M, e);

xv = a*(cos(E) - e);
yv = a*(sqrt(1 - e^2)*sin(E));

v = atan2(yv, xv);
r = hypot(xv, yv);

x = r*(cos(N)*cos(v+w) - sin(N)*sin(v+w)*cos(i));
y = r*(sin(N)*cos(v+w) + cos(N)*sin(v+w)*cos(i));
z = r*(sin(v+w)*sin(i));

rVec = [x; y; z];
lon  = atan2(y, x);
end

function [x, y] = orbit_curve(elem, npts)
N = deg2rad(mod(elem.N,360));
i = deg2rad(mod(elem.i,360));
w = deg2rad(mod(elem.w,360));
a = elem.a; e = elem.e;

v = linspace(0, 2*pi, npts); % true anomaly
r = a*(1 - e^2) ./ (1 + e*cos(v));

x = r .* (cos(N).*cos(v+w) - sin(N).*sin(v+w).*cos(i));
y = r .* (sin(N).*cos(v+w) + cos(N).*sin(v+w).*cos(i));
end

function label_point(rVec, lon, name)
x = rVec(1); y = rVec(2);
scale = 1.06;
text(x*scale, y*scale, sprintf('%s\n', name), 'FontSize', 10);
end

function print_summary(name, rVec, lon)
r = hypot(rVec(1), rVec(2));
fprintf('%-8s: r=%.4f AU, lon=%.2f deg\n', name, r, wrap360(rad2deg(lon)));
end

function ang = wrap360(ang)
ang = mod(ang, 360);
if ang < 0
    ang = ang + 360;
end
end

%% Orbital elements (Paul Schlyter)
function elem = elements_mercury(d)
elem.N = 48.3313 + 3.24587E-5*d;
elem.i = 7.0047  + 5.00E-8*d;
elem.w = 29.1241 + 1.01444E-5*d;
elem.a = 0.387098;
elem.e = 0.205635 + 5.59E-10*d;
elem.M = 168.6562 + 4.0923344368*d;
end

function elem = elements_venus(d)
elem.N = 76.6799 + 2.46590E-5*d;
elem.i = 3.3946  + 2.75E-8*d;
elem.w = 54.8910 + 1.38374E-5*d;
elem.a = 0.723330;
elem.e = 0.006773 - 1.302E-9*d;
elem.M = 48.0052 + 1.6021302244*d;
end

function elem = elements_mars(d)
elem.N = 49.5574 + 2.11081E-5*d;
elem.i = 1.8497  - 1.78E-8*d;
elem.w = 286.5016 + 2.92961E-5*d;
elem.a = 1.523688;
elem.e = 0.093405 + 2.516E-9*d;
elem.M = 18.6021 + 0.5240207766*d;
end

function elem = elements_jupiter(d)
elem.N = 100.4542 + 2.76854E-5*d;
elem.i = 1.3030   - 1.557E-7*d;
elem.w = 273.8777 + 1.64505E-5*d;
elem.a = 5.20256;
elem.e = 0.048498 + 4.469E-9*d;
elem.M = 19.8950  + 0.0830853001*d;
end

function elem = elements_saturn(d)
elem.N = 113.6634 + 2.38980E-5*d;
elem.i = 2.4886   - 1.081E-7*d;
elem.w = 339.3939 + 2.97661E-5*d;
elem.a = 9.55475;
elem.e = 0.055546 - 9.499E-9*d;
elem.M = 316.9670 + 0.0334442282*d;
end

function elem = elements_uranus(d)
elem.N = 74.0005 + 1.3978E-5*d;
elem.i = 0.7733 + 1.9E-8*d;
elem.w = 96.6612 + 3.0565E-5*d;
elem.a = 19.18171 - 1.55E-8*d;
elem.e = 0.047318 + 7.45E-9*d;
elem.M = 142.5905 + 0.011725806*d;
end

function elem = elements_neptune(d)
elem.N = 131.7806 + 3.0173E-5*d;
elem.i = 1.7700 - 2.55E-7*d;
elem.w = 272.8461 - 6.027E-6*d;
elem.a = 30.05826 + 3.313E-8*d;
elem.e = 0.008606 + 2.15E-9*d;
elem.M = 260.2471 + 0.005995147*d;
end

function elem = elements_sun_for_geocentric(d)
elem.N = 0.0;
elem.i = 0.0;
elem.w = 282.9404 + 4.70935E-5*d;
elem.a = 1.000000;
elem.e = 0.016709 - 1.151E-9*d;
elem.M = 356.0470 + 0.9856002585*d;
end