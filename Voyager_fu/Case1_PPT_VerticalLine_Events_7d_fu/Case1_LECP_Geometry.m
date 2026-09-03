function geometry = Case1_LECP_Geometry(cfg)
%Case1_LECP_Geometry Documented nominal Voyager 1 LECP low-energy aperture.
%   Applies only to the P1 low-energy aperture used in Case 1.
%   The passband is read independently from original CDF metadata.
%   Columns are numbered S1...S8. Look vectors point OUT of the aperture;
%   particle velocity vectors point into it. They differ by a minus sign.
%
%   Sources:
%   https://voyager.ftecs.com/Handbook/DataFileDesc/SEDR/tabg1.html
%   https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/fk/vg1_v02.tf
%   https://pds-ppi.igpp.ucla.edu/data/VG1-S-LECP-4-SUMM-SECTOR-15MIN-V1.0/CATALOG/VG1_LECP_INST.CAT

%% user-approved nominal mounting, with no as-built offset calibration
if nargin < 1, cfg = Case1_Config; end
if ~cfg.NominalLECPGeometryApproved
    error('VoyagerCase1:GeometryApproval', ...
        'Official nominal LECP mounting has not been approved.');
end
% S/C clock=200 deg, cone=90 deg from SEDR Table G-1. The official FK
% gives azimuth_SC=-55 deg-clock, HGA=-Z_SC. This gives the signed axis:
axisSC = [-sind(15); cosd(15); 0];
hgaSC = [0; 0; -1];
tangentSC = cross(axisSC, hgaSC);

%% sector zero and handedness
% The S8/S1 boundary is along HGA. Increasing sector number is a positive
% right-handed rotation about axisSC. This is an angular numbering rule,
% not an assumption that the motor always scans in one temporal direction.
% Independent check: at Canopus lock S3 looks north of the ecliptic and
% retrograde, as stated by PDS; reversing the sign violates both checks.
theta = 22.5 + (0:7)*45;
lookSC = hgaSC*cosd(theta) + tangentSC*sind(theta);

%% return both direction conventions explicitly
geometry = struct;
geometry.Sector = 1:8;
geometry.CenterFromHGA_deg = theta;
geometry.AxisSC = axisSC;
geometry.HGASC = hgaSC;
geometry.LookSC = lookSC;
geometry.ParticleSC = -lookSC;
geometry.ActiveSectors = 1:7;
geometry.BlockedSector = 8;
geometry.Channel = 'P1 low-energy aperture; use original CDF energy bounds';
geometry.Model = 'Official nominal mounting; as-built offsets unavailable';
geometry.SourceURL = [ ...
    "https://voyager.ftecs.com/Handbook/DataFileDesc/SEDR/tabg1.html"; ...
    "https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/fk/vg1_v02.tf"; ...
    "https://pds-ppi.igpp.ucla.edu/data/VG1-S-LECP-4-SUMM-SECTOR-15MIN-V1.0/CATALOG/VG1_LECP_INST.CAT"];
end
