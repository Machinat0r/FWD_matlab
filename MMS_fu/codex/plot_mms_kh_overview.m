function result = plot_mms_kh_overview(eventRow,windowInfo,sc,outputRoot)
%PLOT_MMS_KH_OVERVIEW Nine-panel MMS KH overview derived from Overview_download.m.
%
% Panels: B_GSM, E_GSM, ne/ni, Vi_GSM, Ve_GSM, particle-moment J,
% J.E and J.(E+Ve x B), electron omni flux, ion omni flux.
%
% The original script is left untouched.  This parameterized derivative
% fixes its overwritten E+Ve x B calculation and uses Ni for the ion
% current rather than imposing Ne=Ni.

if nargin < 4 || isempty(outputRoot)
    outputRoot = 'C:\Users\Administrator\Documents\KH';
end
validateattributes(sc,{'numeric'},{'scalar','integer','>=',1,'<=',4});

plotStart = windowInfo.PlotStartUTC;
plotEnd = windowInfo.PlotEndUTC;
tt = sprintf('%s/%s',isoText(plotStart),isoText(plotEnd));
tint = irf.tint(tt);

warnings = strings(0,1);
B = []; E = []; Ne = []; Ni = []; Vi = []; Ve = []; J = []; energy = [];
eSpec = []; iSpec = [];

try
    Bts = mms.get_data('B_gsm_brst',tint,sc);
    B = irf.ts2mat(Bts);
    if isempty(B), warnings(end+1) = "B: no data returned"; end
catch ME
    warnings(end+1) = "B: " + string(ME.message);
end
try
    Ets = mms.get_data('E_gse_edp_brst_l2',tint,sc);
    E = irf.ts2mat(irf_gse2gsm(Ets));
    if isempty(E), warnings(end+1) = "E: no data returned"; end
catch ME
    warnings(end+1) = "E: " + string(ME.message);
end
try
    Nets = mms.get_data('Ne_fpi_brst_l2',tint,sc);
    Ne = irf.ts2mat(Nets);
    if isempty(Ne), warnings(end+1) = "Ne: no data returned"; end
catch ME
    warnings(end+1) = "Ne: " + string(ME.message);
end
try
    Nits = mms.get_data('Ni_fpi_brst_l2',tint,sc);
    Ni = irf.ts2mat(Nits);
    if isempty(Ni), warnings(end+1) = "Ni: no data returned"; end
catch ME
    warnings(end+1) = "Ni: " + string(ME.message);
end
try
    Vits = mms.get_data('Vi_gse_fpi_brst_l2',tint,sc);
    Vi = irf.ts2mat(irf_gse2gsm(Vits));
    if isempty(Vi), warnings(end+1) = "Vi: no data returned"; end
catch ME
    warnings(end+1) = "Vi: " + string(ME.message);
end
try
    Vets = mms.get_data('Ve_gse_fpi_brst_l2',tint,sc);
    Ve = irf.ts2mat(irf_gse2gsm(Vets));
    if isempty(Ve), warnings(end+1) = "Ve: no data returned"; end
catch ME
    warnings(end+1) = "Ve: " + string(ME.message);
end

% Particle-moment current and two energy-conversion measures.  Do not use
% irf_resamp here: it extrapolates and can interpolate across burst gaps.
% All derived quantities are evaluated only where every source has support
% inside the same contiguous data segment.
if all(~cellfun(@isempty,{Ne,Ni,Vi,Ve}))
    try
        NiR = gapAwareResamp(Ni,Ne);
        ViR = gapAwareResamp(Vi,Ne);
        VeR = gapAwareResamp(Ve,Ne);
        qe = 1.602176634e-19;
        % cm^-3 * km/s * e, converted to nA/m^2 => factor e*1e18.
        jValues = qe*1e18*(NiR(:,2).*ViR(:,2:4)-Ne(:,2).*VeR(:,2:4));
        validJ = all(isfinite([Ne(:,2),NiR(:,2),ViR(:,2:4),VeR(:,2:4)]),2);
        jValues(~validJ,:) = NaN;
        if any(validJ)
            J = [Ne(:,1),jValues];
        else
            J = [];
            warnings(end+1) = "J: no common contiguous Ne/Ni/Vi/Ve support";
        end
    catch ME
        warnings(end+1) = "J: " + string(ME.message);
        J = [];
    end
end

if ~isempty(J) && ~isempty(E) && ~isempty(B)
    try
            ER = gapAwareResamp(E,J);
            BR = gapAwareResamp(B,J);
            VeJ = gapAwareResamp(Ve,J);
            Eprime = ER(:,2:4) + 1e-3*cross(VeJ(:,2:4),BR(:,2:4));
            % nA/m^2 * mV/m = pW/m^3.
            energyValues = [sum(J(:,2:4).*ER(:,2:4),2), ...
                sum(J(:,2:4).*Eprime,2)];
            validEnergy = all(isfinite([J(:,2:4),ER(:,2:4),BR(:,2:4),VeJ(:,2:4)]),2);
            energyValues(~validEnergy,:) = NaN;
            if any(validEnergy)
                energy = [J(:,1),energyValues];
            else
                energy = [];
                warnings(end+1) = "energy: no common contiguous J/E/B/Ve support";
            end
    catch ME
        warnings(end+1) = "energy: " + string(ME.message);
        energy = [];
    end
end

try
    v = mms.db_get_variable(sprintf('mms%d_fpi_brst_l2_des-moms',sc), ...
        sprintf('mms%d_des_energyspectr_omni_brst',sc),tint);
    eSpec = variableToSpec(v,'electron omni energy flux');
catch ME
    warnings(end+1) = "electron flux: " + string(ME.message);
end
try
    v = mms.db_get_variable(sprintf('mms%d_fpi_brst_l2_dis-moms',sc), ...
        sprintf('mms%d_dis_energyspectr_omni_brst',sc),tint);
    iSpec = variableToSpec(v,'ion omni energy flux');
catch ME
    warnings(end+1) = "ion flux: " + string(ME.message);
end

nPanels = 9;
tRef = epochUnix(plotStart);
fig = figure('Visible','off','Color','w','Units','pixels','Position',[20 20 1500 1800]);
set(fig,'UserData',struct('t_start_epoch',tRef));
cleanupFig = onCleanup(@() closeValid(fig));
ax = gobjects(nPanels,1);

% 1. Magnetic field.
ax(1) = irf_subplot(nPanels,1,-1);
if ~isempty(B)
    Bp = breakTimeGaps(reduceForPlot(B,6000));
    Bmag = sqrt(sum(Bp(:,2:4).^2,2));
    tPlot = Bp(:,1)-tRef;
    plot(ax(1),tPlot,Bp(:,2),'b',tPlot,Bp(:,3),'g',tPlot,Bp(:,4),'r',tPlot,Bmag,'k','LineWidth',0.7);
    legend(ax(1),{'B_x','B_y','B_z','|B|'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(1),'B unavailable'); end
ylabel(ax(1),'B [nT]');

% 2. Electric field.
ax(2) = irf_subplot(nPanels,1,-2);
if ~isempty(E)
    Ep = reduceForPlot(E,6000);
    Ep = breakTimeGaps(Ep);
    tPlot = Ep(:,1)-tRef;
    plot(ax(2),tPlot,Ep(:,2),'b',tPlot,Ep(:,3),'g',tPlot,Ep(:,4),'r','LineWidth',0.6);
    legend(ax(2),{'E_x','E_y','E_z'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(2),'E unavailable'); end
ylabel(ax(2),'E [mV/m]');

% 3. Densities.
ax(3) = irf_subplot(nPanels,1,-3);
hold(ax(3),'on');
densLegend = strings(0,1);
if ~isempty(Ne)
    Nep = breakTimeGaps(reduceForPlot(Ne,6000));
    plot(ax(3),Nep(:,1)-tRef,Nep(:,2),'b','LineWidth',0.8);
    densLegend(end+1) = "n_e";
end
if ~isempty(Ni)
    Nip = breakTimeGaps(reduceForPlot(Ni,6000));
    plot(ax(3),Nip(:,1)-tRef,Nip(:,2),'g','LineWidth',0.8);
    densLegend(end+1) = "n_i";
end
hold(ax(3),'off');
if isempty(Ne) && isempty(Ni), blankPanel(ax(3),'density unavailable');
else, legend(ax(3),cellstr(densLegend),'Location','eastoutside','FontSize',7); end
ylabel(ax(3),'n [cm^{-3}]');

% 4. Ion velocity.
ax(4) = irf_subplot(nPanels,1,-4);
if ~isempty(Vi)
    Vip = breakTimeGaps(reduceForPlot(Vi,6000));
    tPlot = Vip(:,1)-tRef;
    plot(ax(4),tPlot,Vip(:,2),'b',tPlot,Vip(:,3),'g',tPlot,Vip(:,4),'r','LineWidth',0.7);
    legend(ax(4),{'V_{ix}','V_{iy}','V_{iz}'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(4),'ion velocity unavailable'); end
ylabel(ax(4),'V_i [km/s]');

% 5. Electron velocity.
ax(5) = irf_subplot(nPanels,1,-5);
if ~isempty(Ve)
    Vep = breakTimeGaps(reduceForPlot(Ve,6000));
    tPlot = Vep(:,1)-tRef;
    plot(ax(5),tPlot,Vep(:,2),'b',tPlot,Vep(:,3),'g',tPlot,Vep(:,4),'r','LineWidth',0.6);
    legend(ax(5),{'V_{ex}','V_{ey}','V_{ez}'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(5),'electron velocity unavailable'); end
ylabel(ax(5),'V_e [km/s]');

% 6. Current density.
ax(6) = irf_subplot(nPanels,1,-6);
if ~isempty(J)
    Jp = breakTimeGaps(reduceForPlot(J,6000));
    tPlot = Jp(:,1)-tRef;
    plot(ax(6),tPlot,Jp(:,2),'b',tPlot,Jp(:,3),'g',tPlot,Jp(:,4),'r','LineWidth',0.6);
    legend(ax(6),{'J_x','J_y','J_z'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(6),'current unavailable'); end
ylabel(ax(6),'J [nA/m^2]');

% 7. Energy conversion.
ax(7) = irf_subplot(nPanels,1,-7);
if ~isempty(energy)
    energyP = breakTimeGaps(reduceForPlot(energy,6000));
    tPlot = energyP(:,1)-tRef;
    plot(ax(7),tPlot,energyP(:,2),'k',tPlot,energyP(:,3),'m','LineWidth',0.7);
    yline(ax(7),0,'Color',[0.4 0.4 0.4],'LineStyle','--');
    legend(ax(7),{'J\cdotE','J\cdot(E+V_e\timesB)'},'Location','eastoutside','FontSize',7);
else, blankPanel(ax(7),'energy conversion unavailable'); end
ylabel(ax(7),'pW/m^3');

% 8-9. Omni-directional FPI energy flux.
ax(8) = irf_subplot(nPanels,1,-8);
if ~isempty(eSpec)
    irf_spectrogram(ax(8),breakSpecTimeGaps(reduceSpecForPlot(eSpec,3000)));
    set(ax(8),'YScale','log'); colormap(ax(8),turbo);
else, blankPanel(ax(8),'electron flux unavailable'); end
ylabel(ax(8),{'E_e','[eV]'});

ax(9) = irf_subplot(nPanels,1,-9);
if ~isempty(iSpec)
    irf_spectrogram(ax(9),breakSpecTimeGaps(reduceSpecForPlot(iSpec,3000)));
    set(ax(9),'YScale','log'); colormap(ax(9),turbo);
else, blankPanel(ax(9),'ion flux unavailable'); end
ylabel(ax(9),{'E_i','[eV]'});

for k = 1:nPanels
    grid(ax(k),'on');
    set(ax(k),'FontSize',8,'Box','on');
end
try
    % Use the current IRFU syntax.  irf_subplot leaves the Y limits at
    % [0 1] until an explicit data-aware Y zoom is requested.
    irf_zoom(ax,'x',tint);
    irf_zoom(ax(1:7),'y');
    irf_plot_axis_align(ax);
catch
    linkaxes(ax,'x');
    xBounds = [epochUnix(plotStart),epochUnix(plotEnd)]-tRef;
    for k = 1:nPanels, xlim(ax(k),xBounds); end
    for k = 1:7, ylim(ax(k),'auto'); end
end
set(ax(1:end-1),'XTickLabel',[]);

titleLines = {
    sprintf('%s / %s / MMS%d / %s / %s',eventRow.EventID,eventRow.CatalogID,sc,eventRow.Tier,eventRow.Confidence)
    sprintf('event %s to %s | plotted burst %s to %s | %s', ...
        isoText(eventRow.StartUTC),isoText(eventRow.EndUTC),isoText(plotStart),isoText(plotEnd),windowInfo.Mode)
    sprintf('Sources: %s',eventRow.Source)
    sprintf('DOI: %s',eventRow.DOI)};
sgtitle(fig,titleLines,'Interpreter','none','FontSize',9,'FontWeight','bold');

eventDir = fullfile(outputRoot,'figures',char(eventRow.EventID + "_" + string(eventRow.StartUTC,'yyyyMMdd_HHmmss')));
if ~exist(eventDir,'dir'), mkdir(eventDir); end
base = sprintf('%s_%s_MMS%d',eventRow.EventID,char(string(plotStart,'yyyyMMdd_HHmmss')),sc);
pngPath = fullfile(eventDir,[base '.png']);
exportgraphics(fig,pngPath,'Resolution',200,'BackgroundColor','white');

availablePanels = [hasFiniteRows(B,4),hasFiniteRows(E,4), ...
    hasFiniteRows(Ne,2) || hasFiniteRows(Ni,2),hasFiniteRows(Vi,4), ...
    hasFiniteRows(Ve,4),hasFiniteRows(J,4),hasFiniteRows(energy,3), ...
    hasFiniteSpec(eSpec),hasFiniteSpec(iSpec)];
if all(availablePanels)
    figureStatus = "ok";
elseif any(availablePanels)
    figureStatus = "partial";
else
    figureStatus = "no_data";
end
result = table(string(eventRow.EventID),sc,plotStart,plotEnd,string(pngPath), ...
    figureStatus,strjoin(warnings,' | '), ...
    'VariableNames',{'EventID','Spacecraft','PlotStartUTC','PlotEndUTC','OutputFile','Status','Warnings'});
clear cleanupFig
closeValid(fig);
end

function spec = variableToSpec(v,labelText)
spec = struct;
spec.t = irf_time(v.DEPEND_0.data,'ttns>epoch');
energy = double(v.DEPEND_1.data);
if ~isvector(energy), energy = energy(1,:); end
spec.f = energy(:)';
data = double(v.data);
data = squeeze(data);
if size(data,1) ~= numel(spec.t) && size(data,2) == numel(spec.t), data = data'; end
spec.p = data;
spec.f_label = '';
spec.p_label = {' ',labelText};
end

function out = reduceForPlot(in,maxPoints)
if size(in,1) <= maxPoints, out = in; return; end
idx = unique(round(linspace(1,size(in,1),maxPoints)));
out = in(idx,:);
end

function out = gapAwareResamp(source,target)
% Linear interpolation without extrapolation and without crossing gaps.
out = [target(:,1),nan(size(target,1),size(source,2)-1)];
if isempty(source) || isempty(target), return; end
source = source(isfinite(source(:,1)),:);
if isempty(source), return; end
[~,ia] = unique(source(:,1),'stable');
source = source(ia,:);
source = sortrows(source,1);
if size(source,1) == 1
    match = abs(target(:,1)-source(1,1)) <= 1e-6;
    out(match,2:end) = repmat(source(1,2:end),sum(match),1);
    return
end
dt = diff(double(source(:,1)));
baseDt = median(dt(isfinite(dt) & dt > 0),'omitnan');
if isempty(baseDt) || ~isfinite(baseDt), return; end
gapLimit = max(10*baseDt,1);
segmentStart = [1;find(dt > gapLimit)+1];
segmentEnd = [find(dt > gapLimit);size(source,1)];
for k = 1:numel(segmentStart)
    rows = segmentStart(k):segmentEnd(k);
    if numel(rows) < 2, continue; end
    targetRows = target(:,1) >= source(rows(1),1) & target(:,1) <= source(rows(end),1);
    if any(targetRows)
        out(targetRows,2:end) = interp1(source(rows,1),source(rows,2:end), ...
            target(targetRows,1),'linear',NaN);
    end
end
end

function out = reduceSpecForPlot(in,maxTimes)
out = in;
if numel(in.t) <= maxTimes, return; end
idx = unique(round(linspace(1,numel(in.t),maxTimes)));
out.t = in.t(idx);
out.p = in.p(idx,:);
if ~isvector(in.f) && size(in.f,1) == numel(in.t)
    out.f = in.f(idx,:);
end
end

function out = breakSpecTimeGaps(in)
% Blank both samples bordering a discontinuity so a spectrogram does not
% paint a single stretched colour cell across missing burst files.
out = in;
if isempty(in) || numel(in.t) < 3, return; end
t = double(in.t(:));
dt = diff(t);
baseDt = median(dt(isfinite(dt) & dt > 0),'omitnan');
if isempty(baseDt) || ~isfinite(baseDt), return; end
gapRows = find(dt > max(10*baseDt,1));
blankRows = unique([gapRows;gapRows+1]);
out.p(blankRows,:) = NaN;
end

function out = breakTimeGaps(in)
% Insert NaNs at discontinuities so direct MATLAB plot calls do not draw
% artificial straight lines across missing burst files.
out = in;
if isempty(in) || size(in,1) < 3, return; end
dt = diff(double(in(:,1)));
baseDt = median(dt(isfinite(dt) & dt > 0),'omitnan');
if isempty(baseDt) || ~isfinite(baseDt), return; end
gapRows = find(dt > max(10*baseDt,1)) + 1;
out(gapRows,2:end) = NaN;
end

function tf = hasFiniteRows(in,minColumns)
tf = ~isempty(in) && size(in,2) >= minColumns && ...
    any(all(isfinite(in(:,2:minColumns)),2));
end

function tf = hasFiniteSpec(spec)
tf = ~isempty(spec) && isfield(spec,'p') && any(isfinite(spec.p(:)));
end

function blankPanel(ax,msg)
cla(ax);
text(ax,0.5,0.5,msg,'Units','normalized','HorizontalAlignment','center', ...
    'Color',[0.55 0.1 0.1],'FontWeight','bold','Interpreter','none');
end

function s = isoText(t)
s = char(string(t,'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function e = epochUnix(t)
e = posixtime(t);
end

function closeValid(fig)
if ~isempty(fig) && isgraphics(fig), close(fig); end
end
