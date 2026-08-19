function result = kh_fote_event(eventID,startUtc,endUtc,varargin)
%KH_FOTE_EVENT Four-spacecraft FOTE/FOTE-V analysis for one MMS KH event.
%
% The implementation follows the affine four-point expansion used by the
% FOTE code in FWD_matlab/FOTE.  FOTE-V uses exactly the same expansion with
% the vector field changed from B to Vi (or Ve).  All four spacecraft are
% placed on one gap-aware time grid.  No smoothing is applied by default.
%
% Example:
%   r = kh_fote_event('KH005','2015-10-01T18:01:24Z', ...
%       '2015-10-01T18:09:00Z');

% Name-value options:
%   DataRoot          local MMS database
%   OutputRoot        directory for PNG, CSV and MAT output
%   IrfuRoot          irfu-matlab root
%   VelocityField     'Vi' (default) or 'Ve'
%   SmoothSamples     7 (default); one gap-aware movmedian pass
%   QualityPercent    FOTE/FOTE-V residual threshold, default 40
%   MaxNullDistanceL  maximum normalized null distance shown, default 4
%   ReliableDistanceL near-tetrahedron threshold, default 2
%   FigureVisible     'off' (default) or 'on'


p = inputParser;
p.addRequired('eventID',@(x)ischar(x) || isstring(x));
p.addRequired('startUtc',@(x)ischar(x) || isstring(x));
p.addRequired('endUtc',@(x)ischar(x) || isstring(x));
p.addParameter('DataRoot','Z:\SPART-WORK\Data\MMS',@(x)ischar(x) || isstring(x));
p.addParameter('OutputRoot',fullfile(pwd,'outputs'),@(x)ischar(x) || isstring(x));
p.addParameter('IrfuRoot','C:\Users\Administrator\Documents\irfu-matlab-master',@(x)ischar(x) || isstring(x));
p.addParameter('VelocityField','Vi',@(x)any(strcmpi(string(x),["Vi","Ve"])));
p.addParameter('SmoothSamples',7,@(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('QualityPercent',40,@(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('MaxNullDistanceL',4,@(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ReliableDistanceL',2,@(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('FigureVisible','off',@(x)any(strcmpi(string(x),["on","off"])));
p.parse(eventID,startUtc,endUtc,varargin{:});
opt = p.Results;
opt.eventID = char(string(eventID));
opt.startUtc = char(string(startUtc));
opt.endUtc = char(string(endUtc));
opt.DataRoot = char(string(opt.DataRoot));
opt.OutputRoot = char(string(opt.OutputRoot));
opt.IrfuRoot = char(string(opt.IrfuRoot));
opt.VelocityField = char(string(opt.VelocityField));
opt.FigureVisible = char(string(opt.FigureVisible));
opt.SmoothSamples = round(opt.SmoothSamples);

setupIrfu(opt.IrfuRoot);
if ~exist(opt.OutputRoot,'dir'), mkdir(opt.OutputRoot); end
setupMmsDatabase(opt.DataRoot);
tint = irf.tint(sprintf('%s/%s',opt.startUtc,opt.endUtc));

fprintf('Loading %s (%s to %s), velocity=%s ...\n', ...
    opt.eventID,opt.startUtc,opt.endUtc,opt.VelocityField);
raw = loadMmsFourSpacecraft(tint,opt.VelocityField);
aligned = alignAtMomentCadence(raw);
fprintf('Common gap-safe samples: %d (median cadence %.3f s)\n', ...
    numel(aligned.t),median(diff(aligned.t),'omitnan'));

if opt.SmoothSamples > 1
    aligned.B = onePassMedian(aligned.B,opt.SmoothSamples,aligned.t);
    aligned.V = onePassMedian(aligned.V,opt.SmoothSamples,aligned.t);
    fprintf('Applied one movmedian pass, %d samples.\n',opt.SmoothSamples);
else
    fprintf('Smoothing disabled.\n');
end

mag = analyzeAffineField(aligned.R,aligned.B);
vel = analyzeAffineField(aligned.R,aligned.V);

fluxField = aligned.V;
for sc = 1:4
    fluxField(:,:,sc) = aligned.V(:,:,sc).*aligned.N(:,sc);
end
flux = analyzeAffineField(aligned.R,fluxField);

mag.eta = safeRatioPercent(abs(mag.divergence),mag.curlMagnitude);
mag.good = all(isfinite(mag.eigenvalues),2) & isfinite(mag.eta) & ...
    isfinite(mag.xi) & mag.eta <= opt.QualityPercent & ...
    mag.xi <= opt.QualityPercent;
mag.near = mag.good & mag.distanceL <= opt.ReliableDistanceL;

% For FOTE-V, div(V) is physical.  The four-point continuity residual
% alpha = |div(nV)|/|curl(nV)| is therefore used instead of forcing div(V)=0.
vel.alpha = safeRatioPercent(abs(flux.divergence),flux.curlMagnitude);
vel.good = all(isfinite(vel.eigenvalues),2) & isfinite(vel.alpha) & ...
    vel.alpha <= opt.QualityPercent;
vel.near = vel.good & vel.distanceL <= opt.ReliableDistanceL;

agreement = compareScrewSense(mag,vel);
agreement(~(mag.good & vel.good)) = NaN;

base = sprintf('%s_%s_%s_FOTE_FOTEV_%s',safeName(opt.eventID), ...
    compactUtc(opt.startUtc),compactUtc(opt.endUtc),upper(opt.VelocityField));
pdfPath = fullfile(opt.OutputRoot,[base '.pdf']);
csvPath = fullfile(opt.OutputRoot,[base '_timeseries.csv']);
matPath = fullfile(opt.OutputRoot,[base '.mat']);

fig = plotEvent(aligned,mag,vel,opt);
exportgraphics(fig,pdfPath,'ContentType','vector','BackgroundColor','white');
if strcmpi(opt.FigureVisible,'off'), close(fig); end

timeUtc = datetime(aligned.t,'ConvertFrom','posixtime','TimeZone','UTC');
outTable = table(timeUtc,mag.distanceKm,mag.distanceL,string(mag.type), ...
    mag.eta,mag.xi,mag.good,mag.near,vel.distanceKm,vel.distanceL, ...
    string(vel.type),vel.alpha,vel.good,vel.near,agreement, ...
    'VariableNames',{'TimeUTC','B_NullDistance_km','B_NullDistance_L', ...
    'B_Type','B_Eta_percent','B_Xi_percent','B_Good','B_Near', ...
    'V_NullDistance_km','V_NullDistance_L','V_Type','V_Alpha_percent', ...
    'V_Good','V_Near','ScrewAgreement'});
writetable(outTable,csvPath);
save(matPath,'aligned','mag','vel','flux','agreement','opt','-v7.3');

summary = buildSummary(mag,vel,agreement);
result = struct('EventID',opt.eventID,'StartUTC',opt.startUtc,'EndUTC',opt.endUtc, ...
    'VelocityField',opt.VelocityField,'Figure',pdfPath,'TimeSeries',csvPath, ...
    'MatFile',matPath,'Summary',summary);
fprintf('Preview figure: %s\n',pdfPath);
fprintf('B quality40=%d, V quality40=%d, screw same=%d, opposite=%d\n', ...
    summary.BQuality40,summary.VQuality40,summary.ScrewSame,summary.ScrewOpposite);
end


function setupMmsDatabase(dataRoot)
% Session-only registration avoids the persistent datastore write performed
% by mms.db_init.  This is equivalent for reads and avoids an R2026a
% home-session crash observed on this workstation.
global MMS_DB;
MMS_DB = mms_db;
localDb = mms_local_file_db([char(string(dataRoot)) filesep]);
MMS_DB.add_db(localDb);
end


function setupIrfu(root)
% Add the minimum paths explicitly.  This avoids IRFU startup network checks.
if exist('mms.get_data','file') == 2 && exist('irf.tint','file') == 2
    return
end
dirs = {root,fullfile(root,'irf'),fullfile(root,'plots'), ...
    fullfile(root,'plots','mms'),fullfile(root,'mission','cluster'), ...
    fullfile(root,'mission','mms'),fullfile(root,'contrib','nasa_cdf_patch')};
for k = 1:numel(dirs)
    if exist(dirs{k},'dir'), addpath(dirs{k}); end
end
leapFile = fullfile(root,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt');
if exist(leapFile,'file'), setenv('CDF_LEAPSECONDSTABLE',leapFile); end
if exist('mms.get_data','file') ~= 2 || exist('irf.tint','file') ~= 2
    error('IRFU/MMS functions were not found below %s.',root);
end
end


function raw = loadMmsFourSpacecraft(tint,velocityField)
raw = struct;
raw.B = cell(4,1); raw.V = cell(4,1); raw.N = cell(4,1);
for sc = 1:4
    bts = mms.get_data('B_gse_brst_l2',tint,sc);
    if strcmpi(velocityField,'Ve')
        vts = mms.get_data('Ve_gse_fpi_brst_l2',tint,sc);
        nts = mms.get_data('Ne_fpi_brst_l2',tint,sc);
    else
        vts = mms.get_data('Vi_gse_fpi_brst_l2',tint,sc);
        nts = mms.get_data('Ni_fpi_brst_l2',tint,sc);
    end
    raw.B{sc} = tseriesMatrix(bts,3,sprintf('MMS%d B',sc));
    raw.V{sc} = tseriesMatrix(vts,3,sprintf('MMS%d %s',sc,velocityField));
    raw.N{sc} = tseriesMatrix(nts,1,sprintf('MMS%d density',sc));
end

pos = mms.get_data('R_gse',tint);
tp = double(pos.time.epochUnix(:));
raw.R = cell(4,1);
for sc = 1:4
    field = sprintf('gseR%d',sc);
    values = double(pos.(field));
    raw.R{sc} = [tp,values(:,1:3)];
end
end


function out = tseriesMatrix(ts,nComponents,label)
if isempty(ts) || isempty(ts.time)
    error('%s returned no data.',label);
end
t = double(ts.time.epochUnix(:));
values = squeeze(double(ts.data));
if isvector(values), values = values(:); end
if size(values,1) ~= numel(t) && size(values,2) == numel(t), values = values.'; end
n = min(numel(t),size(values,1));
if n == 0 || size(values,2) < nComponents
    error('%s has an unexpected data shape.',label);
end
out = [t(1:n),values(1:n,1:nComponents)];
end


function aligned = alignAtMomentCadence(raw)
% MMS1 moment cadence is the comparison grid.  Linear interpolation is
% allowed only within each contiguous source segment; there is no extrapolation.
t = double(raw.V{1}(:,1));
t = unique(t(isfinite(t)),'stable');
n = numel(t);
if n < 20, error('The reference moment series contains fewer than 20 samples.'); end

B = nan(n,3,4); V = nan(n,3,4); R = nan(n,3,4); N = nan(n,4);
for sc = 1:4
    B(:,:,sc) = gapAwareResamp(raw.B{sc},t);
    V(:,:,sc) = gapAwareResamp(raw.V{sc},t);
    R(:,:,sc) = gapAwareResamp(raw.R{sc},t);
    N(:,sc) = gapAwareResamp(raw.N{sc},t);
end
valid = isfinite(t);
for sc = 1:4
    valid = valid & all(isfinite(B(:,:,sc)),2) & all(isfinite(V(:,:,sc)),2) & ...
        all(isfinite(R(:,:,sc)),2) & isfinite(N(:,sc));
end
if nnz(valid) < 20
    error('Fewer than 20 samples have common contiguous four-spacecraft support.');
end
aligned = struct('t',t(valid),'B',B(valid,:,:),'V',V(valid,:,:), ...
    'R',R(valid,:,:),'N',N(valid,:));
end


function values = gapAwareResamp(source,targetTime)
values = nan(numel(targetTime),size(source,2)-1);
source = source(isfinite(source(:,1)),:);
if isempty(source), return; end
[~,ia] = unique(source(:,1),'stable');
source = sortrows(source(ia,:),1);
if size(source,1) == 1
    hit = abs(targetTime-source(1,1)) <= 1e-6;
    values(hit,:) = repmat(source(1,2:end),nnz(hit),1);
    return
end
dt = diff(source(:,1));
baseDt = median(dt(isfinite(dt) & dt>0),'omitnan');
if ~isfinite(baseDt), return; end
gapLimit = max(10*baseDt,1);
starts = [1;find(dt>gapLimit)+1];
ends = [find(dt>gapLimit);size(source,1)];
for k = 1:numel(starts)
    rows = starts(k):ends(k);
    if numel(rows) < 2, continue; end
    targetRows = targetTime>=source(rows(1),1) & targetTime<=source(rows(end),1);
    if any(targetRows)
        values(targetRows,:) = interp1(source(rows,1),source(rows,2:end), ...
            targetTime(targetRows),'linear',NaN);
    end
end
end


function out = onePassMedian(in,window,t)
out = in;
dt = diff(double(t));
baseDt = median(dt(isfinite(dt) & dt>0),'omitnan');
if ~isfinite(baseDt), return; end
split = find(dt>max(10*baseDt,1))+1;
starts = [1;split];
stops = [split-1;size(in,1)];
for k = 1:numel(starts)
    rows = starts(k):stops(k);
    for sc = 1:size(in,3)
        out(rows,:,sc) = movmedian(in(rows,:,sc),window,1,'omitnan','Endpoints','shrink');
    end
end
end


function a = analyzeAffineField(R,F)
n = size(F,1);
a = struct;
a.gradient = nan(3,3,n);
a.eigenvalues = complex(nan(n,3));
a.divergence = nan(n,1);
a.curlMagnitude = nan(n,1);
a.xi = nan(n,1);
a.distanceKm = nan(n,1);
a.distanceL = nan(n,1);
a.scaleL = nan(n,1);
a.gradientCondition = nan(n,1);
a.geometryCondition = nan(n,1);
a.type = strings(n,1);
a.screwSense = nan(n,1); % -1 screw-in (As), +1 screw-out (Bs)
a.finite = false(n,1);

for i = 1:n
    r = nan(4,3); f = nan(4,3);
    for sc = 1:4
        r(sc,:) = R(i,:,sc);
        f(sc,:) = F(i,:,sc);
    end
    if ~all(isfinite([r(:);f(:)])), continue; end
    r0 = mean(r,1); f0 = mean(f,1);
    x = r-r0; y = f-f0;
    if rank(x) < 3, continue; end
    G = (x\y).'; % G(i,j) = dF_i/dx_j, matching FOTE reshape convention
    if ~all(isfinite(G(:))), continue; end
    a.gradient(:,:,i) = G;
    a.geometryCondition(i) = cond(x);
    a.gradientCondition(i) = cond(G);
    ev = eig(G).';
    a.eigenvalues(i,:) = ev;
    a.divergence(i) = trace(G);
    curlVector = [G(3,2)-G(2,3),G(1,3)-G(3,1),G(2,1)-G(1,2)];
    a.curlMagnitude(i) = norm(curlVector);
    maxEig = max(abs(ev));
    if maxEig > 0, a.xi(i) = 100*abs(sum(ev))/maxEig; end

    pairs = [1 2;1 3;1 4;2 3;2 4;3 4];
    sep = nan(6,1);
    for j = 1:6, sep(j) = norm(r(pairs(j,1),:)-r(pairs(j,2),:)); end
    L = mean(sep,'omitnan');
    a.scaleL(i) = L;
    if rcond(G) > 1e-12
        rNull = r0-(G\f0.').';
        d = sqrt(sum((r-rNull).^2,2));
        a.distanceKm(i) = min(d);
        a.distanceL(i) = a.distanceKm(i)/L;
    end
    [a.type(i),a.screwSense(i)] = classifyGradient(ev);
    a.finite(i) = isfinite(a.distanceL(i)) && isfinite(a.gradientCondition(i));
end
end


function [type,sense] = classifyGradient(ev)
type = "degenerate"; sense = NaN;
scale = max(abs(ev));
if ~isfinite(scale) || scale == 0, type = "X"; return; end
imaginary = max(abs(imag(ev))) > 1e-6*scale;
realPart = real(ev);
zeroTol = 1e-8*scale;
nPos = nnz(realPart>zeroTol);
nNeg = nnz(realPart<-zeroTol);
if imaginary
    if nPos == 1 && nNeg == 2
        type = "As"; sense = -1;
    elseif nPos == 2 && nNeg == 1
        type = "Bs"; sense = 1;
    elseif nPos == 3
        type = "S+"; sense = 2;
    elseif nNeg == 3
        type = "S-"; sense = -2;
    else
        type = "O";
    end
else
    if min(abs(realPart))/scale < 1e-6
        type = "X";
    elseif nPos == 1 && nNeg == 2
        type = "A";
    elseif nPos == 2 && nNeg == 1
        type = "B";
    elseif nPos == 3
        type = "S+";
    elseif nNeg == 3
        type = "S-";
    end
end
end


function pct = safeRatioPercent(numerator,denominator)
pct = nan(size(numerator));
valid = isfinite(numerator) & isfinite(denominator) & denominator>0;
pct(valid) = 100*numerator(valid)./denominator(valid);
end


function agreement = compareScrewSense(mag,vel)
agreement = nan(size(mag.screwSense));
spiral = abs(mag.screwSense)==1 & abs(vel.screwSense)==1;
agreement(spiral & mag.screwSense==vel.screwSense) = 1;
agreement(spiral & mag.screwSense~=vel.screwSense) = -1;
end


function fig = plotEvent(d,mag,vel,opt)
t0 = d.t(1);
x = d.t-t0;
bMean = mean(d.B,3,'omitnan');
vMean = mean(d.V,3,'omitnan');
fig = figure('Visible',opt.FigureVisible,'Color','w','Units','pixels', ...
    'Position',[30 30 1540 1260]);
tl = tiledlayout(fig,5,1,'TileSpacing','compact','Padding','compact');
ax = gobjects(5,1);

ax(1) = nexttile(tl); hold(ax(1),'on');
bMag = sqrt(sum(bMean.^2,2));
plot(ax(1),x,bMean(:,1),'Color',[0.1 0.35 0.85],'LineWidth',0.7,'DisplayName','B_x');
plot(ax(1),x,bMean(:,2),'Color',[0.1 0.60 0.25],'LineWidth',0.7,'DisplayName','B_y');
plot(ax(1),x,bMean(:,3),'Color',[0.85 0.18 0.18],'LineWidth',0.7,'DisplayName','B_z');
plot(ax(1),x,bMag,'k','LineWidth',0.9,'DisplayName','|B|');
ylabel(ax(1),'B_{GSE} [nT]'); legend(ax(1),'Location','eastoutside','NumColumns',4);

ax(2) = nexttile(tl); hold(ax(2),'on');
plot(ax(2),x,vMean(:,1),'Color',[0.1 0.35 0.85],'LineWidth',0.7,'DisplayName','V_x');
plot(ax(2),x,vMean(:,2),'Color',[0.1 0.60 0.25],'LineWidth',0.7,'DisplayName','V_y');
plot(ax(2),x,vMean(:,3),'Color',[0.85 0.18 0.18],'LineWidth',0.7,'DisplayName','V_z');
ylabel(ax(2),sprintf('%s_{GSE} [km/s]',opt.VelocityField));
legend(ax(2),'Location','eastoutside','NumColumns',3);

ax(3) = nexttile(tl); hold(ax(3),'on');
maxError = 100;
bError = max([mag.eta mag.xi],[],2);
plot(ax(3),x,min(mag.eta,maxError),'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.55,'DisplayName','\eta_B');
plot(ax(3),x,min(mag.xi,maxError),'Color',[0.48 0.29 0.65], ...
    'LineWidth',0.55,'DisplayName','\xi_B');
plotTypeMarkers(ax(3),x,bError,mag.type,mag.good,maxError);
yline(ax(3),opt.QualityPercent,'--','40%','Color',[0 0 0], ...
    'HandleVisibility','off');
ylabel(ax(3),'FOTE B error [%]');
ylim(ax(3),[0 maxError]);

ax(4) = nexttile(tl); hold(ax(4),'on');
plot(ax(4),x,min(vel.alpha,maxError),'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.55,'DisplayName','\alpha_V');
plotTypeMarkers(ax(4),x,vel.alpha,vel.type,vel.good,maxError);
yline(ax(4),opt.QualityPercent,'--','40%','Color',[0 0 0], ...
    'HandleVisibility','off');
ylabel(ax(4),sprintf('FOTE-V %s error [%%]',opt.VelocityField));
ylim(ax(4),[0 maxError]);

ax(5) = nexttile(tl); hold(ax(5),'on');
bOut = mag.good & mag.type=="Bs";
bIn = mag.good & mag.type=="As";
vOut = vel.good & vel.type=="Bs";
vIn = vel.good & vel.type=="As";
scatter(ax(5),x(bOut),3*ones(nnz(bOut),1),28,[0.10 0.30 0.85],'>','filled');
scatter(ax(5),x(bIn),2*ones(nnz(bIn),1),28,[0.85 0.15 0.15],'^','filled');
scatter(ax(5),x(vOut),ones(nnz(vOut),1),28,[0.10 0.30 0.85],'>','filled');
scatter(ax(5),x(vIn),zeros(nnz(vIn),1),28,[0.85 0.15 0.15],'^','filled');
yticks(ax(5),0:3);
yticklabels(ax(5),{'V_i screw-in','V_i screw-out','B screw-in','B screw-out'});
ylim(ax(5),[-0.6 3.6]); ylabel(ax(5),'error < 40%');

for k = 1:numel(ax)
    grid(ax(k),'on'); box(ax(k),'on'); set(ax(k),'FontSize',9,'Layer','top');
    xlim(ax(k),[0 max(x)]);
end
set(ax(1:4),'XTickLabel',[]);
formatTimeAxis(ax(5),t0,max(x));
linkaxes(ax,'x');

sgtitle(tl,opt.eventID,'FontWeight','bold','FontSize',14,'Interpreter','none');
end


function plotTypeMarkers(ax,x,y,type,good,maxY)
spec = {
    "A",  '^',[0.85 0.15 0.15],false,'A';
    "B",  '>',[0.10 0.30 0.85],false,'B';
    "As", '^',[0.85 0.15 0.15],true, 'A_s screw-in';
    "Bs", '>',[0.10 0.30 0.85],true, 'B_s screw-out';
    "X",  'x',[0.15 0.15 0.15],false,'X';
    "O",  'o',[0.15 0.15 0.15],false,'O';
    "S+", 's',[0.75 0.15 0.70],true, 'source';
    "S-", 's',[0.00 0.60 0.65],true, 'sink'};
for k = 1:size(spec,1)
    idx = good & type==spec{k,1} & isfinite(y) & y<=maxY;
    if ~any(idx), continue; end
    h = scatter(ax,x(idx),y(idx),9,spec{k,3},spec{k,2}, ...
        'LineWidth',0.55,'MarkerFaceAlpha',0.65,'MarkerEdgeAlpha',0.65, ...
        'DisplayName',spec{k,5});
    if spec{k,4}, h.MarkerFaceColor = spec{k,3}; end
end
if any(good & isfinite(y) & y<=maxY)
    legend(ax,'show','Location','eastoutside','NumColumns',1,'FontSize',8);
end
end


function formatTimeAxis(ax,t0,durationSeconds)
if durationSeconds <= 0, return; end
ticks = linspace(0,durationSeconds,7);
xticks(ax,ticks);
dt = datetime(t0+ticks,'ConvertFrom','posixtime','TimeZone','UTC');
xticklabels(ax,cellstr(string(dt,'HH:mm:ss')));
xlabel(ax,sprintf('UTC on %s',char(string(datetime(t0,'ConvertFrom','posixtime','TimeZone','UTC'),'yyyy-MM-dd'))));
end


function s = buildSummary(mag,vel,agreement)
s = struct;
s.BQuality40 = nnz(mag.good); s.VQuality40 = nnz(vel.good);
s.BAsQuality40 = nnz(mag.good & mag.type=="As");
s.BBsQuality40 = nnz(mag.good & mag.type=="Bs");
s.VAsQuality40 = nnz(vel.good & vel.type=="As");
s.VBsQuality40 = nnz(vel.good & vel.type=="Bs");
s.ScrewSame = nnz(agreement==1);
s.ScrewOpposite = nnz(agreement==-1);
end


function out = safeName(in)
out = regexprep(char(string(in)),'[^A-Za-z0-9_-]','_');
end


function out = compactUtc(in)
out = regexprep(char(string(in)),'[^0-9]','');
if numel(out) > 14, out = out(1:14); end
end
