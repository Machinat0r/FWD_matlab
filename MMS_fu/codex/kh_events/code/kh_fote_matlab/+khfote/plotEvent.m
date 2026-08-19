function fig = plotEvent(raw,d,mag,vel,opt,tStart,tStop)
%KH FOTE 绘图模块：生成单事件7-panel统计图。
%
% Panel定义
%   1. 四星平均原始磁场；2. 四星平均原始速度；
%   3. 四星平均时间平均磁场；4. 四星平均时间平均速度；
%   5. 磁场FOTE eta/xi与零点类型；6. FOTE-V alpha与零点类型；
%   7. 误差、特征值稳定度及可选持续时间均合格的screw方向。
%
% 图中不包含零点距离panel；事件号是唯一标题。

fig = figure('Visible',opt.FigureVisible,'Color','w','Units','pixels', ...
    'Position',[30 20 1600 1750]);
layout = tiledlayout(fig,7,1,'TileSpacing','compact','Padding','compact');
ax = gobjects(7,1);

%% Panel 1-4：原始场和时间平均场
ax(1) = nexttile(layout);
plotVectorPanel(ax(1),raw.t,raw.B,tStart,'B','Raw');
ax(2) = nexttile(layout);
plotVectorPanel(ax(2),raw.t,raw.V,tStart,opt.VelocityField,'Raw');
meanLabel = sprintf('%g s mean',opt.SmoothSeconds);
ax(3) = nexttile(layout);
plotVectorPanel(ax(3),d.t,d.B,tStart,'B',meanLabel);
ax(4) = nexttile(layout);
plotVectorPanel(ax(4),d.t,d.V,tStart,opt.VelocityField,meanLabel);

x = d.t-tStart;
breaks = gapBreakRows(d.t);
maxError = 100;

%% Panel 5：磁场FOTE误差与类型
ax(5) = nexttile(layout);
hold(ax(5),'on');
bError = max([mag.eta mag.xi],[],2);
plot(ax(5),x,breakValues(min(mag.eta,maxError),breaks), ...
    'Color',[.35 .35 .35],'LineWidth',.55,'DisplayName','\eta_B');
plot(ax(5),x,breakValues(min(mag.xi,maxError),breaks), ...
    'Color',[.48 .29 .65],'LineWidth',.55,'DisplayName','\xi_B');
plotTypeMarkers(ax(5),x,bError,mag.type,mag.markerGood,maxError);
yline(ax(5),opt.QualityPercent,'k--','HandleVisibility','off');
ylabel(ax(5),'FOTE B error [%]');
ylim(ax(5),[0 maxError]);

%% Panel 6：FOTE-V连续性误差与类型
ax(6) = nexttile(layout);
hold(ax(6),'on');
plot(ax(6),x,breakValues(min(vel.alpha,maxError),breaks), ...
    'Color',[.35 .35 .35],'LineWidth',.55,'DisplayName','\alpha_V');
plotTypeMarkers(ax(6),x,vel.alpha,vel.type,vel.markerGood,maxError);
yline(ax(6),opt.QualityPercent,'k--','HandleVisibility','off');
ylabel(ax(6),sprintf('FOTE-V %s error [%%]',opt.VelocityField));
ylim(ax(6),[0 maxError]);

%% Panel 7：最终合格的screw-in/screw-out
ax(7) = nexttile(layout);
hold(ax(7),'on');
bOut = mag.sustainedGood & mag.type=="Bs";
bIn = mag.sustainedGood & mag.type=="As";
vOut = vel.sustainedGood & vel.type=="Bs";
vIn = vel.sustainedGood & vel.type=="As";
scatter(ax(7),x(bOut),3*ones(nnz(bOut),1),12,[.10 .30 .85],'>','filled');
scatter(ax(7),x(bIn),2*ones(nnz(bIn),1),12,[.85 .15 .15],'^','filled');
scatter(ax(7),x(vOut),ones(nnz(vOut),1),12,[.10 .30 .85],'>','filled');
scatter(ax(7),x(vIn),zeros(nnz(vIn),1),12,[.85 .15 .15],'^','filled');
yticks(ax(7),0:3);
yticklabels(ax(7), ...
    {sprintf('%s screw-in',opt.VelocityField), ...
    sprintf('%s screw-out',opt.VelocityField), ...
    'B screw-in','B screw-out'});
ylim(ax(7),[-.6 3.6]);
if opt.MinQualityDurationSeconds>0
    durationText = sprintf('run >= %g s',opt.MinQualityDurationSeconds);
else
    durationText = 'pointwise';
end
ylabel(ax(7),sprintf('error <= %g%%, M_\\lambda >= %g, %s', ...
    opt.QualityPercent,opt.EigenvalueStabilityThreshold,durationText));

%% 公共坐标轴格式
for k = 1:7
    grid(ax(k),'on');
    box(ax(k),'on');
    set(ax(k),'FontSize',9,'Layer','top');
    xlim(ax(k),[0 tStop-tStart]);
end
set(ax(1:6),'XTickLabel',[]);
formatTimeAxis(ax(7),tStart,tStop-tStart);
linkaxes(ax,'x');
sgtitle(layout,opt.eventID,'FontWeight','bold','FontSize',14, ...
    'Interpreter','none');
end


function plotVectorPanel(ax,t,field,tStart,fieldName,prefix)
% 仅在四星均有真实样本时绘制四星算术平均；任一卫星缺数即显示断点。
% 磁场panel同时绘制|B|。
hold(ax,'on');
x = t-tStart;
meanField = mean(field,3);
meanField = breakValues(meanField,gapBreakRows(t));
colors = [.10 .35 .85;.10 .60 .25;.85 .18 .18];
if strcmpi(fieldName,'B')
    labels = {'B_x','B_y','B_z'};
else
    labels = {'V_{ix}','V_{iy}','V_{iz}'};
end
for component = 1:3
    plot(ax,x,meanField(:,component),'Color',colors(component,:), ...
        'LineWidth',.65,'DisplayName',labels{component});
end
if strcmpi(fieldName,'B')
    plot(ax,x,sqrt(sum(meanField.^2,2)),'k','LineWidth',.85, ...
        'DisplayName','|B|');
    ylabel(ax,{prefix,'B_{GSE} [nT]'});
else
    ylabel(ax,{prefix,sprintf('%s_{GSE} [km/s]',fieldName)});
end
legend(ax,'Location','eastoutside','NumColumns',1);
end


function rows = gapBreakRows(t)
% 用NaN断开超过正常采样间隔10倍或1 s的线段。
if numel(t)<3, rows = []; return; end
dt = diff(t);
baseDt = median(dt(isfinite(dt) & dt>0),'omitnan');
if ~isfinite(baseDt)
    rows = [];
else
    rows = find(dt>max(10*baseDt,1))+1;
end
end


function out = breakValues(in,rows)
out = in;
out(rows,:) = NaN;
end


function plotTypeMarkers(ax,x,y,type,good,maxY)
% 将原FOTE类型映射为固定marker和颜色，便于不同事件之间比较。
spec = {
    "A", '^',[.85 .15 .15],false,'A';
    "B", '>',[.10 .30 .85],false,'B';
    "As",'^',[.85 .15 .15],true,'A_s screw-in';
    "Bs",'>',[.10 .30 .85],true,'B_s screw-out';
    "X", 'x',[.15 .15 .15],false,'X';
    "O", 'o',[.15 .15 .15],false,'O';
    "S+",'s',[.75 .15 .70],true,'source';
    "S-",'s',[.00 .60 .65],true,'sink'};
for k = 1:size(spec,1)
    idx = good & type==spec{k,1} & isfinite(y) & y<=maxY;
    if ~any(idx), continue; end
    h = scatter(ax,x(idx),y(idx),9,spec{k,3},spec{k,2}, ...
        'LineWidth',.55,'MarkerFaceAlpha',.65, ...
        'MarkerEdgeAlpha',.65,'DisplayName',spec{k,5});
    if spec{k,4}, h.MarkerFaceColor = spec{k,3}; end
end
legend(ax,'show','Location','eastoutside','NumColumns',1,'FontSize',8);
end


function formatTimeAxis(ax,tStart,duration)
ticks = linspace(0,duration,8);
xticks(ax,ticks);
dt = datetime(tStart+ticks,'ConvertFrom','posixtime','TimeZone','UTC');
xticklabels(ax,cellstr(string(dt,'HH:mm:ss')));
dateText = char(string(datetime(tStart,'ConvertFrom','posixtime', ...
    'TimeZone','UTC'),'yyyy-MM-dd'));
xlabel(ax,sprintf('UTC on %s',dateText));
end
