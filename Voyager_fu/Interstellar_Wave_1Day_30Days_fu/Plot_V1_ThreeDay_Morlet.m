function [outputs,audit] = Plot_V1_ThreeDay_Morlet(windows,spec,cfg)
% Three-day |B| and scalar Morlet power; no processing annotations in the figure.
visibility='off'; if cfg.Visible, visibility='on'; end
fig=figure('Color','w','Visible',visibility,'Position',[50 50 1650 850]);
ax(1)=axes(fig,'Position',[0.08 0.57 0.78 0.34]);
plot(ax(1),datenum(windows.EpochUTC),windows.BMean_nT,'k','LineWidth',0.65);
ylabel(ax(1),{'|B| 3-day mean','(nT)'});
titleText='Voyager 1 | 2012-08-25 to 2021-12-16';
if isfield(cfg,'DisplayTitle'), titleText=cfg.DisplayTitle; end
title(ax(1),titleText,'FontSize',15,'FontWeight','bold');
ax(2)=axes(fig,'Position',[0.08 0.13 0.78 0.34]);
t=datenum(spec.TimeUTC(:)); f=spec.FrequencyHz(:);
te=[t-cfg.SampleDays/2;t(end)+cfg.SampleDays/2];
fe=[f(1);sqrt(f(1:end-1).*f(2:end));f(end)];
c=nan(numel(fe),numel(te)); v=spec.Power.'; v(v<=0)=NaN;
c(1:end-1,1:end-1)=log10(v);
periodEdgesDays=1./fe/86400; % Coordinate change only; PSD stays per Hz.
surface(ax(2),te.',periodEdgesDays,zeros(size(c)),c,'EdgeColor','none','FaceColor','flat');
view(ax(2),2); colormap(ax(2),jet(256));
validLog=log10(spec.Power(isfinite(spec.Power)&spec.Power>0));
assert(~isempty(validLog),'No usable wavelet coefficients.');
limits=cfg.PSDColorLimits;
if limits(1)==limits(2), limits=limits+[-1 1]; end
clim(ax(2),limits);
set(ax(2),'YScale','log','YDir','reverse','YLim',cfg.DisplayPeriodDays, ...
    'YTick',[6 10 15 20 30],'YTickLabel',{'6','10','15','20','30'});
ylabel(ax(2),'周期（天）','FontName','Microsoft YaHei'); xlabel(ax(2),'UTC');
cb=colorbar(ax(2),'Position',[0.885 0.13 0.012 0.34]);
cb.Label.String={'log_{10} P_{|B|}','(nT^2 Hz^{-1})'}; cb.FontSize=12;
cb.Ticks=[-1 0 1 2 3];
years=(year(cfg.StartUTC)+1):year(cfg.StopUTC);
if isempty(years), years=year(cfg.StartUTC); end
ticks=datenum(datetime(years,1,1));
for k=1:2
    set(ax(k),'FontSize',12,'Box','on','TickDir','out','Layer','top', ...
        'XLim',datenum([cfg.StartUTC cfg.StopUTC]),'XTick',ticks,'XTickLabel',[]);
    text(ax(k),0.01,0.91,sprintf('(%c)','a'+k-1),'Units','normalized', ...
        'FontWeight','bold','FontSize',12);
end
set(ax(2),'XTickLabel',string(years));
linkaxes(ax,'x');
base=fullfile(cfg.OutputRoot,'V1_three_day_B_Morlet_6day_30days');
exportgraphics(fig,[base '.png'],'Resolution',180);
exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',180);
savefig(fig,[base '.fig']);
outputs=string({[base '.png'];[base '.pdf'];[base '.fig']});
audit=struct('ColorLimitsLog10',limits,'PowerCells','native three-day power, flat display', ...
    'VerticalAxis','Period (days)','YDirection','reverse','PeriodLimitsDays',cfg.DisplayPeriodDays,'AboveNyquist','Outside displayed range','MissingWindows','Linear interpolation of empty three-day means after original averaging','MethodAnnotationsAdded',false);
if ~cfg.Visible, close(fig); end
end




