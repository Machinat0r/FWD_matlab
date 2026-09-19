function [outputs,audit] = Plot_V1_Daily_Fourier(daily,spec,cfg)
% Daily |B| and short-time Fourier power; no processing annotations in the figure.
visibility='off'; if cfg.Visible, visibility='on'; end
fig=figure('Color','w','Visible',visibility,'Position',[50 50 1650 850]);
ax(1)=axes(fig,'Position',[0.08 0.57 0.78 0.34]);
plot(ax(1),datenum(daily.EpochUTC),daily.BMean_nT,'k','LineWidth',0.65);
ylabel(ax(1),{'|B| daily mean','(nT)'});
titleText='Voyager 1 | 2012-08-25 to 2021-12-16';
if isfield(cfg,'DisplayTitle'), titleText=cfg.DisplayTitle; end
title(ax(1),titleText,'FontSize',15,'FontWeight','bold');
ax(2)=axes(fig,'Position',[0.08 0.13 0.78 0.34]);
t=datenum(spec.TimeUTC(:)); f=spec.FrequencyHz(:);
te=[t-cfg.StepDays/2;t(end)+cfg.StepDays/2];
fe=[cfg.ComputableFrequencyHz(1);(f(1:end-1)+f(2:end))/2;cfg.ComputableFrequencyHz(2)];
c=nan(numel(fe),numel(te)); v=spec.Power.'; v(v<=0)=NaN;
c(1:end-1,1:end-1)=log10(v);
periodEdgesDays=1./fe/86400; % Coordinate change only; PSD stays per Hz.
surface(ax(2),te.',periodEdgesDays,zeros(size(c)),c,'EdgeColor','none','FaceColor','flat');
view(ax(2),2); colormap(ax(2),jet(256));
validLog=log10(spec.Power(isfinite(spec.Power)&spec.Power>0));
assert(~isempty(validLog),'No usable Fourier power.');
limits=cfg.PSDColorLimits;
if limits(1)==limits(2), limits=limits+[-1 1]; end
clim(ax(2),limits);
periodTicks=unique([2 3 5 10 20 30 45 60 90 cfg.DisplayPeriodDays]);
periodTicks=periodTicks(periodTicks>=cfg.DisplayPeriodDays(1) & periodTicks<=cfg.DisplayPeriodDays(2));
set(ax(2),'YScale','log','YDir','reverse','YLim',cfg.DisplayPeriodDays, ...
    'YTick',periodTicks,'YTickLabel',string(periodTicks));
ylabel(ax(2),'周期（天）','FontName','Microsoft YaHei'); xlabel(ax(2),'UTC');
cb=colorbar(ax(2),'Position',[0.885 0.13 0.012 0.34]);
cb.Label.String={'log_{10} P_{|B|}','(nT^2 Hz^{-1})'}; cb.FontSize=12;
cb.Ticks=unique([limits(1),ceil(limits(1)):floor(limits(2)),limits(2)]);
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
base=fullfile(cfg.OutputRoot,sprintf('V1_daily_B_Fourier_2day_%gdays',cfg.DisplayPeriodDays(2)));
exportgraphics(fig,[base '.png'],'Resolution',180);
exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',180);
savefig(fig,[base '.fig']);
outputs=string({[base '.png'];[base '.pdf'];[base '.fig']});
audit=struct('ColorLimitsLog10',limits,'PowerCells','native FFT bins, flat display, one-day hop', ...
    'VerticalAxis','Period (days)','YDirection','reverse','PeriodLimitsDays',cfg.DisplayPeriodDays,'AboveNyquist','Outside displayed range','MissingDays','User-authorized linear interpolation in the daily input','MethodAnnotationsAdded',false);
if ~cfg.Visible, close(fig); end
end





