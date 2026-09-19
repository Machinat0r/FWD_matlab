function outputs = Plot_V1_Wave_1Day_30Days(raw,spec,cfg)
% Stacked field and frequency panels; all explanatory text stays in README.
visibility='off'; if cfg.Visible, visibility='on'; end
fig=figure('Color','w','Visible',visibility,'Position',[40 40 1600 1400]);
ax=gobjects(6,1); left=0.085; width=0.755; height=0.125; gap=0.025;
for k=1:6
    ax(k)=axes(fig,'Position',[left 0.835-(k-1)*(height+gap) width height]);
    set(ax(k),'FontSize',12,'Box','on','TickDir','out','Layer','top');
end
t=datenum(raw.EpochUTC);
% Insert NaN only as a plotting break when a source Epoch is entirely absent.
breaks=find(seconds(diff(raw.EpochUTC))>3600)+1;
plotT=t; values=[raw.B_RTN_nT raw.Bmag_nT];
for k=flip(breaks(:).')
    plotT=[plotT(1:k-1);(t(k-1)+t(k))/2;plotT(k:end)]; %#ok<AGROW>
    values=[values(1:k-1,:);nan(1,4);values(k:end,:)]; %#ok<AGROW>
end
colors=[0 0.32 0.8;0 0.6 0.2;0.85 0.12 0.10;0 0 0];
hold(ax(1),'on');
for k=1:4, plot(ax(1),plotT,values(:,k),'Color',colors(k,:),'LineWidth',0.45); end
ylabel(ax(1),'B (nT)');
legend(ax(1),{'B_R','B_T','B_N','|B|'},'Orientation','horizontal','Location','northwest','Box','off');
plotTitle='Voyager 1 | 2012-08-25 to 2021-12-16'; if isfield(cfg,'DisplayTitle'), plotTitle=cfg.DisplayTitle; end
title(ax(1),plotTitle,'FontWeight','bold');
power={spec.TracePSD,spec.MagnitudePSD,spec.ComponentPSD(:,:,1), ...
    spec.ComponentPSD(:,:,2),spec.ComponentPSD(:,:,3)};
labels={'P_{B_R}+P_{B_T}+P_{B_N}','P_{|B|}','P_{B_R}','P_{B_T}','P_{B_N}'};
finiteLog=[];
for k=1:5
    v=power{k}; v=v(isfinite(v)&v>0); finiteLog=[finiteLog;log10(v(:))]; %#ok<AGROW>
end
assert(~isempty(finiteLog),'No finite power to plot.');
limits=[max(floor(min(finiteLog)),ceil(max(finiteLog))-6),ceil(max(finiteLog))];
% Only the displayed color scale is limited to six decades; all PSD values are retained.
if limits(1)==limits(2), limits=limits+[-1 1]; end
for k=1:5
    drawSpectrum(ax(k+1),spec.TimeUTC,spec.FrequencyHz,power{k});
    clim(ax(k+1),limits); colormap(ax(k+1),jet(256));
    cb=colorbar(ax(k+1),'Position',[0.86 ax(k+1).Position(2) 0.012 height]);
    cb.Label.String={['log_{10} ' labels{k}],'(nT^2 Hz^{-1})'}; cb.FontSize=11;
    ylabel(ax(k+1),'f (Hz)');
end
tickYears=(year(cfg.StartUTC)+1):year(cfg.StopUTC);
if isempty(tickYears), tickYears=year(cfg.StartUTC); end
ticks=datenum(datetime(tickYears,1,1));
for k=1:6
    xlim(ax(k),datenum([cfg.StartUTC cfg.StopUTC]));
    set(ax(k),'XTick',ticks,'XTickLabel',[]);
    text(ax(k),0.008,0.10,sprintf('(%c)','a'+k-1),'Units','normalized','FontWeight','bold','FontSize',12);
end
set(ax(6),'XTickLabel',string(tickYears)); xlabel(ax(6),'UTC');
linkaxes(ax,'x');
base=fullfile(cfg.OutputRoot,'V1_magnetic_wave_1day_30days');
exportgraphics(fig,[base '.png'],'Resolution',180);
exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',180);
savefig(fig,[base '.fig']);
outputs=string({[base '.png'];[base '.pdf'];[base '.fig']});

%% Separate sampling diagnostic; no extra annotations in the science figure
diagnostic=figure('Color','w','Visible',visibility,'Position',[80 80 1500 390]);
h=axes(diagnostic,'Position',[0.09 0.20 0.79 0.67]);
drawSpectrum(h,spec.TimeUTC,spec.FrequencyHz,spec.SamplingWindow);
clim(h,[-6 0]); colormap(h,parula(256));
cb=colorbar(h); cb.Label.String='log_{10} |W(f)|^2';
ylabel(h,'f (Hz)'); xlabel(h,'UTC');
xlim(h,datenum([cfg.StartUTC cfg.StopUTC]));
set(h,'XTick',ticks,'XTickLabel',string(tickYears),'FontSize',12,'TickDir','out');
title(h,[plotTitle ' | Sampling window']);
file=fullfile(cfg.OutputRoot,'V1_sampling_window.png');
exportgraphics(diagnostic,file,'Resolution',180); outputs(end+1,1)=string(file);
if ~cfg.Visible, close(fig); close(diagnostic); end
end

function drawSpectrum(ax,time,f,power)
% Flat cells use their computed value; no display interpolation.
t=datenum(time(:));
te=[t-0.5;t(end)+0.5]; % StepDays is explicitly 1 in the formal entry.
fe=[f(1);sqrt(f(1:end-1).*f(2:end));f(end)];
c=nan(numel(fe),numel(te)); v=power.'; v(v<=0)=NaN;
c(1:end-1,1:end-1)=log10(v);
surface(ax,te.',fe,zeros(size(c)),c,'EdgeColor','none','FaceColor','flat');
view(ax,2); set(ax,'YScale','log','YLim',[f(1) f(end)], ...
    'YTick',[1/(30*86400),1e-6,3e-6,1/86400], ...
    'YTickLabel',{'3.86\times10^{-7}','10^{-6}','3\times10^{-6}','1.16\times10^{-5}'},'Box','on','Layer','top');
end



