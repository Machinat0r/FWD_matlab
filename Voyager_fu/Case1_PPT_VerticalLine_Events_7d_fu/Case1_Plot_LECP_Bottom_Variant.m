function output = Case1_Plot_LECP_Bottom_Variant(ax, output, opts)
%Case1_Plot_LECP_Bottom_Variant Baseline-difference PAD or sector heat map.
%   Input comes from the current raw-CDF pipeline in this same run.
%   Keep all S1--S7 values independent. No angular or sector averaging.
%   BaselineUTC is an explicitly approved half-open interval.
%   Modified: 2026-09-08

%% original sector values and display selection
suffix = '1d';
halfWidth = 0.5;
if strcmp(opts.PADCadence, 'hour'), suffix = '1h'; halfWidth = 1/48; end
flux = output{:, cellstr(compose(['Flux_S%d_', suffix], 1:7))};
time = output.EpochUTC;
fluxValid = isfinite(flux) & flux > 0;
allSeven = all(fluxValid, 2);
mode = opts.LECPBottomMode;
audit = output.Properties.UserData;
audit.SectorMergeApplied = false;
audit.BottomMode = mode;
audit.SectorOrderTopToBottom = [5 6 7 8 1 2 3 4];
audit.S8Policy = 'Excluded blocked sector; row 8 stays blank.';
audit.DisplayTimeHalfWidthDays = halfWidth;
audit.DisplayWidthMeaning = 'Nominal cadence glyph, not original exposure support.';
output.SectorUsable = allSeven;
values = nan(size(flux));

%% approved sector-specific temporal baseline subtraction
if strcmp(mode, 'pad_difference')
    interval = opts.BaselineUTC;
    interval.TimeZone = 'UTC';
    inBaseline = time >= interval(1) & time < interval(2);
    assert(any(inBaseline), 'No retained CDF records fall in the baseline interval.');
    selected = flux(inBaseline, :);
    selected(~fluxValid(inBaseline, :)) = NaN;
    baseline = mean(selected, 1, 'omitnan');
    counts = sum(isfinite(selected), 1);
    assert(all(counts > 0), 'A sector has no finite positive baseline flux.');
    values = flux - baseline;
    usable = logical(output.PADUsable);
    pitch = output{:, cellstr(compose('PA_S%d_deg', 1:7))};
    audit.BaselineUTC = interval;
    audit.BaselineRowIndices = find(inBaseline);
    audit.BaselineEpochUTC = time(inBaseline);
    audit.BaselineFluxBySector = selected;
    audit.BaselineMeanBySector = baseline;
    audit.BaselineCountBySector = counts;
    audit.BaselinePolicy = ['Arithmetic mean of retained daily/hourly flux records ', ...
        'within the approved interval, independently for each sector; ', ...
        'missing/nonpositive source flux excluded. Subtract in linear flux units.'];
    audit.DifferencePolicy = 'Signed J_s(t)-baseline_s; negatives and zeros retained.';
    audit.PitchAngleCalculated = true;
    audit.MagneticFieldRequired = true;
    for sector = 1:7
        output.(sprintf('BaselineFlux_S%d', sector)) = ...
            repmat(baseline(sector), height(output), 1);
        output.(sprintf('DifferenceFlux_S%d', sector)) = values(:, sector);
    end
    colorValues = values(repmat(usable, 1, 7) & isfinite(values));
    assert(~isempty(colorValues), 'No drawable baseline-difference PAD.');
    bound = prctile(abs(colorValues), opts.ColorPercentiles(2));
    if bound == 0, bound = 1; end
    if ~isempty(opts.DifferenceColorLimits)
        bound = opts.DifferenceColorLimits(2);
    end
    % Keep the previous daily color limits while averaging each sector.
    % A complete UTC grid prevents crossing gaps by record index.
    audit.PADDisplayAverageDays = opts.PADDisplayAverageDays;
    if opts.PADDisplayAverageDays == 3
        utcDay = dateshift(time,'start','day');
        assert(numel(unique(utcDay)) == numel(utcDay), ...
            'Three-day averaging requires at most one input record per UTC day.');
        dayGrid = (min(utcDay):days(1):max(utcDay)).';
        [~, gridIndex] = ismember(utcDay,dayGrid);
        gridValues = nan(numel(dayGrid),7);
        gridValues(gridIndex,:) = values;
        averagedGrid = movmean(gridValues,3,1,'omitnan');
        sampleCounts = movsum(double(isfinite(gridValues)),3,1);
        values = averagedGrid(gridIndex,:);
        for sector = 1:7
            output.(sprintf('DifferenceFlux3d_S%d',sector)) = values(:,sector);
            output.(sprintf('AverageDayCount_S%d',sector)) = ...
                sampleCounts(gridIndex,sector);
        end
        audit.AverageWindowStartUTC = utcDay-days(1);
        audit.AverageWindowEndUTCExclusive = utcDay+days(2);
        audit.AveragePolicy = ['Centered UTC day-1/day/day+1 arithmetic mean ', ...
            'of signed differences independently per sector; finite days only, ', ...
            'including partial windows at plot endpoints. Original center-day ', ...
            'PAD mask and PA retained; missing center days are not filled.'];
        audit.ColorRangeReference = 'Unaveraged daily differences; unchanged color scale.';
        colorValues = values(repmat(usable,1,7) & isfinite(values));
    end
    audit.NegativeDisplayAsMissing = logical(opts.DifferenceNegativeAsMissing);
    audit.NegativeCellsHidden = 0;
    if opts.DifferenceNegativeAsMissing
        negative = values < 0;
        audit.NegativeCellsHidden = nnz(negative & repmat(usable,1,7));
        audit.NegativePolicy = ['After temporal averaging, negative differences ', ...
            'become NaN separately per sector. Zero remains valid. Signed daily ', ...
            'and averaged differences retained in audit. Original seven-sector ', ...
            'PA geometry and center-day mask unchanged; no gap expansion.'];
        for sector = 1:7
            output.(sprintf('NegativeDifferenceMasked_S%d',sector)) = negative(:,sector);
        end
        values(negative) = NaN;
        colorValues = values(repmat(usable,1,7) & isfinite(values));
    end
    limits = [0 bound];
    if ~isempty(opts.DifferenceColorLimits)
        audit.ColorRangeReference = 'Explicit common limits for baseline comparisons.';
    end
    audit.ColorLimitPolicy = 'Linear difference scale from zero to the existing upper bound; negative values use the lowest color, upper tails saturate for display only.';
    if opts.DifferenceNegativeAsMissing
        audit.ColorLimitPolicy = 'Linear scale from zero; negative differences are blank; upper tails saturate for display only.';
    end
    audit.AbsoluteDifferenceColorPercentile = opts.ColorPercentiles(2);
    audit.ColorSaturatedCellCount = nnz(colorValues<0 | colorValues>bound);
    audit.NegativeCellsShownAtLowestColor = nnz(colorValues<0);
    audit.ColorMapName = 'turbo';
    cmap = turbo(256);
    label = {'\DeltaJ', '(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'};
    note = sprintf('Sector baseline: %s to %s UTC; J - baseline; values below 0 use lowest color', ...
        datestr(interval(1), 'dd-mmm-yyyy'), ...
        datestr(interval(2)-seconds(1), 'dd-mmm-yyyy'));
    if opts.PADDisplayAverageDays == 3
        note = ['3-day centered sector mean; ', note];
        label{1} = '\langle\DeltaJ\rangle_{3d}';
    end
    if opts.DifferenceNegativeAsMissing
        note = strrep(note,'values below 0 use lowest color','negative differences blank');
    end
else
    %% sector-number map: no magnetic field or attitude requirement
    usable = allSeven;
    if strcmp(opts.SectorColorMode, 'normalized')
        denominator = max(flux, [], 2);
        denominator(~usable) = NaN;
        values = flux ./ denominator;
        output.NormalizingMaximumFlux = denominator;
        limits = [0 1];
        label = 'J_s / max(J_1,...,J_7)';
        note = 'Each UTC bin normalized by its own largest S1-S7 flux; S8 excluded';
    else
        values = flux;
        z = log10(flux(repmat(usable, 1, 7)));
        assert(~isempty(z), 'No complete seven-sector flux record.');
        limits = [min(z), max(z)];
        if limits(1) == limits(2), limits = limits + [-0.5 0.5]; end
        label = 'log_{10} J';
        note = 'Original sector flux; S8 excluded';
    end
    cmap = turbo(256);
    audit.PitchAngleCalculated = false;
    audit.MagneticFieldRequired = false;
    audit.SectorColorMode = opts.SectorColorMode;
    audit.NormalizationPolicy = ...
        'Normalized mode uses the same-time S1-S7 maximum; requires all seven positive finite fluxes.';
end
output.BottomPanelUsable = usable;
audit.ColorLimits = limits;
audit.ColorLabel = label;
audit.ZeroAndNegativeDifferenceRetained = strcmp(mode, 'pad_difference');

%% seven separate cells at every usable epoch
orderTop = [5 6 7 8 1 2 3 4];
for sector = 1:7
    output.(sprintf('DisplayValue_S%d', sector)) = nan(height(output), 1);
    output.(sprintf('DisplayLowerEdge_S%d_deg', sector)) = nan(height(output), 1);
    output.(sprintf('DisplayUpperEdge_S%d_deg', sector)) = nan(height(output), 1);
end
for row = find(usable).'
    if strcmp(mode, 'pad_difference')
        [centers, sectors] = sort(pitch(row, :));
        middle = (centers(1:end-1) + centers(2:end))/2;
        edges = [max(0, 2*centers(1)-middle(1)), middle, ...
            min(180, 2*centers(end)-middle(end))];
    else
        sectors = 1:7;
    end
    x = datenum(time(row)); %#ok<DATNM>
    for k = 1:7
        sector = sectors(k);
        value = values(row, sector);
        if strcmp(mode, 'pad_difference')
            low = edges(k); high = edges(k+1);
            colorValue = value;
            output.(sprintf('DisplayPA_S%d_deg', sector))(row) = pitch(row, sector);
            output.(sprintf('DisplayFlux_S%d', sector))(row) = value;
            output.(sprintf('DisplayOrder_S%d', sector))(row) = k;
            output.(sprintf('DisplayLowerEdge_S%d_deg', sector))(row) = low;
            output.(sprintf('DisplayUpperEdge_S%d_deg', sector))(row) = high;
        else
            position = find(orderTop == sector);
            low = position-0.5; high = position+0.5;
            colorValue = value;
            if strcmp(opts.SectorColorMode, 'absolute'), colorValue = log10(value); end
        end
        output.(sprintf('DisplayValue_S%d', sector))(row) = value;
        if ~isfinite(colorValue), continue, end
        surface(ax, [x-halfWidth x+halfWidth; x-halfWidth x+halfWidth], ...
            [low low; high high], zeros(2), ones(2)*colorValue, ...
            'FaceColor', 'flat', 'EdgeColor', 'none', 'HandleVisibility', 'off', ...
            'UserData', struct('Sector',sector,'Value',value,'BottomMode',mode));
    end
end
view(ax, 2);
if strcmp(mode, 'pad_difference')
    set(ax, 'YDir', 'normal', 'YLim', [0 180], 'YTick', [0 45 90 135 180]);
else
    set(ax, 'YDir', 'reverse', 'YLim', [0.5 8.5], ...
        'YTick', 1:8, 'YTickLabel', string(orderTop));
end
colormap(ax, cmap);
clim(ax, limits);
cb = colorbar(ax, 'Location', 'eastoutside');
cb.Label.String = label;
cb.Label.Interpreter = 'tex';
cb.FontSize = 8;
if strcmp(mode, 'pad_difference')
    cb.Ticks = linspace(limits(1),limits(2),5);
    labels = string(compose('%.3g',cb.Ticks));
    labels(1) = "<= " + labels(1);
    if opts.DifferenceNegativeAsMissing, labels(1) = "0"; end
    cb.TickLabels = labels;
end
% User request: method notes belong in the audit, not inside the figure.
audit.ProcessingNote = note;
audit.ProcessingNoteShownInFigure = false;
output.Properties.UserData = audit;
end
