clear;clc

%% input
InputDir = 'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS';
InputFiles = { ...
    fullfile(InputDir,'SMILE_UVI_L2_AURORA-GEO_20260720T230011(1).cdf'), ...
    fullfile(InputDir,'SMILE_UVI_L2_AURORA-GEO_20260720T230041(1).cdf')};
ProductType = 'AURORA-GEO';
OutputDir = fullfile(InputDir,'SMILE_UVI_readout_20261003');
IrfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';
ReaderFile = fullfile(fileparts(mfilename('fullpath')),'Read_SMILE_UVI.m');

%% load data
AllUVI = cell(1,numel(InputFiles));
AllInfo = cell(1,numel(InputFiles));
AllUTC = cell(1,numel(InputFiles));
for iInput = 1:numel(InputFiles)
    FilePath = InputFiles{iInput};
    FileType = ProductType;
    ic = 1;
    run(ReaderFile);
    assert(isscalar(UVI),'主程序每次读取一个文件。');
    assert(isscalar(Epoch),'当前主程序要求每个文件包含一帧图像。');
    AllUVI{iInput} = Data;
    AllInfo{iInput} = Info;
    AllUTC{iInput} = UTC;
end

%% image data
ImageNames = {'SMILE_UVI_IMAGE','SMILE_UVI_GRID_IMG_GEO','SMILE_UVI_GRID_IMG_AACGM'};
ImageTitles = {'Corrected image','GEO grid image','AACGM grid image'};
PlotImages = cell(numel(InputFiles),3);
PlotUnits = cell(numel(InputFiles),3);
ColorLimits = [inf(3,1),-inf(3,1)];
for iInput = 1:numel(InputFiles)
    for iPanel = 1:3
        PlotName = ImageNames{iPanel};
        RawImage = AllUVI{iInput}.(PlotName);
        assert(ismatrix(RawImage),'图像必须是二维数组：%s',PlotName);
        FillEntries = AllInfo{iInput}.VariableAttributes.FILLVAL;
        FillRow = find(strcmp(FillEntries(:,1),PlotName));
        assert(isscalar(FillRow),'图像FILLVAL属性不唯一或缺失：%s',PlotName);
        FillValue = cast(FillEntries{FillRow,2},'like',RawImage);
        PlotImages{iInput,iPanel} = double(RawImage);
        PlotImages{iInput,iPanel}(RawImage==FillValue) = NaN;
        GoodValues = PlotImages{iInput,iPanel}(isfinite(PlotImages{iInput,iPanel}));
        assert(~isempty(GoodValues),'图像没有可显示的有效像素：%s',PlotName);
        ColorLimits(iPanel,1) = min(ColorLimits(iPanel,1),min(GoodValues));
        ColorLimits(iPanel,2) = max(ColorLimits(iPanel,2),max(GoodValues));
        UnitEntries = AllInfo{iInput}.VariableAttributes.UNITS;
        UnitRow = find(strcmp(UnitEntries(:,1),PlotName));
        assert(isscalar(UnitRow),'图像UNITS属性不唯一或缺失：%s',PlotName);
        PlotUnits{iInput,iPanel} = UnitEntries{UnitRow,2};
        assert(strcmp(PlotUnits{iInput,iPanel},PlotUnits{1,iPanel}), ...
            '同类图像的单位不同，不能共用色标。');
    end
end
ColorLimits(:,2) = ColorLimits(:,2)*3/5;
ColorLimits(1,2) = ColorLimits(1,2)*3/5;
assert(all(ColorLimits(:,2)>ColorLimits(:,1)),'图像色标范围必须非零。');

%% image plot
if ~exist(OutputDir,'dir'), mkdir(OutputDir); end
FigHandles = gobjects(1,numel(InputFiles));
OutputFiles = cell(1,numel(InputFiles));
for iInput = 1:numel(InputFiles)
    [~,InputName] = fileparts(InputFiles{iInput});
    FigHandles(iInput) = figure('Color','w','Position',[80 80 1650 650], ...
        'Name',InputName,'NumberTitle','off');
    Layout = tiledlayout(FigHandles(iInput),1,3,'TileSpacing','compact','Padding','compact');
    for iPanel = 1:3
        Ax = nexttile(Layout);
        PlotImage = PlotImages{iInput,iPanel};
        imagesc(Ax,PlotImage,'AlphaData',isfinite(PlotImage));
        axis(Ax,'image');
        set(Ax,'YDir','reverse','Color',[0.85 0.85 0.85],'FontSize',12);
        colormap(Ax,jet(256));
        clim(Ax,ColorLimits(iPanel,:));
        xlabel(Ax,'Array column');
        ylabel(Ax,'Array row');
        title(Ax,sprintf('%s (%d x %d)',ImageTitles{iPanel}, ...
            size(PlotImage,1),size(PlotImage,2)),'FontWeight','normal');
        CB = colorbar(Ax);
        CB.Label.String = PlotUnits{iInput,iPanel};
    end
    EpochText = char(AllUTC{iInput}.EPOCH);
    StartText = char(AllUTC{iInput}.SMILE_UVI_EXPO_START);
    EndText = char(AllUTC{iInput}.SMILE_UVI_EXPO_END);
    sgtitle(Layout,{['SMILE UVI | EPOCH ' strrep(EpochText(1:end-1),'T',' ') ' UTC'], ...
        ['Exposure: ' strrep(StartText(1:end-1),'T',' ') ' to ' ...
        strrep(EndText(1:end-1),'T',' ') ' UTC']},'FontSize',14,'FontWeight','normal');
    xlabel(Layout,'Array rows increase downward; grey = FILLVAL','FontSize',11);

    %% save figure
    OutputFiles{iInput} = fullfile(OutputDir,[InputName '_three_images.png']);
    exportgraphics(FigHandles(iInput),OutputFiles{iInput},'Resolution',150);
    fprintf('已保存：%s\n',OutputFiles{iInput});
end