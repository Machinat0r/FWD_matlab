function checks=Run_Voyager_Daily_Boundaries(syncArchive)
% Two new daily figures from original CDFs; existing deliverables are read-only.
if nargin<1, syncArchive=false; end
%% Paths and preserved originals
CodeDir=fileparts(mfilename('fullpath')); addpath(CodeDir);
addpath(fullfile(fileparts(CodeDir),'Case1_PPT_VerticalLine_Events_7d_fu'));
OriginalDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Boundary_Markers_1990/';
paths=strings(0,1); before=strings(0,1);
for sc=1:2
    files=dir(fullfile(OriginalDir,sprintf('V%d',sc),'*_5panels.*'));
    for j=1:numel(files)
        path=fullfile(files(j).folder,files(j).name);
        paths(end+1,1)=string(path); before(end+1,1)=string(Case1_File_SHA256(path));
    end
end
%% Read raw CDFs, calculate with existing definitions, add the requested markers
Run_Voyager_Extended_Overviews(1:2,syncArchive,'boundary_markers');
%% Verify new data against existing results; no audit data are science inputs
checks=struct;
for sc=1:2
    original=load(fullfile(OriginalDir,sprintf('V%d',sc),'extended_overviews_audit.mat'),'audit');
    marked=load(fullfile(OutputDir,sprintf('V%d',sc),'extended_overviews_audit.mat'),'audit');
    assert(isequaln(original.audit.all_daily.Windows,marked.audit.all_daily.Windows), ...
        'Daily values differ from the existing figure.');
    newFig=openfig(marked.audit.all_daily.OutputFiles(3),'invisible');
    oldFig=openfig(original.audit.all_daily.OutputFiles(3),'invisible');
    ax=findall(newFig,'Type','axes'); oldAx=findall(oldFig,'Type','axes');
    assert(numel(ax)==5 && numel(oldAx)==5);
    dates=Voyager_Boundary_Dates(sc);
    for k=1:5
        tile=ax(k).Layout.Tile;
        other=oldAx(arrayfun(@(h) h.Layout.Tile==tile,oldAx));
        newLine=findall(ax(k),'Type','line'); oldLine=findall(other,'Type','line');
        assert(isequaln(newLine.XData,oldLine.XData) && isequaln(newLine.YData,oldLine.YData));
        assert(isequal(ax(k).XLim,other.XLim) && isequal(ax(k).YLim,other.YLim));
        assert(strcmp(ax(k).YScale,other.YScale));
        for j=1:2
            line=findall(ax(k),'Tag',char(dates.ShortName(j)));
            assert(numel(line)==1 && strcmp(line.LineStyle,'--') && line.Value==dates.TimeUTC(j));
        end
    end
    close(newFig); close(oldFig);
    checks.(sprintf('V%d',sc))=struct('DailyValuesUnchanged',true,'AxesAndCurvesUnchanged',true, ...
        'TwoMarkersPerPanelVerified',true,'Dates',dates,'OutputFiles',marked.audit.all_daily.OutputFiles);
end
%% Confirm that no existing figure file was overwritten
for k=1:numel(paths), assert(before(k)==string(Case1_File_SHA256(paths(k)))); end
checks.OriginalFilesUnchanged=numel(paths);
checks.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(OutputDir,'validation_audit.mat'),'checks');
writetable(table(paths,before,'VariableNames',{'OriginalPath','UnchangedSHA256'}), ...
    fullfile(OutputDir,'original_figures_unchanged.csv'));
disp(checks);
end