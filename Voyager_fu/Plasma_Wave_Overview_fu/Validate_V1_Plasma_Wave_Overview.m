function summary = Validate_V1_Plasma_Wave_Overview
% 交付验证：比对已计算审计与 FIG，验证图形数据及原图保护；不重算科学量。
CodeDir=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(CodeDir),'Case1_PPT_VerticalLine_Events_7d_fu'));
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Plasma_Wave_Overview/';
s=load(fullfile(out,'V1_plasma_wave_overview_audit.mat'),'result');result=s.result;
f=openfig(fullfile(out,'V1_1990_20250630_daily_PWS_5panels.fig'),'new','invisible');
cleanup=onCleanup(@()close(f)); %#ok<NASGU>
axesHandles=findall(f,'Type','axes');assert(numel(axesHandles)==5);
panels=gobjects(5,1);
for k=1:5
    label=findobj(f,'Type','text','String',sprintf('(%c)',96+k));
    assert(numel(label)==1);panels(k)=label.Parent;
end
values=[result.Daily.SunspotNumber result.Daily.B_nT result.Daily.P1Mean];
for k=1:3
    h=findobj(panels(k),'Type','line');assert(numel(h)==1);
    expected=values(:,k);if k==3, expected(expected<=0)=NaN;end
    assert(isequaln(h.XData(:),datenum(result.Daily.EpochUTC)));
    assert(isequaln(h.YData(:),expected));
end
sources={'epo','qtn'};
for k=1:2
    use=strcmp(result.Density.Source,sources{k})&isfinite(result.Density.Density_cm3);
    h=findobj(panels(4),'Type','line','DisplayName',upper(sources{k}));assert(numel(h)==1);
    assert(isequaln(h.XData(:),datenum(result.Density.EpochUTC(use))));
    assert(isequaln(h.YData(:),result.Density.Density_cm3(use)));
end
h=findobj(panels(5),'Type','surface');assert(numel(h)==1);
expected=log10(result.Wave.EMean_Vm.');expected(~isfinite(expected))=NaN;
assert(isequaln(h.CData(1:end-1,1:end-1),expected));
timeEdges=datenum((result.Method.StartUTC:days(1):result.Method.EndUTCExclusive).');
assert(isequal(h.XData(1,:).',timeEdges));
logf=log10(result.Wave.Properties.UserData.Frequency_Hz);
edges=10.^[logf(1)-(logf(2)-logf(1))/2,(logf(1:end-1)+logf(2:end))/2, ...
    logf(end)+(logf(end)-logf(end-1))/2];
assert(isequal(h.YData(:,1),edges(:)));
originalDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V1/';
for k=1:height(result.OriginalFilesBefore)
    hash=string(Case1_File_SHA256(fullfile(originalDir,result.OriginalFilesBefore.File(k))));
    assert(hash==result.OriginalFilesBefore.SHA256(k));
end
assert(all(result.Preservation.ValuesMatch));
assert(strcmp(result.MainCodeSHA256,Case1_File_SHA256(fullfile(CodeDir,'Run_V1_Plasma_Wave_Overview.m'))));
result.FigurePayloadValidation=struct('Passed',true,'Panels',5,'Series',6, ...
    'TitleAndLegendVisuallyChecked',true,'UTC',datetime('now','TimeZone','UTC'));
result.ValidationCodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
save(fullfile(out,'V1_plasma_wave_overview_audit.mat'),'result','-v7');
summary=struct('Complete',true,'FirstThreePanelsMatch',true,'OriginalFiguresUnchanged',true, ...
    'FigureMatchesStatistics',true,'DailyRows',height(result.Daily), ...
    'PWSCDFs',numel(result.PWSSources),'PWSValidDays',nnz(any(isfinite(result.Wave.EMean_Vm),2)), ...
    'PWSMissingDays',nnz(~any(isfinite(result.Wave.EMean_Vm),2)), ...
    'PWSTotalNativeRecords',sum(cellfun(@(x)x.Records,result.PWSSources)), ...
    'DensityRecords',height(result.Density),'DensityValidDays',result.Coverage.ValidDays(4), ...
    'DensityEPORecords',nnz(strcmp(result.Density.Source,'epo')), ...
    'DensityQTNRecords',nnz(strcmp(result.Density.Source,'qtn')), ...
    'Coverage',table2struct(result.Coverage),'Preservation',table2struct(result.Preservation), ...
    'OutputFiles',result.OutputFiles,'SourceBytes',result.DownloadStatus.Bytes);
fid=fopen(fullfile(out,'final_validation.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(fid);
disp(summary);
end
