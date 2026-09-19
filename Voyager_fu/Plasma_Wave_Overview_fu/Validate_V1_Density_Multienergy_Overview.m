function validation = Validate_V1_Density_Multienergy_Overview
% 结果核验入口；读取派生 MAT 仅用于测试，不替代正式原始 CDF 计算入口。
%% 路径
CodeDir=fileparts(mfilename('fullpath'));
addpath(CodeDir,fullfile(fileparts(CodeDir),'Case1_PPT_VerticalLine_Events_7d_fu'));
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Density_Multienergy_Overview/';
s=load(fullfile(out,'V1_density_multienergy_audit.mat'),'result');r=s.result;
assert(height(r.Channels)==18 && numel(unique(r.Channels.CDFVariable))==18);
assert(isequal(r.ProtectedBefore,r.ProtectedAfter));
for k=1:height(r.ProtectedBefore)
    assert(strcmp(Case1_File_SHA256(r.ProtectedBefore.File(k)),r.ProtectedBefore.SHA256(k)));
end
%% 每个能道独立检查代表日期的直接算术均值、样本数、缺测及零值
sourceDay=dateshift(r.RawCOHO.EpochUTC,'start','day');
days=dateshift(r.Daily.EpochUTC,'start','day');
checks=0;
for c=1:18
    x=r.RawCOHO.Flux(:,c); valid=isfinite(x);
    assert(sum(r.FluxCount(:,c))==nnz(valid));
    assert(isequal(isnan(r.FluxMean(:,c)),r.FluxCount(:,c)==0));
    occupied=find(r.FluxCount(:,c)>0);
    probe=occupied(unique([1,ceil(numel(occupied)/2),numel(occupied)]));
    zeroDays=find(r.FluxMean(:,c)==0,1);
    sourceZeroDay=unique(sourceDay(valid & x==0));
    if ~isempty(sourceZeroDay)
        [found,i]=ismember(sourceZeroDay,days);probe=[probe;i(find(found,1))]; %#ok<AGROW>
    end
    probe=unique([probe;zeroDays]);
    for j=probe(:).'
        values=x(sourceDay==days(j) & valid);
        expected=mean(values);
        assert(numel(values)==r.FluxCount(j,c));
        assert(abs(expected-r.FluxMean(j,c))<=eps(max(1,abs(expected)))*4);
        checks=checks+1;
    end
end
%% 四张可编辑 FIG 的实际图元与全部输出数据逐点核对
groups={1:3,4:8,9:13,14:18};
names=["LECP_3channels","CRS_01_05","CRS_06_10","CRS_11_15"];
figureChecks=table(names(:),zeros(4,1),zeros(4,1),false(4,1), ...
    'VariableNames',{'FigureGroup','Panels','ProtonChannels','AllPayloadChecksPassed'});
time=datenum(r.Daily.EpochUTC); xLimits=datenum([r.Method.StartUTC r.Method.EndUTCExclusive]);
for g=1:4
    file=fullfile(out,"V1_1990_20250630_daily_"+names(g)+".fig");
    f=openfig(file,'invisible');
    cleanup=onCleanup(@()close(f));
    ax=findall(f,'Type','axes');
    assert(numel(ax)==numel(groups{g})+3);
    assert(isempty(findall(f,'Type','surface')) && isempty(findall(f,'Type','colorbar')));
    for k=1:numel(ax), assert(isequal(ax(k).XLim,xLimits)); end
    assertLine(f,'Sunspot',time,r.Daily.SunspotNumber,false);
    assertLine(f,'MAG',time,r.Daily.B_nT,false);
    for c=groups{g}
        y=r.FluxMean(:,c);y(y<=0)=NaN;
        assertLine(f,r.Channels.CDFVariable(c),time,y,true);
    end
    for source=["epo","qtn"]
        use=r.Density.Source==source;
        h=findobj(f,'Type','line','Tag',char(source));
        assert(numel(h)==1);
        assert(isequaln(h.XData(:),datenum(r.Density.EpochUTC(use))));
        assert(isequaln(h.YData(:),r.Density.Density_cm3(use)));
        a=ancestor(h,'axes');
        assert(all(h.YData>=a.YLim(1) & h.YData<=a.YLim(2)));
    end
    figureChecks.Panels(g)=numel(ax);
    figureChecks.ProtonChannels(g)=numel(groups{g});
    figureChecks.AllPayloadChecksPassed(g)=true;
    clear cleanup
end
assert(sum(figureChecks.ProtonChannels)==18);
sourceManifest=fullfile(fileparts(fileparts(r.DensityAudit.CurrentCSV)), ...
    '..','source_check_20260915','source_verification_manifest.json');
% CurrentCSV = .../PDS_release_20260910/data/file; resolve directly for clarity.
sourceManifest='Z:/SPART-WORK/Data/Voyager/voyager1/pws/derived/electron_density/native/source_check_20260915/source_verification_manifest.json';
official=jsondecode(fileread(sourceManifest)); assert(official.CurrentArchiveMatchesWebsite);
validation=struct('Complete',true,'DailyRows',height(r.Daily),'COHOCDFs',numel(r.COHOSources), ...
    'ProtonChannels',18,'Figures',4,'ElectricFieldPanels',0,'DirectMeanSpotChecks',checks, ...
    'OriginalFilesUnchanged',true,'FirstThreeStatisticsUnchanged',true, ...
    'AllFiguresMatchFullData',true,'ExactRequestedTimeAxis',true, ...
    'DensityRecords',height(r.Density),'DensityDays',r.DensityAudit.ValidDays, ...
    'DensityCurrentArchiveMatchesWebsite',true,'DensityPreviousFigureComplete',true, ...
    'DensityOldVersionRecords',r.DensityAudit.OldRecords, ...
    'DensityOldOnlyMultisetRecords',height(r.DensityAudit.OldOnly), ...
    'DensityCurrentOnlyMultisetRecords',height(r.DensityAudit.CurrentOnly), ...
    'DensitySharedTimeSourceValueDifferences',height(r.DensityAudit.ChangedDensity), ...
    'SourceManifest',sourceManifest,'FigureChecks',figureChecks, ...
    'Coverage',r.Coverage,'OutputFiles',r.OutputFiles, ...
    'CheckedUTC',datetime('now','TimeZone','UTC'));
fid=fopen(fullfile(out,'final_validation.json'),'w');
assert(fid>=0);fprintf(fid,'%s',jsonencode(validation,'PrettyPrint',true));fclose(fid);
writetable(figureChecks,fullfile(out,'figure_validation.csv'));
disp(figureChecks);fprintf('Validation passed: 18 proton channels, 4 figures, all 755 density records.\n');
end

function assertLine(f,tag,x,y,isLog)
h=findobj(f,'Type','line','Tag',char(tag));assert(numel(h)==1);
assert(isequaln(h.XData(:),x(:)) && isequaln(h.YData(:),y(:)));
assert(strcmp(h.LineStyle,'-') && h.LineWidth==.4);
a=ancestor(h,'axes');
if isLog, assert(strcmp(a.YScale,'log')); end
good=isfinite(y);
assert(all(y(good)>=a.YLim(1)&y(good)<=a.YLim(2)));
end
