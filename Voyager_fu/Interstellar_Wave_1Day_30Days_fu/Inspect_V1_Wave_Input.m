function report = Inspect_V1_Wave_Input
% Inspect original hourly CDF sampling and gaps before spectral processing.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
base = 'Z:/SPART-WORK/Data/Voyager/voyager1/coho/1hr/l2/merged_mag_plasma';
out = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days';
if ~isfolder(out), mkdir(out); end
report = table;
for yy=2012:2021
    files=dir(fullfile(base,num2str(yy),'**','*.cdf'));
    t=[]; b=[];
    for k=1:numel(files)
        q=Voyager_Read_CDF_Product(fullfile(files(k).folder,files(k).name),'coho');
        t=[t;posixtime(q.Epoch)]; %#ok<AGROW>
        b=[b;double(q.BR(:)) double(q.BT(:)) double(q.BN(:))]; %#ok<AGROW>
        if yy==2012 && k==1, disp(q.variable_meta); end
    end
    [t,ix]=sort(t); b=b(ix,:);
    gridTime=posixtime((datetime(yy,1,1,'TimeZone','UTC'):hours(1):datetime(yy+1,1,1,'TimeZone','UTC')-hours(1)).');
    good=false(size(gridTime)); [present,where]=ismember(t,gridTime);
    good(where(present))=all(isfinite(b(present,:)),2);
    edges=diff([false;good;false]); starts=find(edges==1); stops=find(edges==-1)-1;
    lengths=stops-starts+1;
    gapEdges=diff([false;~good;false]); gs=find(gapEdges==1); ge=find(gapEdges==-1)-1;
    gaps=ge-gs+1;
    row=table(yy,numel(t),nnz(good),max(lengths),median(lengths),max(gaps),median(diff(t)), ...
        'VariableNames',{'Year','Records','ValidVectors','MaxContinuousHours','MedianContinuousHours','MaxMissingHours','MedianCadenceSeconds'});
    report=[report;row]; %#ok<AGROW>
    disp(row);
end
writetable(report,fullfile(out,'input_sampling_audit.csv'));
end

