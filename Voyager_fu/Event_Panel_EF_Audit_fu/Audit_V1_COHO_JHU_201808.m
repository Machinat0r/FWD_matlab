function report = Audit_V1_COHO_JHU_201808
%Audit_V1_COHO_JHU_201808 Compare original hourly values without reprocessing.
% Reads the public JHU ASCII and COHO CDF. No averages or science filters.
addpath(genpath('C:/Users/Administrator/Documents/irfu-matlab-master'));
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Case1_PPT_VerticalLine_Events_7d_fu');
root = 'Z:/SPART-WORK/Data/Voyager';
asciiFile = fullfile(root,'voyager1/lecp/1h/calibrated_ascii/jhuapl_proton_flux/2018/v1_2018_prot_flux_1h.txt');
cdfFile = fullfile(root,'voyager1/coho/1hr/l2/merged_mag_plasma/2018/08/voyager1_coho1hr_merged_mag_plasma_20180801_v01.cdf');
out = fullfile(root,'voyager1/coho/1hr/derived/source_audit/2018/panel_e_20180823');
if ~isfolder(out), mkdir(out); end
fidIn = fopen(asciiFile,'r'); assert(fidIn>=0);
values = textscan(fidIn,repmat('%f',1,11),'HeaderLines',2,'CollectOutput',true);
fclose(fidIn); x = values{1};
assert(size(x,2) == 11 && all(isfinite(x(:,1:3)),'all'));
% The published 2018-named file also contains records from 2019-01-01.
% Retain them and match dates using the actual year/doy/hour columns.
t = datetime(x(:,1),ones(size(x,1),1),ones(size(x,1),1),'TimeZone','UTC') + days(x(:,2)-1) + hours(x(:,3));
assert(numel(unique(t)) == numel(t));
c = Voyager_Read_CDF_Product(cdfFile,'coho');
[found,index] = ismember(c.Epoch,t);
comparison = table((1:numel(c.Epoch))',c.Epoch,found,'VariableNames',{'COHO_CDFRecord','EpochUTC','JHU_RecordPresent'});
comparison.JHU_TextLine = nan(numel(c.Epoch),1);
comparison.JHU_TextLine(found) = index(found)+2;
for k = 1:3
    name = sprintf('protonFlux%d_LECP',k);
    val = nan(numel(c.Epoch),1); sigma = val;
    val(found) = x(index(found),4+2*k);
    sigma(found) = x(index(found),5+2*k);
    comparison.(sprintf('COHO_P%d',k)) = c.(name);
    comparison.(sprintf('JHU_P%d',k)) = val;
    comparison.(sprintf('JHU_P%d_Sigma',k)) = sigma;
    comparison.(sprintf('P%d_EqualAfterCDFSingleStorage',k)) = ...
        c.(name) == double(single(val));
    comparison.(sprintf('P%d_AbsoluteDifference',k)) = abs(c.(name)-val);
end
first = datetime(2018,8,20,'TimeZone','UTC');
last = datetime(2018,8,27,'TimeZone','UTC');
event = comparison(comparison.EpochUTC>=first & comparison.EpochUTC<last,:);
report = struct('ASCIIFile',asciiFile,'CDFFile',cdfFile,'ASCIIAnnualRecords',size(x,1), ...
    'ASCIIYears',unique(x(:,1))','CDFAugustRecords',height(comparison),'AugustJHURecordsPresent',nnz(found), ...
    'EventCDFRecords',height(event),'EventJHURecordsPresent',nnz(event.JHU_RecordPresent));
for k = 1:3
    values = comparison.(sprintf('COHO_P%d',k));
    equal = comparison.(sprintf('P%d_EqualAfterCDFSingleStorage',k));
    valid = isfinite(values);
    eventValues = event.(sprintf('COHO_P%d',k));
    eventEqual = event.(sprintf('P%d_EqualAfterCDFSingleStorage',k));
    report.(sprintf('AugustP%dValid',k)) = nnz(valid);
    report.(sprintf('AugustP%dExactMatches',k)) = nnz(valid & equal);
    report.(sprintf('EventP%dValid',k)) = nnz(isfinite(eventValues));
    report.(sprintf('EventP%dExactMatches',k)) = nnz(isfinite(eventValues) & eventEqual);
    report.(sprintf('AugustP%dCDFValueWithoutASCII',k)) = nnz(valid & ~found);
end
writeCSV(comparison,fullfile(out,'COHO_vs_JHU_201808_pointwise.csv'));
writeCSV(event,fullfile(out,'COHO_vs_JHU_event_20180820_26.csv'));
save(fullfile(out,'COHO_vs_JHU_201808_pointwise.mat'),'comparison','event','report','asciiFile','cdfFile');
fid = fopen(fullfile(out,'comparison_summary.json'),'w','n','UTF-8');
assert(fid>=0); cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
disp(report);
disp(event(ismember(event.EpochUTC,[datetime(2018,8,23,10,0,0,'TimeZone','UTC');datetime(2018,8,24,8,0,0,'TimeZone','UTC')]),:));
end
function writeCSV(t,file)
t.EpochUTC.Format = 'yyyy-MM-dd''T''HH:mm:ss''Z''';
t.EpochUTC = string(t.EpochUTC);
writetable(t,file,'Encoding','UTF-8');
end


