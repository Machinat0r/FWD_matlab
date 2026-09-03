function result = Case1_L1_Output_Manifest(dayFile,hourFile,figure4File)
%Case1_L1_Output_Manifest Verify all overwritten image artifacts and hashes.
day = load(dayFile,'result'); hour = load(hourFile,'result');
four = load(figure4File,'result');
dr = day.result.Run.ReportV1; hr = hour.result.Run.ReportV1;
assert(all(strcmp(string(dr.Status),'ok')),'A daily production run failed.');
assert(all(strcmp(string(hr.Status),'ok')),'An hourly production run failed.');
files = [string(dr.FigureFile);string(hr.FigureFile); ...
    string(hr.PeakPADFigureFile);string(four.result.FigureFiles)];
group = [repmat("event_daily",height(dr),1); ...
    repmat("event_hourly",height(hr),1);repmat("event_peak5",height(hr),1); ...
    repmat("Figure4",numel(four.result.FigureFiles),1)];
assert(numel(unique(files)) == numel(files),'Output file overlap.');
sha = strings(size(files)); width = zeros(size(files)); heightPixels = width;
bytes = width; modified = NaT(size(files));
for ii = 1:numel(files)
    assert(isfile(files(ii)),'Expected image missing.');
    info = imfinfo(files(ii)); stat = dir(files(ii));
    assert(strcmpi(info.Format,'png') && info.Width>0 && info.Height>0);
    width(ii) = info.Width; heightPixels(ii) = info.Height;
    bytes(ii) = stat.bytes; modified(ii) = datetime(stat.datenum,'ConvertFrom','datenum');
    sha(ii) = string(Case1_File_SHA256(char(files(ii))));
end
result = struct('CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'DayRunFile',string(dayFile),'HourRunFile',string(hourFile), ...
    'Figure4RunFile',string(figure4File));
result.Images = table(group,files,sha,width,heightPixels,bytes,modified, ...
    'VariableNames',{'Group','ImageFile','SHA256','Width','Height','Bytes','ModifiedLocalTime'});
cfg = Case1_Config;
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','l1_fallback');
result.AuditFile = string(fullfile(folder,'final_output_manifest.mat'));
save(result.AuditFile,'result','-v7.3');
disp(groupsummary(result.Images,'Group'));
end
