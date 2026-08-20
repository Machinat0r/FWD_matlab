% Download Voyager 1 LECP hourly sectored P1 flux for event 2020-E01.
% NASA CDAWeb is accessed programmatically; no browser, interpolation, or
% additional averaging is used.

programRoot = fileparts(mfilename('fullpath'));
monthlyRoot = fullfile(programRoot, 'Voyager_Interstellar_Monthly');
downloadRoot = fullfile(monthlyRoot, 'Voyager1_LECP_Sectored_Hourly');
databaseRoot = fullfile(monthlyRoot, 'Voyager1_Selected_Event_Data', ...
    'voyager1', 'lecp', '1hr', 'l2', 'sectored_flux', '2020', '06-09');

cdfName = 'voyager1_lecp_p1_sectored_hourly_20200622_20200908.cdf';
convenienceFile = fullfile(downloadRoot, cdfName);
databaseFile = fullfile(databaseRoot, cdfName);
metadataFile = fullfile(downloadRoot, ...
    'voyager1_2020e01_lecp_hourly_metadata.json');

if ~isfolder(downloadRoot), mkdir(downloadRoot); end
if ~isfolder(databaseRoot), mkdir(databaseRoot); end

dataset = 'VOYAGER-1_LECP_LEV-2-HOURLY-AVG';
startUTC = '20200622T000000Z';
stopUTC = '20200909T000000Z';
variables = ['FHDU_SectoredFluxes,', ...
    'FHDU_SectoredFluxUncertainties'];
requestURL = sprintf([ ...
    'https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/', ...
    'sp_phys/datasets/%s/data/%s,%s/%s?format=cdf'], ...
    dataset, startUTC, stopUTC, variables);

options = weboptions('Timeout', 180, 'ContentType', 'text');
xmlText = webread(requestURL, options);
if isstring(xmlText), xmlText = char(xmlText); end
urlToken = regexp(xmlText, ...
    '<(?:[A-Za-z0-9_]+:)?Name>([^<]+)</(?:[A-Za-z0-9_]+:)?Name>', ...
    'tokens', 'once');
if isempty(urlToken)
    error('VoyagerDownload:CDAWebResponse', ...
        'CDAWeb returned no generated CDF URL.');
end
cdfURL = strrep(urlToken{1}, '&amp;', '&');

websave(convenienceFile, cdfURL, weboptions('Timeout', 300));
fileInfo = dir(convenienceFile);
if isempty(fileInfo) || fileInfo.bytes < 1000
    error('VoyagerDownload:SmallCDF', ...
        'Downloaded CDF is missing or unexpectedly small.');
end
copyfile(convenienceFile, databaseFile, 'f');

metadata = struct;
metadata.dataset = dataset;
metadata.time_interval = sprintf('%s/%s', startUTC, stopUTC);
metadata.variables = split(string(variables), ',');
metadata.request = requestURL;
metadata.generated_cdf = cdfURL;
metadata.convenience_file = convenienceFile;
metadata.database_file = databaseFile;
metadata.bytes = fileInfo.bytes;
metadata.processing = [ ...
    'NASA hourly averages; no additional averaging or interpolation'];
fid = fopen(metadataFile, 'w');
assert(fid >= 0, 'Unable to create metadata file: %s', metadataFile);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(metadata, 'PrettyPrint', true), 'char');
clear cleanup

fprintf('Downloaded %d bytes\n%s\n', fileInfo.bytes, convenienceFile);
