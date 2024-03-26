function cdawebFilesDownload(Date,ic,varargin)
%------written by Wending Fu, Jul.13 2022 in Beijing------------
% Example：cdawebCrawl('2017-08-23/2017-08-24',[1:4],'inst','fgm','drm','brst')
% Input variable
% Date:'2017-08-23/2017-08-23';
% ic:[1:4]; %Satellite id
% inst:'fgm'; %Instrument（optional）
% drm:'srvy'; %mode（optional）
% dpt:''; %description（optional)
% Noticing!! 
% If an instrument does NOT have this descriptor, ALL data for that instrument will not be downloaded !!!
% see also SDCDataMove,SDCPlot,SDCFilesDownload
end
function [Index,Filename] = NameList(ic,Instrument,Mode,Level,Description,Date)
end
function FilesDownload(Index,Filename,OutputDir)
if ~isfolder(OutputFiles_dir)
    mkdir(OutputFiles_dir);
end

url = ['https://cdaweb.gsfc.nasa.gov', Index, '/', Filename];

end