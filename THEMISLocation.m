%read THEMIS position data and save into mRth.mat file
function [varargout] = THEMISLocation(thIds,years,dataDir)

% dataDir = 'D:\THEMIS\'; 
% dataDir = '/data/themis';
% thIds = 'abcde';

for Id=1:length(thIds)
  thId=thIds(Id);
  R = [];
  for year=years(1:end)
    fullPath = sprintf('%sth%s%sor%sssc%s%02d',...
      dataDir,thId,filesep,filesep,filesep,year);
    if ~exist(fullPath,'dir'), continue, end
    files = dir(sprintf('%s%sth%s_or_ssc_*_v*.cdf',fullPath,filesep,thId));
   
    if ~isempty(files)
      for iFile=1:length(files)
        fileToRead = [fullPath,filesep,files(iFile).name];
        fprintf('reading %s\n',fileToRead)
        tmpData = spdfcdfread(fileToRead,'CombineRecords',true,'Variable','XYZ_GSM');
%         tmpData = tmpData*6371.2; % comvert to kilometers
        tmpEpoch = spdfcdfread(fileToRead,'CombineRecords',true,...
          'KeepEpochAsIs',true,'Variable','Epoch');
        tmpEpoch = irf_time(tmpEpoch,'cdfepoch>epoch');
        R = [R; tmpEpoch tmpData];
        clear tmpData tmpEpoch
      end
    end
  end
  % remove repeating points at month boundary
  ii = diff(R(:,1))==0; R(ii,:) = [];
  eval(['Rth' thId '=R;'])
  eval(['varargout{Id} = Rth' thId ';'])
end
end