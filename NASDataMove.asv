%------written by Wending Fu, May.2021 in Beijing------------
% see also SDCDataMove
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DataDir = '/Volumes/172.17.190.41/Data/Du/MMS_fgm_srvy/2015-09-01To2022-10-20'; %'/Volumes/172.17.190.41/Data/Du/MMS_fgm_srvy';
CDFList = dir([DataDir,'/**/*.cdf']);
% CDFList = regexp(ls([DataDir,'/**/*.cdf']), '\s+', 'split');
i = 0;
while i <= length(CDFList)
    SDCDataMove([CDFList(i).folder,'/'],ParentDir)
    i = i+1;
end