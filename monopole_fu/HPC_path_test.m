clear;clc;
global ParentDir 
% global OutputDir
ParentDir = '/gs/home/by2430102/MMS/'; 
DownloadDir = ParentDir;
TempDir = [DownloadDir,'temp/']; mkdir(TempDir);
 
ic = 1:4; iic = 1:4;
load(['/gs/home/by2430102/fwd_matlab_patch/monopole_fu/NameTags.mat']);

% Date = '2019-09-24/2019-09-25';
Date = '2017-01-01T05:16:00.000Z/2017-01-01T06:28:00.000Z';
splitDate = regexp(Date,'/','split');
startDate = splitDate{1};
endDate   = splitDate{2};
startDT = datetime(startDate, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
endDT   = datetime(endDate,   'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
NameTagsDT = datetime(NameTags, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
mask = (NameTagsDT >= startDT) & (NameTagsDT <= endDT);
NameTags = NameTags(mask);

%% output directory & DB init
splitDate{1} = erase(splitDate{1},{':','.'});
splitDate{2} = erase(splitDate{2},{':','.'});
OutputDir = [ParentDir,'Monopole_Search/',splitDate{1},'To',splitDate{2},'/'];
clear mask NameTagsDT
mkdir(OutputDir)
mms.db_init('local_file_db',ParentDir);
%%
msg = [NameTags{1}(2:end-2),'的数据导入2出现问题'];
fid = fopen([OutputDir,'errorlog.txt'], 'a', 'n', 'UTF-8');
if fid == -1
    warning('无法打开 %s 进行写入', [OutputDir,'errorlog.txt']);
else
    fprintf(fid, '%s\n', msg);
    fclose(fid);
end