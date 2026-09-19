function [density,audit] = V2_Read_Official_PWS_Density(DataDir)
% 直接读取网站原生密度 CSV；不重算频率、不拟合、不插值。
%% 路径、PDS 标签与校验
base=fullfile(DataDir,'pws','derived','electron_density','native','PDS_release_20260910','data');
file=fullfile(base,'vg2-vlism-density-2019-2025.csv');
labelFile=fullfile(base,'vg2-vlism-density-2019-2025.lblx');
label=fileread(labelFile);
token=regexp(label,'<md5_checksum>([^<]+)</md5_checksum>','tokens','once');
fid=fopen(file,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
digest=java.security.MessageDigest.getInstance('MD5');digest.update(typecast(bytes,'int8'));
md5=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
assert(strcmpi(md5,token{1}),'PDS CSV differs from official MD5.');
options=detectImportOptions(file,'VariableNamingRule','preserve','TextType','string');
options=setvartype(options,options.VariableNames{1},'string');raw=readtable(file,options);
n=regexp(label,'<records>(\d+)</records>','tokens','once');
assert(height(raw)==str2double(n{1}) && width(raw)==22);
assert(strcmp(raw.Properties.VariableNames{12},'N_e (cm^-3)'));
%% 保留官方每一行的原始时刻与已发布电子密度
t=datetime(raw{:,1},'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
density=table(t,raw{:,12},raw{:,15},raw{:,16},lower(string(raw{:,18})),(1:height(raw)).', ...
    'VariableNames',{'EpochUTC','ElectronDensity_cm3','Minimum_cm3','Maximum_cm3','Source','OriginalCSVRow'});
density=sortrows(density,'EpochUTC');
assert(all(isfinite(density.ElectronDensity_cm3)&density.ElectronDensity_cm3>0));
start=regexp(label,'<start_date_time>([^<]+)</start_date_time>','tokens','once');
stop=regexp(label,'<stop_date_time>([^<]+)</stop_date_time>','tokens','once');
assert(min(t)==datetime(start{1},'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC'));
assert(max(t)==datetime(stop{1},'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC'));
% 原始 N_e 与原始频率均保存；即使二者有差异也不自行改写官方 N_e。
audit=struct('File',file,'LabelFile',labelFile,'LabelText',label,'RawTable',raw, ...
    'URL','https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/vg2-vlism-density-2019-2025.csv', ...
    'SHA256',Case1_File_SHA256(file),'MD5',md5,'RecordCount',height(density), ...
    'Method','Official N_e values and native SCET plotted as points. No averaging, frequency-to-density recalculation, interpolation, old-version merging or electron/proton conversion. CSV is the native published density product; other variables still read original CDF.');
end
