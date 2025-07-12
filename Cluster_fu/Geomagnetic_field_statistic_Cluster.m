%------written by Wending Fu, Jul.2025 in Beijing------------
clear;clc;
cd /Volumes/SPART-NAS/Data/Cluster/
ParentDir = '/Volumes/SPART-NAS/Data/Cluster/';
ic=1:4;

% TT = '2002-03-17\2003-01-01';
TT = '2002-01-01\2002-01-02';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'\']; mkdir(TaskDir)
Datelist = datenum(Datelist,'yyyy-mm-dd');
Datelist = datestr(Datelist(1):Datelist(2),'yyyy-mm-dd');

mms.db_init('local_file_db',ParentDir);
%%
units = irf_units;
% parfor_progress(length(NameTags)-1);
for tempDate = 1:size(Datelist,1)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;clear B1 B2 B3 B4 R1 R2 R3 R4;

Tsta = [Datelist(tempDate,:) 'T00:00:000Z'];
Tend = [Datelist(tempDate+1,:) 'T00:00:000Z'];
tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH
disp(['当前日期:' Datelist(tempDate,:)])

try 
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
catch
    caa_download(tint,'CL_SP_AUX');
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
end

R = c_caa_var_get('sc_r_xyz_gse__CL_SP_AUX','mat');
TT = irf_time(R(:,1),'epoch>utc');
d = datetime(TT, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS''Z''', ...
    'TimeZone','UTC');
YY = year(d); MM = month(d); DD = day(d);
hh = hour(d); mm = minute(d); ss = floor(second(d)); 

data = [YY, MM, DD, hh, mm, ss, R(:,2:4)./units.RE*1e3,ones(1440,1), -400*ones(1440,1),zeros(1440,1),zeros(1440,1)];
save('/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/Cluster_fu/R.mat','data')

tint = [R(1,1) R(end,1)];
Tsta = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(1))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];
Tend = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(2))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];

% if min(sqrt(R(:,2).^2+R(:,3).^2+R(:,4).^2)) < 20000+6371
%     continue
% end
errorflag = 1;
try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
catch
    try
    %    Magnetic fields
    c_eval("caa_download(tint,'C*_CP_FGM_FULL')",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    catch
        errorflag = 666;
    end
end

try
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
catch 
    try
    c_eval("caa_download(tint,'C*_CP_AUX_POSGSE_1M')",ic);  % position & velocity for each sc
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
    catch
    errorflag = 666;
    end
end
  
if errorflag == 666
    writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],[TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
end

try
% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
c_eval('B? = irf_gse2gsm(B?_gse);',ic);
% c_eval('B?=irf_abs(B?_gsm);',ic);
% c_eval('B? = irf_resamp(B?,B1);',2:4);
c_eval('B? = irf_resamp(B?,B1);',2:4);
catch
        writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
end
end