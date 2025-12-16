%----------written by Wending Fu, Dec.2025 in Beijing------------
clear;clc;close all
cd '/Volumes/SPART-WORK/Data/Cluster/'
ParentDir = '/Volumes/SPART-WORK/Data/Cluster/';
ic=1:4;

% TT = '2002-03-17\2003-01-01';
TT = '2001-02-07\2008-01-01';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'/']; mkdir(TaskDir)
Datelist = datenum(Datelist,'yyyy-mm-dd');
Datelist = datestr(Datelist(1):Datelist(2),'yyyy-mm-dd');

mms.db_init('local_file_db',ParentDir);
%%
units = irf_units;

% dir path
func_path = mfilename("fullpath");
func_path = strrep(func_path, '\', '/');
func_dir = strfind(func_path,'fwd_matlab_patch') - 1;

%%
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
tint = [R(1,1) R(end,1)];
Tsta = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(1))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];
Tend = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(2))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];
tint=[iso2epoch(Tsta) iso2epoch(Tend)];

% if min(sqrt(R(:,2).^2+R(:,3).^2+R(:,4).^2)) > 7*units.RE/1e3
%     continue
% end

try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
catch
    writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],[TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end

try
%% load R & B
% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
c_eval('B? = irf_abs(irf_gse2gsm(B?_gse));',ic);
Rc_gse = c_caa_var_get('sc_r_xyz_gse__CL_SP_AUX','mat');
c_eval("dR? = c_caa_var_get('sc_dr?_xyz_gse__CL_SP_AUX','mat');",ic);
c_eval('R?_gse = [Rc_gse(:,1) Rc_gse(:,2:4) + dR?(:,2:4)];',ic);
c_eval('R? = irf_gse2gsm(R?_gse);',ic);
% c_eval('R? = irf_resamp(R?,B1);',ic)
c_eval('R? = irf_resamp(R?,B1);',ic);
c_eval('R? = irf_abs(R?);',ic);
%%
c_eval('B? = irf_resamp(B?,B1);',2:4);
noise = estimate_noise_strength_4sc(B1(:,2:4), B2(:,2:4), B3(:,2:4), B4(:,2:4));
save([TaskDir, num2str(tempDate),'.mat'], 'noise')

c_eval('B?(find(R?(:,5)<10*units.RE/1e3),:)=[];',ic);
c_eval('[Bmax?, id?] = max(B?(:,5));',ic);

if Bmax1 > 80 || Bmax2 > 80 || Bmax3 > 80 || Bmax4 > 80
    writematrix(['Strong B at: ',datestr(datenum(1970,1,1,0,0,0)+B1(id1,1)/86400,'yyyymmdd HH:MM:SS.FFF'), 'maxB = ', num2str(Bmax1)],...
                    [TaskDir,'StrongB.txt'],'WriteMode','append','Encoding','UTF-8')
end

catch
        writematrix(['Data reading failure at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
end
end