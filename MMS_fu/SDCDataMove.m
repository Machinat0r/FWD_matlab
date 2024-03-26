function SDCDataMove(sourceDir,destDir)
%------written by Wending Fu, May.2021 in Beijing------------
% move data downloaded from MMS website
% script description:
% sourceDir: the path of the old directory you restored the datas
% destDir: the path of the directory you want to restore the datas
% sourceDir='D:\MMS\20170611\';
% destDir='D:\MMS\';
% the above two should be changed according to your special case
% mms1_edp_brst_l2_dce_20151025110614_v0.1.0
% mms1_fgm_brst_l2_20151025110614_v4.18.0
% mms2_fpi_brst_l2_des-moms_20151025110614_v2.1.0
% mms1_edp_brst_l2_hmfe_20151025110614_v0.0.8
% tha_l0_454_20080106.pkt
% tha_l1_mom_20090305_v01.cdf
% tha_l1_spin_20080218_v01.cdf
% tha_l2_fbk_20080214_v01.cdf
% see also SDCFilenames,SDCPlot,SDCFilesDownload

cd(sourceDir)
filenames0=dir;
j=0;
for ii=1:length(filenames0)
% for ii=2850:2870
    if filenames0(ii).name(1) ~= '.'
    % if  ~strcmp(filenames0(ii).name,'.') && ~strcmp(filenames0(ii).name,'..') &&...
    %         ~strcmp(filenames0(ii).name,'.DS_Store')
        filename1=filenames0(ii).name;
      
                 num=find(filename1=='_');
                 if ~isfolder([destDir,filename1(1:4)])
                    mkdir(destDir,filename1(1:4))
                 end
           if  filename1(num(2)+1:num(3)-1)=='brst'
                 destDir1=[destDir,filename1(1:4),'/'];     
                  for numi=1:length(num)-2
                         destDir_mms=[destDir1,filename1(num(numi)+1:num(numi+1)-1),'/'];
                         if ~isfolder([destDir1,filename1(num(numi)+1:num(numi+1)-1)])
                            mkdir(destDir1,filename1(num(numi)+1:num(numi+1)-1))
                         end
                         destDir1=destDir_mms;                        
                  end
                  datefile1=filename1(num(end-1)+1:num(end)-1);
%                   destDir_mms=[destDir_mms,datefile1(1:4),'\',datefile1(5:6),'\',datefile1(7:8),'\'];
                  if ~isfolder([destDir1,datefile1(1:4)])
                      mkdir(destDir1,datefile1(1:4))
                  end
                  destDir_mms=[destDir_mms,datefile1(1:4),'/'];
                  if ~isfolder([destDir_mms,datefile1(5:6)])
                      mkdir(destDir_mms,datefile1(5:6))
                  end
                  destDir_mms=[destDir_mms,datefile1(5:6),'/'];
                  if ~isfolder([destDir_mms,datefile1(7:8)])
                      mkdir(destDir_mms,datefile1(7:8))
                  end
                  destDir_mms=[destDir_mms,datefile1(7:8),'/'];          
            else 
                 destDir1=[destDir,filename1(1:4),'/'];
                  for numi=1:length(num)-2
                         destDir_mms=[destDir1,filename1(num(numi)+1:num(numi+1)-1),'/'];
                         if ~isfolder([destDir1,filename1(num(numi)+1:num(numi+1)-1)])
                            mkdir(destDir1,filename1(num(numi)+1:num(numi+1)-1))
                         end
                         destDir1=destDir_mms;
                  end
                  datefile1=filename1(num(end-1)+1:num(end)-1);
%                  destDir_mms=[destDir_mms,datefile1(1:4),'\',datefile1(5:6),'\'];
                  if ~isfolder([destDir1,datefile1(1:4)])
                      mkdir(destDir1,datefile1(1:4))
                  end
                  destDir_mms=[destDir_mms,datefile1(1:4),'/'];
                  if ~isfolder([destDir_mms,datefile1(5:6)])
                      mkdir(destDir_mms,datefile1(5:6))
                  end
                  destDir_mms=[destDir_mms,datefile1(5:6),'/'];
            end
         movefile([sourceDir,filename1],destDir_mms,'f');
         j=j+1;
    end
end
Identification(mfilename('fullpath'));
end