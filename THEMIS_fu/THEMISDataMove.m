function THEMISDataMove(sourceDir,destDir)
%%
% created by fwd on 2021.12.27
% Different with SDCDataMove.m,this program MOVE files instead of COPY
% example: THEMISDataMove('D:\THEMIS\temp\','D:\THEMIS\')
%%
cd(sourceDir)
filenames0=dir;
j=0;
for ii=1:length(filenames0)
    flag = 0;
    if  ~strcmp(filenames0(ii).name,'.') && ~strcmp(filenames0(ii).name,'..')
        filename1=filenames0(ii).name;
        num=find(filename1=='_');
        
        % classify satellits
        if ~isfolder([destDir,filename1(1:3)])
            mkdir(destDir,filename1(1:3))
        end
        destDir1 = [destDir,filename1(1:3),'/'];
        
        % classify data level
        if ismember(filename1(num(1)+1:num(2)-1), {'l1';'l2';'or'})
            destDir2 = [destDir1,filename1(num(1)+1:num(2)-1),'/'];
            if ~isfolder(destDir2)
                mkdir(destDir2)
            end
        else 
            flag = 1;
            disp(['(+_+)? Unrecognized File:',filename1])
        end
                                   
        % classify instruments
        if ismember(filename1(num(2)+1:num(3)-1), {'ssc'; 'efi'; 'esa'; 'fgm'; 'mom'; 'eff';'scm'})
            destDir3 = [destDir2,filename1(num(2)+1:num(3)-1),'/'];
            if ~isfolder(destDir3)
                mkdir(destDir3)
            end
        else
            flag = 1;
            disp(['(+_+)? Unrecognized File:',filename1])
        end
        
        % classify year
        destDir_THEMIS = [destDir3,filename1(num(3)+1:num(3)+4),'/'];
        if ~isfolder(destDir_THEMIS)
            mkdir(destDir_THEMIS)
        end
        
        % move file
        if flag == 0
            movefile([sourceDir,filename1],destDir_THEMIS,'f');
        end
        j=j+1;
    end
end
disp(['(ง •_•)ง Successfully Classify:', num2str(j),'/',num2str(length(filenames0)-2)])
end