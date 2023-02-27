clear
clc
path='C:\Users\fwd\Desktop\';% to acess h5 file path
File=dir('C:\Users\fwd\Desktop\*.h5');  % to obtain all of *.h5 file
allname={File.name}';        % to save h5 filename into an array
file_count=size(allname,1);  %to obtain number of *.h5 file 
dataX=zeros(file_count,1);   %  to build an array for saving PSD_X
dataY=zeros(file_count,1);   % to build an array for saving PSD_Y
dataZ=zeros(file_count,1);   % to build an array for saving PSD_Z
dataLogLat=zeros(file_count,1);   % to build an array for saving LonLat
dataLogFreULF=zeros(file_count,1);   % to build an array for saving LonLat
for k=1:file_count
    file_path=strcat(path,allname(k,1));    %to read each h5 file name
    dataa{k,1}=h5read(file_path{1,1},'/PSD_X');   % to use subset by h5read  function
dataa{k,2}=h5read(file_path{1,1},'/PSD_Y');   % to use subset by h5read  function
dataa{k,3}=h5read(file_path{1,1},'/PSD_Z');   % to use subset by h5read  function
dataa{k,4}=h5read(file_path{1,1},'/LonLat');   % to use subset by h5read  function
dataa{k,5}=h5read(file_path{1,1},'/FreULF');   % to use subset by h5read  function
end