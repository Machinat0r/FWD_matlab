function [index,varargout]=monopole_index_wrong(B1,B2,B3,B4)
%%Example:
%index=monopole_index(B1,B2,B3,B4);
%B? is in any coordinates, 1st colomn is time, 2nd to 4th colomn is Bx,By,Bz
%
%------written by Wending Fu, Mar.09.2022 in Beijing------------
%% Check TSeries
idB = {B1,B2,B3,B4};

if isa(B1,'TSeries') 
for i=1:4
if isa(idB{i},'TSeries')
    idB{i} =  [idB{i}.time.epochUnix double(idB{i}.data)];
end
end
end
%% Check inputs
nCol=size(idB{1,1},2);
nRow=size(idB{1,1},1);
if nCol<3
	error('Time tag must be included in each input vector. Please do so and try again.')
elseif nCol==3 
    if size(idB{1,2},1)==nRow && size(idB{1,3},1)==nRow && size(idB{1,4},1)==nRow
	isTimeSpecified = false;
    else
    error('The vectors needs to be the same size or have the time vector in them. See usage: help c_4_poincare_index') 
    end
else
	isTimeSpecified = true;
end
if nargin==0
	help c_4_poincare_index;
	return;
elseif nargin < 4
	error('Too few input values. See usage: help c_4_poincare_index')
elseif nargin>4
	error('Too many input values. See usage: help c_4_poincare_index')
end

%% Treat the case when the time column is specified
if isTimeSpecified
	%Each vector contains time tags so all vectors needs to be resampled to
	%establish synchronisation between the S/C's
	idB{1,2}=irf_resamp(idB{1,2},idB{1,1});
	idB{1,3}=irf_resamp(idB{1,3},idB{1,1});
	idB{1,4}=irf_resamp(idB{1,4},idB{1,1});
	
	%Renaming of vectors for simplicty sake and removing time for calculations
	%using B?
	time=idB{1,1}(:,1); %time= first column of SC1
	idB{1,1}=idB{1,1}(:,2:4); %vec1= 2 to 4th column (Bx,By,Bz) for s/c 1 (if that's the one you placed in position 1)
	idB{1,2}=idB{1,2}(:,2:4); %vec2= 2 to 4th column (Bx,By,Bz)
	idB{1,3}=idB{1,3}(:,2:4); %vec3= 2 to 4th column (Bx,By,Bz)
	idB{1,4}=idB{1,4}(:,2:4); %vec4= 2 to 4th column (Bx,By,Bz)
end

%%
lx=size(idB{1,1},1); %lx becomes the number of rows in B1

%Create the zero matrix that will be used to map each point from xyz to
%By,Bx,Bz space

Map_sc1= zeros(lx,3);
Map_sc2= zeros(lx,3);
Map_sc3= zeros(lx,3);
Map_sc4= zeros(lx,3);
% map the points from x,y,z to magnetic three-dimensional field space (Bx,By,Bz instead of
% x,y,z) (ref. Greene, J.M. 1990)
% A null point in configuration space (xyz satellites) corresponds to the
% origin in M space.

%Mapping of the B fields for each s/c but as a unit vector (length=1)
%Calculate the length of each vector
vl1= sqrt(dot(idB{1,1},idB{1,1},2)); %norm (length of vector)=sqrt(dot(vec1,vec1,2)).
%dot(vec1,vec1,2) treats each row as a vector in the matrix so A=dot(vec1,vec1,2)
% would give A(1,:) dot product of vec1(1,:) and vec1(1,:)
vl2= sqrt(dot(idB{1,2},idB{1,2},2));
vl3= sqrt(dot(idB{1,3},idB{1,3},2));
vl4= sqrt(dot(idB{1,4},idB{1,4},2));

%% Unit vectors is used in solid angle so we need divide the vector with
%its norm (magnitude of the vector/length)
for i=1:3   %Each column on Map_sc1 is given unit vector (first column in 3D divided by the length of vector in the
	% direction of each s/c)
	Map_sc1(:,i)=idB{1,1}(:,i)./vl1; %Needs to use the for loop so that matrix dimensions agrees
	Map_sc2(:,i)=idB{1,2}(:,i)./vl2;
	Map_sc3(:,i)=idB{1,3}(:,i)./vl3;
	Map_sc4(:,i)=idB{1,4}(:,i)./vl4;
end

%% Calculate the solid angles
c_eval('flag? = zeros(length(Map_sc?),1);',1:4);
for i=1:4
    ids = 1:4;
    ids(i)=[];
    c_eval(['V = abs(dot(cross(Map_sc',num2str(ids(1)),',Map_sc',num2str(ids(2)),',2),Map_sc',num2str(ids(3)),',2));'])
%     c_eval(['V11 = abs(dot(cross(Map_sc',num2str(ids(1)),',Map_sc',num2str(ids(2)),',2),Map_sc',num2str(i),',2));'])
%     c_eval(['V12 = abs(dot(cross(Map_sc',num2str(ids(1)),',Map_sc',num2str(ids(3)),',2),Map_sc',num2str(i),',2));'])
%     c_eval(['V13 = abs(dot(cross(Map_sc',num2str(ids(2)),',Map_sc',num2str(ids(3)),',2),Map_sc',num2str(i),',2));'])
    
    c_eval('temp_Map_sc? = -Map_sc?;',ids)
    c_eval(['V1 = abs(dot(cross(temp_Map_sc',num2str(ids(1)),',temp_Map_sc',num2str(ids(2)),',2),Map_sc',num2str(i),',2));'])
    c_eval(['V2 = abs(dot(cross(temp_Map_sc',num2str(ids(1)),',temp_Map_sc',num2str(ids(3)),',2),Map_sc',num2str(i),',2));'])
    c_eval(['V3 = abs(dot(cross(temp_Map_sc',num2str(ids(2)),',temp_Map_sc',num2str(ids(3)),',2),Map_sc',num2str(i),',2));'])

    c_eval(['flag',num2str(i),'(V1<V & V2<V & V3<V)=1;']);
end
index = zeros(length(flag1),1);
index(flag1+flag2+flag3+flag4==4) = 1;

if isTimeSpecified
	index=[time index];
end
Identification(mfilename('fullpath'));

c_eval('varargout{?} = flag?;')
end