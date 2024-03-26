function plotPolyhedron(x,y,z,varargin)
% 
%------written by Wending Fu, May.11.2022 in Beijing------------
FaceColor = varargin{1};
FaceAlpha = varargin{2};
flag = '2D';
if max(x)-min(x) == 0
    Proj = x(1);
    Plane1 = y;
    Plane2 = z;
    sign = 'x';
elseif max(y)-min(y) == 0
    Proj = y(1);
    Plane1 = x;
    Plane2 = z;
    sign = 'y';
elseif max(z)-min(z) == 0
    Proj = z(1);
    Plane1 = x;
    Plane2 = y;
    sign = 'z';
else
    flag = '3D';
end
switch flag
    case '3D'
        Poly3D(x,y,z,FaceColor,FaceAlpha)
    case '2D'
        Poly2D(Plane1,Plane2,Proj,sign,FaceColor,FaceAlpha)
end
end
%% 3D case
function Poly3D(x,y,z,FaceColor,FaceAlpha)
shp=alphaShape(x,y,z);
shp.Alpha = 10*shp.Alpha;
% shp.Alpha=20*max([max(x)-min(x),max(y)-min(y),max(z)-min(z)])/2;
[elements,nodes]=boundaryFacets(shp);
normalvec=zeros(size(elements));
for k=1:size(elements,1)
    normalvec(k,:)=cross(nodes(elements(k,2),:)-nodes(elements(k,1),:),nodes(elements(k,3),:)-nodes(elements(k,1),:));
end
faceNear=nchoosek(1:size(elements,1),2);
hold on
trisurf(elements,nodes(:,1),nodes(:,2),nodes(:,3),'FaceColor',FaceColor,'FaceAlpha',FaceAlpha,'EdgeColor','none');
for k=1:size(faceNear,1)
    isEdge=ismember(elements(faceNear(k,1),:),elements(faceNear(k,2),:));
    if sum(isEdge)>1 && vecnorm(cross(normalvec(faceNear(k,1),:),normalvec(faceNear(k,2),:)),1)/(vecnorm(normalvec(faceNear(k,1),:),1)*vecnorm(normalvec(faceNear(k,2),:),1))>1e-6
%         a = plot3(nodes(elements(faceNear(k,1),isEdge),1),nodes(elements(faceNear(k,1),isEdge),2),nodes(elements(faceNear(k,1),isEdge),3),'Color','m');
%         alpha(a,0.2)
    end
end
hold on
end
%% 2D case
function Poly2D(P1,P2,Proj,sign,FaceColor,FaceAlpha)
shp=alphaShape(P1,P2);
shp.Alpha=10*shp.Alpha;
[elements,nodes]=boundaryFacets(shp);
L = size(elements,1);
for i = 1:size(nodes,1)
    tempelements = [elements(1:L,1:2),ones(L,1)*i];
    try
        elements = [elements;tempelements];
    catch
        elements = tempelements;
    end
end

switch sign
    case 'x'
        trisurf(elements,ones(size(nodes,1),1)*Proj,nodes(:,1),nodes(:,2),'FaceColor',FaceColor,'FaceAlpha',FaceAlpha,'EdgeColor','none');
    case 'y'
        trisurf(elements,nodes(:,1),ones(size(nodes,1),1)*Proj,nodes(:,2),'FaceColor',FaceColor,'FaceAlpha',FaceAlpha,'EdgeColor','none');
    case 'z'
        trisurf(elements,nodes(:,1),nodes(:,2),ones(size(nodes,1),1)*Proj,'FaceColor',FaceColor,'FaceAlpha',FaceAlpha,'EdgeColor','none');
end
hold on
end