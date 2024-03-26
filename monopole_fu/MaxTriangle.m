function Loc = MaxTriangle(LocX,LocY,LocZ)
%------written by Wending Fu, May 16.2022 in Beijing------------
id = nchoosek(1:length(LocX),3);
Area = ones(size(id,1),1);
for i = 1:size(id,1)
    Area(i) = TriangleArea(LocX(id(i,:)),LocY(id(i,:)),LocZ(id(i,:)));
end
idx = find(Area == max(Area));
LocX = LocX(id(idx,:));
LocY = LocY(id(idx,:));
LocZ = LocZ(id(idx,:));
Loc = {LocX,LocY,LocZ};
end
function Area = TriangleArea(X,Y,Z)
    a = sqrt((X(2)-X(1))^2+(Y(2)-Y(1))^2+(Z(2)-Z(1))^2);
    b = sqrt((X(3)-X(2))^2+(Y(3)-Y(2))^2+(Z(3)-Z(2))^2);
    c = sqrt((X(1)-X(3))^2+(Y(1)-Y(3))^2+(Z(1)-Z(3))^2);
    p = (a+b+c)/2;
    Area = sqrt(p*(p-a)*(p-b)*(p-c));
end