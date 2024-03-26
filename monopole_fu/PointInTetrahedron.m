function ToF = PointInTetrahedron(Loc,R1,R2,R3,R4)
% function ToF = PointInTetrahedron(Loc,R1,R2,R3,R4)
% This function is to determine the Point whether in Tetrahedron or not
% Loc is the location of the Point
% R1,R2,R3,R4 is the location of vertices of tetrahedron, respectively
% The returned value ToF is a logical value of true or faluse
%------written by Wending Fu, Nov.10.2022 in Beijing------------

if length(R1) ==4
    c_eval('R? = R?(:,2:4);')
end
ToF = SameSide(Loc,R1,R2,R3,R4) && SameSide(Loc,R2,R3,R4,R1) ...
    && SameSide(Loc,R3,R4,R1,R2) && SameSide(Loc,R4,R1,R2,R3);
end

% This function is to determine the Point whether at the sameside of four
% faces or not
function ToF = SameSide(Loc,R1,R2,R3,R4)
    normal = cross(R2-R1,R3-R1);
    normal = normal/(sqrt(sum(normal.*normal,2)));
    V4 = R4-R1;V4 = V4/(sqrt(sum(V4.*V4,2)));
    P = Loc-R1;P = P/(sqrt(sum(P.*P,2)));
    dotV4 = dot(normal,V4);
    dotP = dot(normal,P);
    if abs(dotP) > 1e-4
    ToF = sign(dotV4) == sign(dotP);
    else
        ToF = true;
    end
end