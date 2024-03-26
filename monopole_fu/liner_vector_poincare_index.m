function [PI,varargout]=liner_vector_poincare_index(vec1,vec2,vec3,vec4)
%Example:
%indices=liner_vector_poincare_index(B1,B2,B3,B4);
%B? is in any coordinates, 1st colomn is time, 2nd to 4th colomn is Bx,By,Bz
%
%Citation: DOI: 10.1109/VISUAL.2004.107
%see also c_fgm_poincare_index,c_4_poincare_index
%------written by Wending Fu, Mar.09.2022 in Beijing------------

n=size(vec1,2);
vec = cell(4,1);
if n>3
    c_eval('vec? = irf_resamp(vec?,vec1);',2:4)
    time=vec1(:,1);
    c_eval('vec{?} = vec?(:,2:4);',1:4)
else
    c_eval('vec{?} = vec?;',1:4)
end
c_eval('v? = sqrt(dot(vec{?},vec{?},2));',1:4)
c_eval('vec{?} = vec{?}./v?;',1:4)

for i = 1:4
    ic = 1:4;ic(i)=[];ic = [ic ic(1)];
    for ii = 1:3
    c_eval(['alpha',num2str(ii),' = acos(dot(vec{',num2str(ic(ii)),'},vec{',num2str(ic(ii+1)),'},2));'])
    end
    s = 0.5*(alpha1+alpha2+alpha3);
    eval(['area' num2str(i) '=4*atan(sqrt(tan(s/2).*tan((s-alpha1)/2).*tan((s-alpha2)/2).*tan((s-alpha3)/2)));'])
end
PI = (area1+area2+area3+area4)/4/pi;
%add time tags
if n>3
    PI=[time,PI];
else
    PI=PI;
end
c_eval('varargout{?} = area?;')
Identification(mfilename('fullpath'));
end