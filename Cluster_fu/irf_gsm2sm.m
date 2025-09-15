function [y]=irf_gsm2sm(x,flag_gsm2sm)
%------written by Wending Fu, Sept.2025 in Beijing------------
%function [y]=irf_gse2gsm(x,flag);
%
% x - [t,Xx,Xy,Xz]  column vector, time in isdat epoch
% y - [t,Yx,Yy,Yz]  column vector
% flag - if flag not given do GSM->SM, if flag=-1 do SM->GSM any other flag values do GSM->SM
%
% Uses: IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION
%

%function [y,tilt]=irf_gsm2sm(x,flag); NOT IMPLEMENTED, NEEDED?

if nargin==0, help irf_gsm2sm;return;end
conv='gsm>sm';
if nargin==2 && flag_gsm2sm==-1
  conv='sm>gsm';
end
if isempty(x)
  irf_log('fcal','empty input variable');
  y=x;
  return;
end % if empty input, empty output

y=irf.geocentric_coordinate_transformation(x,conv);
