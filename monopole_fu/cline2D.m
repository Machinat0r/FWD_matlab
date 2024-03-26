% ------written by Wending Fu, Nov.24.2022 in Beijing------------
% This function plots a 3D line (x,y,z) encoded with scalar color data (c)
% using the specified colormap (default=jet);
% The "othercolor" package is supported
% minc & maxc is the boundary of colorbar
% SYNTAX: h=cline(x,y,c,minc,maxc,colormap,linewidth);
function h=cline2D(x,y,c,minc,maxc,cmap,linewidth)
if  nargin<7
  fprintf('Insufficient input arguments\n');
  help cline
  return;
elseif nargin==7
    if length(linewidth) == 1
        linewidth = repmat(linewidth,size(c));
    end
end
%% colormap
lw = linewidth; lw(lw<0) = 0;
try 
cmap = colormap(cmap);
% cmap = cmap(1:240,:);
catch
cmap = othercolor(cmap);
% cc1 = [linspace(0.45,0.45,256);linspace(0.44,0.7,256);linspace(0.44,0.85,256)];
% % cc2 = [linspace(0.2,0.5,224);linspace(0.5,0.6,224);linspace(0.8,1,224)];
% % cc = [cc1,cc2];
% cc = cc1;
% cmap = cc';
cmap = cmap(32:256,:);
end
yy=linspace(minc,maxc,size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',c);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;

%%
% Lot line segment with appropriate color for each data pair...
for i=1:length(x)-1
    h(i)=line([x(i) x(i+1)],[y(i) y(i+1)],'color',[cm(:,i)],'LineWidth',lw(i));
end
return