%% Plot overview figure with focus on elecyric field；Ohm law; 

ic=1

npanels = 7;
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultLineLineWidth', 0.75);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 18; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
cmap = 'jet';
h = irf_plot(npanels);
ic_B =1:4;
ic_J =1:4;
iisub = 0;
cmap = colormap('jet');
zoomy = [];
tint = irf.tint('2017-05-22T10:42:12.00Z/2017-05-22T10:42:14.00Z');
% tint = irf.tint('2016-11-21T07:22:20.000/2016-11-21T07:22:24.000');


if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.99 0.98],'fontsize',12);
  grid(hca,'off');
  hca.YLim = [-10 10];
  % hca.YTick = [-20 -10 0 10 20];
end
if 0 % E x ohm
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ex ohm');
  set(hca,'ColorOrder',mms_colors('xyza'))
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% c_eval('irf_plot(hca,{mvaE?.z, mvaE?t.z},''comp'');',ic);hold(hca,'on');
  % c_eval('irf_plot(hca,{gsmE?.x, -gsmVexB?,  gsmJxB?, },''comp'');',ic);hold(hca,'on');
  plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
  hca.YLabel.String = {'E_{X} (mV/m)'};
  irf_legend(hca,{'E  ',' -VexB ',' JxB/ne', '-\nabla\cdotP_{e}/ne'},[0.05 0.1],'fontsize',10); 
%   ' -\nabla \cdot P_{e}/q_{e}n'
  grid(hca,'off');
  hca.YLim = [-25 25];hca.YTick =  [-50 0 50];
end

if 0 % E x ohm
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E l');
%   set(hca,'ColorOrder',mms_colors('xyza'))
 set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[0.65 0.65 0.65]]);
  c_eval('irf_plot(hca,{gsmE?.x, -gsmVixB?.x, gsmOhmJxB?.x, -divPe?nel, mvares.j?l},''comp'');',ic);hold(hca,'on');
  % c_eval('irf_plot(hca,{gsmE?.x, -mvaVixB?.x, mvaOhmJxB?.x, -divPe?nel},''comp'');',ic);hold(hca,'on');

%   c_eval('irf_plot(hca,{E?par.x},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{mvaE?.x, -mvaVixB?.x, mvaOhmJxB?.x},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{mvaE?lsmooth, -mvaVixB?.x, mvaOhmJxB?.x},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{mvaE?resample.x, -mvaVixB?.x, mvaOhmJxB?.x},''comp'');',ic);hold(hca,'on');
 
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
  hca.YLabel.String = {'E_{x} (mV/m)'};
%   irf_legend(hca,{'E','-VixB','JxB/ne',' VexB'},[0.05 0.9],'fontsize',10); 
  irf_legend(hca,{'E','-VixB','JxB/ne',' -\nabla \cdot P_{e}/ne','res.j'},[0.99 0.98],'fontsize',10); 
  grid(hca,'off');
   hca.YLim = [-20 20]; %hca.YTick = [-20 0 20];
end
if 0 % % E y ohm
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E m');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[0.65 0.65 0.65]]);
%  set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE?.y, -gsmVixB?.y, gsmOhmJxB?.y, -divPe?nem, mvares.j?m},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{E?par.x},''comp'');',ic);hold(hca,'on');

%   c_eval('irf_plot(hca,{mvaE?.y, -mvaVixB?.y, mvaOhmJxB?.y},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{mvaE?msmooth, -mvaVixB?.y, mvaOhmJxB?.y},''comp'');',ic);hold(hca,'on');
   % c_eval('irf_plot(hca,{mvaE?resample.y, -mvaVixB?.y, mvaOhmJxB?.y},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
  hca.YLabel.String = {'E_{y} (mV/m)'};
%   irf_legend(hca,{'E','-VixB','JxB/ne',' VexB'},[0.05 0.9],'fontsize',10); 
  irf_legend(hca,{'E','-VixB','JxB/ne',' -\nabla \cdot P_{e}/ne','res.j'},[0.99 0.98],'fontsize',10); 
  grid(hca,'off');
   hca.YLim = [-20 20];%hca.YTick = [-20  0 20 ];
end
if 0 % E z ohm
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E n');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[0.65 0.65 0.65]]);
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{mvaE?.z, -mvaVixB?.z, mvaOhmJxB?.z, mvaVexB?.z},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{mvaE?.z, -mvaVixB?.z, mvaOhmJxB?.z, divPe?ne},''comp'');',ic);hold(hca,'on');
% c_eval('irf_plot(hca,{mvaE?nsmooth, -mvaVixB?.z, mvaOhmJxB?.z, -divPe?ne},''comp'');',ic);hold(hca,'on');
c_eval('irf_plot(hca,{gsmE?.z, -gsmVixB?.z, gsmOhmJxB?.z, -divPe?nen, mvares.j?n},''comp'');',ic);hold(hca,'on');
% c_eval('irf_plot(hca,{E?par.x},''comp'');',ic);hold(hca,'on');

plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
hca.YLabel.String = {'E_{z} (mV/m)'};
%   irf_legend(hca,{'E','-VixB','JxB/ne',' VexB'},[0.05 0.9],'fontsize',10); 
  irf_legend(hca,{'E','-VixB','JxB/ne',' -\nabla \cdot P_{e}/ne','res.j'},[0.99 0.98],'fontsize',10); 
%   ' -\nabla \cdot P_{e}/q_{e}n'
  grid(hca,'off');
   hca.YLim = [-20 20];%hca.YTick = [-20 0 20];
end

if 0 % J abs
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi abs');
  set(hca,'ColorOrder',mms_colors('1   '))
  c_eval('irf_plot(hca,{irf_abs(gsmJ?).tlim(tint)},''comp'');',ic)
%   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'|J_{mom}|','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1   '))
  % irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  grid(hca,'off');
  hca.YLim = [0 100];  
end
if 0 % Jfpi xyz
    iisub = iisub + 1;
    zoomy = [zoomy iisub];
    hca = irf_panel('J fpi');
    set(hca,'ColorOrder',mms_colors('xyza'))
    c_eval('irf_plot(hca,{gsmJ?.x,gsmJ?.y,gsmJ?.z},''comp'');',ic)
    % c_eval('irf_plot(hca,{gsmJcurl.x,gsmJcurl.y,gsmJcurl.z},''comp'');',ic)

    plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
    hca.YLabel.String = {'J_{mom gsm}','(nA/m^2)'};
    grid(hca,'off');
    % irf_legend(hca,{'mms1','mms2','Jcurl'},[0.05 0.9],'fontsize',10);
      irf_legend(hca,{'x','y','z'},[0.99 0.98],'fontsize',12);
    % irf_legend(hca,{'J_{\perp 1}','J_{\perp 2}','J_{||}'},[0.08 0.98],'fontsize',12);
    hca.YLim = [-25 60];
    %hca.YLim = [-1100 1100];
end
if 0 % E.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi');
  set(hca,'ColorOrder',mms_colors('1   '))
  c_eval('irf_plot(hca,{EdotJ?.tlim(tint)},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E\cdotJ (pW/m^3)'};
  
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-200 700];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end

if 0 % E.J fpi 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi abs');
  set(hca,'ColorOrder',mms_colors('1   '))
  c_eval('irf_plot(hca,{abs(EdotJ?.tlim(tint))},''comp'');',ic);hold(hca,'on');
%   c_eval('irf_plot(hca,{EdotJ?.tlim(tint)},''comp'');',ic);hold(hca,'on');
%   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
  grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'|E.J|','(pW/m^3)'};
  
%   irf_legend(hca,{'El.Jl','Em.Jm','En.Jn','E.J'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-200 700];  
end

clear EdotJ_inter
c_eval('aaa? = irf.ts2mat(EdotJ?.tlim(tint));',ic);
c_eval('EdotJ_inter?(:,1) = aaa?(:,1);',ic);
c_eval('EdotJ_inter?(:,2) = cumsum(aaa?(:,2));',ic);

if 0 % E.J sum 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi sum');
  set(hca,'ColorOrder',mms_colors('12  '))
  % c_eval('irf_plot(hca,{},''comp'');',ic);hold(hca,'on');
  c_eval("irf_plot([EdotJ_inter?(:,1) EdotJ_inter?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
  plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
  grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'\Sigma E\cdotJ','(pW/m^3)'};
%   irf_legend(hca,{'El.Jl','Em.Jm','En.Jn','E.J'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 5000];  
end

if 1 % E.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{EdotJ?.tlim(tint), gsmEdotJ?x, gsmEdotJ?y, gsmEdotJ?z},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E\cdotJ (pW/m^3)'};
  irf_legend(hca,{'E\cdotJ','E\cdotJ_x','E\cdotJ_y',' E\cdotJ_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 800];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end

if 1 % -VixB dot J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('-VixB.J fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{EdotJ?.tlim(tint), gsmEdotJ?x, gsmEdotJ?y, gsmEdotJ?z},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E\cdotJ (pW/m^3)'};
  irf_legend(hca,{'E\cdotJ','E\cdotJ_x','E\cdotJ_y',' E\cdotJ_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 800];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end

if 0 % -VixB 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('-VixB fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{-gsmVixB?.x, -gsmVixB?.y, -gsmVixB?.z},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'-VixB (mV/m)'};
  irf_legend(hca,{'-VixB_x','-VixB_y','-VixB_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [0 5];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end

if 1 % JxB/ne 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('JxB/ne fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{gsmOhmJxB?.x, gsmOhmJxB?.y, gsmOhmJxB?.z},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
   hca.YLabel.String = {'JxB/ne (mV/m)'};
  irf_legend(hca,{'JxB/ne_x','JxB/ne_y','JxB/ne_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-20 20];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end




if 1 % delta Pe 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('-\nabla \cdot P_{e}/ne fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{-divPe?nel, -divPe?nem, -divPe?nen},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
   hca.YLabel.String = {'-\nabla \cdot P_{e}/ne (mV/m)'};
  irf_legend(hca,{'-\nabla \cdot P_{e}/ne_x','-\nabla \cdot P_{e}/ne_y','-\nabla \cdot P_{e}/ne_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-10 10];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end

if 1 % res.j 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('res.j fpi');
  set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  c_eval('irf_plot(hca,{mvares.j?l, mvares.j?m, mvares.j?n},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--','Linewidth',0.5);hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
   hca.YLabel.String = {'res.j (mV/m)'};
  irf_legend(hca,{'res.j_x','res.j_y','res.j_z'},[0.99 0.98],'fontsize',10); 
%   irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-20 20];  
  % irf_legend(hca,{'MMS1 ',' MMS2 '},[0.05 0.1],'fontsize',10);
end


if 1 % E.J sum 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi sum');
  set(hca,'ColorOrder',mms_colors('12  '))
  % c_eval('irf_plot(hca,{},''comp'');',ic);hold(hca,'on');
  c_eval("irf_plot([EdotJ_inter?(:,1) EdotJ_inter?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
  plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
  grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'\Sigma E\cdotJ','(pW/m^3)'};
%   irf_legend(hca,{'El.Jl','Em.Jm','En.Jn','E.J'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 5000];  
end
%% Plot setting
set(gca,"XTickLabelRotation",0)
irf_adjust_panel_position
legends = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'};
nInd = 1;
for ii = [1:npanels]
    irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
    nInd = nInd + 1;
end

irf_zoom(h(1:npanels),'x',tint)
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
if 0
    hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
    for ii = 1:numel(hmark);
        hmark(ii).FaceAlpha = 0.5;
    end
end

for ii = 1:npanels;
    set(h(ii),'ticklength',[0.005 0.005]);
end
set(gcf,'color','w');
set(gcf,'render','painters');