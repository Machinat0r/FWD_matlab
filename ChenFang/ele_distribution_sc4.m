%% Plot S4 overview figure with focus on electrons, including single time electron distributions
% npanels = 9;
% % tint = irf.tint('2017-07-20T11:44:00.00Z/2017-07-20T11:44:05.00Z');
% tint = irf.tint('2017-07-17T07:53:05.00Z/2017-07-17T07:53:07.00Z');
% tintZoom = tint;
% 
% cmap = 'jet';
% [h1,h2] = initialize_combined_plot(npanels,3,4,0.45,'vertical');
% ic = 1:4;
% iisub = 0;
% cmap = colormap('jet');
% 
% isub = 0;
% zoomy = [];
% 
% 
% if 1 % B  
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('B');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
%   c_eval('irf_plot(hca,{gsmB1.z,gsmB2.z,gsmB3.z,gsmB4.z,},''comp'');',ic);
% % grid off;
%   hca.YLabel.String = {'B','(nT)'};
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]); 
%   irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.98 0.9],'fontsize',12);  
% end
% if 1 % ne
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('n');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]])
%   c_eval('irf_plot(hca,{ne1,ne2,ne3,ne4},''comp'');',ic)
%   hca.YLabel.String = {'n_e','(cm^{-3})'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);    
%   
% end
% if 0 % E
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('E');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{bdryE?.x,bdryE?.y,bdryE?.z},''comp'');',ic)
%   hca.YLabel.String = {'E','(mV/m)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
% end
% if 1 % E perp1
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('E perp1');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]])
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmE1perp.x,gsmE2perp.x,gsmE3perp.x,gsmE4perp.x},''comp'');',ic)
%   hca.YLabel.String = {'E_{\perp1}','(mV/m)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
% 
% end
% if 1 % E perp2
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('E perp2');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]])
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmE1perp.y, gsmE2perp.y, gsmE3perp.y, gsmE4perp.y},''comp'');',ic)
%   hca.YLabel.String = {'E_{\perp2}','(mV/m)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% 
% %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
% end
% if 1 % E par
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('E par');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmE1par,gsmE2par,gsmE3par,gsmE4par},''comp'');',ic)
%   hca.YLabel.String = {'E_{||}','(mV/m)'};  
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% 
% end
% if 0 % J  
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('J fpi');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
%   hca.YLabel.String = {'J','(nA/m^2)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
% end
% if 1 % Jn   
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Jn');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmJ1perp.x,gsmJ2perp.x,gsmJ3perp.x,gsmJ4perp.x},''comp'');',ic)
%   hca.YLabel.String = {'Jn','(nA/m^2)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
% end
% if 1 % J m
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Jm');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmJ1perp.y,gsmJ2perp.y,gsmJ3perp.y,gsmJ4perp.y},''comp'');',ic)
%   hca.YLabel.String = {'Jm','(nA/m^2)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
% end
% if 1 % J l
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('J par');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{gsmJ1par,gsmJ2par,gsmJ3par,gsmJ4par},''comp'');',ic)
%   hca.YLabel.String = {'J_{||}','(nA/m^2)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
% %   set(hca,'ColorOrder',mms_colors('xyza'))
% %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
% end
% if 0 % Vi  
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Vi');
%   set(hca,'ColorOrder',mms_colors('xyza'))
% %   c_eval('irf_plot(hca,{bdryVi?.x.tlim(tint),bdryVi?.y.tlim(tint),bdryVi?.z.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_i','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
% end
% if 0 % Ve  
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Ve');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);     
% end
% if 0 % ve perp
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Ve');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
% end
% if 1 % Ve par
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Ve par');
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]]);
%   c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
%    c_eval('irf_plot(hca,{gsmVe1par,gsmVe2par,gsmVe3par,gsmVe4par},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'ve_{||}','(km/s)'};
%   set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1];[0 0 0]])
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % Ve perp par
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Ve perp par');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint),bdryVe?par.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'v_{\perp}','n_{\perp}','B_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % Q
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('sqrtQ');
%   set(hca,'ColorOrder',mms_colors('1234'))
% % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   c_eval('irf_plot(hca,{Q1,Q2,Q3,Q4},''comp'');',ic)
% %   hca.YLabel.String = {'$$\sqrt{Q}$$'};
%   ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
%   set(hca,'ColorOrder',mms_colors('12'))    
% end
% % ic=1;
% if 0 % e DEF omni 64
%     
%   isub = isub + 1;
%   hca = irf_panel('e DEF omni 64');  
%   c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
%   set(hca,'yscale','log');
%   set(hca,'ytick',[1e1 1e2 1e3 1e4]);
% %   hold(hca,'on')
% %   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
% %   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
% %   hold(hca,'off')
%   hca.YLabel.String = {'E_e','(eV)'};   
%   colormap(hca,cmap) 
% end
% if 0 % ePDist pa 64
%   isub = isub + 1;
%   hca = irf_panel('e PA e64 deflux lowe');  
%   eint = [100 40000];  
%   c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
%   hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 0 % Te par perp Ti/Tref
%   isub = isub + 1;
%   zoomy = [zoomy isub];
%   hca = irf_panel('Te');
%   set(hca,'ColorOrder',mms_colors('123'))
%   refTi = 10;
%   c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
%   hca.YLabel.String = {'T','(eV)'};
%   set(hca,'ColorOrder',mms_colors('123'))
%   irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
%   %hca.YLim = [10 400];  
%   %irf_zoom(hca,'y')
% end
% 
% 
% legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
% nInd = 1;
% for ii = 1:npanels
%   irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
%   nInd = nInd + 1;
%   h1(ii).FontSize = 12;
% end
% 
% %irf_zoom(h(1:iisub),'x',fastTint)
% irf_zoom(h1,'x',tintZoom)
% irf_zoom(h1(zoomy),'y')
% irf_plot_axis_align
% h1(1).Title.String = irf_ssub('MMS ?',ic);
% if 0
% hmark = irf_pl_mark(h1(1:6),tintBCS, 'yellow');
% for ii = 1:numel(hmark);
%   hmark(ii).FaceAlpha = 0.5;
% end
% end
%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 7;
% tint = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
tint = irf.tint('2017-07-17T07:53:15.40Z/2017-07-17T07:53:16.50Z');
tintZoom = tint;
cmap = 'jet';
[h1,h2] = initialize_combined_plot(npanels,3,3,0.4,'vertical');
ic = 1:4;
iisub = 0;
cmap = colormap('jet');

isub = 0;
zoomy = [];

if 0 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z,gsmB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?.x,bdryE?.y,bdryE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?perp.x,gsmE?perp.y,gsmE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{gsmE1par,gsmE2par,gsmE3par,gsmE4par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % J  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmJ?.x.tlim(tint),gsmJ?.y.tlim(tint),gsmJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % E.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{EdotJ?.tlim(tint),EdotJ?par.tlim(tint),EdotJ?perp.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{EdotJ?.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E.J','(nW/m^3)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{bdryVi?.x.tlim(tint),bdryVi?.y.tlim(tint),bdryVi?.z.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gsmVi?.x.tlim(tint),gsmVi?.y.tlim(tint),gsmVi?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);     
end
if 0 % ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
end
if 0 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % Ve perp par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint),bdryVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp}','n_{\perp}','B_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % e DEF omni 64
  isub = isub + 1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%   hold(hca,'on')
%   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [5000 40000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 64 
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux 1 ');  
   eint = [100,400]; 
%    eint = [5000,270000]; 
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};  
  hca.YLabel.String = {['0.1-0.4'],'keV'};
  hca.YTick = [45 90 135]; 
%   irf_colormap('standard')
  colormap(hca,cmap)
%   hca.CLim = [6.9 7.1];
end
if 0 % ePDist pa 64 
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux 2 ');  
   eint = [400,1500];  
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};  
  hca.YLabel.String = {['0.4-1.5'],'keV'};
  hca.YTick = [45 90 135]; 
%   irf_colormap('standard')
  colormap(hca,cmap)
%   hca.CLim = [6.9 7.1];
end
if 0 % ePDist pa 64 
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux 3 ');  
   eint = [1500,270000];  
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};  
  hca.YLabel.String = {['1.5-27'],'keV'};
  hca.YTick = [45 90 135]; 
%   irf_colormap('standard')
  colormap(hca,cmap)
%   hca.CLim = [6.9 7.1];
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % Q
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('sqrtQ');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{Q?},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmB1.z,gsmB2.z,gsmB3.z,gsmB4.z,},''comp'');',ic);
% grid off;
  hca.YLabel.String = {'B','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]); 
  irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE1par,gsmE2par,gsmE3par,gsmE4par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);

end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('ne');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{ne1,ne2,ne3,ne4},''comp'');',ic);
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);    
  
end
if 0 % ni
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('ni');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{ni1,ni2,ni3,ni4},''comp'');',ic)
  hca.YLabel.String = {'n_i','(cm^{-3})'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);    
  
end
if 1 % E perp1
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp1');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE1perp.x,gsmE2perp.x,gsmE3perp.x,gsmE4perp.x},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp1}','(mV/m)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  hca.YLim = [-20 20];
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  

end
if 1 % E perp2
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp2');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE1perp.y, gsmE2perp.y, gsmE3perp.y, gsmE4perp.y},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp2}','(mV/m)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
hca.YLim = [-20 20];
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % Jn   
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Jn');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmJ1perp.x,gsmJ2perp.x,gsmJ3perp.x,gsmJ4perp.x},''comp'');',ic)
  hca.YLabel.String = {'Jn','(nA/m^2)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % J m
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Jm');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmJ1perp.y,gsmJ2perp.y,gsmJ3perp.y,gsmJ4perp.y},''comp'');',ic)
  hca.YLabel.String = {'Jm','(nA/m^2)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % J l
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmJ1par,gsmJ2par,gsmJ3par,gsmJ4par},''comp'');',ic)
  hca.YLabel.String = {'J_{||}','(nA/m^2)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmVe1par,gsmVe2par,gsmVe3par,gsmVe4par},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'ve_{||}','(km/s)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmVi?par.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmVi1par,gsmVi2par,gsmVi3par,gsmVi4par},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'vi_{||}','(km/s)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %hca.YLim = [-1100 1100];  
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h1(ii).FontSize = 10;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h1,'x',tintZoom)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h1(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
set(gcf,'color','w');
%% Plot single time particle distributions, 1 sc, 4 projections,
ic=1:4;
% tintDist = irf.tint('2017-07-20T11:44:01.00Z/2017-07-20T11:44:05.00Z');
% tintDist = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
tintDist = irf.tint('2017-07-12T11:54:00.000Z/2017-07-12T11:55:00.000Z');
c_eval('dist? = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('dist_scm? = ePDist?;',ic)
 c_eval('scpot? = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)

c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim =60*1e3;
elevlim = 15;
strCMap = 'jet';
% projclim = [-1 2];  
projclim_int = [-13 -10];
  

c_eval('times = ePDist1.time;',ic)
[tind,~] = times.tlim(tintDist);

doReducedF = 1;
for it_ =1:1:numel(tind);
  it = tind(it_);
  time = times(it);
  it
  if exist('hmark'); delete(hmark); end
  c_eval('hmark_tmp = irf_pl_mark(h1(?),time.epochUnix,''green''); hmark(?) = hmark_tmp;',1:npanels)
  
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  
%   c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
%   par = hatB;
%   newy = [0 1 0]; 
%   perp1 = cross(newy,par);
%   perp2 = cross(par,perp1);  
  

  timeUTC = time.utc;      
  isub = 1;

  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
%   vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
%   
%   if 1 % Perpendicular plane, slice
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;perp2;par]; vlabels = {'v_{Bx[0 1 0]}','v_{Bx([0 1 0]xB)}','v_{B}'};
%    % mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
%  mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);   
% %mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);
%     hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
%   end 
%   if 1 % B plane 1, slice
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;par;-perp2]; vlabels = {'v_{Bx[0 1 0]}','v_{B}','-v_{Bx([0 1 0]xB)}'};
% %     mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%  mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);     
% %mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
%     hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
%   end 
%   if 1 % B plane 2, slice
%     hca = h2(isub); isub = isub + 1;
%     xyz = [perp2;par;perp1]; vlabels = {'v_{Bx([0 1 0]xB)}','v_{B}','v_{Bx[0 1 0]}'};
% %     mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%  mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);   
% %mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
%     hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
%   end
  vint = [-Inf Inf];
  % mms1
  if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{Bx(ExB)}','v_{B}'};
    mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('mms1');
  end 
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{B}','-v_{Bx(ExB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{ExB}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  %mms2
  if 1 % Perpendicular plane, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{Bx(ExB)}','v_{B}'};
     mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
     hca.Title.String = sprintf('mms2');
  end 
  if 1 % B plane 1, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{B}','-v_{Bx(ExB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 1 % B plane 2, integrated, smaller vint
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{ExB}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 
  %mms3
  if 0 % Perpendicular plane, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{Bx(ExB)}','v_{B}'};
     mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
     hca.Title.String = sprintf('mms3');
  end 
  if 0 % B plane 1, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{B}','-v_{Bx(ExB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 0 % B plane 2, integrated, smaller vint
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{ExB}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end
  %mms4
  if 0 % Perpendicular plane, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{Bx(ExB)}','v_{B}'};
     mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
     hca.Title.String = sprintf('mms4');
  end 
  if 0 % B plane 1, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{B}','-v_{Bx(ExB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 0 % B plane 2, integrated, smaller vint
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{ExB}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzcleint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 0 % Pitchangle distribution
    hca = h2(isub); isub = isub + 1;
    plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
    hca.YScale = 'log'; hca.XScale = 'log';
    hca.YLabel.String = ['f_e (' ePitch.units ')'];
    hca.XLabel.String = 'E (eV)';
    hca.XLim = [50 30000];
    legend(hca,{'0','90','180'})
    hca.YTick = 10.^[-4:2];
    hca.YLim = [1e-4 1e2];
  end
  if 0 % invisible
    hca = h2(isub); isub = isub + 1;
    hca.Visible = 'off';
  end
  if 0 % invisible
    hca = h2(isub); isub = isub + 1;
    hca.Visible = 'off';
  end
  % ExB plane, with markings
  if 0
    vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
    mms.plot_int_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    hca.Title.String = '';
  end
  
  h2(1).Title.String = {timeUTC(1:23),h2(1).Title.String};
  for ii = 1:2
    colormap(h2(ii),strCMap)
  end
  %cn.print(['e_proj_in_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath 'proj_int_mms1/'])
  set(gcf,'render','painters');
figname=['20170717_el03e5_' num2str(it)];
print(gcf, '-dpdf', [figname '.pdf']);
end

