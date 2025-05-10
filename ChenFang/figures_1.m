

%  % Plot overview figure with focus on ions and electrons, including hpca
% npanels =7;
% 
% set(0,'DefaultAxesFontSize',10);
% 
% set(0,'DefaultLineLineWidth', 0.5);
% fn=figure(1);clf;
% set(gcf,'PaperUnits','centimeters')
% xSize = 13; ySize = 16; coef=floor(min(800/xSize,800/ySize));
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% 
% cmap = 'jet';
% h = irf_plot(npanels);
% ic =1;
% iisub = 0;
% cmap = colormap('jet');
% zoomy = [];
% tint = irf.tint('2017-07-17T07:52:45.00Z/2017-07-17T07:53:15.00Z');
% % tint = irf.tint('2017-08-06T00:08:55.00Z/2017-08-06T00:09:05.00Z');
% 
% if 1 % B
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('B');
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]])
%   %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
% %   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
%    c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   hca.YLabel.String = {'B','(nT)'};
%    
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]])
%   irf_legend(hca,{'B_{x}   ','B_{y}   ','B_{z}   ','|B|'},[0.1 0.95],'fontsize',12);
%   grid(hca,'off');
% end
% 
% if 0 % FPI: ne HPCA: O H He 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('n e i');
%  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{ne?,nH?,nHe?,nO?},''comp'');',ic)
%   hca.YLabel.String = {'n','(cm^{-3})'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'e','H','He','O'},[0.98 0.9],'fontsize',12);  
% end
% if 0 % ne ni 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('n');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
%   hca.YLabel.String = {'n','(cm^{-3})'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])  
%   irf_legend(hca,{'n_e','n_i'},[0.05 0.98],'fontsize',8);
% end
% if 1 % ne 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('n');
%   set(hca,'ColorOrder',[[0 0 0];[0 1 0]])
%   c_eval('irf_plot(hca,{ne?},''comp'');',ic)
%   hca.YLabel.String = {'n_e','(cm^{-3})'}; 
% grid(hca,'off');
% %   irf_legend(hca,{'n_e'},[0.05 0.95],'fontsize',8);
% end
% if 1 % Vi 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Vi');
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]])
%   c_eval('irf_plot(hca,{gsmVi?.x.tlim(tint),gsmVi?.y.tlim(tint),gsmVi?.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_i','(km/s)'};
%  
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]])
%   irf_legend(hca,{'Vi_{x}   ','Vi_{y}   ','Vi_{z}'},[0.1 0.98],'fontsize',12);     
% grid(hca,'off');
% end
% if 0 % Ve 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Ve');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
% end
% if 0 % ve perp
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Ve');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',12);   
% end
% if 0 % Ve par
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Ve par');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_{e,||}','(km/s)'};
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % Ve perp 1
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Ve perp');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmVe?perp.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % E
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
%   hca.YLabel.String = {'E','(mV/m)'};
%    set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
%   irf_zoom(hca,'y')
% end
% if 0 % E perp
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E perp');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmE?perp.x,gsmE?perp.y,gsmE?perp.z},''comp'');',ic)
%   hca.YLabel.String = {'E_{\perp}','(mV/m)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',12);  
% end
% if 1 % Te par perp
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Te');
%   set(hca,'ColorOrder',[[1 0 0];[0 0 1]])
% %   refTi = 10;
% %   c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic);
%    c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic);
%   hca.YLabel.String = {'T','(eV)'};
% grid(hca,'off');
%   set(hca,'ColorOrder',[[1 0 0];[0 0 1]])
% %   irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
%    irf_legend(hca,{'Te_{||}     ','Te_{\perp}'},[0.1 0.95],'fontsize',12);
%   %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 4000];
%   %hca.YTick
%   irf_zoom(hca,'y')
% 
% end
% if 0 % J 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('J fpi');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmJ?.x.tlim(tint),gsmJ?.y.tlim(tint),gsmJ?.z.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'J','(nA/m^2)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x','y','z'},[0.1 0.98],'fontsize',12);  
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % E.J 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E.J fpi');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{EdotJ?par.tlim(tint),EdotJ?perp.tlim(tint),EdotJ?.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'E.J','(nW/m^3)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'E.Jpar','E.Jperp','E.J'},[0.1 0.9],'fontsize',12);  
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % E'.J 
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E`.J fpi');
% %   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{RedotJ?.tlim(tint)});',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'E`.J','(nW/m^3)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % Te par perp
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Te');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   refTi = 10;
%   c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
%   hca.YLabel.String = {'T','(eV)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
%   %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 4000];
%   %hca.YTick
%   irf_zoom(hca,'y')
% end
% if 0 % Q
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('sqrtQ');
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0,0,0]])
%   c_eval('irf_plot(hca,{Q?},''comp'');',ic)
% %   hca.YLabel.String = {'$$\sqrt{Q}$$'};
%   ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
%   set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0,0,0]])    
% end
% if 0 % E par
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E par');
%   set(hca,'ColorOrder',[[0 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
%   hca.YLabel.String = {'E_{||}','(mV/m)'};  
%   grid off;
% end
% if 0 % Ve perp par
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('Ve perp par');
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint),gsmVe?par.tlim(tint)},''comp'');',ic)
%   %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
%   hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
%   irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.98],'fontsize',12);  
%   %hca.YLim = [-1100 1100];  
% end
% if 0 % gradPe
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('gradPe');
%   set(hca,'ColorOrder',mms_colors('xyz'))
%   irf_plot(hca,{gsmGradPe.x*1e3,gsmGradPe.y*1e3,gsmGradPe.z*1e3},'comp');
%   hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
%   set(hca,'ColorOrder',mms_colors('xyz'))  
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
%   irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
% end
% if 0 % e DEF omni 32
%   iisub = iisub + 1;
%   hca = irf_panel('e DEF omni');  
%   c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
%   hca.YLabel.String = {'E_e','(eV)'};  
%   set(hca,'yscale','log');
%   set(hca,'ytick',[1e1 1e2 1e3 1e4 2e4]);
% end
% 
% if 0 % ePDist pa 64 
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux 1 ');  
%    % eint = [100,400]; 
%    eint = [5000,30000];
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
% %   c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%    try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% end
% 
%  
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% %   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};    
% % hca.YLabel.String = {['5-27'],'keV'}; 
%  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
% %   hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   hca.YTick = [45 90 135]; 
% %   irf_colormap('standard')
%   colormap(hca,cmap)
% hca.CLim = [5.5 7.5];
% end
% if 0 % ePDist pa 64 
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux 2 ');  
%    eint = [400,1500];  
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% %   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};  
%   hca.YLabel.String = {['0.4-1.5'],'keV'};
%   hca.YTick = [45 90 135]; 
% %   irf_colormap('standard')
%   colormap(hca,cmap)
% %   hca.CLim = [6.9 7.1];
% end
% if 0 % ePDist pa 64 
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux 3 ');  
%    eint = [1500,270000];  
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% %   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};  
%   hca.YLabel.String = {['1.5-27'],'keV'};
%   hca.YTick = [45 90 135]; 
% %   irf_colormap('standard')
%   colormap(hca,cmap)
% %   hca.CLim = [6.9 7.1];
% end
% if 0 % ePDist pa 64 2
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux lowe 2');  
%    eint = [43,283];  
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {[ '43-280' ],'eV'};
%   hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 0 % ePDist pa 64 3
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux lowe 3');  
%    eint = [485,4179];  
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% %   hca.YLabel.String = {[num2str(eint(1)) '-' num2str(floor(eint(2)/100))],'eV'};  
% hca.YLabel.String = {['0.48-4.2'],'keV'};
%   hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 0 % ePDist pa 64 4
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux lowe 4');  
%    eint = [5470,27525];  
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {['5.5-27'],'keV'};
%   hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 1 % e DEF omni 64
%   iisub = iisub + 1;
%   hca = irf_panel('e DEF omni 64');  
%   c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
%   set(hca,'yscale','log');
%   set(hca,'ytick',[1e2 1e3 1e4]);
%   set(hca,'ylim',[1e1 3e4]);
%   hold(hca,'on')
%   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
%   hca.YLabel.String = {'E_e','(eV)'};   
%   colormap(hca,cmap) 
%    hca.CLim = [4 8];
%  grid(hca,'off');
% end
% %
% if 0 %low
% 
% %
%  iisub = iisub + 1;
%   hca = irf_panel('e energy low');
% 
% specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
% specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_elow.p=energy_low.data;%data matrix
% specrec_p_elow.f_label='';
% specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
% irf_spectrogram(h(iisub),specrec_p_elow)
% irf_legend(hca,[ '0<E<200' ],[0.99 0.90],'color',0*[1 1 1]) 
% hca.YLabel.String = {'\theta_{low,e}','(\circ)'};   
% hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 0 %mid
%  iisub = iisub + 1;
%   hca = irf_panel('e energy mid');
% specrec_p_emid=struct('t',irf_time(energy_mid.DEPEND_0.data,'ttns>epoch'));
% specrec_p_emid.f=transpose(energy_mid.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_emid.p=energy_mid.data;%data matrix
% specrec_p_emid.f_label='';
% specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
% irf_spectrogram(h(iisub),specrec_p_emid)
% irf_legend(hca,[ '200<E<2000'],[0.99 0.90],'color',0*[1 1 1]) 
% hca.YLabel.String = {'\theta_{mid,e}','(\circ)'};   
%   
% hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% % irf_colormap('standard')
% end
% if 0 %high
%  iisub = iisub + 1;
%   hca = irf_panel('e energy high');
% specrec_p_ehigh=struct('t',irf_time(energy_high.DEPEND_0.data,'ttns>epoch'));
% specrec_p_ehigh.f=transpose(energy_high.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_ehigh.p=energy_high.data;%data matrix
% specrec_p_ehigh.f_label='';
% specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
% irf_spectrogram(h(iisub),specrec_p_ehigh)
% irf_legend(hca,['2000<E<27000' ],[0.99 0.90],'color',0*[1 1 1]) 
% hca.YLabel.String = {'\theta_{high,e}','(\circ)'};   
%   
% hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% % irf_colormap('standard')
% end
% 
% 
% if 0 % iPDist pa 64
%   iisub = iisub + 1;
%   hca = irf_panel('i PA e64 deflux lowe');  
%   eint = [5000 40000];  
%   try
%     c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% end
% 
%   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
%   hca.YTick = [45 90 135];   
%   colormap(hca,cmap)
% end
% if 0 % e DEF omni 64
%   iisub = iisub + 1;
%   hca = irf_panel('e DEF omni 64');  
%   c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
%   set(hca,'yscale','log');
%   set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%   set(hca,'ylim',[1e1 1e4]);
%    hca.YLabel.String = {'E_e','(eV)'}; 
% end
% if 1 % i DEF omni 64
%   iisub = iisub + 1;
%   hca = irf_panel('i DEF omni 64');  
%   c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
%   set(hca,'yscale','log');
%   set(hca,'ytick',[1e2 1e3 1e4]);
%   set(hca,'ylim',[1e1 3e4]);
% %   hold(hca,'on')
% %   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
% %   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
% %   hold(hca,'off')
%   hca.YLabel.String = {'E_i','(eV)'};   
%   colormap(hca,cmap); 
%   hca.CLim = [3 6];
%    grid off;
% end
% if 0 % ePDist pa 64 
%   iisub = iisub + 1;
%   hca = irf_panel('e PA e64 deflux 1 ');  
%    % eint = [100,400]; 
%    eint = [5000,30000];
% %   eint = [1.5*max(scPot1.data) 40000];  
%   try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   end
% %   c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
%    try
%     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
%   catch
%     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% end
% 
%  
%   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% %   hca.YLabel.String = {[num2str(eint(1)),'-',num2str(eint(2))],'eV'};    
% % hca.YLabel.String = {['5-27'],'keV'}; 
%  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
% %   hca.YLabel.String = {'Pitchangle','(\circ)'};   
%   hca.YTick = [45 90 135]; 
% %   irf_colormap('standard')
%   colormap(hca,cmap)
% hca.CLim = [6 7.5];
% end
% if 1 % E par
%   iisub = iisub + 1;
%   zoomy = [zoomy iisub];
%   hca = irf_panel('E par');
%   set(hca,'ColorOrder',[[0 0 0];[0 1 0];[0 0 1]])
%   c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
%   hca.YLabel.String = {'E_{||}','(mV/m)'}; 
%    grid off;
% end
% 
% 
%  legends = {'c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
%  nInd = 1;
%  for ii = [1:npanels]  
%    irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
%   nInd = nInd + 1;
%  end
% 
% irf_zoom(h(1:npanels),'x',tint)
% irf_zoom(h(zoomy),'y')
% %irf_zoom(h([1:3 6 10]),'y')
% %hca = irf_panel('Te'); irf_zoom(hca,'y')
% %hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
% %hca = irf_panel('n'); hca.YLim = [0 11];
% %hca = irf_panel('J zoom'); hca.YLim = [-500 800];
% %hca = irf_panel('Ve'); hca.YLim = [-900 500];
% %hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
% %hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
% %hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
% %hca = irf_panel('Ti'); hca.YLim = [300 700];
% 
% irf_plot_axis_align
% h(1).Title.String = irf_ssub('MMS ? GSM',ic);
% % if 0
% % hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
% % for ii = 1:numel(hmark);
% %   hmark(ii).FaceAlpha = 0.5;
% % end
% 
% for ii=1:7;
% irf_pl_mark(h(ii),[iso2epoch('2017-07-17T07:53:05.50Z')],'k');
% irf_pl_mark(h(ii),[iso2epoch('2017-07-17T07:53:06.10Z')],'k');
% end
% 
% % for ii = 1:npanels;
% %   h(ii).FontSize = 12;
% % end
% %irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
% set(gcf,'color','w');
%  set(gcf,'render','painters');
% figname=['tu1'];
% print(gcf, '-dpdf', [figname '.pdf']);

%% e PAD 32 energy

 grid off;ic=1;
%tint = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
tint = irf.tint('2018-09-11T08:05:38.00Z/2018-09-11T08:05:42.00Z')
cmap = colormap('jet');
energy=ePitch1.depend{1,1}(1,:);
E=zeros(32,2);
for ii=1:32
    E(ii,:)=[energy(ii)-1,energy(ii)+1];  
end

npanels = 14
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
%
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
  grid off;
end
if 1 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % ePDist pa 64 1
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 1');  
   eint=E(1,:);
   %eint=[E(1,1),E(2,2)];
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 2
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 2');  
    eint=E(2,:); 
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 3
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 3');  
   eint=E(3,:); 
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};      
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 4
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 4');  
   eint=E(4,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 5
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 5');  
   eint=E(5,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 6
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 6');  
   eint=E(6,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 7
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 7');  
   eint=E(7,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 8
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 8');  
   eint=E(8,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
 end
if 1 % ePDist pa 64 9
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 9');  
   eint=E(9,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 10
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 10');  
   eint=E(10,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 11
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 11');  
   eint=E(11,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 12
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 12');  
   eint=E(12,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 13
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 13');  
   eint=E(13,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? e',ic);
set(gcf,'color','w');
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['2017 ele_pitch1'];
print(gcf, '-dpng', [figname '.png']);
%

npanels = 10;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
if 0 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

if 1 % ePDist pa 64 14
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 14');  
   eint=E(14,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 15
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 15');  
   eint=E(15,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 16
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 16');  
   eint=E(16,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 17
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 17');  
   eint=E(17,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 18
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 18');  
   eint=E(18,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 19
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 19');  
   eint=E(19,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 20
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 20');  
   eint=E(20,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 21
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 21');  
   eint=E(21,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 22
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 22');  
   eint=E(22,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? e',ic);
set(gcf,'color','w');
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['2017 ele_pitch2'];
print(gcf, '-dpng', [figname '.png']);

figure;
npanels = 11;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(3);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
   set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

if 1 % ePDist pa 64 23
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 23');  
   eint=E(23,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 24
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 24');  
   eint=E(24,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 25
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 25');  
   eint=E(25,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 26
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 26');  
   eint=E(26,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 27
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 27');  
   eint=E(27,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 28
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 28');  
   eint=E(28,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 29
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 29');  
   eint=E(29,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 30
  hca = irf_panel('e PA e64 deflux lowe 30');  
   eint=E(30,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 31
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 31');  
   eint=E(31,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 32
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 32');  
   eint=E(32,:);
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
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? e',ic);
set(gcf,'color','w');
set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['2017 ele_pitch3'];
print(gcf, '-dpng', [figname '.png']);
%% i PAD 32 energy
ic=4;
% tint = irf.tint('2017-07-25T22:10:00.00Z/2017-07-25T22:10:40.00Z');
cmap = colormap('jet');
energy=iPitch1.depend{1,1}(1,:);
E=zeros(32,2);
for ii=1:32
    E(ii,:)=[energy(ii)-1,energy(ii)+1];  
end


npanels = 14;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
%
if 0 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % iPDist pa 64 1
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 1');  
   eint=E(1,:);
   %eint=[E(1,1),E(2,2)];
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 2
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 2');  
    eint=E(2,:); 
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 3
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 3');  
   eint=E(3,:); 
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};      
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 4
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 4');  
   eint=E(4,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 5
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 5');  
   eint=E(5,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 6
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 6');  
   eint=E(6,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 7
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 7');  
   eint=E(7,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 8
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 8');  
   eint=E(8,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
 end
if 1 % iPDist pa 64 9
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 9');  
   eint=E(9,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 10
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 10');  
   eint=E(10,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 11
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 11');  
   eint=E(11,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 12
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 12');  
   eint=E(12,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 13
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 13');  
   eint=E(13,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? i',ic);
set(gcf,'color','w');
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['20170717 ion_pitch1'];
print(gcf, '-dpdf', [figname '.pdf']);%

figure;
npanels = 10;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
if 0 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % iPDist pa 64 14
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 14');  
   eint=E(14,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 15
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 15');  
   eint=E(15,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 16
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 16');  
   eint=E(16,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 17
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 17');  
   eint=E(17,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 18
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 18');  
   eint=E(18,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 19
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 19');  
   eint=E(19,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 20
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 20');  
   eint=E(20,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 21
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 21');  
   eint=E(21,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 22
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 22');  
   eint=E(22,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? i',ic);
set(gcf,'color','w');
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['20170717 ion_pitch2'];
print(gcf, '-dpdf', [figname '.pdf']);%
figure;
npanels = 11;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(3);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;

zoomy = [];
if 0 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 %E
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]])
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % iPDist pa 64 23
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 23');  
   eint=E(23,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 24
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 24');  
   eint=E(24,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 25
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 25');  
   eint=E(25,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 26
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 26');  
   eint=E(26,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 27
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 27');  
   eint=E(27,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 28
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 28');  
   eint=E(28,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 29
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 29');  
   eint=E(29,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 30
  hca = irf_panel('i PA e64 deflux lowe 30');  
   eint=E(30,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 31
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 31');  
   eint=E(31,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % iPDist pa 64 32
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe 32');  
   eint=E(32,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,iPDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? i',ic);
set(gcf,'color','w');
 set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['20170717 ion_pitch3'];
print(gcf, '-dpdf', [figname '.pd']);%

