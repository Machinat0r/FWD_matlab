%% e PAD 32 energy
tint = irf.tint('2021-07-21T13:19:55.00Z/2021-07-21T13:20:15.00Z');
ic=1;

cmap = colormap('jet');
energy=ePitch1.depend{1,1}(1,:);
E=zeros(32,2);
for ii=1:32
    E(ii,:)=[energy(ii)-1,energy(ii)+1];  
end

% % % %% 6.5-165 eV
% % % npanels = 15;
% % % set(0,'DefaultAxesFontSize',8);
% % % set(0,'DefaultLineLineWidth', 0.5);
% % % fn=figure(1);clf;
% % % set(gcf,'PaperUnits','centimeters')
% % % xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
% % % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % % set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% % % set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% % % % cmap = 'jet';
% % % h = irf_plot(npanels);
% % % iisub = 0;
% % % 
% % % zoomy = [];
% % % %
% % % if 1 % B
% % %   iisub = iisub + 1;
% % %   zoomy = [zoomy iisub];
% % %   hca = irf_panel('B');
% % % % % %   set(hca,'ColorOrder',mms_colors('xyza'))
% % %   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
% % % %   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
% % % %    c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
% % %   hca.YLabel.String = {'B','(nT)'};
% % % % % %   set(hca,'ColorOrder',mms_colors('xyza'))
% % %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
% % % end
% % % if 1 % ne ni 
% % %   iisub = iisub + 1;
% % %   zoomy = [zoomy iisub];
% % %   hca = irf_panel('n');
% % % % % %   set(hca,'ColorOrder',mms_colors('12'))
% % %   c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
% % %   hca.YLabel.String = {'n','(cm^{-3})'};
% % % % % %   set(hca,'ColorOrder',mms_colors('12'))  
% % %   irf_legend(hca,{'n_e','n_i'},[0.08 0.98],'fontsize',12);
% % %   grid(hca,'off');
% % % %     hca.YLim = [0 0.65];
% % % end
% % % if 1 % ePDist pa 64 1
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 1');  
% % %    eint=E(1,:);
% % %    %eint=[E(1,1),E(2,2)];
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 2
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 2');  
% % %     eint=E(2,:); 
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 3
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 3');  
% % %    eint=E(3,:); 
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};      
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 4
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 4');  
% % %    eint=E(4,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 5
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 5');  
% % %    eint=E(5,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 6
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 6');  
% % %    eint=E(6,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 7
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 7');  
% % %    eint=E(7,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 8
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 8');  
% % %    eint=E(8,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % %  end
% % % if 1 % ePDist pa 64 9
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 9');  
% % %    eint=E(9,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 10
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 10');  
% % %    eint=E(10,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 11
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 11');  
% % %    eint=E(11,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 12
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 12');  
% % %    eint=E(12,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 13
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 13');  
% % %    eint=E(13,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % %
% % % legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
% % % nInd = 1;
% % % for ii = [1:npanels]  
% % %   irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
% % %   nInd = nInd + 1;
% % % end
% % % 
% % % %irf_zoom(h(1:iisub),'x',fastTint)
% % % irf_zoom(h(1:npanels),'x',tint)
% % % irf_zoom(h(zoomy),'y')
% % % 
% % % irf_plot_axis_align
% % % h(1).Title.String = irf_ssub('MMS ? GSM',ic);
% % % set(gcf,'color','w');
% % % 
% % % %% 216-1863 eV
% % % npanels = 11;
% % % set(0,'DefaultAxesFontSize',8);
% % % set(0,'DefaultLineLineWidth', 0.5);
% % % fn=figure(2);clf;
% % % set(gcf,'PaperUnits','centimeters')
% % % xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
% % % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % % set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% % % set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% % % % cmap = 'jet';
% % % h = irf_plot(npanels);
% % % iisub = 0;
% % % 
% % % zoomy = [];
% % % if 1 % B
% % %   iisub = iisub + 1;
% % %   zoomy = [zoomy iisub];
% % %   hca = irf_panel('B');
% % % % % %   set(hca,'ColorOrder',mms_colors('xyza'))
% % %   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
% % % %   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
% % % %    c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
% % %   hca.YLabel.String = {'B','(nT)'};
% % % % % %   set(hca,'ColorOrder',mms_colors('xyza'))
% % %   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
% % % end
% % % if 1 % ne ni 
% % %   iisub = iisub + 1;
% % %   zoomy = [zoomy iisub];
% % %   hca = irf_panel('n');
% % % % % %   set(hca,'ColorOrder',mms_colors('12'))
% % %   c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
% % %   hca.YLabel.String = {'n','(cm^{-3})'};
% % % % % %   set(hca,'ColorOrder',mms_colors('12'))  
% % %   irf_legend(hca,{'n_e','n_i'},[0.08 0.98],'fontsize',12);
% % %   grid(hca,'off');
% % % %     hca.YLim = [0 0.65];
% % % end
% % % if 1 % ePDist pa 64 14
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 14');  
% % %    eint=E(14,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 15
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 15');  
% % %    eint=E(15,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 16
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 16');  
% % %    eint=E(16,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 17
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 17');  
% % %    eint=E(17,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 18
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 18');  
% % %    eint=E(18,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 19
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 19');  
% % %    eint=E(19,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 20
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 20');  
% % %    eint=E(20,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 21
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 21');  
% % %    eint=E(21,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 22
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 22');  
% % %    eint=E(22,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
% % %   hca.YTick = [45 90 135];   
% % %   colormap(hca,cmap)
% % % end
% % % %
% % % legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)'};
% % % nInd = 1;
% % % for ii = [1:npanels]  
% % %   irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
% % %   nInd = nInd + 1;
% % % end
% % % 
% % % %irf_zoom(h(1:iisub),'x',fastTint)
% % % irf_zoom(h(1:npanels),'x',tint)
% % % irf_zoom(h(zoomy),'y')
% % % 
% % % irf_plot_axis_align
% % % h(1).Title.String = irf_ssub('MMS ? GSM',ic);
% % % set(gcf,'color','w');

%% 2438-27525 eV
npanels = 12;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
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
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
%    c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne ni 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
% % %   set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
% % %   set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.08 0.98],'fontsize',12);
  grid(hca,'off');
%     hca.YLim = [0 0.65];
end
% % % 
% % % if 1 % ePDist pa 64 20-----1087.3
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 20');  
% % %    eint=E(20,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
% % %   hca.YTick = [45 90 135];   
% % %   caxis(hca,[6.5 6.9])
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 21-----1423.3
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 21');  
% % %    eint=E(21,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};     
% % %   hca.YTick = [45 90 135];   
% % %   caxis(hca,[6.6 7.1])
% % %   colormap(hca,cmap)
% % % end
% % % if 1 % ePDist pa 64 22-----1863.2
% % %   iisub = iisub + 1;
% % %   hca = irf_panel('e PA e64 deflux lowe 22');  
% % %    eint=E(22,:);
% % % %   eint = [1.5*max(scPot1.data) 40000];  
% % %   try
% % %     c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   catch
% % %     c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
% % %   end
% % %   %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
% % %   %hca.YLabel.String = {'Pitchangle','(\circ)'};   
% % %   %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
% % % %   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
% % %   hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};  
% % %   hca.YTick = [45 90 135];   
% % %   caxis(hca,[6.8 7.2])
% % %   colormap(hca,cmap)
% % % end


if 1 % ePDist pa 64 23-----2438.9
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
%   caxis(hca,[7 7.35])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 24-----3192.6
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
%   caxis(hca,[7 7.5])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 25-----4179.2
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
%   caxis(hca,[7 7.7])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 26-----5470.7
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
%   caxis(hca,[7 7.7])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 27-----7161.3
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
%   caxis(hca,[6.7 7.7])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 28-----9374.3
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
%   caxis(hca,[6.6 7.6])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 29-----12271.1
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 29');  
   eint=E(29,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  caxis(hca,[6.3 7.5])
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   

  colormap(hca,cmap)
end
if 1 % ePDist pa 64 30-----16063.2
  hca = irf_panel('e PA e64 deflux lowe 30');  
   eint=E(30,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
%   caxis(hca,[6.1 7.2])
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'}; 
  hca.YTick = [45 90 135];   

  colormap(hca,cmap)
end
if 1 % ePDist pa 64 31-----21027.1
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
  caxis(hca,[6.05 7])
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 32-----27525
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 32');  
   eint=E(32,:);
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  caxis(hca,[5.5 7])
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
%   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {[num2str((eint(1)+eint(2))/2,'%.1f')],'eV'};   
  hca.YTick = [45 90 135];   

  colormap(hca,cmap)
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
% % % h(1).Title.String = irf_ssub('MMS ? GSM',ic);
% % % set(gcf,'color','w');
irf_adjust_panel_position

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname=['C:\Users\fwd\Desktop\Ti~mor~\M\Formation of the rolling-pin distribution of suprathermal electrons behind dipolarization fronts\Figure2\2k+'];
print(gcf, '-dpdf', [figname '.pdf']);
