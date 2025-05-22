%% Compact overview
npanels = 6;

set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');
zoomy = [];

if 0 % e DEF omni
  iisub = iisub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 1 % Te par perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Te');
%   set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
%   set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
%   set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
%   set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % J
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmJ?.x.tlim(tint),gsmJ?.y.tlim(tint),gsmJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x.tlim(tint),gsmVi?.y.tlim(tint),gsmVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % ePDist pa 64
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?perp.x.tlim(tint),gsmE?perp.y.tlim(tint),gsmE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))  
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')
%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
%hca = irf_panel('n'); hca.YLim = [0 11];
%hca = irf_panel('J zoom'); hca.YLim = [-500 800];
%hca = irf_panel('Ve'); hca.YLim = [-900 500];
%hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
%hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
%hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
%hca = irf_panel('Ti'); hca.YLim = [300 700];
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h(ii).FontSize = 12;
end
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
% for ipanel = 1:npanels
%    h(ipanel).Position(1) = h(ipanel).Position(1) + 0.07;
%    h(ipanel).Position(3) = h(ipanel).Position(3) - 0.06; 
% end

%% MMS orbit
% tint=irf.tint('2017-07-20T11:44:00Z/2017-07-20T19:44:05Z');
% mms.db_init('local_file_db','G:\data\mms_db\data');
% h = mms.mms4_pl_conf(tint);
flag = 3;
set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
figname=['location'];
% print(gcf, '-dpdf', [figname '.pdf']);

%% Plot overview figure with focus on ions and electrons, including hpca
switch flag
case 1

ic=1;

npanels = 14;

set(0,'DefaultAxesFontSize',10);
set(0,'DefaultLineLineWidth', 01);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 80; ySize = 85; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (20-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

cmap = 'jet';
h = irf_plot(npanels);
% ic =2;
iisub = 0;
cmap = colormap('jet');
zoomy = [];
% tint = irf.tint('2017-08-23T15:38:30.00Z/2017-08-23T15:39:15.00Z');
tint=irf.tint('2018-09-08T14:51:20.00Z/2018-09-08T14:52:00.00Z');
%

% tint1 = irf.tint('2017-06-22T02:33:11.50Z/2017-06-22T02:33:23.00Z');
%   time1=irf_time('2017-06-22T02:33:11.50Z','utc2epoch')-irf_time('2017-06-22T02:33:09.00Z','utc2epoch');
%   time2=irf_time('2017-06-22T02:33:23.00Z','utc2epoch')-irf_time('2017-06-22T02:33:09.00Z','utc2epoch');

%   tint = irf.tint('2017-06-22T03:01:26.00Z/2017-06-22T03:01:43.00Z');
%   tint1 = irf.tint('2017-06-22T03:01:30.50Z/2017-06-22T03:01:43.00Z');
%   time1=irf_time('2017-06-22T03:01:30.50Z','utc2epoch')-irf_time('2017-06-22T03:01:26.00Z','utc2epoch');
%   time2=irf_time('2017-06-22T03:01:43.00Z','utc2epoch')-irf_time('2017-06-22T03:01:26.00Z','utc2epoch');
 
 
%   tint1 = irf.tint('2017-06-22T02:33:12.00Z/2017-06-22T02:33:20.00Z');
% tint1 = irf.tint('2017-06-22T03:01:31.00Z/2017-06-22T03:01:40.00Z');

if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
%    c_eval('irf_plot(hca,{B?.x.tlim(tint),B?.y.tlim(tint),B?.z.tlim(tint),B?.abs.tlim(tint)},''comp'');',ic);
hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  hca.YLabel.String = {'B','(nT)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',12);
%     hca.YLim = [-38 35];
end
if 0 % Bz
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Bz');
  set(hca,'ColorOrder',mms_colors('1234'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.z.tlim(tint)},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  hca.YLabel.String = {'Bz','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',12);
 hca.YLim = [0 16];
end
if 0 % FPI: ne HPCA: O H He 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n e i');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{ne?,nH?,nHe?,nO?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234')); 
  irf_legend(hca,{'e','H','He','O'},[0.98 0.9],'fontsize',12);  
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
if 0 % ni-nc 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('ns');
  set(hca,'ColorOrder',mms_colors('12'))
   c_eval('NC? = N2_psd?+N3_psd?;',ic)
  c_eval(' Ns?=ni?.tlim(tint1)-NC?.tlim(tint1);',ic)
  c_eval('irf_plot(hca,{ni?,NC?.tlim(tint1)},''comp'');',ic);hold(hca,'on');
  c_eval('irf_plot(hca,{Ns?},''k-'');',ic);hold(hca,'off');
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_i','n_c'},[0.98 0.9],'fontsize',12);
  grid(hca,'off');
end
if 0 % ni h c select 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('ni ');
  set(hca,'ColorOrder',mms_colors('1234'))
%   c_eval('NC? = N2_psd?+N3_psd?;',ic)
  c_eval('irf_plot(hca,{ne?,ni?,N1_psd?,N2_psd?},''comp'');',ic)
%     c_eval('irf_plot(hca,{N2_psd?,N3_psd?},''comp'');',ic)
  hca.YLabel.String = {'ni ','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  grid(hca,'off');
%   irf_legend(hca,[num2str(energyrange(1),'%.0f') '<E<' num2str(energyrange(2),'%.0f')],[0.98 0.9],'fontsize',12);
irf_legend(hca,{'n_e','n_i',' ni_H',' ni_C'},[0.08 0.98],'fontsize',12);
% irf_legend(hca,{' ni_C1',' ni_C2'},[0.08 0.98],'fontsize',12);
end
if 0 % ni 1 2 select 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('ni 1 2');
  set(hca,'ColorOrder',mms_colors('123'))
%   c_eval('irf_plot(hca,{ni?,N1_psd?,N2_psd?},''comp'');',ic)
    c_eval('irf_plot(hca,{N2_psd?,N3_psd?},''comp'');',ic)
  hca.YLabel.String = {'ni1/2 ','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('123'))  
%   irf_legend(hca,[num2str(energyrange(1),'%.0f') '<E<' num2str(energyrange(2),'%.0f')],[0.98 0.9],'fontsize',12);
% irf_legend(hca,{'ni',' ni_H',' ni_C'},[0.08 0.98],'fontsize',12);
irf_legend(hca,{' ni_C1',' ni_C2'},[0.08 0.98],'fontsize',12);
end
if 1 % Vi 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  
  c_eval('irf_plot(hca,{gsmVi?.x.tlim(tint),gsmVi?.y.tlim(tint),gsmVi?.z.tlim(tint)},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',12);   
%    hca.YLim = [-800 1000];
end
if 0 % Vi perp12// 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('ViB ');
  set(hca,'ColorOrder',mms_colors('xyza'))
  
  c_eval('irf_plot(hca,{gsmVi?fac.x.tlim(tint),gsmVi?fac.y.tlim(tint),gsmVi?fac.z.tlim(tint)},''comp'');',ic);
  hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'Vi_{\perp 1}','Vi_{\perp 2}','Vi_{||}'},[0.08 0.98],'fontsize',12);   
%    hca.YLim = [-150 500];
end
if 0 % Vi cold 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vic ');
  set(hca,'ColorOrder',mms_colors('xyza'))
  
  c_eval('irf_plot(hca,{gsmV2_psd?fac.x.tlim(tint),gsmV2_psd?fac.y.tlim(tint),gsmV2_psd?fac.z.tlim(tint)},''comp'');',ic);
  hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_iC','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'VC_{\perp 1}','VC_{\perp 2}','VC_{||}'},[0.08 0.98],'fontsize',12);   
%    hca.YLim = [-150 500];
end
if 0 % Vi hot 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vih ');
  set(hca,'ColorOrder',mms_colors('xyza'))
  
  c_eval('irf_plot(hca,{gsmV1_psd?fac.x.tlim(tint),gsmV1_psd?fac.y.tlim(tint),gsmV1_psd?fac.z.tlim(tint)},''comp'');',ic);hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_iH','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'VH_{\perp 1}','VH_{\perp 2}','VH_{||}'},[0.08 0.98],'fontsize',12);   
%    hca.YLim = [-150 500];
end
if 0 % Vix perp select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi\perpx');
  set(hca,'ColorOrder',mms_colors('1234'))

  c_eval('irf_plot(hca,{gsmVi?perp.x.tlim(tint),gsmV1_psd?perp.x.tlim(tint1),gsmV2_psd?perp.x.tlim(tint1),gsmV3_psd?perp.x.tlim(tint1)},''comp'');',ic);
%    c_eval('irf_plot(hca,{V2_psd?.x.tlim(tint),V3_psd?.x.tlim(tint),gsmVExB?.x.tlim(tint)},''comp'');',ic);
%    c_eval('irf_plot(hca,{VOp?.x.tlim(tint)},''k*'');',ic);
%     c_eval('irf_plot(hca,{VHepp?.x.tlim(tint)},''r*'');',ic);
%   legend(hca,'xuyin')
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i\perpx ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'VT','VH','VC1','VC2'},[0.08 0.98],'fontsize',10);   
%   legend(hca,{'VT','VH','VC','VExB','Ohpca'},'fontsize',6);  
%   legend(hca,'boxoff');  
end
if 0 % Viy perp select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Viy');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{gsmVi?perp.y.tlim(tint),gsmV1_psd?perp.y.tlim(tint1),gsmV2_psd?perp.y.tlim(tint1),gsmV3_psd?perp.y.tlim(tint1)},''comp'');',ic);
%    c_eval('irf_plot(hca,{V2_psd?.y.tlim(tint),V3_psd?.y.tlim(tint),gsmVExB?.y.tlim(tint)},''comp'');',ic);
%    c_eval('irf_plot(hca,{VOp?.y.tlim(tint)},''k*'');',ic);
%   c_eval('irf_plot(hca,{VHepp?.y.tlim(tint)},''r*'');',ic);
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i\perpy ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
%     irf_legend(hca,{'VT','VH','VC','VExB'},[0.98 0.95],'fontsize',10);      
end
if 0 % Viz perp select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Viz');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{gsmVi?perp.z.tlim(tint),gsmV1_psd?perp.z.tlim(tint1),gsmV2_psd?perp.z.tlim(tint1),gsmV3_psd?perp.z.tlim(tint1)},''comp'');',ic);
%    c_eval('irf_plot(hca,{V2_psd?.z.tlim(tint),V3_psd?.z.tlim(tint),gsmVExB?.z.tlim(tint)},''comp'');',ic);
%    c_eval('irf_plot(hca,{VOp?.z.tlim(tint)},''k*'');',ic);
%   c_eval('irf_plot(hca,{VHepp?.z.tlim(tint)},''r*'');',ic);
  %   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i\perpz ','(km/s)'};
 set(hca,'ylim',[-300 350]);
  set(hca,'ColorOrder',mms_colors('1234'))
  
%   irf_legend(hca,{'VT','VH','VC','VExB'},[0.98 0.95],'fontsize',10);     
end
if 0 % Vi par select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi par select');
  set(hca,'ColorOrder',mms_colors('1234'))
%   c_eval('irf_plot(hca,{gsmVi?par.tlim(tint),gsmV1_psd?par.tlim(tint1),gsmV2_psd?par.tlim(tint1),gsmV3_psd?par.tlim(tint1)},''comp'');',ic);
  hold(hca,'on');
 Y=linspace(time1,time2,10000)';
 X1=interp1(irf_time(double(gsmV2_psd1par.tlim(tint1).time.epoch),'ttns2epoch')-irf_time('2017-06-22T03:01:26.00Z','utc2epoch'),gsmV2_psd1par.tlim(tint1).data,Y);
 X2=interp1(irf_time(double(gsmV3_psd1par.tlim(tint1).time.epoch),'ttns2epoch')-irf_time('2017-06-22T03:01:26.00Z','utc2epoch'),gsmV3_psd1par.tlim(tint1).data,Y);
%   area(linspace(0,1797.7,size(C1,1))', C1,'facecolor','r');
% area(hca,linspace(6.55,17.89,size(gsmV2_psd1par.tlim(tint1).data,1))',gsmV2_psd1par.tlim(tint1).data,'facecolor','g');hold(hca,'on');
% area(hca,linspace(6.55,17.89,size(gsmV3_psd1par.tlim(tint1).data,1))',gsmV3_psd1par.tlim(tint1).data,'facecolor','b');hold(hca,'on');
area(hca,Y',X1,'facecolor','g');hold(hca,'on');
area(hca,Y',X2,'facecolor','b');hold(hca,'on');

c_eval('irf_plot(hca,{gsmVi?par.tlim(tint),gsmV1_psd?par.tlim(tint1),gsmV2_psd?par.tlim(tint1),gsmV3_psd?par.tlim(tint1)},''comp'');',ic);
  hold(hca,'on');
M1=nanmean(gsmV1_psd1par.tlim(tint1).data)
M2=nanmean(gsmV2_psd1par.tlim(tint1).data)
M3=nanmean(gsmV3_psd1par.tlim(tint1).data)
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
     grid(hca,'off');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  irf_legend(hca,{'VT','VH','VC1','VC2'},[0.05 0.95],'fontsize',10);  
  irf_legend(hca,{'VC1=',num2str(M2(1),3),'VC2=',num2str(M3(1),3)},[0.05 0.05],'color','k','fontsize',10); 
  hca.YLim = [-380 300];  
end
if 0 % Vi par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % VC par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('VC par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmV2_psd??par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{C,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % VC t
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('VCt par');
  set(hca,'ColorOrder',mms_colors('yx'))
  c_eval('irf_plot(hca,{V2_psd?.abs.tlim(tint),V3_psd?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{C}','(km/s)'};
   irf_legend(hca,{'VC1','VC2'},[0.08 0.95],'fontsize',10);
  %hca.YLim = [-1100 1100];  
end
if 0 % VC-VExB
    c_eval('VC?=V2_psd?-gsmVExB?.resample(V2_psd?);',ic);
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('VC-');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{VC?.x.tlim(tint),VC?.y.tlim(tint),VC?.z.tlim(tint),VC?.abs.tlim(tint)},''comp'');',ic);hold(hca,'on');
  %c_eval('irf_plot(hca,{VOp?.x.tlim(tint)},''m'');',ic);
%   c_eval('irf_plot(hca,{VOp?.x.tlim(tint)},''k*'');',ic);
%   legend(hca,'xuyin')
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'VC- ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',10);   
%   legend(hca,{'VT','VH','VC','VExB','Ohpca'},'fontsize',6);  
%   legend(hca,'boxoff');  
end
if 1 % Ve 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',12);  
  grid(hca,'off');
end
if 0 % Ve EXB x 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vex');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVExB?.x.tlim(tint)},''comp'');',ic);hold(hca,'on');
%   legend(hca,'xuyin')
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_Ex ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'Ve','VExB'},[0.08 0.98],'fontsize',10);   
%   legend(hca,{'VT','VH','VC','VExB','Ohpca'},'fontsize',6);  
%   legend(hca,'boxoff');  
end
if 0 % Ve EXB y 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Viy');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVe?.y.tlim(tint),gsmVExB?.y.tlim(tint)},''comp'');',ic);hold(hca,'on');
%   legend(hca,'xuyin')
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_ey ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
%   irf_legend(hca,{'VT','VH','VC','VExB'},[0.08 0.98],'fontsize',10);   
%   legend(hca,{'VT','VH','VC','VExB','Ohpca'},'fontsize',6);  
%   legend(hca,'boxoff');  
end
if 0 % Ve EXB z 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vez');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVe?.z.tlim(tint),gsmVExB?.z.tlim(tint)},''comp'');',ic);hold(hca,'on');
%   legend(hca,'xuyin')
%   set(hca,'marker','o');
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_iz ','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'VT','VH','VC','VExB'},[0.08 0.98],'fontsize',10);   
%   legend(hca,{'VT','VH','VC','VExB','Ohpca'},'fontsize',6);  
%   legend(hca,'boxoff');  
end
if 0 % ve perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
end
if 0 % Ve par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 1 % E perppar
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?fac.x.tlim(tint),gsmE?fac.y.tlim(tint),gsmE?fac.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.08 0.98],'fontsize',12);
  irf_zoom(hca,'y')
  grid(hca,'off');
end
if 0 % E xyz
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Exyz');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x.tlim(tint),gsmE?.y.tlim(tint),gsmE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.08 0.98],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?perp.x,gsmE?perp.y,gsmE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{gsmJ?.x.tlim(tint),gsmJ?.y.tlim(tint),gsmJ?.z.tlim(tint)},''comp'');',ic)
 c_eval('irf_plot(hca,{gsmJ?fac.x.tlim(tint),gsmJ?fac.y.tlim(tint),gsmJ?fac.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
%   irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12); 
irf_legend(hca,{'J_{\perp 1}','J_{\perp 2}','J_{||}'},[0.08 0.98],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % J curlometer
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J cur');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmJcurl.x.tlim(tint),gsmJcurl.y.tlim(tint),gsmJcurl.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hold(hca,'on');
   plot(hca,[-1000,1000],[0,0],'k--');hold(hca,'off');
   hca.YLabel.String = {'J cur','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.95],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % E.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi');
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{EdotJ?par.tlim(tint),EdotJ?perp.tlim(tint),EdotJ?.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E.J','(nW/m^3)'};
% % %   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E.Jpar','E.Jperp','E.J'},[0.08 0.98],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
  grid(hca,'off');
end
if 0 % E'.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E`.J fpi');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{RedotJ?.tlim(tint)});',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E`.J','(nW/m^3)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Te par perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Te');
% % %   set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'Te','(eV)'};
% % %   set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.08 0.98],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
  grid(hca,'off');
end
if 0 % Ti
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ti');
%   set(hca,'ColorOrder',mms_colors('12'))
%   c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
% irf_plot(hca,[Ti1_para(:,1) (Ti1_para(:,2)+2*Ti1_perp(:,2))/3], 'color','k', 'Linewidth',0.75);
% irf_plot(hca,[Ti_perp(:,1) Ti_perp(:,2)], 'color','k', 'Linewidth',0.75); hold(hca,'on');
% irf_plot(hca,[Ti_para(:,1) Ti_para(:,2)], 'color','r', 'Linewidth',0.75); hold(hca,'off');
c_eval('irf_plot(hca,{Ti?},''color'',''k'',''linewidth'',1.5);',ic);hold(hca,'on'); 
  hca.YLabel.String = {'Ti','(eV)'};
%   set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'T_{i,\perp}','T_{i,||}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % Ti 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('T_{i}');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{Ti?,facTi?.xx.tlim(tint),(facTi?.yy+facTi?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_{i}','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'Ti','T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
   grid(hca,'off');
end
if 0 % Ti H select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('T_{i}H');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{T1_para_psd?,T1_perp_psd?},''comp'');',ic)
  hca.YLabel.String = {'T_{i}C','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||,H}',' T_{\perp,H}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % Ti C select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('T_{i}C');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{T2_para_psd?,T2_perp_psd?},''comp'');',ic)
  hca.YLabel.String = {'T_{i}C','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||,C}',' T_{\perp,C}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 600];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % Ti  para select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('T_{i,||}');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),T2_para_psd?},''comp'');',ic)
  hca.YLabel.String = {'T_{i,||}','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||,T}','T_{||,C}'},[0.98 0.90],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % Ti  perp select
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('T_{i,\perp}');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{(facTe?.yy+facTe?.zz)/2,T2_perp_psd?},''comp'');',ic)
  hca.YLabel.String = {'T_{i,\perp}','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{\perp,T}','T_{\perp,C}'},[0.98 0.90],'fontsize',10);

  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % B curvature 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Curv B');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmCurvB},''comp'');',ic)
  hca.YLabel.String = {'Curv B'};
  set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % curvBradius 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('R_c');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{curvBradius},''comp'');',ic)
  hca.YLabel.String = {'R_c','(km)'};
  set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 5000];
  %hca.YTick
  set(hca,'yscale','log');
       set(hca,'ylim',[5e2 1e5]);
  set(hca,'ytick',[1e3 1e4 1e5]);
  grid(hca,'off');
%   irf_zoom(hca,'y')
end
if 0 % p gyroradius
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('p gyroradius');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{rp?},''comp'');',ic)
  hca.YLabel.String = {'Rp','(km)'};
  set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 5000];
  %hca.YTick

  irf_zoom(hca,'y')
end
if 0 % k^2 Rc/Rg 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('k^2');
  set(hca,'ColorOrder',mms_colors('123'))
   c_eval('k?=curvBradius/rp?;',ic)
  c_eval('irf_plot(hca,{k?},''comp'');',ic);  hold(hca,'on'); 
 c_eval(' n=length(k?);',ic)
 c_eval(' L1=[k?.time.epoch 10*ones(n,1)];',ic);
 c_eval(' L2=[k?.time.epoch 25*ones(n,1)];',ic);
  c_eval('irf_plot(hca,{L1},''k--'');',ic);  hold(hca,'on'); 
 c_eval('irf_plot(hca,{L2},''k--'');',ic);  hold(hca,'off');
 grid(hca,'off');
  hca.YLabel.String = {'k^2'};
  set(hca,'ColorOrder',mms_colors('123'))
%   irf_legend(hca,{'T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 5000];
  %hca.YTick
  set(hca,'yscale','log');
       set(hca,'ylim',[5e-1 1e2]);
  set(hca,'ytick',[ 1e0 1e1 ]);
%   irf_zoom(hca,'y')
end
if 0 % p inertial length
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('p inertial length');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{Lp?},''comp'');',ic)
  hca.YLabel.String = {'Lp','(km)'};
  set(hca,'ColorOrder',mms_colors('12'))
%   irf_legend(hca,{'T_{||}',' T_{\perp}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
%   hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % P 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('P');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{Pt?,Pi?,Pd?,PB?},''comp'');',ic)
  hca.YLabel.String = {'P','nPa'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'P_t',' P_i',' P_d','P_B'},[0.08 0.2],'fontsize',10);
%   irf_legend(hca,{'P_B',' P_{i,||}',' P_{i,\perp}',' P_d'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % P hot
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('PH');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{PtH?,PiH?,PdH?},''comp'');',ic)
  hca.YLabel.String = {'PH','nPa'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'P_{t,H}','P_{i,H}',' P_{d,H}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % P cold
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('PC');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{PtC?,PiC?,PdC?},''comp'');',ic)
  hca.YLabel.String = {'PC','nPa'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'P_{t,C}','P_{i,C}',' P_{d,C}'},[0.08 0.2],'fontsize',10);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [0 5000];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % Q
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('sqrtQ');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{Q?},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 1 % Le
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Le');
%   t=gseR1.time.epoch;
  
% % %   set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{Le?},''comp'');',ic);hold(hca,'on');
%    plot(hca,[-100 100], [39.1519 39.0810],'r--');hold(hca,'off');
  hca.YLabel.String = {'Le (km)'};
  hca.YLim = [0 27];
% % %   set(hca,'ColorOrder',mms_colors('12'))   
   grid(hca,'off');
end
if 0 % Ve perp par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint),gsmVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % gradPe
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gsmGradPe.x*1e3,gsmGradPe.y*1e3,gsmGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % e DEF omni 32
  iisub = iisub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist pa 64 
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux 1 ');  
%    eint = [100,400]; 
   eint = [5000,270000]; 
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
  hca.YLabel.String = {['5-27'],'keV'};
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
if 1 % ePDist pa 64 2
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 2');  
   eint = [43,283];  
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 3
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 3');  
   eint = [300,2000];  
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
 
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
%   hca.YLabel.String = {[num2str(eint(1)) '-' num2str(floor(eint(2)/100))],'eV'};  
hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64 4
  iisub = iisub + 1;
  hca = irf_panel('e PA e64 deflux lowe 4');  
   eint = [3000,27525];  
%   eint = [1.5*max(scPot1.data) 40000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end



%
if 0 %low

%
 iisub = iisub + 1;
  hca = irf_panel('e energy low');

specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
specrec_p_elow.p=energy_low.data;%data matrix
specrec_p_elow.f_label='';
specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
irf_spectrogram(h(iisub),specrec_p_elow)
irf_legend(hca,[ '0<E<200' ],[0.99 0.90],'color',0*[1 1 1]) 
hca.YLabel.String = {'\theta_{low,e}','(\circ)'};   
hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 %mid
 iisub = iisub + 1;
  hca = irf_panel('e energy mid');
specrec_p_emid=struct('t',irf_time(energy_mid.DEPEND_0.data,'ttns>epoch'));
specrec_p_emid.f=transpose(energy_mid.DEPEND_1.data(1,1:30));%energy levels
specrec_p_emid.p=energy_mid.data;%data matrix
specrec_p_emid.f_label='';
specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
irf_spectrogram(h(iisub),specrec_p_emid)
irf_legend(hca,[ '200<E<2000'],[0.99 0.90],'color',0*[1 1 1]) 
hca.YLabel.String = {'\theta_{mid,e}','(\circ)'};   
  
hca.YTick = [45 90 135];   
  colormap(hca,cmap)
% irf_colormap('standard')
end
if 0 %high
 iisub = iisub + 1;
  hca = irf_panel('e energy high');
specrec_p_ehigh=struct('t',irf_time(energy_high.DEPEND_0.data,'ttns>epoch'));
specrec_p_ehigh.f=transpose(energy_high.DEPEND_1.data(1,1:30));%energy levels
specrec_p_ehigh.p=energy_high.data;%data matrix
specrec_p_ehigh.f_label='';
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
irf_spectrogram(h(iisub),specrec_p_ehigh)
irf_legend(hca,['2000<E<27000' ],[0.99 0.90],'color',0*[1 1 1]) 
hca.YLabel.String = {'\theta_{high,e}','(\circ)'};   
  
hca.YTick = [45 90 135];   
  colormap(hca,cmap)
% irf_colormap('standard')
end
if 1 % e DEF omni 64
  iisub = iisub + 1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
%   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
% irf_colormap('standard')
%  hca.CLim = [4 8];
end

if 0 % iPDist pa 64 all
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux ');  
  eint = [10 28000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
end

  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

if 0 % iPDist pa 64 1
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe1');  
%   eint = [600,1300];  
  eint = [600,3000];
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
end

  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
   hca.CLim = [4.5 6];
end
if 0 % iPDist pa 64 2
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe2');  
  eint = [1600,1800];  
%   eint = [1600,3000];
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
end

  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % iPDist pa 64 3
  iisub = iisub + 1;
  hca = irf_panel('i PA e64 deflux lowe3');  
  eint = [2000,3000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
end

  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 %directional spectrograms
  iisub = iisub + 1;
  PAlim=[0 15];
  hca = irf_panel('Ei parallel');
%   irf_spectrogram(h(hca),PDistiVI2.pitchangles(Bfast,PAlim).deflux.specrec,'log');
   c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,PAlim).deflux.specrec,''log'');',ic)
%    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).palim(PAlim).eflux.specrec,''log'');',ic)
set(hca,'yscale','log');
   set(hca,'ylim',[1e2 3e4]);
  set(hca,'ytick',[1e2 1e3 1e4]); 
irf_legend(hca,'Ions parallel to B',[0.01,0.9],'FontSize',10);
irf_legend(hca,{[num2str(PAlim(1),'%.0f') '<\theta_B<' num2str(PAlim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
hold(hca,'on');
c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',1);',ic);hold(hca,'on'); 
Line1=irf.ts_scalar(tint,repmat(energyrange2(2), tint.length,1));
Line2=irf.ts_scalar(tint,repmat(energyrange3(2), tint.length,1));
irf_plot(hca,Line1,'k--','LineWidth',1);
irf_plot(hca,Line2,'k--','LineWidth',1);
 hca.YLabel.String = {'i+ par','(eV)'};  
end
if 0 %directional spectrograms
  iisub = iisub + 1;
  PAlim=[30 60];
  hca = irf_panel('Ei 1');
%   irf_spectrogram(h(hca),PDistiVI2.pitchangles(Bfast,PAlim).deflux.specrec,'log');
   c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,PAlim).deflux.specrec,''log'');',ic)
%    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,18).palim(PAlim).deflux.specrec,''log'');',ic)
set(hca,'yscale','log');
   set(hca,'ylim',[5e2 5e3]);
  set(hca,'ytick',[1e2 1e3 1e4]);

% irf_legend(hca,'Ions perpendicular to B',[0.01,0.9],'FontSize',10);
irf_legend(hca,{[num2str(PAlim(1),'%.0f') '<\theta_B<' num2str(PAlim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
hold(hca,'on');
% c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',1);',ic);hold(hca,'on'); 
ylabel(hca,'i+ (eV)');
Line1=irf.ts_scalar(tint,repmat(energyrange2(2), tint.length,1));
Line2=irf.ts_scalar(tint,repmat(energyrange3(2), tint.length,1));
irf_plot(hca,Line1,'k--','LineWidth',1);
irf_plot(hca,Line2,'k--','LineWidth',1);
 hca.YLabel.String = {'i+ per','(eV)'}; 
   hca.CLim = [5 6];
end
if 0 %directional spectrograms
  iisub = iisub + 1;
  PAlim=[120 150];
  hca = irf_panel('Ei 2');
%   irf_spectrogram(h(hca),PDistiVI2.pitchangles(Bfast,PAlim).deflux.specrec,'log');
   c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,PAlim).deflux.specrec,''log'');',ic)
%    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,18).palim(PAlim).deflux.specrec,''log'');',ic)
set(hca,'yscale','log');
   set(hca,'ylim',[5e2 5e3]);
  set(hca,'ytick',[1e2 1e3 1e4]);

irf_legend(hca,'Ions anti-parallel to B',[0.01,0.9],'FontSize',10);
irf_legend(hca,{[num2str(PAlim(1),'%.0f') '<\theta_B<' num2str(PAlim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
hold(hca,'on');
% c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',1);',ic);hold(hca,'on'); 
ylabel(hca,'i+ (eV)');
Line1=irf.ts_scalar(tint,repmat(energyrange2(2), tint.length,1));
Line2=irf.ts_scalar(tint,repmat(energyrange3(2), tint.length,1));
irf_plot(hca,Line1,'k--','LineWidth',1);
irf_plot(hca,Line2,'k--','LineWidth',1);
 hca.YLabel.String = {'i+ anti par','(eV)'};  
  hca.CLim = [5 6];
end
if 0 %directional spectrograms
  iisub = iisub + 1;
  PAlim=[30 60];
  hca = irf_panel('Ei 1');
%   irf_spectrogram(h(hca),PDistiVI2.pitchangles(Bfast,PAlim).deflux.specrec,'log');
   c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,PAlim).deflux.specrec,''log'');',ic)
%    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).palim(PAlim).eflux.specrec,''log'');',ic)
set(hca,'yscale','log');
   set(hca,'ylim',[1e2 3e4]);
  set(hca,'ytick',[1e2 1e3 1e4]); 
irf_legend(hca,'Ions parallel to B',[0.01,0.9],'FontSize',10);
irf_legend(hca,{[num2str(PAlim(1),'%.0f') '<\theta_B<' num2str(PAlim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
hold(hca,'on');
c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',1);',ic);hold(hca,'on'); 
Line1=irf.ts_scalar(tint,repmat(energyrange2(2), tint.length,1));
Line2=irf.ts_scalar(tint,repmat(energyrange3(2), tint.length,1));
irf_plot(hca,Line1,'k--','LineWidth',1);
irf_plot(hca,Line2,'k--','LineWidth',1);
 hca.YLabel.String = {'i+ par','(eV)'};  
end
if 0 % i DEF omni 64
  iisub = iisub + 1;
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic);
 hold(hca,'on'); 
%   c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',2);',ic);hold(hca,'on'); 
  c_eval('irf_plot(hca,{Ti?},''color'',''k'',''linewidth'',1.5);',ic);hold(hca,'on'); 
%    c_eval(' n=length(iPDist?.time.epoch);',ic)
%  c_eval(' Lc=[iPDist?.time.epoch energyrange(2)*ones(n,1)];',ic);
%   c_eval('irf_plot(hca,{Lc},''w--'');',ic);  hold(hca,'off');
%   plot([-1000,1000],[energyrange2(2),energyrange2(2)],'k','linewidth',2);hold(hca,'on');
%   plot(hca,[-1000,1000],[energyrange3(2),energyrange3(2)],'k','linewidth',2);hold(hca,'off');
  set(hca,'yscale','log');
   set(hca,'ylim',[1e1 3e4]);
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);



%   c_eval('irf_plot([iPDist?.time.epoch energyrange(2)],''k'');',ic)  
%   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};   
  colormap(hca,cmap) 
%   hca.CLim = [5 5.8];
  hca.CLim = [4.3 6];
end
if 1 % i DEF omni cold
  iisub = iisub + 1;
  hca = irf_panel('i DEF omni cold');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic);
 hold(hca,'on'); 
%   c_eval('irf_plot(hca,{gsmVExB?.abs.tlim(tint)},''color'',''w'',''linewidth'',2);',ic);hold(hca,'on'); 
%   c_eval('irf_plot(hca,{Ti?},''color'',''k'',''linewidth'',1.5);',ic);hold(hca,'on'); 
%    c_eval(' n=length(iPDist?.time.epoch);',ic)
%  c_eval(' Lc=[iPDist?.time.epoch energyrange(2)*ones(n,1)];',ic);
%   c_eval('irf_plot(hca,{Lc},''w--'');',ic);  hold(hca,'off');
%   plot([-1000,1000],[energyrange2(2),energyrange2(2)],'k','linewidth',2.5);hold(hca,'on');
%   plot(hca,[-1000,1000],[energyrange3(2),energyrange3(2)],'k','linewidth',2.5);hold(hca,'off');
  set(hca,'yscale','log');
   set(hca,'ylim',[1e1 3e4]);
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);



%   c_eval('irf_plot([iPDist?.time.epoch energyrange(2)],''k'');',ic)  
%   c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};   
  colormap(hca,cmap) 
%   hca.CLim = [5 5.7];
%   hca.CLim = [4.7 6];
end
if 0 % fi_red entire energy channels set to 0
 iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('fi_red');
  irf_spectrogram(hca,if1D_rmE12.specrec);
  set(hca,'ylim',[-2000 2000]);
   set(hca,'ytick',[-1500 -1000 -500 0  500 1000 1500]);
  hca.YLabel.String = 'v_i|| (km/s)'; 

% hca.CLim = [-2.8 -0.5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
% irf_zoom(h(zoomy),'y')

%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
%hca = irf_panel('n'); hca.YLim = [0 11];
%hca = irf_panel('J zoom'); hca.YLim = [-500 800];
%hca = irf_panel('Ve'); hca.YLim = [-900 500];
% hca = irf_panel('Vix'); hca.YLim = [-100 500];
% hca = irf_panel('Viy'); hca.YLim = [-100 400];
% hca = irf_panel('Viz'); hca.YLim = [-300 350];
%hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
%hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
%hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
%hca = irf_panel('Ti'); hca.YLim = [300 700];

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
% for ii = 1:npanels;
%   h(ii).FontSize = 12;
% end
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
set(gcf,'color','w');
% export_fig figure1.pdf


%% e PAD 32 energy
case 2
ic=1;

% tint = irf.tint('2017-07-25T22:10:00.00Z/2017-07-25T22:10:40.00Z');
cmap = colormap('jet');
energy=ePitch1.depend{1,1}(1,:);
E=zeros(32,2);
for ii=1:32
    E(ii,:)=[energy(ii)-1,energy(ii)+1];  
end

%
npanels = 15;
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
%
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
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint)
irf_zoom(h(zoomy),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

%
npanels = 11;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
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
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

%
npanels = 12;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(3);clf;
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
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

%% i PAD 32 energy
case 3
ic=1;
% tint = irf.tint('2017-07-25T22:10:00.00Z/2017-07-25T22:10:40.00Z');
cmap = colormap('jet');
energy=iPitch1.depend{1,1}(1,:);
E=zeros(32,2);
for ii=1:32
    E(ii,:)=[energy(ii)-1,energy(ii)+1];  
end

%
npanels = 14;
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
%
if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

%
npanels = 10;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
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
  set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

%
npanels = 11;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(3);clf;
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
  set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
%   set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint),gsmB?.abs.tlim(tint)},''comp'');',ic)
%   c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
%    c_eval('irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},''comp'');',ic)
   c_eval("irf_plot(hca,{gsmB?.x.tlim(tint),gsmB?.y.tlim(tint),gsmB?.z.tlim(tint)},'comp');",ic)
  hca.YLabel.String = {'B','(nT)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
set(hca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
h(1).Title.String = irf_ssub('MMS ? GSM',ic);
set(gcf,'color','w');

end

