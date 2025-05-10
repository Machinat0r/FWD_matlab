% clc
% clear
% ic = 1;
% % tint = irf.tint('2017-07-18T13:03:53.00Z/2017-07-18T13:05:30.00Z');
% tint = irf.tint('2017-07-11T22:34:01.00Z/2017-07-11T22:34:03.00Z');
% mms.db_init('local_file_db','G:\data\mms_db\data');
% db_info = datastore('mms_db');  
% disp('Loading skymaps...')
% c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
% c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
% c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
% c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
% c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 12;
tint = irf.tint('2017-07-06T17:47:06.00Z/2017-07-06T17:47:14.00Z');
 %tint = irf.tint('2017-07-12T11:54:40.00Z/2017-07-12T11:54:43.00Z');
tintZoom = tint;

cmap = 'jet';
[h1,h2] = initialize_combined_plot(npanels,3,4,0.4,'vertical');
ic = 2;
iisub = 0;
cmap = colormap('jet');

isub = 0;
zoomy = [];

if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z,gsmB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
%   set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
%   set(hca,'ColorOrder',mms_colors('12'))    
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
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE?perp.x,gsmE?perp.y,gsmE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % J  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmJ?.x.tlim(tint),gsmJ?.y.tlim(tint),gsmJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 1 % E.J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E.J fpi');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{EdotJ?.tlim(tint),EdotJ?par.tlim(tint),EdotJ?perp.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E.J','(nW/m^3)'};
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  irf_legend(hca,{'E.J','E.Jpar','E.Jperp'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
%   set(hca,'ColorOrder',mms_colors('xyza'))
%   c_eval('irf_plot(hca,{bdryVi?.x.tlim(tint),bdryVi?.y.tlim(tint),bdryVi?.z.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gsmVi?.x.tlim(tint),gsmVi?.y.tlim(tint),gsmVi?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);     
end
if 0 % ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmVe?perp.x.tlim(tint),gsmVe?perp.y.tlim(tint),gsmVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
end
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % Ve perp par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
%   set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint),bdryVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
%   set(hca,'ColorOrder',mms_colors('xyza'))
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
if 1 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

if 1 % Q
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('sqrtQ1');
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{Q1},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
set(hca,'Ylim',[0 0.1], 'ytick',[0 0.02 0.04 0.06 0.08 0.1]);
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
end
if 1 % Q
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('sqrtQ2');
set(hca,'ColorOrder',[[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{Q2},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
set(hca,'Ylim',[0 0.1], 'ytick',[0 0.02 0.04 0.06 0.08 0.1]);
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
end
if 1 % Q
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('sqrtQ3');
set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{Q3},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
set(hca,'Ylim',[0 0.1], 'ytick',[0 0.02 0.04 0.06 0.08 0.1]);
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
end
if 1 % Q
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('sqrtQ4');
set(hca,'ColorOrder',[[0 0 1]]);
  c_eval('irf_plot(hca,{Q4},''comp'');',ic)
%   hca.YLabel.String = {'$$\sqrt{Q}$$'};
  ylabel(hca,'$$\sqrt{Q}$$','Interpreter','latex');
set(hca,'Ylim',[0 0.1], 'ytick',[0 0.02 0.04 0.06 0.08 0.1]); 
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
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

%% Plot single time particle distributions, 1 sc, 4 projections,
% %
% npanels = 12;
% tintZoom = tint;
% 
% cmap = 'jet';
% [h1,h2] = initialize_combined_plot(npanels,3,3,0.4,'vertical');
% iisub = 0;
% cmap = colormap('jet');
% 
% isub = 0;
% zoomy = [];
ic=1:4;
% tintDist = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
tintDist = irf.tint('2017-07-06T17:47:06.00Z/2017-07-06T17:47:14.00Z');

% tintDist = irf.tint('2017-07-18T13:04:53.00Z/2017-07-18T13:04:58.00Z');
c_eval('dist? = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('dist_scm? = ePDist?;',ic)
c_eval('scpot? = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim =20*1e3;
elevlim = 20;
strCMap = 'jet';
% projclim = [-1.5 0.6]; %07.191
projclim = [-3 2]; 
%projclim_int =[-11.6 -9.5];
projclim_int = [-11.4 -8.6];  

c_eval('times = ePDist?.time;',ic)
[tind,~] = times.tlim(tintDist);

doReducedF = 1;
for it_ =1:1:numel(tind);
  it = tind(it_);
  time = times(it);
  it_
  if exist('hmark'); delete(hmark); end
  c_eval('hmark_tmp = irf_pl_mark(h1(?),time.epochUnix,''green''); hmark(?) = hmark_tmp;',1:npanels)
  
%   c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
%   c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
%   c_eval('hatExB = cross(hatE,hatB);',ic)
%   par = hatB;
%   perp1 = hatExB;
%   perp2 = cross(par,perp1);  
  
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  par = hatB;
  newy = [1 0 0]; 
  perp2 = cross(par,newy);
  perp1 = cross(perp2,par);  
  

  timeUTC = time.utc;      
  isub = 1;

%   vectors = {hatExB,'(bxe)xb'; hatE,'bxe'; hatB,'b'};
%   vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  if 0% Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,distccm1,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot1,'vlabel',vlabels);
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
  end 
%   if 1 % B plane 1, slice
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;par;-perp2]; vlabels = {'v_{BxE}','v_{B}','-v_{Bx(ExB)}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
%     hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
%   end 
%   if 1 % B plane 2, slice
%     hca = h2(isub); isub = isub + 1;
%     xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{BxE}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
%     hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
%   end

% 
  if 0 % Perpendicular plane, integrated
     vint = [-inf inf];
      hca = h2(isub); isub = isub + 1; 
    hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    tt=0:pi/100:2*pi;
    xx=8*sin(tt);yy=8*cos(tt);
    plot(xx,yy,'color','k','linewidth',0.75);
    hold (gca,'off');
  end 
 %% MMS1
   vint = [-inf inf];
  if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
%     hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
    [FF1]=mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%     tt=0:pi/100:2*pi;
%     xx=8*sin(tt);yy=8*cos(tt);
%     plot(xx,yy,'color','k','linewidth',0.75);
%     hold (gca,'off');
agyrotropy10{it_,1}=time;
agyrotropy10{it_,2}=FF1;
  end

     vint = [-inf inf];
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp2;par;perp1]; vlabels = {'v_{bxe}','v_{b}','v_{(bxe)xb}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF1]=mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
agyrotropy11{it_,1}=time;
agyrotropy11{it_,2}=FF1;
    
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [par;perp1;perp2]; vlabels = {'v_{b}','v_{(bxe)xb}','v_{bxe}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF1]= mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    agyrotropy12{it_,1}=time;
agyrotropy12{it_,2}=FF1;
  end 
 %% MMS2
   if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
%     hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
    [FF2]=mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%     tt=0:pi/100:2*pi;
%     xx=8*sin(tt);yy=8*cos(tt);
%     plot(xx,yy,'color','k','linewidth',0.75);
%     hold (gca,'off');
agyrotropy20{it_,1}=time;
agyrotropy20{it_,2}=FF2;
  end

     vint = [-inf inf];
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp2;par;perp1]; vlabels = {'v_{bxe}','v_{b}','v_{(bxe)xb}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF2]=mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
agyrotropy21{it_,1}=time;
agyrotropy21{it_,2}=FF2;
    
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [par;perp1;perp2]; vlabels = {'v_{b}','v_{(bxe)xb}','v_{bxe}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF2]= mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    agyrotropy22{it_,1}=time;
agyrotropy22{it_,2}=FF2;
  end 
 %% MMS3 
  if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
%     hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
    [FF3]=mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%     tt=0:pi/100:2*pi;
%     xx=8*sin(tt);yy=8*cos(tt);
%     plot(xx,yy,'color','k','linewidth',0.75);
%     hold (gca,'off');
agyrotropy30{it_,1}=time;
agyrotropy30{it_,2}=FF3;
  end

     vint = [-inf inf];
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp2;par;perp1]; vlabels = {'v_{bxe}','v_{b}','v_{(bxe)xb}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF3]=mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
agyrotropy31{it_,1}=time;
agyrotropy31{it_,2}=FF3;
    
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [par;perp1;perp2]; vlabels = {'v_{b}','v_{(bxe)xb}','v_{bxe}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF3]= mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    agyrotropy32{it_,1}=time;
agyrotropy32{it_,2}=FF3;
  end  
  %% MMS4
   if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
%     hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
    [FF4]=mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%     tt=0:pi/100:2*pi;
%     xx=8*sin(tt);yy=8*cos(tt);
%     plot(xx,yy,'color','k','linewidth',0.75);
%     hold (gca,'off');
agyrotropy40{it_,1}=time;
agyrotropy40{it_,2}=FF4;
  end

     vint = [-inf inf];
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp2;par;perp1]; vlabels = {'v_{bxe}','v_{b}','v_{(bxe)xb}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF4]=mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
agyrotropy41{it_,1}=time;
agyrotropy41{it_,2}=FF4;
    
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [par;perp1;perp2]; vlabels = {'v_{b}','v_{(bxe)xb}','v_{bxe}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
     [FF4]= mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    agyrotropy42{it_,1}=time;
agyrotropy42{it_,2}=FF4;
  end 
%   if 1 % B plane 1, integrated
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;par;-perp2]; vlabels = {'v_{BxE}','v_{B}','-v_{Bx(ExB)}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
%     hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%   end  
%   if 1 % B plane 2, integrated
%     hca = h2(isub); isub = isub + 1;
%     xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{BxE}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
%     hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%   end  
%   vint = 40000*[-1 1];
%   if 1 % Perpendicular plane, integrated, smaller vint
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;perp2;par]; vlabels = {'v_{(bxe)xb}','v_{bxe}','v_{b}'};
%     mms.plot_int_projection(hca,dist_scm1,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);
%     hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%   end 
%   if 1 % B plane 1, integrated, smaller vint
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;par;-perp2]; vlabels = {'v_{BxE}','v_{B}','-v_{Bx(ExB)}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);
%     hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%   end  
%   if 1 % B plane 2, integrated, smaller vint
%     hca = h2(isub); isub = isub + 1;
%     xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(ExB)}','v_{B}','v_{BxE}'};
%     %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%     mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);        
%     hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
%   end  
%   
%   if 0 % Pitchangle distribution
%     hca = h2(isub); isub = isub + 1;
%     plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
%     hca.YScale = 'log'; hca.XScale = 'log';
%     hca.YLabel.String = ['f_e (' ePitch.units ')'];
%     hca.XLabel.String = 'E (eV)';
%     hca.XLim = [50 30000];
%     legend(hca,{'0','90','180'})
%     hca.YTick = 10.^[-4:2];
%     hca.YLim = [1e-4 1e2];
%   end
%   if 0 % invisible
%     hca = h2(isub); isub = isub + 1;
%     hca.Visible = 'off';
%   end
%   if 0 % invisible
%     hca = h2(isub); isub = isub + 1;
%     hca.Visible = 'off';
%   end
%   % ExB plane, with markings
%   if 0
%     vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
%     hca = h2(isub); isub = isub + 1; 
%     xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
%     hca.Title.String = '';
%   end
   if 0 % Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,distscm2,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot2,'vlabel',vlabels);
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
   end 
   vint = [-inf inf];
 if 0 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    tt=0:pi/100:2*pi;
    xx=8*sin(tt);yy=8*cos(tt);
    plot(xx,yy,'color','k','linewidth',0.75);
    hold (gca,'off');
  end 
 
  if 0 % Perpendicular plane, integrated, smaller vint
    vint = 40000*[-1 1];hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm2,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 



 if 0 % Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,distscm3,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot3,'vlabel',vlabels);
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
 end 
   vint = [-inf inf];
 if 0 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    tt=0:pi/100:2*pi;
    xx=8*sin(tt);yy=8*cos(tt);
    plot(xx,yy,'color','k','linewidth',0.75);
    hold (gca,'off');
  end 
 
  if 0 % Perpendicular plane, integrated, smaller vint
    vint = 40000*[-1 1];hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm3,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 





 if 0 % Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
%     mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,dist4,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot4,'vlabel',vlabels);
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
 end 
   vint = [-inf inf];
 if 0 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    hold (gca,'on');
    xyz = [perp1;perp2;par];vlabels = {'v_{bxe}','v_{(bxe)xb}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
    tt=0:pi/100:2*pi;
    xx=8*sin(tt);yy=8*cos(tt);
    plot(xx,yy,'color','k','linewidth',0.75);
    hold (gca,'off');
  end 
 
  if 0 % Perpendicular plane, integrated, smaller vint
   vint = 40000*[-1 1]; hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{bxe}','v_{(bxe)xb}}','v_{b}'};
    mms.plot_int_projection(hca,dist_scm4,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int_2,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 



  h2(1).Title.String = {timeUTC(1:23),h2(1).Title.String};
  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  %cn.print(['e_proj_in_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath 'proj_int_mms1/'])
  set(gcf,'render','painters');
figname=['20170706e1' num2str(it)];
print(gcf, '-dpng', [figname '.png']);
end
 af = ePitch1.depend(1);
af{1};
 af=af{1};
energy_pad = af(1:2,:);
save('energy_pad6','energy_pad')
save('agyrotropy610.mat','agyrotropy10')
save('agyrotropy620.mat','agyrotropy20')
save('agyrotropy630.mat','agyrotropy30')
save('agyrotropy640.mat','agyrotropy40')
