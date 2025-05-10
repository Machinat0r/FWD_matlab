clc
clear
ic = 1:4;

% tint = irf.tint('2017-06-11T17:55:17.40Z/2017-06-11T17:55:18.40Z');
% tint = irf.tint('2017-07-17T07:53:05.50Z/2017-07-17T07:53:06.10Z');
% tint = irf.tint('2017-07-17T07:53:15.50Z/2017-07-17T07:53:16.00Z');
tint = irf.tint('2017-07-17T07:53:15.60Z/2017-07-17T07:53:16.00Z');

%% Load datastore
mms.db_init('local_file_db','D:\Works\mms_db\data');
db_info = datastore('mms_db');   
%% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);
%c_eval('tic; gseB?scm = mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint); toc',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval(['gsmB?scm=irf_gse2gsm(gseB?scm);'],ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval(['gsmE?=irf_gse2gsm(gseE?);'],ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
c_eval('[gsmE?par,gsmE?perp] = irf_dec_parperp(gsmB?,gsmE?); gsmE?par.name = ''E par''; gsmE?perp.name = ''E perp'';',ic)
c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
% V E, 144 * [0.70 -0.69 -0.19]
%V_E=132 * [0.92 0.37 -0.09];
 V_E=185 * [0.92 0.37 -0.13];
% V_E=256 * [0.87 0.49 -0.06];
c_eval('vec?=repmat(V_E,size(gsmE?par.data,1),1);',1:4);
c_eval('Vvec?=irf.ts_vec_xyz(gsmE?par.time,vec?);',1:4);
c_eval('[vvec?par,vvec?perp] = irf_dec_parperp(gsmB?,Vvec?);',ic)

%% spacecraft position
% % % Load spacecraft position
% disp('Loading spacecraft position...')
% R = mms.get_data('R_gse',tint);

% if size(R.gseR1,2) == 4
%   c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
%   c_eval(['gsmR?=irf_gse2gsm(gseR?);'],ic);
% else
%   c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
%   c_eval(['gsmR?=irf_gse2gsm(gseR?);'],ic);
% end

%%
% vel=2530;
c_eval('E?para=gsmE?par;',1:4);
% c_eval('vvec?=repmat(vel,size(E?para,1),1);',1:4);
% Tevec=[ones(size(Epar3,1),1) repmat(1/Teav,size(Epar3,1),1)];
c_eval('Phipar?t=irf_integrate(E?para);',1:4)
c_eval('Phipar?=Phipar?t.data.*vvec?par.data;',1:4)
c_eval('Phipar?=irf.ts_scalar(E?para.time,Phipar?);',1:4);
%%
figure('name','Electron/ion hole');
npanels=3;
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 23; ySize = 15; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

cmap = 'jet';
h = irf_plot(npanels);
iisub = 0;
cmap = colormap('jet');
zoomy = [];
isub=1;
if 1 % Bz
   hca=h(isub);isub=isub+1;
  zoomy = [zoomy isub];
 
 % set(hca,'ColorOrder',mms_colors('1234'))
  set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmB1.z,gsmB2.z,gsmB3.z,gsmB4.z},''comp'');')
%   c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
% grid off;
  hca.YLabel.String = {'Bz','(nT)'};
%  set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
 irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.05 0.9],'fontsize',10); 
  grid(hca,'off');
% irf_legend(hca,{'Bx','By','Bz',''},[0.05 0.9],'fontsize',10);  
end
if 0 % Bx  
   hca=h(isub);isub=isub+1;
  zoomy = [zoomy isub];
 
%   set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmB1.x,gsmB2.x,gsmB3.x,gsmB4.x},''comp'');')
%   c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
% grid off;
  hca.YLabel.String = {'Bx','(nT)'};
  %set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.05 0.1],'fontsize',10); 
% irf_legend(hca,{'Bx','By','Bz',''},[0.05 0.9],'fontsize',10);  
end
if 0 % By  
   hca=h(isub);isub=isub+1;
  zoomy = [zoomy isub];
 
%   set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmB1.y,gsmB2.y,gsmB3.y,gsmB4.y},''comp'');')
%   c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
% grid off;
  hca.YLabel.String = {'By','(nT)'};
%   set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
%   irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.05 0.9],'fontsize',10); 
% irf_legend(hca,{'Bx','By','Bz',''},[0.05 0.9],'fontsize',10);  
end
if 1 % Electric field

hca=h(isub);isub=isub+1;
  zoomy = [zoomy isub];
  
 % set(hca,'ColorOrder',mms_colors('1234'))
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{gsmE1par,gsmE2par,gsmE3par,gsmE4par},''comp'');')
   hca.YLim = [-10 10];
  hca.YLabel.String = {'E_{||}','(mv/m)'};  
   grid(hca,'off');
%     set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
%     irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.02 0.05]); 
%     set(hca,'ColorOrder',[0 0 0]);       

%  
end

if 1 % Potential 
    hca=h(isub); isub=isub+1;
    zoomy = [zoomy isub];
%     c_eval('Pt?=gsmE?par.time.epoch;',1:4);
%     c_eval('Ppar?=[Pt? Phipar?];',1:4);
set(hca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
  c_eval('irf_plot(hca,{Phipar1,Phipar2,Phipar3,Phipar4},''comp'');')   
   hca.YLim = [-10 50];
     hca.YLabel.String = {'Potential_{||}','(V)'};  
%     irf_plot(hca,Phipar1,'k'); hold on;    
%     irf_plot(hca,Phipar2,'r'); hold on;
%     irf_plot(hca,Phipar3,'g'); hold on;    
%     irf_plot(hca,Phipar4,'b'); hold on;
%     irf_zoom('x',[cn_toepoch(t1) cn_toepoch(t2)]);
%     set(hca,'ColorOrder',[[0.8 0.5 0.2];[0.4 0.2 0.8]]);
%     irf_legend(hca,{'C3','C4'},[0.02 0.05])
      grid(hca,'off');
end

  
%title(['Electrostatic potential and magnetic field fluctuations  \newline (n_e=',num2str(peaNeav,'%.3f'),' cm^{-3} B_0=',num2str(B0,'%.1f'),' nT)'])
% irf_zoom(h,'x',[Epar3(1,1) Epar3(end,1)])


irf_zoom(h(1:npanels),'x',tint)
% irf_zoom(h(zoomy-1),'y')
% 
%     irf_zoom(hca,'y');
% irf_plot_axis_align
% h(1).Title.String = irf_ssub('MMS ? GSM',ic);
% 
% 
% 
% set(gcf,'color','w');