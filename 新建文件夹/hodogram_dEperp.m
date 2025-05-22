clear
clc

mms.db_init('local_file_db','/Volumes/172.17.190.41/Data/MMS/')


% Tsta='2017-08-12T05:23:30.00Z';   
% Tend='2017-08-12T05:24:00.00Z';

% Tsta='2017-08-12T05:24:00.00Z';   
% Tend='2017-08-12T05:24:10.00Z';

% Tsta='2017-08-12T05:24:10.00Z';   
% Tend='2017-08-12T05:24:20.00Z';

%  Tsta='2017-08-12T05:23:50.00Z';   
%  Tend='2017-08-12T05:24:10.00Z';

% tsta1='2017-06-11T01:59:36Z';   
% tend1='2017-06-11T01:59:44Z';
% tsta='2017-06-11T01:59:36.10Z';   
% tend='2017-06-11T01:59:44.10Z';
tsta='2018-07-06T12:38:20.000Z';   
tend='2018-07-06T12:39:10.000Z';

tint=irf.tint(tsta,tend);
% tint1=irf.tint(tsta1,tend1);
%% Load Data 

N=[0.51,-0.78,0.38]
q=[0 0.95 0.31]
L=cross(N,q)
M=cross(N,L)




% 
% for ic=1:4
% c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tint);',ic);
% c_eval('B?=irf.ts2mat(Bxyz?);',ic);
% % c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
% % c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,5);',ic);
% % c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);
% % c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);
% c_eval('Exyz?_gse=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
% c_eval('Exyz?_gsm=irf_gse2gsm(Exyz?_gse);',ic);
% c_eval('E?=irf.ts2mat(Exyz?_gsm);',ic);
% c_eval('Exyz? = irf.ts_vec_xyz(Exyz?_gse.time,Exyz?_gsm.data);',ic);


for ic=1
% load B


c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);
% lvbo
c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
c_eval('Bbf? = B?_ts.filt(0.086,0.087,dfB?,3);',ic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);

c_eval(['B?_srvy=mms.get_data(''B_gsm_srvy'',tint,?);'],ic);
c_eval(['B?_srvy=irf.ts2mat(B?_srvy);'],ic);
% c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);
c_eval(['Blmn?=irf_convert_fac(Bbf?,B?_srvy,[1,0,0]);'],ic);
% c_eval('Blmn? = B?;',ic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?=irf.ts2mat(E?_ts);'],ic);
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
c_eval('Ebf? = E?_gsm.filt(0.086,0.087,dfE?,3);',ic);
c_eval(['Ebf?=irf.ts2mat(Ebf?);'],ic);

% c_eval('Elmn?=irf_newxyz(Ebf1,L,M,N);',ic);
c_eval(['Elmn?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);


% c_eval('magB? = Bxyz?.abs;',ic);
% c_eval('magBexb? = magB?.resample(Exyz?_gsm.abs);',ic);
% 
% c_eval('Bt?=irf.ts2mat(magB?);',ic);

end


%% Rotate E and B into field-aligned coordinates
% for ic=1:4
%     c_eval('Exyzfac? = irf_convert_fac(Exyz?,Bxyz?,[1 0 0]);',ic);
%     c_eval('Bscmfac? = irf_convert_fac(Bscm?,Bxyz?,[1 0 0]);',ic);
%     c_eval('Rfac? = irf_convert_fac(R?,Bxyz?,[1 0 0]);',ic);
% 
% %     c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
% %     c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,5);',ic);
% %     c_eval(['Bbf?=irf.ts2mat(Bbf?);'],ic);
% %     c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);
%     
% end
% for ic=1:4
%    c_eval('Rfac?_mat=irf.ts2mat(Rfac?);',ic); 
%    c_eval('R?_mat=irf.ts2mat(R?);',ic); 
% end


%% Bandpass filter E and B waveforms



% 
% for ic=1:4
%      c_eval('dfE? = 1/median(diff(Exyz?.time.epochUnix));',ic);
%      c_eval('dfB? = 1/median(diff(Bscm?.time.epochUnix));',ic);
%      c_eval('Exyzhf? = Exyz?_gsm.filt(0.8,1.2,dfE?,3);',ic);
%      c_eval('Exyzhf?_mat=irf.ts2mat(Exyzhf?);',ic);
%      c_eval('Bscmhf? = Bscm?.filt(0.8,1.2,dfB?,3);',ic);
%      c_eval('Bscmhf?_mat=irf.ts2mat(Bscmhf?);',ic);
% end


xSize = 10; ySize = 10; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])





subplot(1,1,1);


plot(Blmn1(:,2),Blmn1(:,3),'k');hold on;
plot(Blmn1(1,2),Blmn1(1,3),'+');hold on;
plot(Blmn1(end,2),Blmn1(end,3),'o');hold off;
% plot(Elmn1(:,2),Elmn1(:,3),'k');hold on;
% plot(Elmn1(1,2),Elmn1(1,3),'+');hold on;
% plot(Elmn1(end,2),Elmn1(end,3),'o');hold off;
% 
% plot(Exyzfachf4_mat(1,2),Exyzfachf4_mat(1,3),'+');hold on;
% plot(Exyzfachf4_mat(end,2),Exyzfachf4_mat(end,3),'o');hold off;
% set(gca,'Ylim',[-10 10],'Xlim',[-10 10]);
% set(gca,'Ylim',[-0.2 0.2],'Xlim',[-0.2 0.2]);
set(gca,'fontsize',10);
xlabel('\deltaB_{\perp1} [nT]');
ylabel('\deltaB_{\perp2} [nT]');
set(gca,'Xlim',[-2 2], 'xtick',[-0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35]);
set(gca,'Ylim',[-2 2], 'ytick',[-0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35]);
% irf_legend(gca,'MMS4',[0.98 0.98]);

% set(gcf,'render','painters');
% % set(gcf,'paperpositionmode','auto')
% figname=['b_pianzhen'];
% print(gcf, '-dpdf', [figname '.pdf']);