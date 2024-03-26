clear;clc
close all
%% location data
% dataDir = 'E:\THEMIS\'; 
global ParentDir 
ParentDir = 'C:\THEMIS\'; 
TempDir = 'C:\THEMIS\temp\'; 
% dataDir = '/data/themis';
thIds = 'b';
% TT = ['2014-08-16T01:40:00Z'];
TT = '2014-02-18T05:09:00Z';
tmpDay = strrep(TT,'-','');
for ic = thIds
THEMISDownload(tmpDay(1:8),['th',ic],'ssc',TempDir)
end
THEMISDataMove(TempDir,ParentDir)
tint = iso2epoch(TT);

for thId=thIds
  R = [];
  for year=str2num(TT(1:4))
    fullPath = sprintf('%sth%s%sor%sssc%s%02d',...
      ParentDir,thId,filesep,filesep,filesep,year);
    if ~exist(fullPath,'dir'), continue, end
    files = dir(sprintf('%s%sth%s_or_ssc_*_v*.cdf',fullPath,filesep,thId));
   
    if ~isempty(files)
      for iFile=1:length(files)
        fileToRead = [fullPath,filesep,files(iFile).name];
        fprintf('reading %s\n',fileToRead)
        tmpData = spdfcdfread(fileToRead,'CombineRecords',true,'Variable','XYZ_GSE');
        tmpData = tmpData*6371.2; % comvert to kilometers
        tmpEpoch = spdfcdfread(fileToRead,'CombineRecords',true,...
          'KeepEpochAsIs',true,'Variable','Epoch');
        tmpEpoch = irf_time(tmpEpoch,'cdfepoch>epoch');
        R = [R; tmpEpoch tmpData];
        clear tmpData tmpEpoch
      end
    end
  end
  % remove repeating points at month boundary
  ii = find(diff(R(:,1))==0); R(ii,:) = []; 
  eval(['Rth' thId '=R;'])
  fprintf('Rth%s >> mRth.mat\n',thId);
%   if exist('./mRth.mat','file')
%       eval(['save(''mRth'',''Rth' thId ''',''-append'')'])
%   else
%       eval(['save(''mRth'',''Rth' thId ''')'])
%   end
end

%% 计算neutral sheet （AEN）
% Bomni_gsm = irf_get_data(tint+[-3600 3600],'Bxgsm,Bygsm,Bzgsm','omni');
% r = sqrt(Bomni_gsm(2,2)^2+Bomni_gsm(2,3));
%tilt: 磁倾角，由idl种geopack_recalc程序算出
tilt = -16.171933;%01
% tilt = -10.264578;%16
h0 = 12.6/pi; 
% tilt = -13; %倾角
%Y0 = 20.2; % Re
%只画XZ平面，取Y=0
% deltaZ = H0*sin(Blat(2,2));
x = -13:0.1:0; %单位Re
y=0;
dz2NS = -h0 * sind(tilt) * atan(x/5) * (2 * cos(y/6));
% dz2NS = irf_gse2gsm([tint*ones(131,1),x',0*ones(131,1),dz2NS'],-1);
%% Init plot 
c_eval('flag? = find(Rth?(:,1)==iso2epoch(TT));',thIds)
c_eval('Loc? = Rth?(flag?,:);',thIds);

set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 180; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% subplot 1
h(1) = subplot(1,2,1);
if strcmp(thIds,'ae')
plot(Locb(1,2)/6371.2,Locb(1,3)/6371.2,'m.','MarkerSize',15);hold on;
% plot(Loce(1,2)/6371.2,Loce(1,3)/6371.2,'b.','MarkerSize',15);hold on;

plot(Rthb(flagb-400:flagb+1040,2)/6371.2,Rthb(flagb-400:flagb+1040,3)/6371.2,'m--');hold on;
% plot(Rthe(flaga-400:flage+1040,2)/6371.2,Rthe(flaga-400:flage+1040,3)/6371.2,'b--');hold on;
else
%     c_eval("plot(Loc?(1,2)/6371.2,Loc?(1,3)/6371.2,'.','MarkerSize',15);hold on;",thIds);
%     c_eval("plot(Rth?(flag?-4000:flag?+1040,2)/6371.2,Rth?(flag?-4000:flag?+1040,3)/6371.2,'--');hold on;",thIds);
    c_eval("plot(Rth?(:,2)/6371.2,Rth?(:,3)/6371.2,'--');hold on;",thIds);
end
grid on

% plot(Rtha(:,2)/6371.2,Rtha(:,3)/6371.2,'m--');hold on;
% plot(Rthe(:,2)/6371.2,Rthe(:,3)/6371.2,'b--');hold on;

set(gca,'XDir','reverse')%对X方向反转
xlabel('X [R_E]  GSE')
ylabel('Y [R_E]  GSE')

% set(gca,'Ylim',[-6 4], 'ytick',[-5:1:5]);
set(gca,'Xlim',[-60,0]);
% set(gca,'Ylim',[-80,80]);
%% subplot 2
h(2) = subplot(1,2,2);
if strcmp(thIds,'ae')
plot(Locb(1,2)/6371.2,Locb(1,4)/6371.2,'m.','MarkerSize',15);hold on;
% plot(Loce(1,2)/6371.2,Loce(1,4)/6371.2,'b.','MarkerSize',15);hold on;

% plot(Rtha(flaga:flaga+1440,2)/6371.2,Rtha(flaga:flaga+1440,4)/6371.2,'m--');hold on;
% plot(Rthe(flaga:flage+1440,2)/6371.2,Rthe(flaga:flage+1440,4)/6371.2,'b--');hold on;
plot(Rthb(flagb:flagb+1440,2)/6371.2,Rthb(flagb:flagb+1440,4)/6371.2,'m--');hold on;
% plot(Rthe(flaga:flage+1440,2)/6371.2,Rthe(flaga:flage+1440,4)/6371.2,'b--');hold on;

legend(gca,{'THEMIS-A (P5)','THEMIS-E (P4)'});
else
%     c_eval("plot(Loc?(1,2)/6371.2,Loc?(1,4)/6371.2,'.','MarkerSize',15);hold on;",thIds);
%     c_eval("plot(Rth?(flag?-400:flag?+1040,2)/6371.2,Rth?(flag?-400:flag?+1040,4)/6371.2,'--');hold on;",thIds);
    c_eval("plot(Rth?(:,2)/6371.2,Rth?(:,4)/6371.2,'--');hold on;",thIds);
    legend(gca,{'a','b','c','d','e'})
end
% plot(x,dz2NS,'k-.')
grid on

% plot(Rtha(:,2)/6371.2,Rtha(:,3)/6371.2,'m--');hold on;
% plot(Rthe(:,2)/6371.2,Rthe(:,3)/6371.2,'b--');hold on;

set(gca,'XDir','reverse')%对X方向反转
xlabel('X [R_E]  GSE')
ylabel('Z [R_E]  GSE')
% 
% add_magnetopause(h(2),tint);
% add_bowshock(h(2),tint);
% add_Earth(h(2));

% set(gca,'Ylim',[-3.5 1.5]);
set(gca,'Xlim',[-60,0]);



sgtitle(TT)
%% 保存
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
Dir =   'C:\Users\fwd\Desktop\Ti~mor~\M\Electron Fermi acceleration by flow vortex\2-Figures\Figure2\';
figname = [Dir TT(1:10) ' Uequal-ratio'];
% print(gcf, '-dpdf', [figname '.pdf']);   
%% Magnetopaus & Bowshock & Earth function
function [flag_omni,omni]=add_magnetopause(h,tint)
    % flag_omni=1 - using OMNI, flag_omni=0 - using default values
%     tMP=getfield(get(gcf,'userdata'),'t');
    [xMP,yMP,omni]=irf_magnetosphere('mp_shue1998',tint);
    if isempty(xMP)
        flag_omni=0;
        [xMP,yMP,omni]=irf_magnetosphere('mp_shue1998');
    else
        flag_omni=1;
    end
    xMP=[fliplr(xMP) xMP];
    yMP=[fliplr(yMP) -yMP];
    line(xMP,yMP,'parent',h,'linewidth',0.5,'linestyle','-','color','k');
end

function [flag_omni,omni]=add_bowshock(h,tint)
    % flag_omni=1 - using OMNI, flag_omni=0 - using default values
%     t=getfield(get(gcf,'userdata'),'t');
    [xBS,yBS,omni]=irf_magnetosphere('bs',tint);
    if isempty(xBS)
        flag_omni=0;
        [xBS,yBS,omni]=irf_magnetosphere('bs');
    else
        flag_omni=1;
    end
    xBS=[fliplr(xBS) xBS];
    yBS=[fliplr(yBS) -yBS];
    line(xBS,yBS,'parent',h,'linewidth',0.5,'linestyle','-','color','k');
end

function add_Earth(h,flag)
    if nargin == 1
        flag='terminator';
    end
    switch flag
        case 'terminator'
            theta=0:pi/20:pi;
            xEarth=sin(theta);yEarth=cos(theta);
            patch(-xEarth,yEarth,'k','edgecolor','none','parent',h)
            patch(xEarth,yEarth,'w','edgecolor','k','parent',h)
        case 'day'
            theta=0:pi/20:2*pi;
            xEarth=sin(theta);yEarth=cos(theta);
            patch(xEarth,yEarth,'w','edgecolor','k','parent',h)
        otherwise
            error('unknown flag');
    end
end