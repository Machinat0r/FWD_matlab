% continue after PSDwithE111111.m
%%
ePSD1 = [timeomni, PSDomni(:,7:32)];
Energy1 = energye0(7:32);
ePSD2 = [timeomni, PSDomni(:,7:28)];
Energy2 = energye0(7:28);
ePSD3 = [timeomni, PSDomni(:,28:32)];
Energy3 = energye0(28:32);
[Ne1, Te1] = PSDtoNeTe(Energy1, ePSD1);
[Ne2, Te2] = PSDtoNeTe(Energy2, ePSD2);
[Ne3, Te3] = PSDtoNeTe(Energy3, ePSD3);
irf_plot([Ne1(:,1) Ne1(:,2)], 'color','k', 'Linewidth',0.75);hold on;
irf_plot([Ne2(:,1) Ne2(:,2)], 'color','b', 'Linewidth',0.75);hold on;
irf_plot([Ne3(:,1) Ne3(:,2)], 'color','r', 'Linewidth',0.75);hold on;
hold off;grid off;
set(gca,'Ylim',[0 0.15]);
set(gca,"XTickLabelRotation",0)
%% PSDtoNeTe
function [Ne, Te] = PSDtoNeTe(Energy, ePSD)
% temperature and density
Te=ePSD(:,1);
Ne=ePSD(:,1);
E_tot=ePSD(:,1);
for ii=1:size(ePSD,1)
if length(Energy(:,1))==1
PSD_tmp(:,1)=Energy';
else
PSD_tmp(:,1)=Energy(ii,:)';                            % 第一列 能量 eV
end
PSD_tmp(:,2)=ePSD(ii,2:end)';                          % 第二列  PSD s^3 km^-6
PSD_tmp(:,2)=PSD_tmp(:,2)*1e-18;                       % 第二列  PSD s^3 m^-6
PSD_tmp(:,3)=PSD_tmp(:,1)*4.4204e+12.*PSD_tmp(:,2);    % 第三列，电子速率分布密度 4piv2f(v)
PSD_tmp(:,4)=(3.5176e+11*PSD_tmp(:,1)).^0.5;           % 第四列，电子标量速度坐标
PSD_tmp(:,5)=PSD_tmp(:,3).*PSD_tmp(:,1)*1.6022e-19;    % 第五列，能量密度 焦耳
PSD_tmp(:,6)=(PSD_tmp(:,1)*4.4204e+12./4./pi).*PSD_tmp(:,3);    % 第六列v2*4piv2f(v)
Ne(ii,2)=trapz(PSD_tmp(:,4), PSD_tmp(:,3));
E_tot(ii,2)=trapz(PSD_tmp(:,4), PSD_tmp(:,5));         % J
E_tot(ii,3)=trapz(PSD_tmp(:,4), PSD_tmp(:,6));         % 根均方速率的平方（m2/s2）
end
E_tot(:,2)=0.667*E_tot(:,2)./Ne(:,2)/1.3806e-23;
Te(:,2)=E_tot(:,2)/11605;   % eV
Te(:,3)= E_tot(:,3).*(9.1*10^(-31))./(3*1.38*10^(-23))./(11605)./Ne(:,2); % eV
Ne(:,2)=Ne(:,2)/1e6;   % cm^-3
ePSD(ePSD==0)=NaN;
end