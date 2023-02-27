 clear;clc;close all
%--------------------------------------------
%�¶Ȳ�ͬ�����ʿ˺��嵥ɫ���������벨��������
%--------------------------------------------
maxP_x=[];maxP_m=[];%�洢��ֵ��
for T=[1e3,1e5,1e6]% ���÷����¶�(K),500,600,800,1000,1200,1500,2000,2400,3000,4000,5000,6000
l_st = 1e-4;
step = 1e-5;
l= l_st:step:1e3; % ���ò�����Χ�����㲽��

M=plk(l,T);
%M=(c1./(l.^5)./(exp(c2./(l.*T))-1)); % ����ָ���¶ȹ��׷�����
loglog(l,M,'LineWidth',1.4) % ���ƹ��׷���������
maxM=max(M); % �ҳ�ָ���¶������׷�����
[~,M1_index] = min(abs(M(1:1000)-1));
i=find(maxM==M); % �ҷ�ֵ������
hold on

maxLentgh=l_st+step*i; % �ҷ�ֵ����
M1Length = l_st+step*M1_index(1)
stem(maxLentgh,maxM,'--','filled') %��ֱ��ͼ
maxP_x=[maxP_x maxLentgh];
maxP_m=[maxP_m maxM];

maxStr=['(',num2str(T),'k',',',num2str(maxLentgh),'\mum',',',round(num2str(maxM/1e6)),'MW/cm^{2}/\mum',')'];
text(l(i+10),M(i+10),maxStr,'VerticalAlignment','baseline','HorizontalAlignment','left','fontsize',7);
hold on
end

% set(0,'DefaultAxesFontSize',8);
% set(0,'DefaultLineLineWidth', 0.5);
% fn=figure(1);clf;
% set(gcf,'PaperUnits','centimeters')
% xSize = 100; ySize = 50; coef=floor(min(800/xSize,800/ySize));
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',[10 10 xSize*coef ySize*coef])

plot(maxP_x,maxP_m) %����ֵ����

titleStr =('����������ȷֲ�����');
title(titleStr);
% set(gca,'xscale','log');
% set(gca,'yscale','log');
set(gca,'XLim',[1e-4 1e3])
set(gca,'YLim',[1e-3 1e16])
xlabel('\lambda / \mum') % ���������Ƽ���λ
ylabel('M_{b\lambda} / W\cdotcm^{-2}\cdot\mum^{-1}') % ���������Ƽ���λ

%plk.m�ļ�
% �������ʿ˹�ʽ�����������׷���ͨ���ܶȣ��������ȣ�
 
function M = plk(x,T)
% x:����um
% T:�¶�K
% w:���׷���ͨ���ܶȣ�w/cm2/um

c1=3.742.*1e4;
c2=1.4388.*1e4;
M = c1./((x.^5).*((exp(c2./(x.*T)))-1));
end
