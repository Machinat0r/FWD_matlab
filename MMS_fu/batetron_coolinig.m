clc; close all
clear;
%% consequence
load('thetapad.mat','thetapad');%投掷角
load('energypad.mat','energypad');%能道
load('paddist1.mat','paddistsmooth');%PSD

%% source
load('paddist2.mat','paddistsmooth2');%PSD

ene_chan=energypad';%能道
pad=thetapad;%投掷角
% initial PSD
PSD_ini=paddistsmooth2;

% final PSD
PSD_fin=paddistsmooth;

%% for good-looking plot
% Pitang=Pitang_org;  Pitang(1)=0;  Pitang(end)=pi; 
Pitang=pad./180*pi; 
Pitang(1)=0;  
Pitang(end)=pi;
NEe=length(ene_chan); %12个能段
Eele=ene_chan;
NPa=length(pad);%13个投掷角

PSD_i=transpose(log10(PSD_ini));

Eele_par=zeros(NPa,NEe);
Eele_perp=zeros(NPa,NEe);

for ii=1:NEe
    for jj=1:NPa
        Eele_par(jj,ii)=Eele(ii)*cos(Pitang(jj));  %Eele_par=[9x2]
        Eele_perp(jj,ii)=Eele(ii)*sin(Pitang(jj)); %Eele_per=[9x2]
    end
end
%interpolation for both Eele and PAD
%interp2 =  returns the interpolated values on a refined grid 
%   formed by repeatedly dividing the intervals 2 times in each dimension.

Nintp=2;
Eele_par_intp=interp2(Eele_par,Nintp);%(让值增多) 平行方向能量
Eele_perp_intp=interp2(Eele_perp,Nintp);% 垂直方向的电子能量
PSDele_intp=interp2(PSD_i,Nintp);  
 
 


% af=0.96;
% ab=1.15
af=1.2;
ab=1;
% af=0.95;
% ab=1.1
%First, Fermi accleration, velocity becomes 17/16 of the original one
FA_factor=af; FA_factor=FA_factor^2;   %Notice Ea=1/2 mv^2
Eele_par_Fermi=Eele_par_intp*FA_factor;%费米因素乘平行方向的能量
Eele_perp_Fermi=Eele_perp_intp;% 垂直方向的能量

%Then, betatron accleration, suppose the accleration is multiply 2
Eele_par_Fb=Eele_par_Fermi;%平行方向的费米能量为垂直平行方向的betatron能量
Eele_perp_Fb=Eele_perp_Fermi*ab;%垂直方向的能量乘betatron因素为betatron垂直方向的能量

%return to energy and pitch angle
Eele_Fb=sqrt(Eele_par_Fb.^2+Eele_perp_Fb.^2);%投掷角所对应能段
Pitang_Fb=acos(Eele_par_Fb./Eele_Fb)*(180/pi);%投掷角

%% Begin plot
n_subplots=2;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(10);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 16; ySize = 10; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


PSDele_intp=10.^PSDele_intp;

% Eele_want=zeros(1,12);
% Eele_want=[Eele(7) Eele(8) Eele(9) Eele(10)];  % choose energy for comparison
Eele_want=[ Eele(14) Eele(15) ...
           Eele(16) Eele(17) Eele(18) Eele(19) Eele(20) Eele(21)];  % choose energy for comparison
Eele_want(end)=Eele_want(end);
NEele_want=length(Eele_want);

N_extr=18;
Pitang_extr=linspace(0,180,N_extr+1);
Pitang_extplot=linspace(0,180,N_extr);

Eele_Fb_extra=zeros(NEele_want,N_extr);
PSD_Fb_extra=zeros(NEele_want,N_extr);

Eele_temp=10000000000.0;
for mm=1:NEele_want %能段
for ii=1:N_extr  %投掷角
    for jj=1:length(Pitang_Fb(:,1)) %线性插值得到的能道
        for kk=1:length(Pitang_Fb(1,:)) %线性插值得到的投掷角
            if ( (Pitang_Fb(jj,kk)>Pitang_extr(ii)) & (Pitang_Fb(jj,kk)<Pitang_extr(ii+1)) )
                if (Eele_Fb(jj,kk)>Eele_want(mm))
                    if (Eele_Fb(jj,kk)<Eele_temp)
                       Eele_temp=Eele_Fb(jj,kk);
                       PSD_temp=PSDele_intp(jj,kk); 
                        pause=1;
                    end
                end
            end
         end
    end
    
    Eele_Fb_extra(mm,ii)=Eele_temp;
    PSD_Fb_extra(mm,ii)=PSD_temp;
    
    Eele_temp=10000000000.0;
end
end

clr=colormap(jet(NEele_want));
clr=fliplr(clr);
Pitang=Pitang.*(180/pi);

%% Original distribution
h(i_subplot)=irf_subplot(1,n_subplots,-i_subplot);i_subplot=i_subplot+1;
PSD_ini=transpose(PSD_ini);
for i_E=1:NEele_want
    xx=pad;
    yy=PSD_ini(:,(i_E-1)*1+14);  % corresponding energy
    semilogy(xx,yy, '-', 'LineWidth',1.3, 'color',clr(i_E,:)); % y轴对数坐标
    set(gca,'Xtick',[0 45 90 135 180]);
    set(gca,'Xlim',[0 180]);
    set(gca,'ylim',[1e-2 200]);
%     set(gca,'ytick',[1e3 1e4]);
    hold on;
end
E_legend=cell(1,NEele_want);
for i_legend=1:NEele_want
    E_legend{i_legend}=num2str(Eele(-(i_legend-1)*1+21),3);
end
E_legend=wrev(E_legend);

str=['''' E_legend{1} 'eV'','];
for ii=2:(NEele_want)
    str=[str ['''' E_legend{ii} 'eV'',']];
end
eval(['legend(h(1),' str '''' 'Location' '''' ',' '''' 'Best' '''' ')']);
legend('boxoff');
xlabel('Pitch angle [deg]')
title('origin')
hold off;

%% Final distribution together with Fermi and Betatron acelariton
h(i_subplot)=irf_subplot(1,n_subplots,-i_subplot);i_subplot=i_subplot+1;
PSD_fin=transpose(PSD_fin);
for i_E=1:NEele_want
    xx=pad;
    yy=PSD_fin(:,(i_E-1)*1+14);
    semilogy(xx,yy, '-', 'LineWidth',1.3, 'color',clr(i_E,:)); 
    set(gca,'Xtick',[0 45 90 135 180]);
    set(gca,'Xlim',[0 180]);
    set(gca,'ylim',[1e-2 200]);
    hold on;
end
% hold off;
% h(i_subplot)=irf_subplot(1,n_subplots,-i_subplot);i_subplot=i_subplot+1;
for i_E=1:NEele_want
    xx=Pitang_extplot;
    yy=PSD_Fb_extra(i_E,:);
    semilogy(xx,yy, '--', 'LineWidth',1.3, 'color',clr(i_E,:)); 
    set(gca,'Xtick',[0 45 90 135 180]);
    set(gca,'Xlim',[0 180]);
    set(gca,'ylim',[1e-2 200]);
    hold on;
end
xlabel('Pitch angle [deg]');
title('final')
irf_legend(gca,'af =',[0.1 0.1]);
irf_legend(gca,num2str(af),[0.2 0.1]);
irf_legend(gca,'ab =',[0.1 0.05]);
irf_legend(gca,num2str(ab),[0.2 0.05]);
hold off;

%% save figure
  set(gcf,'render','painters');%矢量图
%   set(gcf,'visible','off');
  figname=['MMS3_PSDmodel_nonaver'];
 % print(gcf, '-dpdf', [figname '.pdf'])
%   print(gcf, '-dpng', [figname '.png'])