clc;
clear;
Tsta='2017-11-15T22:20:00.000Z';%检查的时间段
Tend='2017-11-15T22:30:01.000Z';


% Tsta='2017-11-15T18:20:00.000Z';
% Tend='2017-11-15T18:30:00.000Z';

% Tsta='2017-11-15T05:55:00.000Z';
% Tend='2017-11-15T06:00:00.000Z';

% Tsta='2017-11-15T23:50:00.000Z';
% Tend='2017-11-15T23:59:00.000Z';

% Tsta='2017-11-15T22:20:00.000Z';
% Tend='2017-11-15T22:30:00.000Z';

% Tsta='2017-11-15T20:00:00.000Z';
% Tend='2017-11-15T20:10:00.000Z';

Tmia='/';
tintStr=[Tsta,Tmia,Tend];

Xsvyspec_NAME ='D:\maven\OVERVIEW\mvn_swe_l2_svyspec_20171115_v04_r04.cdf';
Xsvyspec_file= dataobj(Xsvyspec_NAME);%读文件
espectrum=get_variable(Xsvyspec_file,'diff_en_fluxes');%单位为eV/eV cm2 s sr，?

Energy=espectrum.DEPEND_1.data;%能道
T_espec=espectrum.DEPEND_0.data;%时间
T_espec=irf_time(T_espec,'ttns>epoch');

espectrum=[T_espec double(espectrum.data)];
espectrum=irf_tlim(espectrum,tintStr);

Energy=double(fliplr(Energy')); %矩阵转置、能道
espectrum(:,2:end)=fliplr(espectrum(:,2:end));%数组左右翻转（大到小）变（小到大）

%% 相空间密度的计算
ePSD=espectrum;
% Energy2=Energy.^2; %能道的平方 单位为eV
% Energy2=repmat(Energy2,size(espectrum,1),1);
% ePSD(:,2:end)=(ePSD(:,2:end)./Energy2).*0.1593;% s^3 km(-6)


%% plot
% TIME1='2017-11-15T21:28:30.07Z'
% TIME1='2017-11-15T18:25:30.07Z'
% TIME1='2017-11-15T11:10:00.07Z'
TIME2='2017-11-15T22:30:00.07Z'%曲线叠加的截止时间点
AA_t1=iso2epoch(Tsta);
t_dif1=abs(ePSD(:,1)-AA_t1);
At1=find(min(min(t_dif1))==t_dif1);
TIME_sure1=ePSD(At1,1);
cc1=find(ePSD(:,1)==TIME_sure1);
AA_t2=iso2epoch(TIME2);
t_dif2=abs(ePSD(:,1)-AA_t2);
At2=find(min(min(t_dif2))==t_dif2);
TIME_sure2=ePSD(At2,1);
cc2=find(ePSD(:,1)==TIME_sure2);

%TIME_sure1=epoch2iso(TIME_sure1);

YY_PSD1=ePSD(cc1:cc2,2:end);



fn=figure;
set(fn,'Position',[200 200 800 400])
%     h(1)=axes('position',[0.08 0.12 0.3 0.80]); 
%     h(1)=axes('position',[0.30 0.10 0.30 0.80]); 
h(1)=axes('position',[0.10 0.10 0.80 0.80]); 
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 
% plot(h(1),log10(Energy),smooth(log10(YY_PSD1),5),'-','LineWidth',0.8,'color','b');
for i =cc1:cc2   %曲线叠加循环
plot(h(1),(Energy),smooth(YY_PSD1(i,1:end),9),'-','LineWidth',0.8,'color','b');

% plot(h(1),(Energy),YY_PSD1,'-','LineWidth',0.8,'color','b');
hold(h(1),'on');
end
% ylabel(h(1),'log 10(f_e (s^3 km^{-6}))');
% ylabel(h(1),'f_e (s^3 km^{-6})');
ylabel(h(1),'J(cm^{-2} s^{-1} sr^{-1} keV^{-1})');
xlabel(h(1),'Ee (eV)')
set(h(1),'yscale','log');
set(h(1),'xscale','log');
% set(h(1),'xlim',[1 6000]);
% set(h(1),'ylim',[0.0001 100000000])
set(h(1),'xtick',[1e1 1e2 1e3])

irf_legend(gca,num2str(Tsta),[0.1 0.12],'fontsize',10);
irf_legend(gca,'-',[0.37 0.12],'fontsize',10);
irf_legend(gca,num2str(Tend),[0.64 0.12],'fontsize',10);


% titlename=[Tsta(1:4) Tsta(6:7) Tsta(9:10) Tsta(12:13) Tsta(15:16) Tsta(18:19) '-' ...
%            Tend(1:4) Tend(6:7) Tend(9:10) Tend(12:13) Tend(15:16) Tend(18:19)];