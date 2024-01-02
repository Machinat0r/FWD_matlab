clc
clear
B_ss_data_path='/Users/wangchaoqi/Downloads/fgm_jno_l3_2018035pc_r1s_v01.sts';
fileID = fopen(B_ss_data_path);
% 定义要查找的字符串
end_object = 'END_OBJECT';
start_data = num2str(2018);
% 初始化变量
last_end_object = 0;
found_start_data = false;
B_ss_data = [];
% 循环读取文件的每一行，直到找到数据开始的位置
while ~found_start_data
    % 读取一行
    line = fgetl(fileID);
    % 查找END_OBJECT在这一行中的所有位置
    k = strfind(line,end_object);
    % 如果找到了，更新最后一次出现的位置
    if ~isempty(k)
        last_end_object = max(k);
    end
    % 检查是否包含2017
    if contains(line,start_data)
        % 找到数据开始的位置
        found_start_data = true;
        % 把这一行转换为数值数组，并添加到数据矩阵中
        B_ss_data = [B_ss_data;str2num(line)];
    end
end

% 跳出说明文字，定位到数据的第一行
fseek(fileID,last_end_object,'bof');

% 继续读取文件的剩余部分，直到文件结束
while ~feof(fileID)
    % 读取一行
    line = fgetl(fileID);
    % 把这一行转换为数值数组，并添加到数据矩阵中
    B_ss_data = [B_ss_data;str2num(line)];
end

% 关闭文件
fclose(fileID);




%% 导入数据及说明文件
%JADE data
load( "/Volumes/ChxDu_2/JEDI_data/2018/35/organized_data.mat")
lbl_text=fileread( "/Volumes/ChxDu_2/JEDI_data/2018/35/JAD_L50_LRS_ELC_ANY_DEF_2018035_V01.LBL");
% load( "/Volumes/ChxDu_2/find_case/2017_073/specEe_high.mat")
% load( "/Volumes/ChxDu_2/find_case/2017_073/specepadhigh_1.mat")
% load( "/Volumes/ChxDu_2/find_case/2017_073/specepadhigh_2.mat")
% load( "/Volumes/ChxDu_2/find_case/2017_073/specepadhigh_3.mat")
%% 处理JADE数据
%获取2017年第78天对应具体日期
date=split(string(datetime(DIM0_UTC(1,:), 'InputFormat', 'uuuu-DDD''T''HH:mm:ss.SSS')));
date=char(date(1));
energy=DIM1_E_DIM1(1,:,1);
energy_eV=energy;%以eV为单位的能道
energy_keV=energy/(1e3);%以keV为单位的能道
%扣取具体时间
for i=1:length(DIM0_UTC)
time2(i,1:13) = strrep(DIM0_UTC(i,:), '2018-035', '');
% MAG_time2(i,1:13) = strrep(MAG_UTC(i,:), '2018-035', '');
MAG_time2(i,1:13) = MAG_UTC(i,9:end);
end

%将字符串进行拼接，以便于进行时间转换
for i=1:length(time2)
    time3(i,1:24)=[date,time2(i,:),'Z'];
    MAG_time3(i,1:24) =[date,MAG_time2(i,:),'Z'];
end

B_data=zeros(length(DIM0_UTC),5);

%组合数据：时间、磁场总值及磁场三分量
for i=1:length(time2)
    B_data(i,1)=iso2epoch(time3(i,:));
    MAG_time_end(i,1)=iso2epoch(MAG_time3(i,:));
     B_data(i,2)=sqrt(MAG_VECTOR_JSSRTP(i,1).^2+MAG_VECTOR_JSSRTP(i,2).^2+MAG_VECTOR_JSSRTP(i,3).^2);
    B_data(i,3:5)=MAG_VECTOR_JSSRTP(i,1:3);
end
time_end=B_data(:,1);
B_data(:,1)=MAG_time_end;

%% 电子能谱
Ee_flux0_1=nanmean(DATA_DIM1,3);%此时的单位是Differential Energy Flux [1/( m^2 sr s     )] 
Ee_flux0_2=nanmean(DATA_SIGMA_DIM1,3);
Ee_flux0_3=nanmean(BACKGROUND_DIM1,3);%此时的单位是1/( m^2 sr s，需要转化为1/s/cm^2/sr/kev,除以能量再除以一个转换因子
Ee_flux0_4=nanmean(BACKGROUND_SIGMA_DIM1,3);
% Ee_flux0= sqrt(Ee_flux0_1.^2);
Ee_flux0= Ee_flux0_1+Ee_flux0_3;
Ee_flux_1=Ee_flux0;
%单位转换
for i=1:64
    Ee_flux_1(:,i)=(Ee_flux0(:,i)/(1e4))./energy_keV(1,i);%从1/(m^2 s sr)转换成1/( cm^2 s sr keV) 
    Ee_flux_2(:,i)=(Ee_flux0(:,i)/(1e4));%从1/(m^2 s sr)转换成1/( cm^2 s sr),用于后面求密度和温度
end

Ee_t_0=time_end;
Ee_flux_low=Ee_flux_1;
%选取对应时间段数据进行绘图
% Tsta='2018-02-04T08:50:00Z'; 
% Tend='2018-02-04T09:40:00Z'; 
Tsta='2018-02-04T09:00:00Z'; 
Tend='2018-02-04T09:50:00Z'; 
T1=iso2epoch(Tsta);
T2=iso2epoch(Tend);
% tint_cut = [iso2epoch(Tsta) iso2epoch(Tend)];
tstart = irf_time(Tsta, 'utc>epoch');
tend = irf_time(Tend, 'utc>epoch');
[~, idx] = min(abs(Ee_t_0(:,1) - T1)); % 找到距离a最近的元素位置
[~, idx1] = min(abs(Ee_t_0(:,1) - T2)); % 找到距离b最近的元素位置
[~, iidx] = min(abs(MAG_time_end(:,1) - T1)); % 找到距离a最近的元素位置
[~, iidx1] = min(abs(MAG_time_end(:,1) - T2)); % 找到距离a最近的元素位置

% TT_ss = linspace(TT1,TT2,size(B_ss_data,1))';%得到的TT矩阵是一个时间矩阵
TT_ss = irf_time([B_ss_data(:,1:6),zeros(size(B_ss_data,1),2)],'doy8>epoch');
B_JSO=[TT_ss,B_ss_data(:,8),B_ss_data(:,9),B_ss_data(:,10)];% nT，磁场的三个分量BX BY BZ，这是一天的磁场数据
B_plot=irf_tlim(B_JSO,tstart, tend);%得到所选时间段内的磁场数据B_plot
%% |B| 算总磁场
Bt=[];
Bt(:,2)=sqrt(B_plot(:,2).^2+B_plot(:,3).^2+B_plot(:,4).^2);
Bt(:,1)=B_plot(:,1);

%画0.1-30keV的电子能谱（1-54列）
Ee_t_low=Ee_t_0(idx:idx1,1);
% Ee_t=time_end(idx:idx1,1);
energycut=[1:54];
Ee_p_low=Ee_flux_low(idx:idx1,energycut);
Ee_f_low=energy_keV(energycut);
Ee_low = struct('t',Ee_t_low,'f',Ee_f_low,'p',Ee_p_low);
specEe_low=struct('t',Ee_low.t);
specEe_low.p=Ee_low.p;
specEe_low.f=Ee_low.f;
specEe_low.f_label='';
specEe_low.p_label={' ','cm-2 s-1 sr-1 keV-1'};

%% 利用低能粒子数据计算密度和温度
Electron_caculate_data=[time_end,Ee_flux_1];%输入通量单位为1/(cm^2 s sr keV)
[Te,Ne,ePSD,Ve]=Juno_neTe_spart_keV(energy_keV,Electron_caculate_data);
[~, idx2] = min(abs(Ne(:,1) - T1)); % 找到距离a最近的元素位置
[~, idx3] =min(abs(Ne(:,1) - T2)); % 找到距离b最近的元素位置
Ne_Te_cut=[idx2:idx3];



%% 算低能部分电子投掷角
%先换单位
for i=1:64
    PA_flux_1_0(:,i,:)=DATA_DIM1(:,i,:)/1e4/energy_keV(1,i);
    PA_flux_2_0(:,i,:)=DATA_SIGMA_DIM1(:,i,:)/1e4/energy_keV(1,i);
    PA_flux_3_0(:,i,:)=BACKGROUND_DIM1(:,i,:)/1e4/energy_keV(1,i);
    PA_flux_4_0(:,i,:)=BACKGROUND_SIGMA_DIM1(:,i,:)/1e4/energy_keV(1,i);
end
x1=PA_flux_1_0;
x2=PA_flux_2_0;
x3=PA_flux_3_0;
x4=PA_flux_4_0;
x5=DIM3_PITCH_ANGLES_DIM1;
x6=energy_keV;
x7=time_end;
x8=tint_cut;
specepadlow_1=Calculate_JADE_PA(x1,x2,x3,x4,x5,[1:44],x7,x8);
specepadlow_2=Calculate_JADE_PA(x1,x2,x3,x4,x5,[45:48],x7,x8);
specepadlow_3=Calculate_JADE_PA(x1,x2,x3,x4,x5,[49:50],x7,x8);
specepadlow_4=Calculate_JADE_PA(x1,x2,x3,x4,x5,[51:54],x7,x8);
energylow_1=nanmean(energy_keV(1,42:44),2);%获取所选能道切片的平均值
energylow_2=nanmean(energy_keV(1,45:48),2);
energylow_3=nanmean(energy_keV(1,49:50),2);
energylow_4=nanmean(energy_keV(1,51:54),2);

load('/Volumes/ChxDu_2/JEDI_data/2018/35/JED_data_180.mat')%这是低精度的数据
load('/Volumes/ChxDu_2/JEDI_data/2018/35/specEe_high.mat')
JED_data_270=JED_data_180;
%% 处理JEDI数据
time_180 = JED_data_180.t;
time_270 = JED_data_270.t;
energy180 = JED_data_180.energy;
energy_180 = JED_data_180.energy;
pitch_angle = JED_data_180.pa;
flux_180_0 = JED_data_180.flux;
flux_270_0 = JED_data_270.flux;
flux_mean_1=squeeze(nanmean(flux_180_0,2));
flux_mean_2=squeeze(nanmean(flux_270_0,2));

f180 = irf_tlim([time_180,flux_mean_1],tint_cut);
f270 = irf_tlim([time_270,flux_mean_2],tint_cut);
f_comb = zeros(size(f180,1),9);
for i0 = 1:size(f180,1)
    for j=2:9
        f_comb(i0,1)=f180(i0,1);
        f_comb(i0,j) = nanmean([f180(i0,j),f270(i0,j)]);
    end
end

energycut1=[3:3];
energycut2=[4:4];
energycut3=[5:5];
energycut4=[6:7];

specepadhigh_1=Calculate_JEDI_PA(JED_data_180,JED_data_270,energycut1,tint_cut);
specepadhigh_2=Calculate_JEDI_PA(JED_data_180,JED_data_270,energycut2,tint_cut);
specepadhigh_3=Calculate_JEDI_PA(JED_data_180,JED_data_270,energycut3,tint_cut);
specepadhigh_4=Calculate_JEDI_PA(JED_data_180,JED_data_270,energycut4,tint_cut);

energyhigh_1=nanmean(energy_180(1,energycut1),2);%能道;
energyhigh_2=nanmean(energy_180(1,energycut2),2);
energyhigh_3=nanmean(energy_180(1,energycut3),2);
energyhigh_4=nanmean(energy_180(1,energycut4),2);


%% plot
n_subplots=10;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters');
xSize = 19; ySize = 29; coef=floor(min(600/xSize,600/ySize));
xLeft = (20-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[20 20 400 800]);
colormap(jet);
% % B plot
% B plot
%第一幅图 总磁场
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B_data(iidx:iidx1,1) B_data(iidx:iidx1,2)], 'color','k', 'Linewidth',1.5); hold on;
% irf_plot([B_lmn(idx:idx1,1) B_lmn(idx:idx1,2)], 'color','b', 'Linewidth',1.5); hold on;
% irf_plot([B_lmn(idx:idx1,1) B_lmn(idx:idx1,3)], 'color','r', 'Linewidth',1.5); hold on;
% irf_plot([B_lmn(idx:idx1,1) B_lmn(idx:idx1,4)], 'color','g', 'Linewidth',1.5); hold on;
grid off;
B=B_data(idx:idx1,2);
y=ylabel({'|B|','(nT)'},'fontsize',8);
set(gca,'Ylim',[0.9*min(B) 1.1*max(B)]);
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',16)
B1=B_data(iidx:iidx1,3);
B2=B_data(iidx:iidx1,4);
B3=B_data(iidx:iidx1,5);

%第二幅图 磁场三分量
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B_data(iidx:iidx1,1) B1], 'color','b', 'Linewidth',1.5); hold on;
irf_plot([B_data(iidx:iidx1,1) B2], 'color','r', 'Linewidth',1.5); hold on;
irf_plot([B_data(iidx:iidx1,1) B3], 'color','g', 'Linewidth',1.5); hold on;
irf_plot([B_data(iidx:iidx1,1) B_data(idx:idx1,2)*0],'k--', 'Linewidth',0.8);hold off;
grid off;
y=ylabel({'BRTP','(nT)'},'fontsize',9);
BBB=[B1 B2 B3];
minB=min(min(BBB));
maxB=max(max(BBB));
set(gca,'Ylim',[1.1*minB 1.1*max(B)]);
set(gca,'ColorOrder',[[0 0 0];[0 1 0];[1 0 0];[0 0 1]]);%这一行不同的代号代表着不同的颜色
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',16)
irf_legend(h(2),'Br',[0.15 0.95],'color','b','fontsize',16)
irf_legend(h(2),'B\theta',[0.25 0.95],'color','r','fontsize',16)
irf_legend(h(2),'B\phi',[0.05 0.95],'color','g','fontsize',16)
set(gca,'ytick',[-20 -10  0  10  20 ]);

% %% B plot
% %第一幅图 总磁场
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([Bt(:,1) Bt(:,2)], 'color','k', 'Linewidth',1.5); hold on;
% irf_plot([Bt(:,1)  Bt(:,2)*0],'k--', 'Linewidth',1);hold off;
% grid off;
% y=ylabel({'|B|','(nT)'},'fontsize',8);
% [maxB,I1]=max(Bt(:,2));
% [minB,I1]=min(Bt(:,2));
% A=min(minB,maxB);
% B=max(minB,maxB);
% set(gca,'Ylim',[A B]);
% set(gca,'ColorOrder',[[0 0 0];[0 1 0];[1 0 0];[0 0 1]]);%这一行不同的代号代表着不同的颜色
% irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)
% 
% %第二幅图 磁场三分量
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([B_plot(:,1) B_plot(:,2)], 'color','b', 'Linewidth',1.5); hold on;
% irf_plot([B_plot(:,1) B_plot(:,3)], 'color','g', 'Linewidth',1.5); hold on;
% irf_plot([B_plot(:,1) B_plot(:,4)], 'color','r', 'Linewidth',1.5); hold on;
% %irf_plot([Bt(:,1) Bt(:,2)], 'color','b', 'Linewidth',1); hold on;
% 
% irf_plot([B_plot(:,1) B_plot(:,2)*0],'k--', 'Linewidth',1);hold off;
% grid off;
% y=ylabel({'B','(nT)'},'fontsize',9);
% BB1=B_plot(:,2);BB2=B_plot(:,3);BB3=B_plot(:,4);
% BBB=[BB1 BB2 BB3];
% minB=min(min(BBB));
% maxB=max(max(BBB));
% C=min(minB,maxB);
% D=max(minB,maxB);
% set(gca,'Ylim',[C D]);
% set(gca,'ColorOrder',[[0 0 0];[0 1 0];[1 0 0];[0 0 1]]);%这一行不同的代号代表着不同的颜色
% irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)
% irf_legend(h(2),'Bx',[0.05 0.95],'color','b','fontsize',12)
% irf_legend(h(2),'By',[0.15 0.95],'color','g','fontsize',12)
% irf_legend(h(2),'Bz',[0.25 0.95],'color','r','fontsize',12)




% Ne h(3)
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Ne(Ne_Te_cut,1) Ne(Ne_Te_cut,2)], 'color','k', 'Linewidth',1.5); hold on;
% irf_plot([Ni(Ni_Ti_cut,1)   Ni(Ni_Ti_cut,2)], 'color','g', 'Linewidth',1.5); hold on;
% irf_plot([Bt(:,1) Bt(:,2)*0+Bg],'k--', 'Linewidth',0.5);hold off;
grid off;
ylabel({'Ne','(cm-3)'},'fontsize',8);
set(gca,'yscale','linear');
set(gca,'ColorOrder',[0 0 0]);
% irf_legend(h(3),{'Ne'},[0.1 0.20],'fontsize',16);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',16)



% % Te h(4)
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Te(Ne_Te_cut,1) Te(Ne_Te_cut,2)/1000], 'color','r', 'Linewidth',1.5); hold on;
grid off;
ylabel({'Te','(keV)'},'fontsize',8);
% set(gca,'yscale','log');
% set(gca,'ylim',[1 1000]);
% set(gca,'ytick',[1e0 1e1 1e2 1e3 1e4]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0]]);
% irf_legend(h(4),{'Te'},[0.1 0.25],'color','r','fontsize',16);  
 irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',16)

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(5) hcb5]=irf_spectrogram(h(5),specEe_low);
ylabel({'Electron','energy(keV)'},'fontsize',8)
set(gca,'yscale','log');
% set(gca,'ytick',[1e2 1e3 ]);
% set(gca,'ylim',[0 1000]);
caxis(h(5),[4.5 6.6]);
poscbarx=get(hcb5,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb5,'pos',poscbarx);
set(hcb5,'fontsize',8);
irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',16)
% 

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(6) hcb6]=irf_spectrogram(h(6),specEe_high);
ylabel({'Electron','energy(keV)'},'fontsize',8)
set(gca,'yscale','log');
set(gca,'ytick',[1e2 1e3 ]);
set(gca,'ylim',[0 1000]);
caxis(h(6),[1 4]);
poscbarx=get(hcb6,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb6,'pos',poscbarx);
set(hcb6,'fontsize',8);
irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',16)

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(7) hcb7]=irf_spectrogram(h(7),specepadhigh_1,'log');
set(h(7),'ytick',[45 90 135]);
set(h(7),'ylim',[0 180]);
caxis(h(7),[3.8 4.4]);
% caxis(h(4),[6 9]);
ylabel({'PA (\circ)'},'fontsize',10)
poscbarx=get(hcb7,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb7,'pos',poscbarx);
set(hcb7,'fontsize',8);
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',16)
irf_legend(h(7),[num2str(energyhigh_1),'keV'],[0.02 0.75],'color','k','fontsize',16)


h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(8) hcb8]=irf_spectrogram(h(8),specepadhigh_2,'log');
set(h(8),'ytick',[45 90 135]);
set(h(8),'ylim',[0 180]);
caxis(h(8),[3 3.6]);
% caxis(h(4),[6 9]);
ylabel({'PA (\circ)'},'fontsize',10)
poscbarx=get(hcb8,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb8,'pos',poscbarx);
set(hcb8,'fontsize',8);
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',16)
irf_legend(h(8),[num2str(energyhigh_2),'keV'],[0.02 0.75],'color','k','fontsize',16)
% % 
% % 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(9) hcb9]=irf_spectrogram(h(9),specepadhigh_3,'log');
set(h(9),'ytick',[45 90 135]);
set(h(9),'ylim',[0 180]);
% caxis(h(7),[5.2 5.6]);
% caxis(h(4),[6 9]);
ylabel({'PA (\circ)'},'fontsize',10)
poscbarx=get(hcb9,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb9,'pos',poscbarx);
set(hcb9,'fontsize',8);
irf_legend(h(9),'(i)',[0.99 0.98],'color','k','fontsize',16)
irf_legend(h(9),[num2str(energyhigh_3),'keV'],[0.02 0.75],'color','k','fontsize',16)
% % 
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(10) hcb10]=irf_spectrogram(h(10),specepadhigh_4,'log');
set(h(10),'ytick',[45 90 135]);
set(h(10),'ylim',[0 180]);
% caxis(h(7),[5.2 5.6]);
% caxis(h(4),[6 9]);
ylabel({'PA (\circ)'},'fontsize',10)
poscbarx=get(hcb10,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb10,'pos',poscbarx);
set(hcb10,'fontsize',8);
irf_legend(h(10),'(j)',[0.99 0.98],'color','k','fontsize',16)
irf_legend(h(10),[num2str(energyhigh_4),'keV'],[0.02 0.75],'color','k','fontsize',16)

tint1=irf.tint(Tsta,Tend);%irf.tint语法：irf.tint(起始时间,终止时间）
tint = irf_time(tint1,'epochtt>epoch');
tintS=tint1;
irf_zoom(tintS,'x',h(1:end));
irf_plot_axis_align;

