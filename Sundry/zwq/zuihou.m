clear,clc
close all
filepath='C:\Matlab\bin\新建文件夹\fwd\Sundry\zwq\GOLD\';
C = {};
 for z=214:245
    A = load([filepath, 'GOLD_L2_TDISK_2020_', num2str(z), '_v03_r01_c01.mat'],'tdisk');
    C = [C, A];
 end
  ttt=0;
    tt=[];I=0;
for z=1:31
   tt=C{1,z}.tdisk ;
   if size(tt,1)==46 && size(tt,2)==52 && size(tt,3)==68;
   ttt=tt+ttt;
   I=I+1;
   end
end
   Tt=ttt/I;
   %将一整月的tdisk进行平均，其中筛掉格式不为46*52*68的数据，否则无法相加
load('C:\Matlab\bin\新建文件夹\fwd\Sundry\zwq\GOLD_L2_TDISK_2020_230_v03_r01_c01.mat');
benjuzhen(1,1,1) = single(nan);
benjuzhen(1,1,2) = single(nan);
for m=1:46
    for n=1:52
        flag = 0;
        sza=solar_zenith_angle(m,n,:)-70;
        Sza=abs(sza);
        [~,index] = sort(Sza);
        
        tmpTime1 = transpose(time_utc(:,m,n,index(1)));
        if isempty(strfind(tmpTime1,'*')) %&& Sza(index(1))<=2
            TimeFlag1 = iso2epoch(tmpTime1(1:20));
        for i = 2:length(index)
            tmpTime2 = transpose(time_utc(:,m,n,index(i)));
            if isempty(strfind(tmpTime2,'*')) && iso2epoch(tmpTime2(1:20)) - TimeFlag1 > 3600*8 %&& abs(Sza(index(2))-Sza(index(1)))<=1;
                time1{:,m,n,1}=tmpTime1;
                time2{:,m,n,1}=tmpTime2;
                T{m,n,:} = Tt(m,n,index([1,i]));
                flag = 1;
%                 T{m,n,1} = tdisk(m,n,index(1));
                break
            end
        end
        end
        
        if flag == 0
            time1{:,m,n,1}=nan;
            time2{:,m,n,1}=nan;
            T{m,n,:} = benjuzhen;
        end
    end
end
T=cell2mat(T);
T1=T(:,:,1);
T2=T(:,:,2);
t=T2-T1;
A=repmat(nanmean(t),46,1);%将一行矩阵复刻为46行
B=t-A;
[~,I]=max(abs(B));
Mt=B(I+[0:size(B,2)-1]*size(B,1));%%%%%%%找出每列中元素最大的并生成矩阵
Mt=max(abs(Mt));
tproxy=t./(2*Mt);%线性化后的温度差，归一化
for a=1:46
    for b=1:52
        time11=time1{1,a,b};
        time22=time2{1,a,b};
        if ~isnan([time11,time22])
        Time1{a,b} = iso2epoch(time11) - iso2epoch('2020-08-17T00:00:00Z');
        Time2{a,b} = iso2epoch(time22) - iso2epoch('2020-08-17T00:00:00Z');
        else
        Time1{a,b} = 0;
        Time2{a,b} = 0;
        end
%         Time1{a,b}=str2num(time11{1,1}(12))*36000+str2num(time11{1,1}(13))*3600+str2num(time11{1,1}(15))*600+str2num(time11{1,1}(16))*60+str2num(time11{1,1}(18))*10+str2num(time11{1,1}(19));
%         Time2{a,b}=str2num(time22{1,1}(12))*36000+str2num(time22{1,1}(13))*3600+str2num(time22{1,1}(15))*600+str2num(time22{1,1}(16))*60+str2num(time22{1,1}(18))*10+str2num(time22{1,1}(19));
    end
end
%把时间转换为秒计算、
ravote=(15)*ones(46,52);%地球自转角速率矩阵化
longitude=double(longitude);
% Time1(cellfun('isempty',Time1)) = {0};
% Time2(cellfun('isempty',Time2)) = {0};%将空元胞数组赋值为0，以便于转换成double类型
for a=1:46
    for b=1:52
        if Time1{a,b}==0
            TIME11(a,b)=NaN;
        else
            TIME11(a,b)=Time1{a,b};
        end
    end
end 
TIME11=TIME11/3600;%秒转换为时
Tlt1=TIME11+longitude./ravote;%转换为地方时
for a=1:46
    for b=1:52
        if Time2{a,b}==0
            TIME22(a,b)=NaN;
        else
            TIME22(a,b)=Time2{a,b};
        end
    end
end  
TIME22=TIME22/3600;%秒转换为时
Tlt2=TIME22+longitude./ravote;%转换为地方时
Dt=Tlt2-Tlt1-12;%将地方时做差 
% Dt=Dt.*(-3.5<Dt);
% Dt=Dt.*(Dt<-0.5);
% Dt(Dt==0)=NaN;%去掉比较离谱的点
% Tlt2(isnan(Dt))=NaN;%   将Dt不合适的的点同步到时间2
% tproxy(isnan(Dt))=NaN;
t2=nanmean(nanmean(Tlt2));%由于t2是常数 所以我取了筛选后的平均值
%开始做拟合参数 fun function f=fun(X,xdata)
%方程为f=2*X1*cos(Ravote*dt/2)*cos(Ravote*t2+Ravote*dt/2-4*xdata-X2)+2*X3*cos(Ravote*dt/2)*cos(Ravote*t2+Ravote*dt/2-3*xdata-X4)
Ravote=15;
for p=1:52;

xdata=(longitude(:,p))';
ydata=(double(tproxy(:,p)))';
dt=-2;
X0=[3;5;6;8];%初始值
% X0=[3.6 30 7.2 30];%初始值

xData = [];yData = [];
for tmp = 1:length(xdata)
if ~isnan(xdata(tmp)) && ~isnan(ydata(tmp))
    xData = [xData xdata(tmp)];
    yData = [yData ydata(tmp)];%将NaN值删掉
end
end
fun = @(X,xdata)2*X(1).*cosd(Ravote.*dt/2).*cosd(Ravote.*t2+Ravote.*dt/2-4.*xdata-X(2))...
    +2.*X(3).*cosd(Ravote*dt/2).*cosd(Ravote.*t2+Ravote.*dt/2-3.*xdata-X(4));
% [x,resnorm]=lsqcurvefit(fun,X0,xData,yData);
if ~isempty(xData)
[x,resnorm]=nlinfit(xData,yData,fun,X0);
else
    x = [nan;nan;nan;nan];resnorm = nan;
end
res{p} = resnorm;
X{p}=x;

% if ~isempty(xData)
%     plot(xData,yData);hold on;
%     plot(xData,fun(x,xData))
% end
end


Q=[];
for c=1:52;

    Q(1,c)=X{c}(1);
    Q(2,c)=X{c}(2);
    Q(3,c)=X{c}(3);
    Q(4,c)=X{c}(4);

end
for c=1:52;
    if Q(1,c)==3;
        Q(1,c)=NaN;Q(2,c)=NaN;Q(3,c)=NaN;Q(4,c)=NaN;
    end
end%将没有拟合的值赋予NaN

%% Figure

figure
subplot(2,2,1),plot(Q(1,:));
for c=1:52;
if abs(Q(4,c))>=24;%将相位大于24的搞成NaN
    Q(4,c)=NaN;
end
end
set(gca,'xtick',linspace(1,53,5));
set(gca,'xticklabel',linspace(-24,24));
xlabel('纬度');
ylabel('K');
title('DE3振幅');

subplot(2,2,2),plot(Q(2,:));
set(gca,'xtick',linspace(1,53,5));
set(gca,'xticklabel',linspace(-70,70,5));
xlabel('纬度');
ylabel('h');
title('DE3相位');

subplot(2,2,3),plot(Q(3,:));
set(gca,'xtick',linspace(1,53,5));
set(gca,'xticklabel',linspace(-70,70,5));
xlabel('纬度');
ylabel('K');
title('DE2振幅');

subplot(2,2,4),plot(Q(4,:));
set(gca,'xtick',linspace(1,53,5));
set(gca,'xticklabel',linspace(-70,70,5));
xlabel('纬度');
ylabel('h');
title('DE2相位');
 