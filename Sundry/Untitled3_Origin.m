clear,clc
load('C:\Matlab\bin\新建文件夹\fwd\Sundry\GOLD_L2_TDISK_2020_300_v03_r01_c01.mat');
for m=1:46
    for n=1:52
        sza=solar_zenith_angle(m,n,:)-70;
        [~,index] = sort(abs(sza));
        T{m,n,:} = tdisk(m,n,index(1:2));
        time(:,m,n,:)=time_utc(:,m,n,index(1:2));
        time1{:,m,n,1}=transpose(time(:,m,n,1));
        time2{:,m,n,1}=transpose(time(:,m,n,2));
    end
end
T=cell2mat(T);
T1=T(:,:,1);
T2=T(:,:,2);
t=T2-T1;
A=repmat(nanmean(t),46,1);
B=t-A;
[~,I]=max(abs(B));
Mt=B(I+[0:size(B,2)-1]*size(B,1));
Mt=repmat(Mt,46,1);
tproxy=t./(2*Mt);%线性化后的温度差
for a=1:46
    for b=1:52
        time11=time1(1,a,b);
        time22=time2(1,a,b);
        Time1{a,b}=str2num(time11{1,1}(12))*36000+str2num(time11{1,1}(13))*3600+str2num(time11{1,1}(15))*600+str2num(time11{1,1}(16))*60+str2num(time11{1,1}(18))*10+str2num(time11{1,1}(19));
        Time2{a,b}=str2num(time22{1,1}(12))*36000+str2num(time22{1,1}(13))*3600+str2num(time22{1,1}(15))*600+str2num(time22{1,1}(16))*60+str2num(time22{1,1}(18))*10+str2num(time22{1,1}(19));
    end
end
%把时间转换为秒计算、
ravote=(7.292e-005)*ones(46,52);%地球自转角速率矩阵化
longitude=double(longitude);
Time1(cellfun('isempty',Time1)) = {0};
Time2(cellfun('isempty',Time2)) = {0};%将空元胞数组赋值为0，以便于转换成double类型
for a=1:46
    for b=1:52
        if Time1{a,b}==0
            TIME11(a,b)=NaN;
        else
            TIME11(a,b)=Time1{a,b};
        end
    end
end     
Tlt1=TIME11+ravote./longitude;%转换为地方时
for a=1:46
    for b=1:52
        if Time2{a,b}==0
            TIME22(a,b)=NaN;
        else
            TIME22(a,b)=Time2{a,b};
        end
    end
end     
Tlt2=TIME22+ravote./longitude;%转换为地方时
Dt=Tlt2-Tlt1-43200;%将地方时做差
longitude(Dt<-14400)=NaN;
Dt(Dt<-14400)=NaN;%去掉比较离谱的点
%开始做拟合参数 fun function f=fun(X,xdata)
%方程为f=2*X1*cos(Ravote*dt/2)*cos(Ravote*t2+Ravote*dt/2-4*xdata-X2)+2*X3*cos(Ravote*dt/2)*cos(Ravote*t2+Ravote*dt/2-3*xdata-X4)
Ravote=ravote(:,1);

xdata=(longitude(:,14))';
ydata=(double(tproxy(:,14)))';
dt=32040;
t2=81133;
X0=[1 1 1 1];
fun = @(X1,X2,X3,X4)(2*X1.*cos(Ravote.*dt/2).*cos(Ravote.*t2+Ravote.*dt/2-4.*xdata'-X2)...
    +2.*X3.*cos(Ravote*dt/2).*cos(Ravote.*t2+Ravote.*dt/2-3.*xdata'-X4));
[X1,X2,X3,X4,resnorm]=lsqcurvefit(fun,X0,xdata,ydata);


 


