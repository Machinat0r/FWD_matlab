clear,clc
close all

%% Load Data
filepath='C:\Matlab\bin\新建文件夹\fwd\Sundry\zwq\GOLD\';
C = {};
 for z=214:245
    A = load([filepath, 'GOLD_L2_TDISK_2020_', num2str(z), '_v03_r01_c01.mat'],'tdisk');
    C = [C, A];%20208月整月进行平均
 end
%%
  ttt=0;
    tt=[];I=0;
for z=1:31
   tt=C{1,z}.tdisk ;
   if size(tt,1)==46 && size(tt,2)==52 && size(tt,3)==68
   ttt=tt+ttt;
   I=I+1;
   end
end
   Tt=ttt/I;
   %将一整月的tdisk进行平均，其中筛掉格式不为46*52*68的数据，否则无法相加
load('C:\Matlab\bin\新建文件夹\fwd\Sundry\zwq\GOLD_L2_TDISK_2020_230_v03_r01_c01.mat');
[~,rank1,rank2,rank3] = size(time_utc);
TimeFlag=cell(rank1,rank2,rank3);
for m=1:rank1
    for n=1:rank2
        for z=1:rank3
            tmpTime1 = transpose(time_utc(:,m,n,z));
        if isempty(strfind(tmpTime1,'*'))
            TimeFlag{m,n,z} = iso2epoch(tmpTime1(1:20));
        else 
           TimeFlag{m,n,z}=NaN;
        end
        end
    end
end
TimeFlag=(cell2mat(TimeFlag)-iso2epoch('2020-08-17T00:00:00Z'))/3600;
%将当天的时间求出来，然后下面进行世界时转成地方时
ravote=15;%地球自转角速度
for m=1:68
   TimeFlag(:,:,m)=TimeFlag(:,:,m)+longitude./ravote; 
end
%timeFlag为地方时

%%
%接下来寻找 早上的时间点是最早可用的时间点，下午点即为最晚可用的点，

T1=zeros(46,52);
T2=zeros(46,52);
Time_id1=zeros(46,52);
Time_id2=zeros(46,52);
sza1=zeros(46,52);
sza2=zeros(46,52);

for m=1:46
for n=1:52 
    [~,index] =sort(TimeFlag(m,n,:));
    flag1 = 0;
    for  p=1:68
        if ~isnan(Tt(m,n,index(p)))
            if flag1 == 0
                T1(m,n) =Tt(m,n,index(p));
                Time_id1(m,n)=TimeFlag(m,n,index(p));
                sza1(m,n)=solar_zenith_angle(m,n,index(p));
                flag1 = 1;
            else
                T2(m,n) =Tt(m,n,index(p));
                Time_id2(m,n)=TimeFlag(m,n,index(p));
                sza2(m,n)=solar_zenith_angle(m,n,index(p));
            end
        end
    end
    
    if isempty(T2(m,n)) && flag1 == 1
        T2(m,n) =T1(m,n);
        Time_id2(m,n)=Time_id1(m,n);
        sza2(m,n)=ssza2(m,n);
    else
        continue
    end    
end
end
%%
  %接下来将早上太晚的点或者晚上太早的点将它们从分析中去掉。同时去掉时间不合适所对应下的天顶角以及温度。
  Time_id1(Time_id11 == 0 && Time_id1 > 8)=nan;
  Time_id2(Time_id2 == 0 && Time_id2 < 16)=nan;
  T1(isnan(Time_id1))=nan;
  sza1(isnan(Time_id1))=nan;
  T2(isnan(Time_id2))=nan;
  sza2(isnan(Time_id2))=nan;
  %接下来将所有的上午点插值到最迟的上午点，并将所有下午点插值到这行像素上最早的下午点。这样就可以得到一个固定的早上和下午时间
Maxmorning=max(Time_id1,[],'all');%最迟的上午的点
Minevening=min(Time_id2,[],'all');%最早的下午的点
% % %  for m=1:46
% % %     for n=1:52 
% % %         TempT=Tt(m,n,:);
% % %         TempTime=TimeFlag(m,n,:);
% % %         TempT2=TempT;
% % %         TempT(isnan(TempT)+isnan(TempTime)~=0)=[];
% % %         TempT=reshape(TempT,length(TempT),1);
% % %         TempTime(isnan(TempT2)+isnan(TempTime)~=0)=[];
% % %         TempTime=reshape(TempTime,length(TempTime),1);
% % %         if length(TempTime)>=2
% % %         Tm(m,n)=interp1(TempTime,TempT,Maxmorning,'spline');
% % %              Te(m,n)=interp1(TempTime,TempT,Minevening,'spline');
% % %         else
% % %             Tm(m,n)=nan;
% % %             Te(m,n)=nan;
% % %         end
% % %     end
% % %  end
% % %       
 for m=1:46
    for n=1:52 
        TempT=Tt(m,n,:);
        TempTime=TimeFlag(m,n,:);
        TempT2=TempT;
        TempT(isnan(TempT)+isnan(TempTime)~=0)=[];
        TempT=reshape(TempT,length(TempT),1);
        TempTime(isnan(TempT2)+isnan(TempTime)~=0)=[];
        TempTime=reshape(TempTime,length(TempTime),1);
        if length(TempTime)>=2
        Tm(m,n)=interp1(TempTime,TempT,Maxmorning,'spline');
             Te(m,n)=interp1(TempTime,TempT,Minevening,'spline');
        else
            Tm(m,n)=nan;
            Te(m,n)=nan;
        end
    end
 end
             