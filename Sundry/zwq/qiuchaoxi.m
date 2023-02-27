clear,clc
close all
filepath='C:\matgold\';
C = {};
 for z=214:245;
    A = load([filepath, 'GOLD_L2_TDISK_2020_', num2str(z), '_v03_r01_c01.mat'],'tdisk');
    C = [C, A];
 end
  ttt=0;
    tt=[];I=0;
for z=1:31;
   tt=C{1,z}.tdisk ;
   if size(tt,1)==46 && size(tt,2)==52 && size(tt,3)==68;
   ttt=tt+ttt;
   I=I+1;
   end
end
   Tt=ttt/I;
   %将一整月的tdisk进行平均，其中筛掉格式不为46*52*68的数据，否则无法相加
load('C:\matgold\GOLD_L2_TDISK_2020_230_v03_r01_c01.mat');
TimeFlag={};
for m=1:46
    for n=1:52
        for z=1:68
       tmpTime1 = transpose(time_utc(:,m,n,z));
       if isempty(strfind(tmpTime1,'*'));
            TimeFlag{m,n,z} = iso2epoch(tmpTime1(1:20));
       else 
           TimeFlag{m,n,z}=NaN;
    end
    end
    end
end
timeFlag=cell2mat(TimeFlag)-iso2epoch('2020-08-17T00:00:00Z');
timeFlag=timeFlag./3600;
%将当天的时间求出来，然后下面进行世界时转成地方时
ravote=15;%地球自转角速度
for m=1:68;
   timeFlag(:,:,m)=timeFlag(:,:,m)+longitude./ravote; 
end

for m=1:46
    for n = 1:52
temp_SZA = solar_zenith_angle(m,n,:);temp_SZA2 = temp_SZA;
temp_Tt = Tt(m,n,:);
temp_SZA(isnan(temp_SZA) + isnan(temp_Tt) ~= 0)=[];
temp_SZA = reshape(temp_SZA,length(temp_SZA),1);
temp_Tt(isnan(temp_SZA2) + isnan(temp_Tt) ~= 0)=[];
temp_Tt = reshape(temp_Tt,length(temp_Tt),1);

        morning=timeFlag(m,n,:)-7;
        [~,index1] = sort(abs(morning));
        evening=timeFlag(m,n,:)-17;
        [~,index2] = sort(abs(evening));
%         if ~isnan(solar_zenith_angle(m,n,:)) || ~isnan(Tt(m,n,:)) ;
        if solar_zenith_angle(m,n,index1)-solar_zenith_angle(m,n,index2)<=1;
        Tm(m,n)=Tt(m,n,index1(1));%Tm为早间温度
        Te(m,n)=Tt(m,n,index2(1));%Te为晚间温度
        else 
            if length(temp_SZA)>=2
             Tm(m,n)=interp1(temp_SZA,temp_Tt,70,'spline');
             Te(m,n)=interp1(temp_SZA,temp_Tt,70,'spline');
            else
             Tm(m,n)=nan;
             Te(m,n)=nan;
            end
        end
    end
end
%上方为 分别求出早上七点半的温度以及下午五点半的温度





   