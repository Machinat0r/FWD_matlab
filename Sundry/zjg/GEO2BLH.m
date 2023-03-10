clear; clc
load('/Users/fwd/Documents/MATLAB/新建文件夹/fwd/Sundry/zjg/GEOdata.mat');
GEOdata = GEOdata * 1e3;
a = 6378137;
f = 1/298.3;
b = a-a*f;
ec = sqrt(a^2-b^2)/a;
% X = GEOdata(500:700,1);
% Y = GEOdata(500:700,2);
% Z = GEOdata(500:700,3);
X = GEOdata(:,1);
Y = GEOdata(:,2);
Z = GEOdata(:,3);

wgs84 = wgs84Ellipsoid("meter");
[Lat,Lon,H] = ecef2geodetic(X,Y,Z,wgs84);
Lat = Lat*180/pi;
Lon = Lon*180/pi;

% % % L = atan(Y./X);
% % % L(X<0) = L(X<0)+pi;
% % % L = L*180/pi;
% % % 
% % % r = sqrt(X.^2+Y.^2);
% % % B1 = atand(Z./r);
% % % 
% % % for i = 1:length(X)
% % %     B2 = 0;
% % %     iter = 0;
% % %     while abs(B1(i)-B2) <= 0.1
% % %         iter = iter+1;
% % %         W1 = sqrt(1-ec*(sin(B1(i)))^2);
% % %         N1 = a/W1;
% % %         B2 = atand((Z(i)+N1*ec*sind(B1(i)))/r(i));
% % %         B1(i) = B2;
% % %         if iter >= 1000
% % %             break
% % %         end
% % %     end
% % % end

scatter(Lon,Lat);