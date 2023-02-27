function plotCircle3(C,R,n,color)
% 
%------written by Wending Fu, May.18.2022 in Beijing------------

theta = linspace(0, 2*pi, 100);
x = C(1)+R*n(2)/sqrt(n(1)^2+n(2)^2).*cos(theta)+...
	n(1)*n(2)/(sqrt(n(1)^2+n(2)^2)*sqrt(n(1)^2+n(2)^2+n(3)^2)).*sin(theta);
y = C(2)-R*n(1)/sqrt(n(1)^2+n(2)^2).*cos(theta)+...
	n(2)*n(3)/(sqrt(n(1)^2+n(2)^2)*sqrt(n(1)^2+n(2)^2+n(3)^2)).*sin(theta);
z = C(3)-R*sqrt(n(1)^2+n(2)^2)/sqrt(n(1)^2+n(2)^2+n(3)^2).*sin(theta);
plot3(x,y,z,'-', 'LineWidth', 0.1,'color',color)
hold on
end