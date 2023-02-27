clear;clc;close all

[X,Y] = meshgrid(-5:.01:5);
% Z = Y.*sin(X) - X.*cos(Y);
Z = sin(Y);
s = surf(X,Y,Z,'FaceAlpha',0.5);
s.EdgeColor = 'none';
grid off