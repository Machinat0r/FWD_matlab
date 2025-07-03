function idx = checkLineIntersect(R1,R2,R3,R4, B1,B2,B3,B4, thresh1, thresh2)
%CHECKLINEINTERSECTALLTIMES：找到每两个矢量之间相距最短的位置，及其离卫星的距离
%
% 输入：
%   R?      n×4 矩阵，[time, x, y, z]；四颗卫星分别为 R1,R2,R3,R4  
%   B?      n×4 矩阵，[time, Bx, By, Bz]；四颗卫星分别为 B1,B2,B3,B4  
%   thresh1 ： 若某时刻的12个距离中大于该阈值倍的卫星间距直接跳过
%   thresh2 :  若某时刻的12个距离均大于该阈值倍的未行间距则直接跳过

    RR12 = irf_abs(R1-R2); RR13 = irf_abs(R1-R3); RR14 = irf_abs(R1-R4); 
    RR23 = irf_abs(R2-R3); RR24 = irf_abs(R2-R4); RR34 = irf_abs(R3-R4); 
    RR_mean = (RR12(:,5) + RR13(:,5) + RR14(:,5) + RR23(:,5) + RR24(:,5) + RR34(:,5))/6;

    [d1,distSat11,distSat12] = findLineClosestDistence(R1,B1,R2,B2);
    [d2,distSat21,distSat22] = findLineClosestDistence(R1,B1,R3,B3);
    [d3,distSat31,distSat32] = findLineClosestDistence(R1,B1,R4,B4);
    [d4,distSat41,distSat42] = findLineClosestDistence(R2,B2,R3,B3);
    [d5,distSat51,distSat52] = findLineClosestDistence(R2,B2,R4,B4);
    [d6,distSat61,distSat62] = findLineClosestDistence(R3,B3,R4,B4);

    distSat = [distSat11, distSat12, distSat21, distSat22, distSat31, distSat32, ...
        distSat41, distSat42, distSat51, distSat52, distSat61, distSat62];
    meandist = mean(distSat,2);

    id1 = any(distSat >= thresh1 * RR_mean, 2);
    id2 = any(distSat <= thresh2 * RR_mean, 2);
    id3 = any(meandist <= thresh2 * RR_mean, 2);
    
    idx = transpose(1:size(B1,1));
    idx1 = idx(id1 == 0 & id2 == 1);
    idx2 = idx(id3 == 1);
    idx = sort(unique([idx1;idx2]));
end

function [d,distSat1,distSat2] = findLineClosestDistence(R1,B1,R2,B2)
%FINDLINECLOSESTAPPROACH  向量化求两条线最近点及距离
%
% 输入（均为 n×4 矩阵）：
%   R1, R2  ：[time, x, y, z]，卫星 1、2 的位置
%   B1, B2  ：[time, Bx, By, Bz]，卫星 1、2 测得的磁场矢量
%
% 输出：
%   t1, t2      ：n×1，最近点对应的参数
%   P1, P2      ：n×3，分别是两条直线上距离最近的点坐标
%   d           ：n×1，两线最近点之间的距离
%   distSat1,distSat2：n×1，|t1| 和 |t2|，即两最近点到各自卫星的距离

    % 1) 位置和方向归一化
    r1 = R1(:,2:4);              % n×3
    r2 = R2(:,2:4);              % n×3
    v1 = B1(:,2:4);  v1 = v1./vecnorm(v1,2,2);  % n×3 单位化
    v2 = B2(:,2:4);  v2 = v2./vecnorm(v2,2,2);  % n×3

    % 2) 中间量
    dp  = r2 - r1;               % n×3
    a   = sum(v1.*v1,2);         % n×1
    b   = sum(v1.*v2,2);         % n×1
    c   = sum(v2.*v2,2);         % n×1
    dP  = sum(v1.*dp,2);         % n×1
    eP  = sum(v2.*dp,2);         % n×1
    denom = a.*c - b.^2;         % n×1

    % 3) 求 t1, t2
    t1 = ( b.*eP - c.*dP ) ./ denom;  % n×1
    t2 = ( a.*eP - b.*dP ) ./ denom;  % n×1

    % 4) 最近点坐标
    P1 = r1 + v1 .* t1;          % n×3  （隐式广播 t1→n×3）
    P2 = r2 + v2 .* t2;          % n×3

    % 5) 两直线最近距离
    d = sqrt( sum( (P2 - P1).^2, 2 ) );  % n×1

    % 6) 离卫星的距离
    distSat1 = abs(t1);
    distSat2 = abs(t2);
end
