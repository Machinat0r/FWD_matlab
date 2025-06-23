function [Q,resQ,locCentroid,resLoc] = CalErrorFast( ...
        R1,R2,R3,R4, ...
        B1,B2,B3,B4, ...
        MultiPower,idx_flag,LocationSkew)

%----------------------------------------------------------------------
% 快速版 CalError —— 对整个矩阵一次性回归（避免逐行回归）
%----------------------------------------------------------------------

% ---------- 常量与初值 ----------
x0  = [idx_flag*1e4*MultiPower+15, LocationSkew*[1 1 1]];
opt = optimoptions('lsqnonlin','Display','off', ...
                   'Jacobian','on','MaxIter',100);

% ---------- 对整个矩阵进行回归：优化一次 ----------
[t,~,flag] = lsqnonlin(@residualAll,x0,[],[],opt);
if flag <= 0
    Q = NaN; resQ = NaN(6,1); resLoc = NaN(6,3); locCentroid = NaN(1,3); return;
end

% ---------- 结果 ----------
resAll      = residualAll(t);       % 18×1
resQ        = reshape(resAll,3,6).'; % 6×3
Q           = mean(resQ(:,1));       % 计算平均 Q
resLoc      = sat2locs(t);          % 求解六个交点位置
locCentroid = mean(resLoc,1);      % 计算位置的重心
Rmat = [R1; R2; R3; R4];          % 4×3

% ---------- 嵌套函数：残差与雅可比矩阵 ------------------------------
    function [F, J] = residualAll(v)
        % v = [Q, x, y, z]
        loc  = v(2:4);                    % 提取位置向量（1×3）
        Br   = [B1(:,2:4); B2(:,2:4); B3(:,2:4); B4(:,2:4)]; % 4×3


        % 计算每个卫星到当前位置的距离向量
        d    = Rmat - loc;                % 4×3
        r2   = sum(d.^2,2);               % 4×1，距离的平方
        r    = sqrt(r2);                  % 4×1，实际距离
        dn   = d ./ r;                    % 单位向量

        % 单极子场模型：B = Q / (4π r²) * r̂
        Bm   = (v(1) ./ (4 * pi * r2)) .* dn; % 4×3

        % 对六个方程对 (B1-B2, B2-B3,...，即每两对卫星) 计算残差
        pair = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];
        F    = reshape(Bm(pair(:,1),:) - Bm(pair(:,2),:) ...
                       - (Br(pair(:,1),:) - Br(pair(:,2),:)),[],1);

        if nargout > 1
            % 计算解析雅可比矩阵（18×4）
            dB_dQ   = (1 ./ (4 * pi * r2)) .* dn;          % dB/dQ
            dB_dxyz = (-2 * v(1) ./ (4 * pi * r2.^2)) .* dn + ...
                      (v(1) ./ (4 * pi * r2)) .* ...
                      ((eye(3) - dn .* dn) ./ r);         % dB/dx,dy,dz

            J = zeros(18, 4);
            for k = 1:6
                i = pair(k, 1); j = pair(k, 2); rows = (k - 1) * 3 + (1:3);
                J(rows, 1)   = (dB_dQ(i,:) - dB_dQ(j,:)).';
                J(rows, 2:4) = squeeze(dB_dxyz(i,:,:) - dB_dxyz(j,:,:));
            end
        end
    end

    function locs = sat2locs(v)
        % 求解六个交点的位置
        pair = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4]; 
        locs = zeros(6, 3);
        for k = 1:6
            i = pair(k, 1); j = pair(k, 2);
            Bi = Bmodel(i); Bj = Bmodel(j);
            t  = dot(Rmat(i,:) - Rmat(j,:), Bi - Bj) / dot(Bi - Bj, Bi - Bj);
            locs(k, :) = Rmat(j, :) + (Bi - Bj) * t;   % 交点位置（近似）
        end
        % 单极子模型的场计算
        function Bm = Bmodel(idx)
            d  = Rmat(idx, :) - v(2:4); 
            r2 = dot(d, d); 
            Bm = v(1) ./ (4 * pi * r2) .* d ./ sqrt(r2);  % 返回 B
        end
    end
end
