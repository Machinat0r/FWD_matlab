#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 14:52:55 2025

------- written by Wending Fu in Beijing ------
"""

import sys
import time
import torch
import numpy as np

__all__ = ['CalError']

def CalError(R1, R2, R3, R4,
             B1, B2, B3, B4,
             max_iter=50, tol=1e-6,
             init_lambda=1.0,
             device='cuda'):
    """
    GPU 加速 Levenberg–Marquardt 最小二乘，类似 MATLAB lsqnonlin。
    对每个事件和 6 对组合分别求解 [Q, x, y, z]。

    输入:
      R1..R4, B1..B4: numpy.ndarray, shape [N,3]
    可选:
      max_iter   : 最大 LM 迭代次数
      tol        : 收敛阈值
      init_lambda: LM 阻尼初始值
      device     : 'cuda' 或 'cpu'

    输出:
      Q_out:   numpy.ndarray, [N,6]
      Loc_out: numpy.ndarray, [N,6,3]
    """
    device = torch.device(device)
    # 转张量
    Rs = [torch.tensor(arr, dtype=torch.float32, device=device) for arr in (R1, R2, R3, R4)]
    Bs = [torch.tensor(arr, dtype=torch.float32, device=device) for arr in (B1, B2, B3, B4)]
    N = R1.shape[0]
    # 6 对组合
    pairs = [(Rs[i], Bs[i], Rs[(i+1)%4], Bs[(i+1)%4]) for i in range(4)] + \
            [(Rs[0], Bs[0], Rs[2], Bs[2]), (Rs[1], Bs[1], Rs[3], Bs[3])]

    Q_out   = np.zeros((N, 6), dtype=np.float32)
    Loc_out = np.zeros((N, 6, 3), dtype=np.float32)

    # 对每对组合和每条记录分别 LM 求解
    for pi, (Ra, Ba, Rb, Bb) in enumerate(pairs):
        for i in range(N):
            # 初始参数 [Q, x, y, z]
            P = torch.tensor([1e5, 0.0, 0.0, 0.0],
                             dtype=torch.float32, device=device, requires_grad=True)
            lamb = init_lambda

            def residual(p):
                Q, x, y, z = p
                loc = torch.stack([x, y, z])  # [3]
                # 计算两组残差
                da = torch.norm(Ra[i] - loc)
                na = (Ra[i] - loc) / da
                fa = Q / (4 * torch.pi * da**2) * na - Ba[i]
                db = torch.norm(Rb[i] - loc)
                nb = (Rb[i] - loc) / db
                fb = Q / (4 * torch.pi * db**2) * nb - Bb[i]
                return torch.cat([fa, fb])  # [6]

            # LM 迭代
            for _ in range(max_iter):
                # 计算残差和雅可比
                r = residual(P)
                cost = (r**2).sum()
                J = torch.autograd.functional.jacobian(residual, P)  # [6,4]
                # 构造解方程
                JTJ = J.T @ J                                    # [4,4]
                JTr = J.T @ r                                    # [4]
                A = JTJ + lamb * torch.eye(4, device=device)
                delta = torch.linalg.solve(A, -JTr)              # [4]

                # 更新尝试
                P_new = P + delta
                r_new = residual(P_new)
                cost_new = (r_new**2).sum()

                # 判断接受与阻尼调整
                if cost_new < cost:
                    # 接受更新，减小阻尼
                    P = P_new.detach().requires_grad_(True)
                    lamb = lamb * 0.1
                else:
                    # 拒绝，增大阻尼
                    lamb = lamb * 10.0

                if abs(cost - cost_new) < tol:
                    break

            # 存储结果
            P_cpu = P.detach().cpu().numpy()
            Q_out[i, pi]      = P_cpu[0]
            Loc_out[i, pi, :] = P_cpu[1:]

    return Q_out, Loc_out


if __name__ == "__main__":
    
    """
    R1 = sys.argv[1]
    R2 = sys.argv[2]
    R3 = sys.argv[3]
    R4 = sys.argv[4]
    
    B1 = sys.argv[5]
    B2 = sys.argv[6]
    B3 = sys.argv[7]
    B4 = sys.argv[8]
    """
    
    R1 = np.array(200*[[0,0,30.000]])
    R2 = np.array(200*[[-20*np.sqrt(2)-10,20,-10]])
    R3 = np.array(200*[[10*np.sqrt(2)+20,10*np.sqrt(6)+10,-20]])
    R4 = np.array(200*[[10*np.sqrt(2)-10,-10*np.sqrt(6)+10,-30]])
    
    B1 = np.array(200*[[-0.652226388383353,2.60890555353341,17.6101124863505]])
    B2 = np.array(200*[[-4.76945737867553,2.91381190113354,-1.57831477978067]])
    B3 = np.array(200*[[2.54039879266100,2.95069672212825,-1.76298754234068]])
    B4 = np.array(200*[[0.993307057965726,-3.31769755861572,-10.4321190547298]])
    
    st_time = time.time()
    [Q, Loc] = CalError(R1, R2, R3, R4,B1, B2, B3, B4)
    ed_time = time.time()
    print(f'{ed_time-st_time}')
    
