#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 14:52:55 2025

------- written by Wending Fu in Beijing ------
"""

import sys
import torch
import numpy as np
from torch import nn
from torch.optim import Adam

__all__ = ['CalError']

def CalError(R1_list, R2_list, R3_list, R4_list,
             B1_list, B2_list, B3_list, B4_list,
             lr=1e-2, steps=500, device='cuda'):

    # 转为 torch 张量并搬到 GPU
    R1 = torch.tensor(R1_list, dtype=torch.float32, device=device)  # [N,M,3]
    R2 = torch.tensor(R2_list, dtype=torch.float32, device=device)
    R3 = torch.tensor(R3_list, dtype=torch.float32, device=device)
    R4 = torch.tensor(R4_list, dtype=torch.float32, device=device)
    B1 = torch.tensor(B1_list, dtype=torch.float32, device=device)
    B2 = torch.tensor(B2_list, dtype=torch.float32, device=device)
    B3 = torch.tensor(B3_list, dtype=torch.float32, device=device)
    B4 = torch.tensor(B4_list, dtype=torch.float32, device=device)

    N, M, _ = R1.shape
    # 参数：Q,x,y,z, shape [N,4]
    params = nn.Parameter(torch.zeros(N, 4, device=device))
    # 初始化 Q 为 1，其它为 0
    params.data[:, 0] = 1.0

    optimizer = Adam([params], lr=lr)

    def compute_loss():
        # params: [N,4] -> Q:[N,1,1], loc:[N,1,3]
        Q = params[:, [0]].unsqueeze(-1)  # [N,1,1]
        loc = params[:, 1:].unsqueeze(1)  # [N,1,3]
        loss = 0.0
        # 定义事件对列表
        pairs = [
            (R1, B1, R2, B2),
            (R2, B2, R3, B3),
            (R3, B3, R4, B4),
            (R4, B4, R1, B1),
            (R1, B1, R3, B3),
            (R2, B2, R4, B4)
        ]
        for Ra, Ba, Rb, Bb in pairs:
            # 计算距离和归一化方向
            da = torch.norm(Ra - loc, dim=2, keepdim=True)  # [N,M,1]
            na = (Ra - loc) / da
            db = torch.norm(Rb - loc, dim=2, keepdim=True)
            nb = (Rb - loc) / db
            # 残差
            fa = Q / (4 * torch.pi * da**2) * na - Ba.unsqueeze(1)  # Ba:[N,M,3]
            fb = Q / (4 * torch.pi * db**2) * nb - Bb.unsqueeze(1)
            loss = loss + (fa**2).sum(dim=(1,2)) + (fb**2).sum(dim=(1,2))  # [N]
        return loss.sum()  # 对所有事件累加

    # 迭代优化
    for _ in range(steps):
        optimizer.zero_grad()
        loss = compute_loss()
        loss.backward()
        optimizer.step()

    # 获取结果
    params_cpu = params.detach().cpu().numpy()  # [N,4]
    Q_out   = params_cpu[:, 0]
    Loc_out = params_cpu[:, 1:]
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

    R1 = np.array(10*[[0,0,3],[0,0,3]])
    R2 = np.array(10*[[-2*np.sqrt(2)-1,2,-1],[-2*np.sqrt(2)-1,2,-1]])
    R3 = np.array(10*[[np.sqrt(2)+2,np.sqrt(6)+1,-2],[np.sqrt(2)+2,np.sqrt(6)+1,-2]])
    R4 = np.array(10*[[np.sqrt(2)-1,-np.sqrt(6)+1,-3],[np.sqrt(2)-1,-np.sqrt(6)+1,-3]])
    
    B1 = np.array([[0,0,8.84194128288308],[0,0,8.84194128288308]])
    B2 = np.array([[-3.49574298101714,1.82620322503796,-0.913101612518978],[-3.49574298101714,1.82620322503796,-0.913101612518978]])
    B3 = np.array([[1.87828065524854,1.89768734028592,-1.10027133390157],[1.87828065524854,1.89768734028592,-1.10027133390157]])
    B4 = np.array([[0.870920849817176,-3.04768108352988,-6.30776678214638],[0.870920849817176,-3.04768108352988,-6.30776678214638]])
    [Q, Loc] = CalError(R1, R2, R3, R4,B1, B2, B3, B4,lr=1e-2, steps=500, device='cuda')
    
    
