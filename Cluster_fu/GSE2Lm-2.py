#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 15:43:41 2025

------- written by Wending Fu in Beijing ------
"""

import numpy as np
import sys
import os
from datetime import datetime

# 如果你还需要写回 mat 文件
from scipy.io import savemat

# 引入 geopack
from geopack import geopack

def calculate_lshell(pos_gsm_re):
    positions = np.asarray(pos_gsm_re, dtype=object)
    n = positions.shape[0]
    lshell = np.full(n, np.nan, dtype=float)
    epoch0 = datetime(1970, 1, 1)

    for i in range(n):
        YY, MM, DD, hh, mm, ss, x, y, z, kp, vx, vy, vz = positions[i]
        dt = datetime(int(YY), int(MM), int(DD), int(hh), int(mm), int(ss))
        ut = (dt - epoch0).total_seconds()

        geopack.recalc(ut, vxgse=vx, vygse=vy, vzgse=vz)

        try:
            xf, yf, zf, xs, ys, zs = geopack.trace(
                x, y, z,
                dir=1,
                rlim=60.0,
                r0=1.0,
                exname='t89',
                parmod=int(kp)
            )
        except Exception:
            continue

        idx = np.where(np.diff(np.sign(zs)) != 0)[0]
        if idx.size == 0:
            continue
        j = idx[0]

        z1, z2 = zs[j], zs[j+1]
        x1, x2 = xs[j], xs[j+1]
        y1, y2 = ys[j], ys[j+1]
        t = -z1 / (z2 - z1)
        xe = x1 + t * (x2 - x1)
        ye = y1 + t * (y2 - y1)

        lshell[i] = np.sqrt(xe**2 + ye**2)

    return lshell

def load_data(matfile):
    """
    尝试用 scipy.io.loadmat 加载；若遇到 v7.3 文件则用 h5py 读取。
    假设变量名为 'data'。
    """
    try:
        from scipy.io import loadmat
        mat = loadmat(matfile)
        return mat['data']
    except NotImplementedError:
        import h5py
        with h5py.File(matfile, 'r') as f:
            # 直接读取整个 dataset
            return f['data'][:]

if __name__ == '__main__':
    '''
    if len(sys.argv) < 2:
        print("Usage: python GSE2Lm.py <path_to_matfile>")
        sys.exit(1)

    matfile = sys.argv[1]'''
    matfile = 'C:/Users/pc/Documents/Timor/fwd_matlab_patch/Cluster_fu/R.mat'
    data = load_data(matfile)
    data = data.transpose

    L = calculate_lshell(data)

    # 输出路径可以根据需要调整
    out_dir = os.path.dirname(matfile)
    out_path = os.path.join(out_dir, 'L.mat')
    savemat(out_path, {'L': L})
    print(f"Saved L-shell array to {out_path}")
