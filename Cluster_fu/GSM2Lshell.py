#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 15:43:41 2025

------- written by Wending Fu in Beijing ------
"""


import numpy as np
from geopack import geopack
import sys
import os
from scipy.io import loadmat, savemat
from datetime import datetime

def calculate_lshell(pos_gsm_re):
    positions = np.asarray(pos_gsm_re, dtype=object)
    n = positions.shape[0]
    lshell = np.full(n, np.nan, dtype=float)
    epoch0 = datetime(1970,1,1)
        
    for i in range(n):
        YY, MM, DD, hh, mm, ss, x, y, z, kp, vx, vy, vz = positions[i]
        direction = np.sign(z)
        time_val = datetime(int(YY), int(MM), int(DD), int(hh), int(mm), int(ss))
        if isinstance(time_val, datetime):
            dt = time_val
        else:
            dt = datetime.utcfromtimestamp(time_val)
        ut = (dt - epoch0).total_seconds()
        
        geopack.recalc(ut, vxgse=vx, vygse=vy, vzgse=vz)
             
        try:
            xf, yf, zf, xs, ys, zs = geopack.trace(
                x, y, z,
                dir=direction,
                rlim=60.0,
                r0=0.1,
                exname = 't89',
                parmod = int(kp)
            )
        except Exception:
            continue
        
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        zs = np.asarray(zs, dtype=float)
        try:
            # 大多数 geopack 版本支持数组输入
            xsm, ysm, zsm = geopack.gsm2sm(xs, ys, zs)
        except Exception:
            # 若你的版本仅支持标量，退化为逐点转换
            xsm = np.empty_like(xs); ysm = np.empty_like(ys); zsm = np.empty_like(zs)
            for k in range(xs.size):
                xsm[k], ysm[k], zsm[k] = geopack.gsm2sm(xs[k], ys[k], zs[k])

        idx = np.where(np.diff(np.sign(zsm)) != 0)[0]
        if idx.size == 0:
            continue
        j = idx[0]

        z1, z2 = zsm[j], zsm[j+1]
        x1, x2 = xsm[j], xsm[j+1]
        y1, y2 = ysm[j], ysm[j+1]
        t = -z1 / (z2 - z1)
        xe = x1 + t * (x2 - x1)
        ye = y1 + t * (y2 - y1)

        lshell[i] = np.sqrt(xe**2 + ye**2)

    return lshell


if __name__ == '__main__':
    
    matfile = sys.argv[1]
    
    #test = [[ 2002.000000000,1,1,1,1,1, 6.000000000,0.0,0.000000000,1.000000000,-400.000000000,0,0 ],
     #       [ 2002,1,1,1,1,1, 10.0,8.0,7.0,1,-400,0,0]]
    mat = loadmat(matfile)
    L = calculate_lshell(mat['data'])
    #L = calculate_lshell(test)
    #print(mat['data'])
    current_dir = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/Cluster_fu/'
    file_name = os.path.join(current_dir, 'L.mat')
    savemat(file_name, {'L': L})
