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
                dir=1,
                rlim=60.0,
                r0=1.0,
                exname = 't89',
                parmod = int(kp)
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
