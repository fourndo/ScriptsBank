# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 18:05:02 2017

@author: DominiqueFournier
"""

import numpy as np


fid = open('GocadFlower.pl', 'wb')

fid.write(('GOCAD PLine 1\n').encode())
fid.write(('HEADER {\n').encode())
fid.write(('name: NerdyRose\n').encode())
fid.write(('}\n').encode())
fid.write(('PROPERTIES theta\n').encode())
fid.write(('PROP_LEGAL\n').encode())
fid.write(('ESIZES 1 1 1\n').encode())

k = 1/10

theta = np.linspace(0, 2/k*np.pi, 1000)


fid.write(('ILINE\n').encode())

for tt in range(theta.shape[0]):
    
    x = 1000*np.cos(theta[tt]*k) * np.cos(theta[tt]) + 315500
    y = 1000*np.cos(theta[tt]*k) * np.sin(theta[tt]) + 6071000
    fid.write(('PVRTX ').encode())
    np.savetxt(fid, np.c_[tt, x, y, theta[tt], theta[tt]], fmt=['%i', '%e', '%e', '%e','%e'], delimiter=' ')

for tt in range(theta.shape[0]-1):
    fid.write(('SEG ').encode())
    np.savetxt(fid, np.c_[tt,tt+1].reshape((1,2)), fmt='%i',delimiter=' ')

np.savetxt(fid, np.c_[tt+1,0].reshape((1,2)), fmt='%i',delimiter=' ')
fid.write(('END').encode())
fid.close()