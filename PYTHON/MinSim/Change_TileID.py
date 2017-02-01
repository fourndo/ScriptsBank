# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 12:18:24 2016

@author: dominiquef
"""

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import griddata
from SimPEG import mkvc, Mesh
import SimPEG.PF as PF
import SimPEG.EM.Static.DC as DC
import SimPEG.EM.FDEM as FDEM
from SimPEG.EM.Static.Utils import drapeTopotoLoc, gen_DCIPsurvey, writeUBC_DCobs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
import os
import glob
import re
import shutil

work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\FWR_TDEM\\OUT_DIR'
meshfile = '..\\..\\Mesh_20m.msh'
out_dir = 'NewGrid\\'
mesh = Mesh.TensorMesh.readUBC(work_dir + '\\' + meshfile)
os.chdir(work_dir)     


# Create data CC grid
offTile = 500
lenTile = 125
gCx = np.asarray(range(72))*62.5 + 500 + mesh.x0[0]
gCy = np.asarray(range(184))*62.5 + 500 + mesh.x0[1]

gCx_old = np.asarray(range(36))*lenTile + 500 + mesh.x0[0]
gCy_old = np.asarray(range(92))*lenTile + 500 + mesh.x0[1]

X,Y = np.meshgrid(gCx,gCy)
X, Y = mkvc(X), mkvc(Y)

X_old,Y_old = np.meshgrid(gCx_old,gCy_old)
X_old, Y_old = mkvc(X_old), mkvc(Y_old)

for file in glob.glob("*Tile*") :
    
    # Get tile number
    tID = int(re.search('\d+',file).group())
    
    name = re.split('_\d*\.',file)
    # Find equivalent in new grid
    indx = np.where(np.logical_and(X_old[tID]==X,Y_old[tID]==Y))[0][0]
    
    shutil.copy(file, out_dir + 'FWR_Tile_' + str(indx) + '.dat')
#    print('Old ID:' + str(tID) + ' New ID: ' + str(indx))