# -*- coding: utf-8 -*-
"""
Run H3DTD forward models on MIRA cluster

Run through a list of Observation, mesh, model files and call H3DTD on the
cluster. The program than copies and rename the pred0 file to final directory.

Created on Sat Nov 12 10:27:25 2016

@author: dominiquef
"""

import os
import glob
import numpy as np
import re

#work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\FWR_TDEM'

#os.chdir(work_dir)

input_dir = './INP_DIR'
out_dir = './OUT_DIR'

os.chdir(input_dir)

# Look at input directory and load in all models
mesh = []
for file in glob.glob("*.msh") :
    mesh.append(file)
    
os.chdir('../')

# Loop through the tiles
for ii in range(len(mesh)):

    tID = int(re.search('\d+',mesh[ii]).group())
    
    fid = open('h3dtd.inp', 'w')    
    fid.write(input_dir + '/Tile_' + str(tID) + '.msh\n')
    fid.write('FILE ' + input_dir + '/Tile_' + str(tID) + '.con\n')
    fid.write('0 \n')
    fid.write(input_dir + '/Trx_Tile_' + str(tID) + '.loc\n')
    fid.write('VTEM.wf ! wave file\n')
    fid.write('times.txt\n')
    
    fid.close()
    
    os.system('mpiexec -n 24 -f $PBS_NODEFILE h3dtd_mumps_sym')

    # Move and rename pred0 and delete inversion files
    os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    os.system('rm recv_h3dtd.txt')
    os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')     