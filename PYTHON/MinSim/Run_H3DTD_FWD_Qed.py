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

# # INPUT # #
input_dir = './INP_DIR/INP_DIR'
out_dir = './OUT_DIR'
nodes = ['g1', 'g2', 'g3']
ppn = ['12', '12', '8']

# # MAIN # #
root_dir = os.getcwd()
os.chdir(input_dir)

# Look at input directory and load in all models
mesh = []
for file in glob.glob("*.msh") :
    mesh.append(file)

# Move back to project folder
os.chdir(root_dir)

# List of all meshed
qlist = range(len(mesh))

# Loop through the tiles and Q a job
for ii in qlist[100:300]:

    # Cycle through the nodes
    idn = np.remainder(ii, len(nodes))

    tID = int(re.search('\d+',mesh[ii]).group())

    # Create directory
    tile_dir = out_dir + '/Tile_' + str(tID)
    os.system('mkdir ' + tile_dir)

    os.chdir(tile_dir)

    # Move files to directory
    fid = open('h3dtd.inp', 'w')
    fid.write(root_dir + os.path.sep + input_dir + os.path.sep + 'Tile_' + str(tID) + '.msh\n')
    fid.write('FILE ' + root_dir + os.path.sep + input_dir + os.path.sep + 'Tile_' + str(tID) + '.con\n')
    fid.write('0 \n')
    fid.write(root_dir + os.path.sep + input_dir + os.path.sep + 'Trx_Tile_' + str(tID) + '.loc\n')
    fid.write(root_dir + os.path.sep + 'VTEM.wf ! wave file\n')
    fid.write(root_dir + os.path.sep + 'times.txt\n')
    fid.close()

    # Write a pbs file and Q
    fid = open('Tile_' + str(tID) + '.pbs', 'w')
    fid.write('#PBS -l nodes=' + nodes[idn] + ':ppn=' + ppn[idn] + ':gold\n')
    fid.write('cd $PBS_O_WORKDIR\n')
    fid.write('cat $PBS_NODEFILE >nodes.txt\n')
    fid.write('mpiexec -n ' + ppn[idn] + ' -f $PBS_NODEFILE h3dtd_mumps_sym\n')
    fid.close()

    os.system('qsub ' + 'Tile_' + str(tID) + '.pbs')

    os.chdir(root_dir)

    # Move and rename pred0 and delete inversion files
    # os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    # os.system('rm recv_h3dtd.txt')
    # os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')
