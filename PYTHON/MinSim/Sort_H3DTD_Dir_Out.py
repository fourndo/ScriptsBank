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
input_dir = './INP_DIR'
out_dir = './OUT_DIR'

# # MAIN # #
root_dir = os.getcwd()

dirs = os.listdir(root_dir + os.path.sep + out_dir)

# Look at input directory and load in all models
for ii in dirs:

    curr_dir = root_dir + os.path.sep + out_dir + os.path.sep + ii

    if os.path.isdir(curr_dir):

        os.chdir(curr_dir)

        if glob.glob("recv_h3dtd.txt"):
            tID = int(re.search('\d+', ii).group())
            os.system('cp recv_h3dtd.txt ' + curr_dir + '/../Grid250/FWR_Tile_'+str(tID)+'.dat')



    # Move and rename pred0 and delete inversion files
    # os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    # os.system('rm recv_h3dtd.txt')
    # os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')
