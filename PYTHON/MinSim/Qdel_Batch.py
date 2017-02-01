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
Qend = 10896
Qstart = 907

for ii in range(Qend-Qstart):

    # Cycle through the nodes
    os.system('qdel ' + str(ii+Qstart+1))


    # Move and rename pred0 and delete inversion files
    # os.system('cp recv_h3dtd.txt '+ out_dir + '/FWR_Tile_'+str(tID)+'.dat')
    # os.system('rm recv_h3dtd.txt')
    # os.system('rm h3dtd.log')
    #os.system('rm h3dtdinv.out')
    #os.system('rm h3dtdinv_stat.txt')
