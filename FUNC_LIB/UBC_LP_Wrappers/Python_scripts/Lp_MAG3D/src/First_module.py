'''
Created on Jul 17, 2013

@author: dominiquef
'''
# MAG_lp_Driver_PyDev_v1.0
# Author: D Fournier
# Last Update: July 16th, 2013
#
# This program invert for "compact" model using a iterative lp-norm on the
# model and gradient. It is an adapted algorithm based on the classic 
# paper of Last & Kubic. The central solver is the UBC - MAGINV3D code. 
import os 
import shutil

         
Home_dir = 'C:\\Users\\dominiquef\\Documents\\Eclipse\\Projects\\Test_proj\\data\\'
UBC_dir = 'C:\\Users\\dominiquef\\Documents\\Eclipse\\Projects\\Test_proj\data\\UBC_Results\\'

from get_UBC_mesh import get_UBC_mesh
from read_input import read_input
from read_UBC_log import read_UBC_log
from numpy import *
from comp_dmdx import comp_dmdx
from comp_dmdy import comp_dmdy
from comp_dmdz import comp_dmdz
from get_Ws3D import get_Ws3D

iter_start = 2

# Target chi factor
chifact = 1.0

# Lp-norm on gradient (l2-norm by default)
p = array([1.8])

# Lq-norm on model
q = array([0])

# Scale phi_cxyy and phi_cm
l = array([1.25])

# Trade-off cooling schedule
cool_beta = 0.5

# Maximum number of iterations
iter_max = 20

Work_dir = Home_dir + 'Workspace\\'
            

if not os.path.exists(Work_dir): os.makedirs(Work_dir)
else: shutil.rmtree(Work_dir); os.makedirs(Work_dir)


for ll in range(1 , 2):
    
    for pp in range(1 , 2):

        for qq in range(1 , 2):
            
            
            file_list = []
            new_dir = Home_dir + 'lp' + str(p) + 'lq' + str(q) + 'lambda' + str(l)
            

            if not os.path.exists(new_dir): os.makedirs(new_dir)
            else: 
                shutil.rmtree(new_dir)
                os.makedirs(new_dir)
            
            shutil.copy(UBC_dir + 'dcinv3d_0' + str(iter_start) + '.con', Work_dir)
            #os.rename(work_dir + 'dcinv3d_0' + str(iter_start) + '.con','dcinv3d_01.con')
            shutil.copy(UBC_dir + 'dcinv3d.log', Work_dir)     
                 
            # Read all the input files and extract parameters
            meshfile, obsfile, topofile, actv_cells = read_input(Home_dir)
             
            # Load mesh and create dimension vectors
            numcell, origin, dX, dY, dZ = get_UBC_mesh(Home_dir + meshfile)

            nX = numcell[0]
            nY = numcell[1]
            nZ = numcell[2] 
               
            mcell=nX*nY*nZ                                      
             
            actv_cells= array(actv_cells)
            #actv_cells.resize((nY,nX,nZ))
            
            
            # Read *.log and extract information
            ndata, phid, Lxyz, beta_in, ref_model = read_UBC_log(Work_dir,iter_start)

            ## Set-up parameters

      

            # Length scales
            Lx=Lxyz[0]         
            Ly=Lxyz[1] 
            Lz=Lxyz[2] 

            # Scaling factors
            alphas = 1/pow(Lx,2)    # Coefficient applied to WstWs in UBC
            alphac = alphas    # Weight on model minimum support
            alphax = 1.0e-0    # Weight on gradx 
            alphay = 1.0e-0    # Weight on grady
            alphaz = 1.0e-0    # Weight on gradz
            
            lbda = l[ll-1]     # Scaling between smooth and compact
            beta = beta_in  # Trade-off parameter

            # Create derivative matrices
            wcx=ones(((nX-1)*nY*nZ,1), dtype=float)
            wcy=ones((nX*(nY-1)*nZ,1), dtype=float)
            wcz=ones((nX*nY*(nZ-1),1), dtype=float) 

            
            Ws = get_Ws3D(mcell,dX,dY,dZ,actv_cells);
            
#            wxdepth, wydepth = get_depthweight(wdepth',nZ,nX,nY,topo_model);
            
            oo=1
                  
            fid = open(Work_dir + 'dcinv3d_0' + str(iter_start) + '.con','r');
            tempo = fid.read()
            tempo = tempo.split('\n')
            model = []
            for tt in range(len(tempo)-1):
                model.append(float(tempo[tt]))
            model = array(model)
            
            
            logmodel = log10(model) - log10(ref_model)            
            logmodel = logmodel * actv_cells
                       
            model = model - ref_model
            model = model * actv_cells
        
            # Compute derivative vectors
            dmdx, Vx = comp_dmdx( logmodel , dX , dY , dZ , actv_cells )
            dmdy, Vy = comp_dmdy( logmodel , dX , dY , dZ , actv_cells )
            dmdz, Vz = comp_dmdz( logmodel , dX , dY , dZ , actv_cells )
            #savetxt('dmdx.dat',dmdx)
            #savetxt('dmdy.dat',dmdy)
            #savetxt('dmdz.dat',dmdz)
            
            