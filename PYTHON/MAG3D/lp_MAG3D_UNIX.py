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
import shutil, re
from numpy import *
from numpy.linalg import *

from get_UBC_mesh import get_UBC_mesh
from read_input import read_input
from read_UBC_log import read_UBC_log
from read_UBC_log_primer import read_UBC_log_primer
from comp_dmdx import comp_dmdx
from comp_dmdy import comp_dmdy
from comp_dmdz import comp_dmdz
from make_depthW import make_depthW
from get_input_files import get_input_files

from get_Ws3D import get_Ws3D

Home_dir, UBC_dir, iter_start, chifact, pvec, qvec, lvec, cool_beta, iter_max = read_input('lp_MAG3D_input.inp')

Home_dir = Home_dir + '/'
UBC_dir = UBC_dir + '/'
Work_dir = Home_dir + 'Workspace/'   


if not os.path.exists(Work_dir): os.makedirs(Work_dir);    
else: shutil.rmtree(Work_dir); os.makedirs(Work_dir); 

if iter_start>10:
    strnum = str(iter_start)
else:
    strnum = '0' + str(iter_start)
    

for ll in range(size(lvec)):
    
    for pp in range(size(pvec)):

        for qq in range(size(qvec)):
            
            p = float(pvec[pp])
            q = float(qvec[qq])
            l = float(lvec[ll])
            
            oo=0
            delta = 1e-6
            file_list = []
            new_dir = Home_dir + 'lp' + str(p) + 'lq' + str(q) + 'lambda' + str(l) + '/'
            

            if not os.path.exists(new_dir): os.makedirs(new_dir)
            else: 
                shutil.rmtree(new_dir)
                os.makedirs(new_dir)
            
            shutil.copy(UBC_dir + 'maginv3d_00' + str(iter_start) + '.sus', Work_dir)
            #os.rename(work_dir + 'maginv3d_0' + str(iter_start) + '.con','maginv3d_01.con')
            shutil.copy(UBC_dir + 'maginv3d.log', Work_dir)     
                 
            # Read all the input files and extract parameters
            meshfile, obsfile, topo_check = get_input_files(Home_dir)
             
            # Load mesh and create dimension vectors
            nX, nY, nZ, origin, dX, dY, dZ = get_UBC_mesh(Home_dir + meshfile)

            print dX.size
            print dY.size
            print dZ.size
            
            mcell=nX*nY*nZ                                      
             
            
            #actv_cells.resize((nY,nX,nZ))
            
            # Read *.log and extract information
            ndata, phid, Lxyz, beta_in, ref_model = read_UBC_log_primer(Work_dir,Home_dir,iter_start,mcell)

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
            
            #lbda = l     # Scaling between smooth and compact
            beta = beta_in  # Trade-off parameter
                       
            #wzx, wzy = make_depthW(layers,nZ,nX,nY,actv_cells);
            
            #model_in = loadtxt(Work_dir + 'maginv3d_00' + str(iter_start) + '.sus')
            #model_in = array(Work_dir + 'maginv3d_00' + str(iter_start) + '.sus')
            fid = open(Work_dir + 'maginv3d_00' + str(iter_start) + '.sus','r')
            fid2 = open(Home_dir + topo_check)
            
            model=zeros((mcell,1), dtype=float)
            actv_cells=zeros((mcell,1), dtype=int)
            limit = mcell-1

            for ii in range (mcell):         
                line = fid.readline()
                #line = line.split('\n')
                model[ii]=float(line)
                line2 = fid2.readline()
                actv_cells[ii]=int(line2)
                
            
            fid.close()
            fid2.close()

            # Recompute weights and adjust beta
            
                        # Create derivative matrices
            wcx=ones(((nX-1)*nY*nZ,1), dtype=float)
            wcy=ones((nX*(nY-1)*nZ,1), dtype=float)
            wcz=ones((nX*nY*(nZ-1),1), dtype=float) 
            
            
            
            Ws = get_Ws3D(mcell,dX,dY,dZ,actv_cells)
            
            ## Inversion iterations
            while (oo<= iter_max ) & (phid >= ndata*chifact):

                model = model - ref_model
                print model.size
                model = model * actv_cells
                print model.size
                # Compute derivative vectors
                dmdx, Vx = comp_dmdx( model , dX , dY , dZ , actv_cells )
                dmdy, Vy = comp_dmdy( model , dX , dY , dZ , actv_cells )
                dmdz, Vz = comp_dmdz( model , dX , dY , dZ , actv_cells )
                print dmdx.size
                phi_xyz = float(alphax * dot( transpose(Vx * dmdx) , (Vx * dmdx) ) + \
                    alphay * dot( transpose(Vy * dmdy) , (Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (Vz * dmdz) ))
                #print '{:e}.'.format(phi_xyz)
                # lp-norm on gradient ||GRAD(m)||p
                wcx = ( (abs(dmdx) ) ** ( 2 - p ) + delta ) ** -1
            
                wcy = ( (abs(dmdy) ) ** ( 2 - p ) + delta ) ** -1
               
                wcz = ( ( abs(dmdz) ) ** ( 2 - p ) + delta ) ** -1
                print wcx.size
                # Compute norm of | GRAD (m) |p and scale on | GRAD (m) |2
                phi_cxyz = float(alphax * dot( transpose( Vx * dmdx) , (wcx * Vx * dmdx) ) + \
                    alphay * dot( transpose( Vy * dmdy) , (wcy * Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (wcz * Vz * dmdz) ))
                
                scale_xyz = phi_xyz / phi_cxyz
     
                # Minimum support on model | m |q
                wcm =  ( ( abs(model) ) ** ( 2 - q ) + delta ) ** -1

                # Compute norm of |Wm| and scale on |GRAD(m)|
                phi_cm = float(dot( transpose(Ws * model) , (wcm * Ws * model)))
                
                scale_cm = (scale_xyz * phi_cxyz) / (alphas * phi_cm)
            
                phi_m = scale_cm * alphas *phi_cm + scale_xyz*phi_cxyz   
                
                
                # Apply scale on weighting vectors... Ready to export
                wcm = (l * scale_cm * wcm)
                #wcx = abs(2.0 - l) * (scale_xyz * wcx)
                #wcy = abs(2.0 - l) * (scale_xyz * wcy)
                #wcz = abs(2.0 - l) * (scale_xyz * wcz)

                # Write weight file
                #w_out = concatenate(((wcm),wcx,wcy,wcz))
                
                fid = open(Work_dir + 'input_w.dat','w')
                
                for ii in range (mcell):
                    fid.write('%12.8e\n'%wcm[ii])
                    
                for ii in range ((nX-1)*nY*nZ):
                    fid.write('%12.8e\n'%wcx[ii])
                                
                for ii in range (nX*(nY-1)*nZ):
                    fid.write('%12.8e\n'%wcy[ii])
                                    
                for ii in range (nX*nY*(nZ-1)):
                    fid.write('%12.8e\n'%wcz[ii]) 
                           
                fid.close();
                #savetxt(Work_dir + 'input_w.dat',w_out)
                mweights = mcell + (nX-1)*nY*nZ + nX*(nY-1)*nZ + nX*nY*(nZ-1)
                savetxt(Home_dir + 'ref_model.dat',ref_model)
                # Adjust trade-off parameter for next iteration
                if (phid < ndata*2) & (oo >0):
                    beta =cool_beta *beta;
                elif oo >= 0: 
                    beta = cool_beta * beta;
                    
                # Write input file
                fid = open(Work_dir + 'lp_ctrl_file.inp','w')
                fid.write('0 !!! restart=1, max # of iterations\n')
                fid.write('2 !!! mode\n')
                fid.write('%12.8e 0.1 !!! mode chifact \n'% beta)
                fid.write('%s%s  ! observations file\n'% (Home_dir,obsfile))
                fid.write('%smaginv3d.mtx  !!! mesh or ncel,aspr\n'% (Home_dir));
                fid.write('%s%s%s%s !!! initial model (NULL = same as reference model)\n'%(Work_dir,'maginv3d_0',strnum,'.sus'))
                fid.write('%sref_model.dat !!! reference model (NULL = calculate)\n'% (Home_dir))
                fid.write('null !!! active cells file (NULL = all cells are active)\n')
                fid.write('0 1 !!! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file  ! bounds\n')
                fid.write('%d %d %d  \n'%(Lx,Ly,Lz))
                fid.write('SMOOTH_MOD\n')
                fid.write('input_w.dat  ! weighting file\n')
                fid.write('0')
                fid.close();
                
                
                oo=oo+1
                
                os.chdir(Work_dir)
                
                #os.system('qsub DClp_3796.pbs')
                os.system('maginv3d_newbounds lp_ctrl_file.inp')
                
                #Extract last phid
                ndata, phid, Lxyz, beta_in = read_UBC_log(Work_dir,Home_dir,1,mcell)
                
                
                if phid ==99999:
                    print('Inversion has stopped at convergence. Iteration: %i\n'%oo)
                else:
                    
                    # Copy and rename result to "new_dir"
                    os.system('cp maginv3d_001.sus '+ new_dir + 'maginv3d_0'+str(oo)+'.sus')
                    os.system('cp maginv3d_001.pre '+ new_dir + 'maginv3d_0'+str(oo)+'.pre')
                    os.system('cp maginv3d.log '+ new_dir + 'maginv3d_0'+str(oo)+'.log')

                    
                    fid = open(Work_dir + 'maginv3d_001.sus','r')

                    for ii in range (mcell):         
                        line = fid.readline()
                        model[ii]=float(line)
       
                    fid.close()

                
                
                