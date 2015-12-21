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


from get_UBC_mesh import get_UBC_mesh
from read_input import read_input
from read_UBC_log import read_UBC_log
from numpy import *
from numpy.linalg import *
from comp_dmdx import comp_dmdx
from comp_dmdy import comp_dmdy
from comp_dmdz import comp_dmdz
from make_depthW import make_depthW
from get_input_files import get_input_files

from get_Ws3D import get_Ws3D

Home_dir, UBC_dir, iter_start, chifact, pvec, qvec, lvec, cool_beta, iter_max, layers = read_input('lp_DC_input.inp')

Home_dir = Home_dir + '\\'
UBC_dir = UBC_dir + '\\'

Work_dir = Home_dir + 'Workspace\\'

if not os.path.exists(Work_dir): os.makedirs(Work_dir)
else: shutil.rmtree(Work_dir); os.makedirs(Work_dir)

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
            delta = 1e-08
            file_list = []
            new_dir = Home_dir + 'lp' + str(p) + 'lq' + str(q) + 'lambda' + str(l) + '/'
            

            if not os.path.exists(new_dir): os.makedirs(new_dir)
            else: 
                shutil.rmtree(new_dir)
                os.makedirs(new_dir)
            
            shutil.copy(UBC_dir + 'dcinv3d_0' + str(iter_start) + '.con', Work_dir)
            #os.rename(work_dir + 'dcinv3d_0' + str(iter_start) + '.con','dcinv3d_01.con')
            shutil.copy(UBC_dir + 'dcinv3d.log', Work_dir)     
                 
            # Read all the input files and extract parameters
            meshfile, obsfile, topofile, actv_cells = get_input_files(Home_dir)
             
            # Load mesh and create dimension vectors
            numcell, origin, dX, dY, dZ = get_UBC_mesh(Home_dir + meshfile)

            nX = numcell[0]
            nY = numcell[1]
            nZ = numcell[2] 
               
            mcell=nX*nY*nZ                                      
             
            actv_cells= array(actv_cells)
            #actv_cells.resize((nY,nX,nZ))
            
            
            # Read *.log and extract information
            ndata, phid, Lxyz, beta_in, ref_model, str_ref_model = read_UBC_log(Work_dir,Home_dir,iter_start,mcell)

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
            
            beta = beta_in  # Trade-off parameter

            # Create derivative matrices
            wcx=ones(((nX-1)*nY*nZ,1), dtype=float)
            wcy=ones((nX*(nY-1)*nZ,1), dtype=float)
            wcz=ones((nX*nY*(nZ-1),1), dtype=float) 
          
            Ws = get_Ws3D(mcell,dX,dY,dZ,actv_cells);
            
            wzx, wzy = make_depthW(layers,nZ,nX,nY,actv_cells);
            
            model = loadtxt(Work_dir + 'dcinv3d_0' + str(iter_start) + '.con')
            logmodel = log10(model)
    ## Inversion iterations
            # Recompute weights and adjust beta
            
            while (oo<= iter_max) & (phid >= ndata*chifact):           
                                      
                for rr in range(len(ref_model)):       
                    model[rr] = model[rr] - ref_model[rr]
                    logmodel[rr] = logmodel[rr] - log10(ref_model[rr]) 
                    
                model = model * actv_cells
                logmodel = logmodel * actv_cells
                
                # Compute derivative vectors
                dmdx, Vx = comp_dmdx( logmodel , dX , dY , dZ , actv_cells )
                dmdy, Vy = comp_dmdy( logmodel , dX , dY , dZ , actv_cells )
                dmdz, Vz = comp_dmdz( logmodel , dX , dY , dZ , actv_cells )
            
                phi_xyz = alphax * dot( transpose(wzx * Vx * dmdx) , (Vx * dmdx) ) + \
                    alphay * dot( transpose(wzy * Vy * dmdy) , (Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (Vz * dmdz) )
                    
                # lp-norm on gradient ||GRAD(m)||p
                wcx = wzx * ( (abs(dmdx) ) ** ( 2 - p ) + delta ) ** -1
            
                wcy = wzy * ( (abs(dmdy) ) ** ( 2 - p ) + delta ) ** -1
               
                wcz =  ( ( abs(dmdz) ) ** ( 2 - p ) + delta ) ** -1
            
                # Compute norm of | GRAD (m) |p and scale on | GRAD (m) |2
                phi_cxyz = alphax * dot( transpose( Vx * dmdx) , (wcx * Vx * dmdx) ) + \
                    alphay * dot( transpose( Vy * dmdy) , (wcy * Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (wcz * Vz * dmdz) )
                    
                scale_xyz = phi_xyz / phi_cxyz
            
            
                # Minimum support on model | m |q
                wcm =  ( ( abs(model) ) ** ( 2 - q ) + delta ) ** -1

            
                # Compute norm of |Wm| and scale on |GRAD(m)|
                phi_cm = dot( transpose(Ws * logmodel) , (wcm * Ws * logmodel))
            
                scale_cm = (scale_xyz * phi_cxyz) / (alphas * phi_cm)
            
                phi_m = scale_cm * alphas *phi_cm + scale_xyz*phi_cxyz   
            
                # Apply scale on weighting vectors... Ready to export
                wcm = (l * scale_cm * wcm)
                wcx = abs(2.0 - l) * (scale_xyz * wcx)
                wcy = abs(2.0 - l) * (scale_xyz * wcy)
                wcz = abs(2.0 - l) * (scale_xyz * wcz)
                
                # Write weight file
                w_out = concatenate((transpose(wcm),wcx,wcy,wcz))
                savetxt(Work_dir + 'input_w.dat',w_out)
                
                # Adjust trade-off parameter for next iteration
                if (phid < ndata*2) & (oo >0):
                    beta =cool_beta *beta;
                elif oo >= 0: 
                    beta = cool_beta * beta;
                    
                # Write input file
                fid = open(Work_dir + 'lp_ctrl_file.inp','w')
                fid.write('1 0 !!! restart=1, max # of iterations\n')
                fid.write('2, %12.8f !!! mode chifact \n'% beta)
                fid.write('%s%s  ! observations file\n'% (Home_dir,obsfile))
                fid.write('%s%s  !!! mesh or ncel,aspr\n'% (Home_dir,meshfile));
                fid.write('%s%s !!! topography (NULL = no topography)\n'% (Home_dir,topofile))
                fid.write('%s%s%s !!! initial model (NULL = same as reference model)\n'%('dcinv3d_',strnum,'.con'))
                fid.write('%s !!! reference model (NULL = calculate)\n'%(str_ref_model))
                fid.write('null !!! active cells file (NULL = all cells are active)\n')
                fid.write('BOUNDS_NONE !!! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file  ! bounds\n')
                fid.write('%d %d %d  \n'%(Lx,Ly,Lz))
                fid.write('daub2 !!! sensitivity wavelet type: Daubechies-4, 2 vanishing moments\n')
                fid.write('1 0.02 !!! relative construction error applied\n')
                fid.write('input_w.dat  ! weighting file\n')
                fid.write('0 !!! store sensitivity matrix in disk\n')
                fid.write('1e-005 !!! tolerance\n')
                fid.write('-1 !!! no limit to use memeory\n')
                fid.close();
                
                
                oo=oo+1
                
                os.chdir(Work_dir)
                # Call UBC Inversion
                fid = open(Work_dir + 'DClp_3796.pbs','w')
                fid.write('#PBS -l nodes=g11:ppn=4:core\n')
                fid.write('# the first line requests cluster node, 12 cores on the gold cluster\n')
                fid.write('\n')
                fid.write('cd $PBS_O_WORKDIR\n')
                fid.write('\n')
                fid.write('\n')
                fid.write('export OMP_NUM_THREADS=12\n')
                fid.write('\n')
                fid.write('dcinv3d lp_ctrl_file.inp\n')
                fid.close();
                
                #os.system('qsub DClp_3796.pbs')
                os.system('dcinv3d lp_ctrl_file.inp')
                
                #Extract last phid
                ndata, phid, Lxyz, beta_in, ref_model, str_ref_model = read_UBC_log(Work_dir,Home_dir,1,mcell)
                
                
                if phid ==99999:
                    print('Inversion has stopped at convergence. Iteration: %i\n'%oo)
                else:
                    
                    strnum = '01'
                    # Copy and rename result to "new_dir"
                    os.system('cp dcinv3d_01.con '+ new_dir + 'dcinv3d_0'+str(oo)+'.con')
                    os.system('cp dcinv3d.pre '+ new_dir + 'dcinv3d_0'+str(oo)+'.pre')
                    os.system('cp dcinv3d.log '+ new_dir + 'dcinv3d_0'+str(oo)+'.log')
                
                    model = loadtxt(Work_dir + 'dcinv3d_' + strnum + '.con')
                
                
                