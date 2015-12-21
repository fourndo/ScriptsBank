# EM_inv_run.py
# Author: D Fournier
# Last Update: July 16th, 2013
#
# This program runs a series of EM inversions on a unix cluster, as specified
# a list of lines.
 
import os 
import shutil

         
from read_input import read_input
from numpy import *

Home_dir, UBC_dir, iter_start, chifact, pvec, qvec, lvec, cool_beta, iter_max = read_input('Python_input.inp')

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
            # Stabilizing term
            delta = 1e-6
            file_list = []
            new_dir = Home_dir + 'lp' + str(p) + 'lq' + str(q) + 'lambda' + str(l)
            

            if not os.path.exists(new_dir): os.makedirs(new_dir)
            else: 
                shutil.rmtree(new_dir)
                os.makedirs(new_dir)
            
            shutil.copy(UBC_dir + 'maginv3d_00' + str(iter_start) + '.sus', Work_dir)
            #os.rename(work_dir + 'maginv3d_0' + str(iter_start) + '.con','maginv3d_01.con')
            shutil.copy(UBC_dir + 'maginv3d.log', Work_dir)     
                 
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
            
            print(Work_dir)
            # Read *.log and extract information
            ndata, phid, Lxyz, beta_in = read_UBC_log_NEWCODE(Work_dir,iter_start)

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
            
            lbda = l     # Scaling between smooth and compact
            beta = beta_in  # Trade-off parameter

            # Create derivative matrices
            wcx=ones(((nX-1)*nY*nZ,1), dtype=float)
            wcy=ones((nX*(nY-1)*nZ,1), dtype=float)
            wcz=ones((nX*nY*(nZ-1),1), dtype=float) 

            
            Ws = get_Ws3D(mcell,dX,dY,dZ,actv_cells);
            
            #wzx, wzy = make_depthW(layers,nZ,nX,nY,actv_cells);
            
            
            
                  
            #fid = open(Work_dir + 'dcinv3d_0' + str(iter_start) + '.con','r');
            #tempo = fid.read()
            #tempo = tempo.split('\n')
            #model = []
            #for tt in range(len(tempo)-1):
            #    model.append(float(tempo[tt]))
            #model = array(model)
            
            model = loadtxt(Work_dir + 'maginv3d_00' + str(iter_start) + '.sus')
            
            
    ## Inversion iterations
            # Recompute weights and adjust beta
            
            while (oo<= iter_max ) & (phid >= ndata*chifact):
                
                #logmodel = log10(model) - log10(ref_model)            
                #logmodel = logmodel * actv_cells
                       
                #model = model - ref_model
                model = model * actv_cells

                # Compute derivative vectors
                dmdx, Vx = comp_dmdx( model , dX , dY , dZ , actv_cells )
                dmdy, Vy = comp_dmdy( model , dX , dY , dZ , actv_cells )
                dmdz, Vz = comp_dmdz( model , dX , dY , dZ , actv_cells )
                
                phi_xyz = alphax * dot( transpose(Vx * dmdx) , (Vx * dmdx) ) + \
                    alphay * dot( transpose(Vy * dmdy) , (Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (Vz * dmdz) )
  
                # lp-norm on gradient ||GRAD(m)||p
                wcx = ( (abs(dmdx) ) ** ( 2.0 - p ) + delta ) ** -1
            
                wcy = ( (abs(dmdy) ) ** ( 2.0 - p ) + delta ) ** -1
               
                wcz = ( ( abs(dmdz) ) ** ( 2.0 - p ) + delta ) ** -1
                
                # Compute norm of | GRAD (m) |p and scale on | GRAD (m) |2
                phi_cxyz = alphax * dot( transpose( Vx * dmdx) , (wcx * Vx * dmdx) ) + \
                    alphay * dot( transpose( Vy * dmdy) , (wcy * Vy * dmdy) ) + \
                    alphaz * dot( transpose(Vz * dmdz) , (wcz * Vz * dmdz) )
                  
                scale_xyz = phi_xyz / phi_cxyz
                        
                # Minimum support on model | m |q
                wcm =  ( ( abs(model) ) ** ( 2 - q ) + delta ) ** -1
                
                # Compute norm of |Wm| and scale on |GRAD(m)|
                phi_cm = dot( transpose(Ws * model) , (wcm * Ws * model))
                
                scale_cm = (scale_xyz * phi_cxyz) / (alphas * phi_cm)
            
                phi_m = scale_cm * alphas *phi_cm + scale_xyz*phi_cxyz   
                
               
                # Apply scale on weighting vectors... Ready to export
                wcm = (l * scale_cm * wcm)
                wcx = abs(2-l) * (scale_xyz * wcx)
                wcy = abs(2-l) * (scale_xyz * wcy)
                wcz = abs(2-l) * (scale_xyz * wcz)
                
                
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
                #fid.write('0 !!! restart=1, max # of iterations\n')
                fid.write('2 !!! mode\n')
                fid.write('%12.8f 0.1 !!! mode chifact \n'% beta)
                fid.write('%s%s  ! observations file\n'% (Home_dir,obsfile))
                fid.write('%smaginv3d.mtx  !!! mesh or ncel,aspr\n'% (Home_dir));
                fid.write('%s%s%s%s !!! initial model (NULL = same as reference model)\n'%(Work_dir,'maginv3d_0',strnum,'.sus'))
                fid.write('VALUE 0 !!! reference model (NULL = calculate)\n')
                fid.write('null !!! active cells file (NULL = all cells are active)\n')
                #fid.write('-1 1 !!! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file  ! bounds\n')
                fid.write('VALUE 0 !!! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file  ! bounds\n')
                fid.write('VALUE 1 !!! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file  ! bounds\n')
                fid.write('%d %d %d  \n'%(Lx,Ly,Lz))
                fid.write('SMOOTH_MOD\n')
                fid.write('input_w.dat  ! weighting file\n')
                #fid.write('0')
                fid.close();
                
                
                oo=oo+1
                
                os.chdir(Work_dir)
                #os.system('ls')
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
                os.system('maginv3d_newbounds lp_ctrl_file.inp')
                
                #Extract last phid
                ndata, phid, Lxyz, beta_in = read_UBC_log_NEWCODE(Work_dir,1)
                
                
                if phid ==99999:
                    print('Inversion has stopped at convergence. Iteration: %i\n'%oo)
                else:
                    
                    strnum = '01'
                    # Copy and rename result to "new_dir"
                    os.system('copy maginv3d_001.sus '+ new_dir + '\maginv3d_0'+str(oo)+'.sus')
                    os.system('copy maginv3d_001.pre '+ new_dir + '\maginv3d_0'+str(oo)+'.pre')
                    os.system('copy maginv3d.log '+ new_dir + '\maginv3d_0'+str(oo)+'.log')
                
                    model = loadtxt(Work_dir + 'maginv3d_0' + strnum + '.sus')
                
                
                