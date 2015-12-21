'''
Created on Jul 18, 2013

@author: dominiquef
'''
def read_UBC_log_NEWCODE(Work_dir,iter_start):
    """ Read UBC log file and extract parameters""" 
    
    import re
  
    fid = open(Work_dir + 'gzinv3d.log','r')
    
    line = fid.readline()
    iteration = 1
    phid = -1
    beta = -1
    while line:
        
        
        if  "# of data:" in line:
            ndata = float(line[22:])

                    
        if  "Iteration:" in line:
            iteration = int(line[21:])
            phid = -1

                    
        if  "multiplier:" in line:
            beta = float(line[21:])
                    
        if  "Le, Ln, Lz:" in line:

            L_scale=re.findall(r'\d+',line)
            Lxyz=[]
            for jj in range(len(L_scale)):
                if float(L_scale[jj])!=0:
                    Lxyz.append(float(L_scale[jj]))
                    
        if  "data misfit:" in line:
            phid = float(line[21:])
                    
        if (iteration==iter_start) & (phid>0) & (beta>0) :                    
            break
                      
        line = fid.readline()
    
    if not line:
        phid = 99999
           
    return ndata, phid, Lxyz, beta
        