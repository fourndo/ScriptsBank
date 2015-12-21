'''
Created on Jul 18, 2013

@author: dominiquef
'''
def read_UBC_log_primer(Work_dir,Home_dir,iter_start,mcell):
    """ Read UBC log file and extract parameters""" 
    
    import re
    from numpy import zeros, ones
    
    fid = open(Work_dir + 'maginv3d.log','r')
    
    line = fid.readline()
    iteration = 1
    phid = -1
    while line:
        
        
        if  "# of data:" in line:
            ndata = float(line[25:])

                    
        if  "Iteration:" in line:
            iteration = int(line[22:])
            phid = -1

                    
        if  "multiplier:" in line:
            beta = float(line[22:])
                    
        if  "Le, Ln, Lz:" in line:

            L_scale=re.findall(r'\d+',line)
            Lxyz=[]
            for jj in range(len(L_scale)):
                if float(L_scale[jj])!=0:
                    Lxyz.append(float(L_scale[jj]))
                    
        if  "data misfit:" in line:
            phid = float(line[22:])
            
        if  "Reference Model:" in line:
            str_ref_model = line[22:].rstrip()
            look_value = re.findall(r'\d+',line)
            if not look_value:
                if "null" in line.lower():
                    ref_model = zeros((mcell,1), dtype=float)
                else:
                    ref_file = Home_dir + str_ref_model
                    fid2 = open(ref_file,'r')
                    ref_model=zeros((mcell,1), dtype=float)
                    limit = mcell-1

                    for ii in range (limit):         
                        line2 = fid2.readline()
                        #line = line.split('\n')
                        ref_model[ii]=float(line2)
                        
                    fid2.close()
            else:
                ref_model = ones((mcell,1), dtype=float) * float(line[22:])
                
                    
        if (iteration==iter_start) & (phid>0):                    
            break
                      
        line = fid.readline()
    
    if not line:
        phid = 99999
           
    return ndata, phid, Lxyz, beta, ref_model
        