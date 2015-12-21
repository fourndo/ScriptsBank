'''
Created on Jul 18, 2013

@author: dominiquef
'''
def read_UBC_log(Work_dir,Home_dir,iter_start,mcell):
    """ Read UBC log file and extract parameters""" 
    
    from numpy import loadtxt   
    import re
    from numpy import zeros, ones
  
    fid = open(Work_dir + 'dcinv3d.log','r')
    
    line = fid.readline()
    iteration = 1
    phid = -1
    while line:
        
        
        if  "# of data:" in line:
            ndata = float(line[26:])

                    
        if  "Iteration:" in line:
            iteration = int(line[15:])
            phid = -1

                    
        if  "tradeoff par:" in line:
            beta = float(line[20:])
                    
        if  "Le, Ln, Lz:" in line:

            L_scale=re.findall(r'\d+',line)
            Lxyz=[]
            for jj in range(len(L_scale)):
                if float(L_scale[jj])!=0:
                    Lxyz.append(float(L_scale[jj]))
                    
        if  "achieved misfit:" in line:
            phid = float(line[20:])
            
        if  "Reference Model:" in line:
            str_ref_model = line[25:].rstrip()
            look_value = re.findall(r'\d+',line)
            if not look_value:
                if "null" in line.lower():
                    ref_model = zeros((mcell,1), dtype=float)
                else:
                    ref_file = Home_dir + str_ref_model
                    ref_model = loadtxt(ref_file)
            else:
                ref_model = ones((mcell,1), dtype=float) * float(line[31:])
                    
        if (iteration==iter_start) & (phid>0):                    
            break
                      
        line = fid.readline()
    
    if not line:
        phid = 99999
           
    return ndata, phid, Lxyz, beta, ref_model, str_ref_model
        