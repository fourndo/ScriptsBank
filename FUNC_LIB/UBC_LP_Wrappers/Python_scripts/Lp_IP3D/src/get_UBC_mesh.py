'''
Created on Jul 17, 2013

@author: dominiquef
'''
def get_UBC_mesh(meshfile):
    """ Read UBC mesh file and extract parameters
         Works for the condenced version (20 * 3) --> [20 20 20] """ 

    fid = open(meshfile,'r')
    
    
    # Go through the log file and extract data and the last achieved misfit
    for ii in range (1, 6): 
                        
        line = fid.readline()
        line = line.split(' ')    
        
            # First line: number of cells in i, j, k 
        if ii == 1:
            
            numcell=[]
                
            for jj in range(len(line)):
                t = int(line[jj])
                numcell.append(t)
                
            
            # Second line: origin coordinate (X,Y,Z)
        elif ii==2:
                
            origin = []
            
            for jj in range(len(line)):
                t = float(line[jj])
                origin.append(t)
                
            
        # Other lines for the dX, dY ,dZ
        elif ii==3:
            
            dX=[]
            
            for jj in range(len(line)-1):
                
                if line[jj].find('*') != -1:
                    
                    ndx = line[jj].split('*')
                    
                    for kk in range(int(ndx[0])):
                        dX.append(ndx[1])
                         
                else:
                    
                    t = float(line[jj])
                    dX.append(t)
                        
        elif ii==4:
            
            dY=[]
            
            for jj in range(len(line)-1):
                
                if line[jj].find('*') != -1:
                    
                    ndy = line[jj].split('*')
                    
                    for kk in range(int(ndy[0])):
                        dY.append(ndy[1])
                         
                else:
                    
                    t = float(line[jj])
                    dY.append(t)
                    
                    
        elif ii==5:
            
            dZ=[]
            
            for jj in range(len(line)-1):

                if line[jj].find('*') != -1:
                    
                    ndz = line[jj].split('*')
                    
                    for kk in range(int(ndz[0])):
                        dZ.append(ndz[1])
                         
                else:
                    
                    t = float(line[jj])
                    dZ.append(t)  
                        
            return numcell, origin, dX, dY, dZ
