'''
Created on Jul 17, 2013

@author: dominiquef
'''
def get_UBC_mesh(meshfile):
    """ Read UBC mesh file and extract parameters
         Works for the condenced version (20 * 3) --> [20 20 20] """ 

    fid = open(meshfile,'r')
    from numpy import zeros
    
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
                
            nX = numcell[0]
            nY = numcell[1]
            nZ = numcell[2] 
            # Second line: origin coordinate (X,Y,Z)
        elif ii==2:
                
            origin = []
            
            for jj in range(len(line)):
                t = float(line[jj])
                origin.append(t)
                
            
        # Other lines for the dX, dY ,dZ
        elif ii==3:
            
            dX=zeros((nX,1), dtype=float)
            count_entry = 0;
            count = 0;
            while (count<nX):

                if line[count_entry].find('*') != -1:
                    
                    ndx = line[count_entry].split('*')
                    
                    for kk in range(int(ndx[0])):
                        dX[count] = (ndx[1])
                        count = count+1
                    count_entry=count_entry+1
                         
                else:
                    
                    t = float(line[count_entry])
                    dX[count]= t
                    count = count+1  
                    count_entry=count_entry+1
                        
        elif ii==4:
            
            dY=zeros((nY,1), dtype=float)
            count_entry = 0;
            count = 0;
            while (count<nY):
                print line[count_entry]
                if line[count_entry].find('*') != -1:
                    
                    ndy = line[count_entry].split('*')
                    
                    for kk in range(int(ndy[0])):
                        dY[count] = (ndy[1])
                        count = count+1
                    count_entry=count_entry+1
                         
                else:
                    
                    t = float(line[count_entry])
                    dY[count]= t
                    count = count+1  
                    count_entry=count_entry+1
                    
        elif ii==5:
            
            dZ=zeros((nZ,1), dtype=float)
            count_entry = 0;
            count = 0;
            while (count<nZ):

                if line[count_entry].find('*') != -1:
                    
                    ndz = line[count_entry].split('*')
                    
                    for kk in range(int(ndz[0])):
                        dZ[count] = (ndz[1])
                        count = count+1
                    count_entry=count_entry+1     
                else:
                    
                    t = float(line[count_entry])
                    dZ[count]= t
                    count = count+1  
                    count_entry=count_entry+1
            fid.close();                
            return nX, nY, nZ, origin, dX, dY, dZ
