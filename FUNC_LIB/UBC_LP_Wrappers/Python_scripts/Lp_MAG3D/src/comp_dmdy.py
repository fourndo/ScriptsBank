'''
Created on Jul 18, 2013

@author: dominiquef
'''
def comp_dmdy( m , dX , dY , dZ , actv_cells ):
    """Computes dm/dy in the other of UBC mesh""" 
   
    from numpy import zeros
    
    nX = len(dX)
    nY = len(dY)
    nZ = len(dZ)
    
    dmdy = zeros((nX*(nY-1)*nZ,1), dtype=float);
    Vy = zeros((nX*(nY-1)*nZ,1), dtype=float);
    count = 0
    skip = 0
    
    for jj in range(nY):
        
        for ii in range(nX):
            
            for kk in range(nZ):
                
                if (jj == nY-1):
                    skip = skip + 1
                
                else:
                                        
                    if (actv_cells[count]==0) | (actv_cells[count + skip + nZ * nX]==0):
                        count = count + 1
                    
                    else:
                        dYmid = ( float(dY[jj]) + float(dY[jj+1]) ) /2
            
                        dmdy[count]= (-m[count+skip] + m[count + skip + nZ * nX] )  / dYmid 

                        Vy[count]= (float(dZ[kk]) * float(dX[ii]) * dYmid )** 0.5
                    
                        count=count+1
                    
    
    return dmdy, Vy