'''
Created on Jul 18, 2013

@author: dominiquef
'''
def comp_dmdx( m , dX , dY , dZ , actv_cells ):
    """Computes dm/dx in the other of UBC mesh""" 
   
    from numpy import zeros
    
    nX = len(dX)
    nY = len(dY)
    nZ = len(dZ)
    
    dmdx = zeros(((nX-1)*nY*nZ,1), dtype=float);
    Vx = zeros(((nX-1)*nY*nZ,1), dtype=float);
    
    count = 0
    skip = 0
    
    for jj in range(nY):
        
        for ii in range(nX):
            
            for kk in range(nZ):
                
                if (ii == nX-1):
                    skip = skip + 1
                
                else:
                                        
                    if (actv_cells[count]==0) | (actv_cells[count+skip+nZ]==0):
                        count = count + 1
                    
                    else:
                        dXmid = ( float(dX[ii]) + float(dX[ii+1]) ) /2
            
                        dmdx[count]= (-m[count+skip] + m[count+skip+nZ] )  / dXmid 
                    
                        Vx[count]= (float(dZ[kk]) * float(dY[jj]) * dXmid )** 0.5
                        count=count+1
                    
    
    return dmdx, Vx