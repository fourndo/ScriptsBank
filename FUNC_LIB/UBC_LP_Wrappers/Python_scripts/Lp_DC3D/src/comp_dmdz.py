'''
Created on Jul 18, 2013

@author: dominiquef
'''
def comp_dmdz( m , dX , dY , dZ , actv_cells ):
    """Computes dm/dy in the other of UBC mesh""" 
   
    from numpy import zeros
    
    nX = len(dX)
    nY = len(dY)
    nZ = len(dZ)
    
    dmdz = zeros((nX*nY*(nZ-1),1), dtype=float);
    Vz = zeros((nX*nY*(nZ-1),1), dtype=float);
    count = 0
    skip = 0
    
    for jj in range(nY):
        
        for ii in range(nX):
            
            for kk in range(nZ):
                
                if (kk == nZ-1):
                    skip = skip + 1
                
                else:
                                        
                    #if actv_cells[jj,ii,kk]==0:
                    #count = count + 1
                    
                    #else:
                    dZmid = ( float(dZ[kk]) + float(dZ[kk+1]) ) /2
            
                    dmdz[count]= (-m[count+skip] + m[count + skip + 1] )  / dZmid 

                    Vz[count]= (float(dY[jj]) * float(dX[ii]) * dZmid )** 0.5
                    
                    count=count+1
                    
    
    return dmdz, Vz