'''
Created on Jul 18, 2013

@author: dominiquef
'''
def comp_wxm(mcell,dX,dY,dZ,actv_cells):
    """Computes dm/dx in the other of UBC mesh""" 
   
    from numpy import zeros
    
    nX = len(dX)
    nY = len(dY)
    nZ = len(dZ)
    
    Ws = zeros((mcell,1), dtype=float);
    
    count = 0
    
    for jj in range(nY):
        
        for ii in range(nX):
            
            for kk in range(nZ):
                
                if actv_cells[jj,ii,kk]==0:
                    count = count + 1
                    
                else:
                    
                    Ws[count]= (float(dZ[kk]) * float(dY[jj]) * float(dX[ii]) )** 0.5
                    count = count + 1
                    
    
    return Ws