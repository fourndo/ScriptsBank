'''
Created on Jul 18, 2013

@author: dominiquef
'''
def get_Ws3D(mcell,dX,dY,dZ,actv_cells):
    
    from numpy import zeros
    
    nX = len(dX)
    nY = len(dY)
    nZ = len(dZ)
    
    #Ws = zeros((mcell,1), dtype=float);
    Ws = []
    count = 0
    
    for jj in range(nY):
        
        for ii in range(nX):
            
            for kk in range(nZ):
                
                if actv_cells[count]==0:
                    count = count + 1
                    Ws.append(0)
                    
                else:
                    
                    Ws.append((float(dZ[kk]) * float(dY[jj]) * float(dX[ii]) ) ** 0.5)
                    count = count + 1
                    
    
    return Ws
                    