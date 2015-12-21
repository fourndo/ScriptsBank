'''
Created on Jul 18, 2013

@author: dominiquef
'''
def make_depthW(layers,nZ,nX,nY,actv_cells):
    """ Make depth weighting matrices in vector form, UBC style"""
    from numpy import zeros 
    from numpy import resize
    from numpy import size
    
    wzx = zeros((nZ*(nX-1)*nY,1), dtype=float)
    wzy = zeros((nZ*nX*(nY-1),1), dtype=float)           
    
    topo = resize(actv_cells,(nY,nX,nZ))
    
    # Build wzx
    count_a = 0
    for jj in range(nY):
        
        for ii in range(nX-1):
            
            count_b = 0
            
            for kk in range(nZ):
                
                
                
                if topo[jj,ii,kk]==1:
                    
                    
                    if count_b<=(size(layers)-1):
                        
                        wzx[count_a] = layers[count_b]
                        count_b = count_b + 1
                    else:
                    
                        wzx[count_a] = 1
                       
                count_a = count_a + 1
                
    # Build wzy
    count_a = 0            
    for jj in range(nY-1):
        
        for ii in range(nX):
            
            count_b = 0
            
            for kk in range(nZ):
                
                if topo[jj,ii,kk]==1:
                    
                    if count_b<=(len(layers)-1):
                        
                        wzy[count_a] = layers[count_b]
                        count_b = count_b + 1
                    else:
                    
                        wzy[count_a] = 1
                        
                count_a = count_a + 1
    
    return wzx, wzy
                    
                
                