'''
Created on Jul 17, 2013

@author: dominiquef
'''
def create_dir(p,q,l,Home_dir):
    """ Create new work directory for inversion using
    the input parameters for name """ 
    
    file_list = []
    new_dir = 'lp' + p + 'lq' + q + 'lambda' + l
    
    import os 
    for directories in os.listdir(Home_dir):
        file_list.append(directories)
        if directories.match(new_dir)==1:
            print 'MATCH'
            
        else:
            print 'CREATE DIR:' + new_dir
            
            return new_dir
    
    
    
    
    
    
    