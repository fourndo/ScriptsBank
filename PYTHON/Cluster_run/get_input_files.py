'''
Created on Jul 17, 2013

@author: dominiquef
'''
def get_input_files(Home_dir):
    
    topo_check = []
    obsfile = []
    meshfile = []
    
    import os 
    for directories in os.listdir(Home_dir):
        if  "msh" in directories.lower():
            
            meshfile = directories
            
            
        elif "topo_model.txt" in directories.lower():
            
            topo_check =  directories;
            
        elif "obs" in directories.lower():
            
            obsfile = directories
            
    

    if not topo_check:
        print "Topocheck file is missing in directory...Please verify"
        return meshfile, obsfile, [], []
        #topo_command = 'topocheck ' + meshfile + ' ' + topofile 
        #chg_dir = 'cd ' + Home_dir
        #print chg_dir
        #os.chdir(Home_dir)
        
        #os.system('topocheck ' + meshfile + ' ' + topofile)
        
        
    if not meshfile:
        print "Mesh file is missing in directory...Please verify !!"
        return [], obsfile, topo_check
    elif not obsfile: 
        print "Observation file is missing in directory...Please verify"
        return meshfile, [], topo_check
    else:   
        return meshfile, obsfile, topo_check