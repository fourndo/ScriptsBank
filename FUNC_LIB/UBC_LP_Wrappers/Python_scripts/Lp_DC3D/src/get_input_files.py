'''
Created on Jul 17, 2013

@author: dominiquef
'''
def get_input_files(Home_dir):
    
    topo_check = []
    topofile = []
    obsfile = []
    meshfile = []
    
    import os 
    for directories in os.listdir(Home_dir):
        if  "msh" in directories.lower():
            
            meshfile = directories
            
        elif  "topo" in directories[-4:]: 
              
            topofile = directories
            
        elif "topo_model.txt" in directories.lower():
            
            topo_check =  directories;
            
        elif "obs" in directories.lower():
            
            obsfile = directories
            

            
    if not topofile:
        print "Topo file is missing in directory...Please verify"
        return meshfile, obsfile, [], []
    
    
    
    else:
        if not topo_check:
            #topo_command = 'topocheck ' + meshfile + ' ' + topofile 
            chg_dir = 'cd ' + Home_dir
            print chg_dir
            os.chdir(Home_dir)
        
            os.system('topocheck ' + meshfile + ' ' + topofile)
        

        # Load topocheck file and create active cell [array]
        fid = open(Home_dir + 'topo_model.txt','r')           
        topo = fid.read()
        topo = topo.split('\n')
        active_cells = []
        for tt in range(len(topo)-1):
            active_cells.append(int(topo[tt]))
        
        if not meshfile:
            print "Mesh file is missing in directory...Please verify !!"
            return [], obsfile, topofile, active_cells
        elif not obsfile: 
            print "Observation file is missing in directory...Please verify"
            return meshfile, [], topofile, active_cells
        else:   
            return meshfile, obsfile, topofile, active_cells