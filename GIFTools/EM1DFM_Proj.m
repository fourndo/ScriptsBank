clear all
dsep = '\';
work_dir = 'C:\Users\DominiqueFournier\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\DIGHEM_TKC_Line10';
data1D = 'DIGHEM_line10.obs';
topomodel = 'nullcell.dat';
mesh3Dfile = 'UBC_mesh_small_TOP.msh';
mesh1Dfile = 'start.con';

% Create project and load core files
proj = GIFproject;
proj.setWorkDir(work_dir);
% mesh3D = mesh3D(proj); mesh3D.readFile([proj.workDir dsep mesh3Dfile]);

% Create 1D mesh
mesh1D = mesh1D(proj); 
mesh1D.readFile([proj.workDir dsep mesh1Dfile]);

% Import 1D model
mod = GIFmodel(proj);
mod.setMesh(mesh1D);
mod.readFile([proj.workDir dsep mesh1Dfile]);

% ALTERNATIVE - Load model and mesh at once
mod = GIFmodel.importGIFMeshModelFiles(proj,[proj.workDir dsep mesh1Dfile],[proj.workDir dsep mesh1Dfile],'cond','test');

% topo = TOPOdata(proj); topo.readFile([proj.workDir dsep topofile]);

% Create EM1DFM inversion object and pair
inv = EM1DFMinversion(proj);
data = EM1DFMinversion.readDobs(inv, [proj.workDir dsep data1D]);
inv.setDobs(data)

% Assign the mesh1D, which will generate a 3D mesh and ijk link with data
inv.setMesh(mesh1D);

% Write input file
EM1DFMinversion.writeInp(inv)

% Write model to file
EM1Dinversion.write1Dmodel(inv,1,1e-3,[proj.workDir dsep 'em1dfm.con'])

% Write data out
inv.writeDobs(inv.getDobs, [proj.workDir dsep 'em1dfm.obs'], inv.getDataOptions)

% Run inversion
EM1DFMinversion.runInversion(inv,true);

% Read in the inverted conductivities in 1D format
invMod = inv.readModel([proj.workDir dsep 'em1dfm_con.mod']);

% Map the 1D models to 3D mesh
mod3D = GIFmodel(inv); mod3D.setMesh(inv.mesh);
inv.setConductivityModel(mod3D, invMod, inv.data2ijk);

% Write out the projected 3D conductivity
inv.getConductivityModel.writeFile([proj.workDir dsep 'Inv3Dmodel.con']);

% Setup interpolation matrix
inv.setInterpolationMatrix( 3, 200, 1, 1);

% Interpolate from soundings to full mesh
m3D = inv.interp1D_to_3D(inv.getItem{end}.getValue);
% Write to file
mod3D.setValue(m3D);
mod3D.writeFile([proj.workDir dsep 'Inv1DInterp3D.con']);


