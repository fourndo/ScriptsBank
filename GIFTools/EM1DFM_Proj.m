clear all
dsep = '\';
work_dir = 'C:\Users\DominiqueFournier\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\DIGHEM_TKC_Line10';
data1D = 'DIGHEM_line10.obs';
topomodel = 'nullcell.dat';
toposurf = '..\TKCtopo.dat';
mesh3Dfile = 'UBC_mesh_small_TOP.msh';
mesh1Dfile = 'start.con';


% Create project and load core files
proj = GIFproject;
proj.setWorkDir(work_dir);

% mesh3D = mesh3D(proj); mesh3D.readFile([proj.workDir dsep mesh3Dfile]);
data3D = 'E3DinvWrite.obs';
data3DObj = E3Dinversion.readDobs(proj,[work_dir dsep data3D]);

% Read in topo
topo = TOPOdata(proj);
topo.readFile([work_dir dsep toposurf]);

% Create 1D mesh
msh1D = mesh1D(proj); 
msh1D.readFile([proj.workDir dsep mesh1Dfile]);

% Import 1D model
mod = GIFmodel(proj);
mod.setMesh(msh1D);
mod.readFile([proj.workDir dsep mesh1Dfile]);

% Write 1D model and mesh two ways
msh1D.writeFile([work_dir '\mesh1Dout.con'], mod);
mod.writeFile([work_dir '\mod1Dout.con']);

% ALTERNATIVE - Load model and mesh at once
mod = GIFmodel.importGIFMeshModelFiles(proj,[proj.workDir dsep mesh1Dfile],[proj.workDir dsep mesh1Dfile],'cond','test');

% topo = TOPOdata(proj); topo.readFile([proj.workDir dsep topofile]);

% Create EM1DFM inversion object and pair
inv = EM1DFMinversion(proj);
data = EM1DFMinversion.readDobs(inv, [proj.workDir dsep data1D]);
inv.setDobs(data)

% Assign the mesh1D, which will generate a 3D mesh and ijk link with data
inv.setMesh(msh1D);
msh3D = inv.getMesh();
msh3D.writeFile([proj.workDir dsep 'Mesh3D.msh'])
% Write input file
EM1DFMinversion.writeInp(inv)

% Write model to file
EM1Dinversion.write1Dmodel(inv,1,1e-3,[proj.workDir dsep 'em1dfm.con'])

% Write data out
inv.writeDobs(inv.getDobs, [proj.workDir dsep 'em1dfm.obs'], inv.getDataOptions)

% Run inversion
EM1DFMinversion.runInversion(inv,true);

% Read the predicted data
EM1DFMinversion.readPred(inv, [proj.workDir dsep 'em1dfm.prd']);

% Read the outfile
inv.readOut();

% Read in the inverted conductivities in 1D format
invMod = inv.readModel([proj.workDir dsep 'em1dfm_con.mod']);

% Map the 1D models to 3D mesh
mod3D = inv.insert1Dto3Dmodel(invMod);%GIFmodel(inv); mod3D.setMesh(inv.mesh);
inv.setConductivityModel(mod3D);

% Write out the projected 3D conductivity
inv.getConductivityModel.writeFile([proj.workDir dsep 'Inv3Dmodel.con']);

% % Setup interpolation matrix
% inv.setInterpolationMatrix(1, 1);
% 
% % Interpolate from soundings to full mesh
% m3D = inv.interpModel1D(inv.getItem{end}.getValue);
% 
% % Write to file
% mod3D.setValue(m3D);
% mod3D.writeFile([proj.workDir dsep 'Inv1DInterp3D.con']);

% 
