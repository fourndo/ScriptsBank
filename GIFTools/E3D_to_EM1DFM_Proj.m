clear all
dsep = '\';
work_dir = 'C:\Users\DominiqueFournier\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\DIGHEM_TKC_Line10';
data3D = 'DIGHEM_E3D.obs';
% topomodel = 'nullcell.dat';
toposurf = '..\TKCtopo.dat';
mesh3Dfile = 'UBC_mesh_small_TOP.msh';
mesh1Dfile = 'start.con';


% Create project and load core files
proj = GIFproject;
proj.setWorkDir(work_dir);

% mesh3D = mesh3D(proj); mesh3D.readFile([proj.workDir dsep mesh3Dfile]);
data3DObj = FEMxyz.readDobs(proj,[work_dir dsep data3D]);

% Read in topo
topo = TOPOdata(proj);
topo.readFile([work_dir dsep toposurf]);

%% WRITE OUT EM1DFM
% OPTION 1: Calculate height above topo and write to file
data3DObj.calculateHeightAboveSurface(topo,'HEIGHT')
% data3DObj.writeEM1DFMobs([work_dir dsep 'E3D_to_EM1DFM.obs']);

% OPTION 2: Convert to a data1D object and write to file
data1D = data3DObj.convertXYZtoSounding();
% data1D.writeDobs([work_dir dsep 'E3D_to_EM1DFM.obs']);

%%
% Create 1D mesh
msh1D = mesh1D(proj); 
dz = [ones(1,10) 1.3.^[1:8] ones(1,20)*10 10*1.3.^[1:8] ones(1,10)*100 0];
msh1D.setdz(dz);
msh1D.setz0(max(topo.getioData(:,'Z')));
% Create a 1D model and write to file
mod = GIFmodel(proj);
mod.setMesh(msh1D);
mod.setValue(1e-2);
% Write out to file
msh1D.writeFile([work_dir '\em1dfm.con']);

inv = EM1DFMinversion(proj);
inv.setDobs(data1D);
inv.setMesh(msh1D);
inv.setTopo(topo);
inv.setInterpNstn(5);

% A mesh 3D has been generated -> write to file
msh3D = inv.getMesh();
msh3D.writeFile([proj.workDir dsep 'Mesh3D.msh'])

% Write input files
% EM1DFMinversion.writeInp(inv)

% Write data out
% inv.writeDobs(inv.getDobs, [proj.workDir dsep 'em1dfm.obs'], inv.getDataOptions)

% Run inversion
% EM1DFMinversion.runInversion(inv,true);
inv.setConductivityModel(ones(inv.mesh.fetchNc,1)*1e-3);

inv.runInversionLC;

% Read the predicted data
dpred = EM1DFMinversion.readPred(inv, [proj.workDir dsep 'em1dfm.prd']);

% Read the outfile
inv.readOut();

% Read in the inverted conductivities in 1D format
invMod = inv.readModel([proj.workDir dsep 'em1dfm_con.mod']);

% Map the 1D models to 3D mesh
m3D = inv.map1Dmodel3D(invMod,[]);

mod3D = GIFmodel(inv,inv.mesh,m3D);
mod3D.writeFile([proj.workDir dsep 'Inv3Dmodel.con']);

% Setup interpolation matrix
% inv.setInterpolationMatrix( 3, 200, 1, 1);

% Interpolate from soundings to full mesh
m3D = inv.interpModel1D(m3D);

% Write to file
mod3D.setValue(m3D);
mod3D.writeFile([proj.workDir dsep 'Inv1DInterp3D.con']);

% 
