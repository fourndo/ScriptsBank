% Create depth weight on octree mesh
clear all 
close all


work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Real_topo\Inv9_FixBeta_1DRef';
meshfile = 'octree_large.msh';
actv = 'active_cells_topo.txt';
mod1 = 'DIGHEM_Inv10.con';
mod2 = 'PARA_octree.con';
wght = 'DOI_octree.con';
distw = 50;

% Create GIFProject and load mesh
proj = GIFproject;

GIFtools(proj);
% Load obsfile
% [trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% load mesh and extract cell center location
Omesh = meshOctree(proj);
Omesh.readFile([work_dir '\' meshfile]);


%% Create active cell
% activecell = GIFmodel(proj);
% activecell.setMesh(Omesh);
% ground = XYZ(:,3)<421;
% 
% activecell.setValue(ground.*1);
% activecell.writeFile([work_dir '\nullcell.dat']);

%% Load active cell topo file
nullcell = GIFmodel(proj);
nullcell.setMesh(Omesh);
nullcell.readFile([work_dir '\' actv]);

mactv = nullcell.value==1;

% Load model A
mconA = GIFmodel(proj);
mconA.setMesh(Omesh);
mconA.readFile([work_dir '\' mod1]);

% Load model B
mconB = GIFmodel(proj);
mconB.setMesh(Omesh);
mconB.readFile([work_dir '\' mod2]);

% Load weight file
mwght = GIFmodel(proj);
mwght.setMesh(Omesh);
mwght.readFile([work_dir '\' wght]);

m1 = mconA.value;
m2 = mconB.value;
w = mwght.value;

mref = m1;
mref(mactv) = (m1(mactv).*(1-w(mactv)) + m2(mactv).*w(mactv)) ./ ( w(mactv) + (1-w(mactv)));

mout = GIFmodel(proj);
mout.setMesh(Omesh);
mout.setValue(mref);

% Save weights
mout.writeFile([work_dir '\VTEM_DIGHEM_DOI.con']);