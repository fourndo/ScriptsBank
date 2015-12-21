% InterpModelDriver
% Original script from Seogi Kang@UBC-GIF
% Takes an octree model with bad topography and transfers it onto
% a different octree mesh
% Last update: June 13th,2014

clear all
close all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\VTEM\3D\DO27\FWR';
data_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Tile_DO27_DO18\Inv1_sqrtWcone_ay3'

Octreemesh      = '..\octree_mesh_big.txt';
Tensormesh      = 'VTEM_DO27.msh';
model           = 'inv_015.con';

p = GIFproject;
GIFtools(p);
%%
%Get octree mesh and model 
meshOct = meshOctree(p);
meshOct.readFile([data_dir '\' Octreemesh]);

%%
modoctree = GIFmodel(p);
modoctree.setMesh(meshOct);
modoctree.readFile([data_dir '\' model]);

val = modoctree.value;
val(val>1e-3) = val(val>1e-3)*10;
modoctree.setValue(val);

%% Get second mesh 
Tmesh = mesh3D(p);
Tmesh.readFile([work_dir '\' Tensormesh]);

% Interpolate onto mesh
modtensor = GIFmodel(p, Tmesh);
modtensor.interpFromModel(modoctree, 1, 0);
modtensor.setName('sigma_homo_tensor')


%% Write to file
modtensor.writeFile([work_dir '\3Dmodel_' model(1:end-4) '.dat']);



