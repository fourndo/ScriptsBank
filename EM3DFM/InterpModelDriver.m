% InterpModelDriver
% Original script from Seogi Kang@UBC-GIF
% Takes an octree model with bad topography and transfers it onto
% a different octree mesh
% Last update: June 13th,2014

clear all
close all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Tile_DO27_DO18\Newmesh';

p = GIFproject;
GIFtools(p);
%%
%Get octree mesh and model 
% meshOct1 = meshOctree(p);
mesh_in = mesh3D(p);
mesh_in.readFile([work_dir '\' 'UBC_mesh_small.msh']);

%%
% model_in = GIFmodel(p);
% model_in.setMesh(meshOct1);
% model_in.readFile([work_dir '\inv_17.con']);

model_in = GIFmodel(p);
model_in.setMesh(mesh_in);
model_in.readFile([work_dir '\Inv_1D_DIGHEM.con']);
% val = model_in.value;
% val(val<2e-4) = 2e-4;
% model_in.setValue(val);

%% Project all column models to top
temp = model_in.value;
nullcell = temp~=1e-8;

nx = length(mesh_in.dx);
ny = length(mesh_in.dy);
nz = length(mesh_in.dz);

temp = reshape(temp,nz,nx,ny);
nullcell = reshape(nullcell,nz,nx,ny);

for ii = 1 : nx
    for jj = 1 : ny
        
        ground = temp(nullcell(:,ii,jj)==1,ii,jj);
        temp(nullcell(:,ii,jj)==0,ii,jj) = ground(1);
        
    end
end
SigtensorProj2top = GIFmodel(p, mesh_in);
SigtensorProj2top.setValue(reshape(temp,nz*nx*ny,1));
SigtensorProj2top.setName('sigma_2_top');

%%
%%
% inactind = Resoctree.value==1e-8;
% val = 1./Resoctree.value;
% val(inactind) = 1/1.38221676e03;
% Sigoctree = GIFmodel(p);
% Sigoctree.setMesh(meshOct1);
% Sigoctree.setValue(Resoctree.value);
% Sigoctree.setName('sigma_homo_octree')
%% Get second octree mesh 
% meshOct2 = meshOctree(p);
% meshOct2.readFile([work_dir '\' 'octree_mesh.txt']);
% 
% XYZ = fetchCenter(meshOct2);
% V = fetchVolume(meshOct2);

%Get topo object


%% Load topo file
% topo = TOPOdata(p);
% topo.importXYZ([work_dir '\' 'CDED_076c05_NAD27.topo']);

%Get active model for octree meshs
% active_octree1 = GIFmodel.createActiveCellsFromTopo(p,meshOct1,topo,1);
% active_octree1.setName('ActiveModel_octree')
% 
% active_octree2 = GIFmodel.createActiveCellsFromTopo(p,meshOct2,topo,1);
% active_octree2.setName('ActiveModel_tensor')
%%

% val(active_octree1.value==0) = 1e-8;
% SigoctreeAir = GIFmodel(p);
% SigoctreeAir.setMesh(meshOct1);
% SigoctreeAir.setValue(val);
% SigoctreeAir.setName('sigma_air_octree')

% Interpolate sigma model from octree to octree mesh
%%
% Siginterp = GIFmodel(p, meshOct2);
% Siginterp.interpFromModel(Sigoctree, 1, 0);
% Siginterp.setName('sigma_homo_tensor')
% %%
% SigtensorAir = GIFmodel(p, meshOct2);
% val = Siginterp.value;
% 
% % Use active cells
% % val(active_octree2.value==0) = 1e-8;
% 
% % Or constant value
% val(XYZ(:,3) > 423) = 1e-8;
% 
% % Add surface conductivity (lake)
% % Or constant value
% val(XYZ(:,3) > 413) = 1e-5;
% 
% % Or constant value
% val(XYZ(:,3) > 423) = 1e-8;
% 
% SigtensorAir.setValue(val);
% SigtensorAir.setName('sigma_air_out')

%%
% meshOct2.writeFile([work_dir '\' 'mesh_octree.msh']);
% SigoctreeAir.writeFile([work_dir '\' 'sig_octree.con']);
% SigtensorAir.writeFile([work_dir '\' 'sig_octree_interp.con']);
% active_octree1.writeFile([work_dir '\' 'active_octree.txt']);
% active_octree2.writeFile([work_dir '\' 'active_octree.txt']);

%% Interpolate onto tensor
% mesh_out = mesh3D(p);
mesh_out = meshOctree(p);
mesh_out.readFile([work_dir '\octree_mesh_big.txt']);

% Interpolate onto mesh
model_out = GIFmodel(p, mesh_out);
model_out.interpFromModel(SigtensorProj2top, 0, 0);
% model_out.setName('sigma_homo_tensor')

XYZ = fetchCenter(mesh_out);
val = model_out.value;

val(XYZ(:,3) > 423) = 1e-8;
model_out.setValue(val);

model_out.writeFile([work_dir '\DIGHEM_1Dmodel_Tensor_5e3_Treshold.con']);


