clear; clc
cd 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Tile_DO27_DO18\Newmesh'
p = GIFproject;
GIFtools(p);

%%
%Get octree mesh and model 
m = meshOctree(p);
m.readFile('octree_mesh_big.txt');

%%
Resoctree = GIFmodel(p);
Resoctree.setMesh(m);
Resoctree.readFile('TKC_parametric.con');

%%
% inactind = Resoctree.value==1e-8;
val = Resoctree.value;
% val(inactind) = 1/1.38221676e03;
Sigoctree = GIFmodel(p);
Sigoctree.setMesh(m);
Sigoctree.setValue(val);
Sigoctree.setName('TKC_parametric_Mike')
%% Get tensor mesh 
mtensor = mesh3D(p);
mtensor.readFile('UBC_mesh_small.msh');

%Get topo object
topo = TOPOdata(p);
topo.importXYZ('CDED_076c05_NAD27.topo');
%%
%Get active model for octreemesh and tensormesh
active_octree = GIFmodel.createActiveCellsFromTopo(p,m,topo,1);
active_octree.setName('ActiveModel_octree')

active_tensor = GIFmodel.createActiveCellsFromTopo(p,mtensor,topo,1);
active_tensor.setName('ActiveModel_tensor')
%%

val(active_octree.value==0) = 1e-8;
SigoctreeAir = GIFmodel(p);
SigoctreeAir.setMesh(m);
SigoctreeAir.setValue(val);
SigoctreeAir.setName('sigma_air_octree')

% Interpolate sigma model from octree to tensor mesh
%%
Sigtensor = GIFmodel(p, mtensor);
Sigtensor.interpFromModel(Sigoctree, 1, 0);
Sigtensor.setName('sigma_homo_tensor')
%%
SigtensorAir = GIFmodel(p, mtensor);
val = Sigtensor.value;
val(active_tensor.value==0) = 1e-8;
SigtensorAir.setValue(val);
SigtensorAir.setName('sigma_air_tensor')

% 
%% Project all column models to top
temp = SigtensorAir.value;
nullcell = active_tensor.value;

nx = length(mtensor.dx);
ny = length(mtensor.dy);
nz = length(mtensor.dz);

temp = reshape(temp,nz,nx,ny);
nullcell = reshape(nullcell,nz,nx,ny);

for ii = 1 : nx
    for jj = 1 : ny
        
        ground = temp(nullcell(:,ii,jj)==1,ii,jj);
        temp(nullcell(:,ii,jj)==0,ii,jj) = ground(1);
        
    end
end
SigtensorProj2top = GIFmodel(p, mtensor);
SigtensorProj2top.setValue(reshape(temp,nz*nx*ny,1));
SigtensorProj2top.setName('sigma_2_top');

% GIFtools(p)
%%
mtensor.writeFile('mesh_tensor.msh');
SigoctreeAir.writeFile('sig_octree.con');
SigtensorAir.writeFile('sig_tensor.con');
SigtensorProj2top.writeFile('sig_tensor_proj2top.con');
active_octree.writeFile('active_octree.txt');
active_tensor.writeFile('active_tensor.txt');


