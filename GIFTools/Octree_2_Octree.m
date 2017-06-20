clear; clc
cd 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Real_topo\Inv11_Round1'
p = GIFproject;
GIFtools(p);

%% Get tensor mesh and model 
m = meshOctree(p);
m.readFile('octree_mesh_big.txt');

ConOct1 = GIFmodel(p);
ConOct1.setMesh(m);
ConOct1.readFile('VTEM_inv_004.con');

val = ConOct1.value;
noval = 1e-8;

%%
% inactind = Resoctree.value==1e-8;

% val(inactind) = 1/1.38221676e03;
% Sigoctree = GIFmodel(p);
% Sigoctree.setMesh(m);
% Sigoctree.setValue(val);
% Sigoctree.setName('TKC_parametric_Mike')
%% Get Octree mesh 
mOctree = meshOctree(p);
mOctree.readFile('octree_large.msh');

% %Get topo object
% topo = TOPOdata(p);
% topo.readFile('..\CDED_Lake_warp.topo');

%% Get active model for octreemesh and tensormesh

% active_octree = ACTIVEmodel.createFromTopo(p,mOctree,topo,1);
% active_octree.setName('ActiveModel_octree');

active_octree = GIFmodel(p, mOctree);
active_octree.readFile('active_cells_topo.txt');
% active_tensor = ACTIVEmodel.createFromTopo(p,m,topo,1);
% active_tensor.setName('ActiveModel_tensor');

%% Project all column models to top
% temp = ConOct1.value;
% % temp(temp==1e-8) = 5e-4; % Change null values outside the VOI
% nullcell = temp~=noval;
% 
% nx = length(m.dx);
% ny = length(m.dy);
% nz = length(m.dz);
% 
% temp = reshape(temp,nz,nx,ny);
% nullcell = reshape(nullcell,nz,nx,ny);
% 
% for ii = 1 : nx
%     for jj = 1 : ny
%         
%         ground = temp(nullcell(:,ii,jj)==1,ii,jj);
%         
%         if ~isempty(ground)
%             temp(nullcell(:,ii,jj)==0,ii,jj) = ground(1);
%         end
%         
%     end
% end
% 
% temp(temp==noval) = 1;%median(temp(temp~=noval));
% sigTensorProj2top = GIFmodel(p, m);
% sigTensorProj2top.setValue(reshape(temp,nz*nx*ny,1));
% sigTensorProj2top.setName('sigma_2_top');

%%
% 
% SigoctreeAir = GIFmodel(p);
% SigoctreeAir.setMesh(m);
% SigoctreeAir.setValue(val);
% SigoctreeAir.setName('sigma_air_octree')

%% Interpolate sigma model from octree to octree
sigOctree = GIFmodel(p, mOctree);
sigOctree.interpFromModel(ConOct1, 0, 0);
sigOctree.setName('sigma_homo_tensor')

val = sigOctree.value;
val(active_octree.value==0) = noval;

sigOctree.setValue(val);
sigOctree.writeFile('VTEM_R2_octree.con');
% active_octree.writeFile('active_octree.txt');
%%
% SigtensorAir = GIFmodel(p, mOctree);
% val = sigOctree.value;
% val(active_tensor.value==0) = 1e-8;
% SigtensorAir.setValue(val);
% SigtensorAir.setName('sigma_air_tensor')

% 


% GIFtools(p)
%%
% mOctree.writeFile('mesh_tensor.msh');
% SigoctreeAir.writeFile('sig_octree.con');
% SigtensorAir.writeFile('sig_tensor.con');
% sigTensorProj2top.writeFile('sig_tensor_proj2top.con');
% active_octree.writeFile('active_octree.txt');
% active_tensor.writeFile('active_tensor.txt');


