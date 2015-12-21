% Create depth weight on octree mesh
clear all 
close all


work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Tile_DO27_DO18\Coarse\Inv3_900';
meshfile = 'octree_large.msh';
activefile = 'active_cells_topo.txt';
obsfile = 'dighem_dwns_25m_900_flr2ppm_TF.dat';
sensfile = 'Speudo_sens.txt';

distw = 50;

% Create GIFProject and load mesh
proj = GIFproject();

% Load obsfile
[trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% load mesh and extract cell center location
Omesh = meshOctree(proj);
Omesh.readFile([work_dir '\' meshfile]);

XYZ = fetchCenter(Omesh);
V = fetchVolume(Omesh);


%% Load active cell topo file
nullcell = GIFmodel(proj);
nullcell.setMesh(Omesh);
nullcell.readFile([work_dir '\' activefile]);
active = nullcell.value==1;
weight = ones(length(active),1);


%% Load speudo sensitivities and/or create conic weights

%% Only use if you have speudo-sensitivty computed
% pseudoJ = GIFmodel(proj);
% pseudoJ.setMesh(Omesh);
% readFile(pseudoJ,[work_dir '\' sensfile]);
% 
% weight = getValue(pseudoJ);
% weight = abs(weight).^-1;
% 
% % Normalize
% weight(active) = (weight(active) / max(weight(active)));
%%

%% Or/else compute inverse distance conic weighting
% Add conic weighting
hweight = zeros(sum(active),1);

for ii = 1 : size(data,1)
    
    dx = data(ii,3) - XYZ(active,2);
    dy = data(ii,2) - XYZ(active,1);
    dz = abs(data(ii,4) - XYZ(active,3));
    
    r = sqrt(dx.^2 + dy.^2);
    
    temp =  1./ (1 + r)  ;
    
    hweight = hweight + temp ;
    
end

hweight = sqrt(hweight / max(hweight));

% Add coninc to weights
% weight(active) = weight(active) .* hweight;
weight(active) = hweight;
% Normalize again
% weight(active) = ( weight(active) / max(weight(active)) );

% Paint air cells
weight(active==0) = 1e-8;
Wr = GIFmodel(proj);
Wr.setMesh(Omesh);

Wr.setValue(weight);
% Save weights
Wr.writeFile([work_dir '\Weighting_sqrt_cone.dat']);
