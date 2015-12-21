% Create depth weight on octree mesh
clear all 
close all


work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Inv38_18stn_square';
meshfile = 'octree_mesh.txt';
activefile = 'active_cells_topo.txt';
% obsfile = 'DIGHEMdwn.txt';
% sensfile = 'Speudo_sens.txt';
distw = 50;

% Create GIFProject and load mesh
proj = GIFproject();

% Load obsfile
% [trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% load mesh and extract cell center location
Omesh = meshOctree(proj);
Omesh.readFile([work_dir '\' meshfile]);

XYZ = fetchCenter(Omesh);
V = fetchVolume(Omesh);

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
nullcell.readFile([work_dir '\' activefile]);
active = nullcell.value==1;
weight = zeros(length(active),1);

%% Create synthetic model
% 
% m = GIFmodel(proj);
% m.setMesh(Omesh);
% 
% model = ones(length(active),1)*1e-8;
% model(active) = 1e-5;
% 
% xmin = 557200;
% xmax = 557400; 
% ymin = 7133500;
% ymax = 7133700;
% zmin = 300;
% zmax = 350;
% 
% selector = (XYZ(:,1) > xmin) & (XYZ(:,1) < xmax) &...
%             (XYZ(:,2) > ymin) & (XYZ(:,2) < ymax) &...
%             (XYZ(:,3) > zmin) & (XYZ(:,3) < zmax);
%         
% model(selector) = 1e-1;
% 
% m.setValue(model);
% % Save weights
% m.writeFile([work_dir '\Synthetic.dat']);

%% Load 1D model and interpolate on 3D mesh
load([work_dir '\m_1D']);
% Specify Z of 1D model
Ztopo = 423;

m = GIFmodel(proj);
m.setMesh(Omesh);

model = ones(length(active),1)*1e-8;
model(active) = 1e-4;

% Get unique cell center values
cellZ = unique(XYZ(active,3));
z1D = Ztopo - cumsum(m_1D(:,2));
interp_m = zeros(length(cellZ),1);

for ii = 1 : length(cellZ)
    
    index = XYZ(:,3) == cellZ(ii);
    
    interp_m(ii) = interp1( z1D ,  m_1D(:,1) , cellZ(ii) );
    
    if isnan(interp_m(ii)) == 1
        
        interp_m(ii) = 1e-4;
        
    end
    
    model(index) = interp_m(ii);
    
end

% xmin = 556900;
% xmax = 557700; 
% ymin = 7133200;
% ymax = 7134000;
% % zmin = 300;
% % zmax = 350;
% 
% selector = (XYZ(:,1) < xmin) | (XYZ(:,1) > xmax) | ...
%             (XYZ(:,2) < ymin) | (XYZ(:,2) > ymax) ;
%         
% model(selector) = 1e-4;
model(active==0) = 1e-8;

m.setValue(model);
% Save weights
m.writeFile([work_dir '\Model_1D_interp.dat']);