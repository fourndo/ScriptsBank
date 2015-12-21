% Create depth weight on octree mesh
clear all 
close all


work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Inv25_noref_conic';
meshfile = 'octree_mesh.txt';
activefile = 'active_cells_topo.txt';
distw = 50;

% Create GIFProject and load mesh
proj = GIFproject();

% load mesh and extract cell center location
Omesh = meshOctree(proj);
Omesh.readFile([work_dir '\' meshfile]);

XYZ = fetchCenter(Omesh);

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
weights = zeros(length(active),1);

%% Create volume weights
% Wr = GIFmodel(proj);
% Wr.setMesh(Omesh);
% 
% % Get volume elements
% V = fetchVolume(Omesh);
% 
% Wr.setValue(sqrt(V));
% Wr.writeFile([work_dir '\Vweight.dat']);

%% Create weight model and populate with "layers"
Wr = GIFmodel(proj);
Wr.setMesh(Omesh);

% Cycle through the active cells and assign weights
zmin = min( XYZ(nullcell.value==1,3) );
zmax = 451;%max( XYZ(nullcell.value==1,3) );

% Create weighting function
% weights(active) = exp((XYZ(active,3) - zmax)/20) * 64 + 1;
% weights(weights>64) = 64;
weights = ones(length(active),1) * 1e-8;
weights(active) = (...
    ((7133602.000000-XYZ(active,2)).^2 +...
    (557299.312500 - XYZ(active,1)).^2).^0.25 ./...
    ( sqrt(2)/2*abs(zmax - XYZ(active,3))) / 2).^(-0.5);
weights(active) = (weights(active) / min(weights(active)));

Wr.setValue(weights);
% Save weights
Wr.writeFile([work_dir '\Cone_weight.dat']);

%% Or create synthetic model
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
% load([work_dir '\m_1D']);
% % Specify Z of 1D model
% Ztopo = 421;
% 
% m = GIFmodel(proj);
% m.setMesh(Omesh);
% 
% model = ones(length(active),1)*1e-8;
% model(active) = 1e-4;
% 
% % Get unique cell center values
% cellZ = unique(XYZ(active,3));
% z1D = Ztopo - cumsum(m_1D(:,2));
% interp_m = zeros(length(cellZ),1);
% 
% for ii = 1 : length(cellZ)
%     
%     index = XYZ(:,3) == cellZ(ii);
%     
%     interp_m(ii) = interp1( z1D ,  m_1D(:,1) , cellZ(ii) );
%     
%     if isnan(interp_m(ii)) == 1
%         
%         interp_m(ii) = 1e-4;
%         
%     end
%     
%     model(index) = interp_m(ii);
%     
% end
% 
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
% model(active==0) = 1e-8;
% 
% m.setValue(model);
% % Save weights
% m.writeFile([work_dir '\Model_1D_interp.dat']);