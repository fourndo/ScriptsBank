% Create depth weight on octree mesh
clear all 
close all


work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Sens_calc';
meshfile = 'octree_mesh.txt';
activefile = 'active_cells_topo.txt';
obsfile = 'DIGHEMdwn_DO27.txt';

% Frequencies
freq = [900 7200 56000];


% Create GIFProject and load mesh
proj = GIFproject();

% load mesh and extract cell center location
Omesh = meshOctree(proj);
Omesh.readFile([work_dir '\' meshfile]);

% Get octree cell center location
XYZ = fetchCenter(Omesh);

% Get volume elements
V = fetchVolume(Omesh);

% Load active cell topo file
nullcell = GIFmodel(proj);
nullcell.setMesh(Omesh);
nullcell.readFile([work_dir '\' activefile]);
active = nullcell.value==1;
weights = zeros(length(active),1);

V = V(active);

%% Select every two points in core region
actvm = XYZ(active,:);
xmin = 5.57138140E+05; xmax = xmin + 300;
ymin = 7.13339100E+06; ymax = ymin + 426; 
zmax = 5.06000000E+02; zmin = zmax - 160;

core = actvm( actvm(:,1) > xmin & actvm(:,2) > ymin & actvm(:,3) > zmin , : );
core = core( core(:,1) < xmax & core(:,2) < ymax & core(:,3) < zmax , : );

% Keep only every 3 points
core = core(1:5:end,:);
V = V(1:5:end);

%% Load data/tx loc
[trx,d] = load_E3D_obs([work_dir '\' obsfile]);
freq = [900 7200 56000];

%% Write new obsfile for all active files
write_e3d_obs([work_dir '\FWR_mesh_J.dat'],freq,trx,core); 

%% RUN FORWARD MODEL

%% LOAD FIELDS AND BUILD WEIGHTS
[d_pred] = load_E3D_pred([work_dir '\e3d_data.txt']);
sigma = 1e-4;

% Cycle through the transmitter location and compute speudo sentivity
for ii = 1 : 29
    
    rx_x = trx(ii,1);
    rx_y = trx(ii,2);
    rx_z = trx(ii,3);
    
    drx_dx = rx_x - XYZ(:,1);
    drx_dy = rx_y - XYZ(:,2);
    drx_dz = rx_z - XYZ(:,3);
    
    R = sqrt( drx_dx.^2 + drx_dy.^2 + drx_dz.^2 );
    

    for jj = 1 : length(freq)

        if rr == 1 
            j = zeros(mcell,1);
        end
            
        j(:) = j(:) + ((4*pi*1e-7) * V ./ (4*pi*R.^3) .* (drx_dy .* sqrt( d_pred(jj,4).^2 + d_pred(jj,5).^2) - drx_dx .* Jy{jj})) .^2; 
        
    end
    
    
end

j = j.^(1/4) ./ dV.^0.5;
j = 10 * j / max(j) ;

