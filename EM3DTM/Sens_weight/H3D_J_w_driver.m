%% Load J for all time channel and build sensitivity weighting
% Need to run H3D_FWD_J_Driver before for each transmitter loop
%
% Written by: D.Fournier
% Last Update: July 12th, 2014
%% INPUT files
clear all
close all

work_dir    = 'C:\LC\Private\dominiquef\Projects\4264_Vedanta_TEM\Forward\Workspace\Tank_Hill\J_analysis\Loop2';

% May repeat the process for multiple transmitter. Start at 1.
% The program will load previous approximated J and add the next tx_ID
% i.e. J_tx1 + J_tx2 + ...
tx_ID = 1;

meshfile    = '..\..\TankHill_mesh_50m_air.msh';
timefile    = '..\..\Times.dat';


rx = load([work_dir '\rx_loc.dat']);


%% Script starts here
% Load mesh
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

dx = xn(2:end) - xn(1:end-1);
dy = yn(2:end) - yn(1:end-1);
dz = zn(1:end-1) - zn(2:end);

x = xn(1) + cumsum(dx) - dx/2;
y = yn(1) + cumsum(dy) - dy/2;
z = zn(1) - cumsum(dz) + dz/2;

nx = length(dx);
ny = length(dy);
nz = length(dz);

mcell = nx*ny*nz;

[dZ,dX,dY] = ndgrid(dz,dx,dy);
dZ = dZ(:);
dX = dX(:);
dY = dY(:);

[Z,X,Y] = ndgrid(z,x,y);
Z = Z(:);
X = X(:);
Y = Y(:);

% Compute volume elements
dV = dZ .* dX .* dY;

%% Load time sensitivity from previous tx (can only be used AFTER tx=1)
if tx_ID == 1;
    
    J = zeros(mcell,1);
    
else
    
    load([work_dir '\..\J']);
    
end

%% Write to file location, time and transmitter parameters
% Read in time channels
times = load([work_dir '\' timefile]);
nt = length(times);

%% Load Jx, Jy and compute weights models and pro
  
cd(work_dir)

file_list = ls;
countx = 1;
county = 1;
for ii = 1:size(file_list,1)-2;

    look_at = strtrim(file_list(ii+2,:));

    if strcmp(look_at(1:2),'Jx')==1 
        
        % Get the time       
        Jx{countx} = load([work_dir '\' look_at]);
        countx = countx +1;
        
    elseif strcmp(look_at(1:2),'Jy')==1
        
        % Get the time       
        Jy{county} = load([work_dir '\' look_at]);
        county = county +1;
        
    end
    
end
% load([work_dir '\Jx']);
% load([work_dir '\Jy']);


% Compute observation location distance and offset
for rr = 1 : size(rx,1)
    
    rx_x = rx(rr,1);
    rx_y = rx(rr,2);
    rx_z = rx(rr,3);
    
    drx_dx = rx_x - X;
    drx_dy = rx_y - Y;
    drx_dz = rx_z - Z;
    
    R = sqrt( drx_dx.^2 + drx_dy.^2 + drx_dz.^2 );
    

    for tt = 1 : nt

        if rr == 1 
            j = zeros(mcell,1);
        end
            
        j(:) = j(:) + ((4*pi*1e-7) * dV ./ (4*pi*R.^3) .* (drx_dy .* Jx{tt} - drx_dx .* Jy{tt})) .^2; 
        
    end
    
    
end

j = j.^(1/4) ./ dV.^0.5;
j = 10 * j / max(j) ;

J = J + j;

if tx_ID == 2
    
    J = 10 * J / max(J);
    J(J<1) = 1;

end
save([work_dir '\..\J'],'J');