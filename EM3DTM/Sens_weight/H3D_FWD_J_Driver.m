%% Create cell center observation for forward model of J
% Written by: D.Fournier
% Last Update: June 26th, 2014
%% INPUT files
clear all
close all

work_dir    = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\Sens_weights\UTEMPhaseII';
meshfile    = 'mesh25_air_loop_6_and_1.msh';
modelfile   = 'Halfspace_1000ohms.con';
timefile    = 'Times.txt';
txfile      = 'Tx_loc.txt';

npadx = 10; npady = 10; npadz = 10; 
air =  1e-8;

%% Script starts here
% Load mesh
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

dx = xn(2:end) - xn(1:end-1);
dy = yn(2:end) - yn(1:end-1);
dz = zn(1:end-1) - zn(2:end);

% Create cell center coordinates
xc = xn(1) + cumsum( dx ) - dx/2; nx = length(xc);
yc = yn(1) + cumsum( dy ) - dy/2; ny = length(yc);
zc = zn(1) - cumsum( dz ) + dz/2; nz = length(zc);

% Create full array of mesh
[Zc,Xc,Yc] = ndgrid(zc,xc,yc);

% Load model
m = load([work_dir '\' modelfile]);

% Padding cells cannot be forward model (too close to edge)
% Only querry cells in core region, everu two cells
iqx = (npadx+1):2:(nx-npadx);
iqy = (npady+1):2:(ny-npady);
iqz = (npadz+1):2:(nz-npadz);

% Create querry location within the mesh
[Qz,Qx,Qy] = ndgrid( zc(iqz) , xc(iqx) , yc(iqy) );


%% Create projection matrix (cells to querry cells)
% Only run the first time. Then only load.
P = get_projector(Xc,Yc,Zc,iqx,iqy,iqz);
save([work_dir '\P'],'P');

load([work_dir '\P']);

%% Write to file location, time and transmitter parameters
% Read in time channels
times = load([work_dir '\' timefile]);
nt = length(times);

% read in transmitter location
tx = load([work_dir '\' txfile]);

% Write to file
write_h3d_obs(work_dir,'Obs_loc.dat',Qx(:),Qy(:),Qz(:),tx,times);

%% RUN FORWARD MODEL ON CLUSTER
return
%% Load predicted and map back to mesh

% Load predicted
data = load([work_dir '\dpred_0.txt']);

njx = length(iqx);
njy = length(iqy);
njz = length(iqz);

% Save model files for Jx and Jy, all time channels seperately.
for jj = 1 : nt
    
    Ex = data(jj:nt:end,5);
    Ey = data(jj:nt:end,6);
%     Ez = data(jj:nt:end,7);
    
    % Get current everywhere and save model
    Jx = abs( ( P * Ex ) .*m );
    Jy = abs( ( P * Ey ) .*m );
%     Jz = (P * Ez).*m;
    
    save([work_dir '\Jx_tc' num2str(times(jj)) '.mod'],'-ascii','Jx');
    save([work_dir '\Jy_tc' num2str(times(jj)) '.mod'],'-ascii','Jy');
%     save([work_dir '\Jz_tc' num2str(times(jj)) '.mod'],'-ascii','Jz');
    
end
