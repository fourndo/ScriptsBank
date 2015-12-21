function [ Wx , Wy, Wz, Vx, Vy, Vz ]=get_GRAD_op3D_center(dx,dy,dz)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1
% clear all
% 
% dx = ones(1,4)*2;
% dy = ones(1,4)*3;
% dz = ones(1,4)*4;

nx = length(dx);
ny = length(dy);
nz = length(dz);

dx = dx(:);
dy = dy(:);
dz = dz(:);

%% FORWARD - Create cell center distance
dxm =  dx(1:end-1)/2 + dx(2:end)/2 ; dxmf = [dxm(1);dxm];
dXmf = spdiags(1./dxmf,0,nx,nx);

dym =  dy(1:end-1)/2 + dy(2:end)/2 ; dymf = [dym(1);dym];
dYmf = spdiags(1./dymf,0,ny,ny);

dzm = dz(1:end-1)/2 + dz(2:end)/2 ;  dzmf = [dzm(1);dzm];
dZmf = spdiags(1./dzmf,0,nz,nz);

% Gradient operators in 1D
d_dx = spdiags (ones(nx+1,1)*[-1,1],[-1,0],nx,nx); d_dx(1,1)=0;
d_dy = spdiags (ones(ny+1,1)*[-1,1],[-1,0],ny,ny); d_dy(1,1)=0;
d_dz = spdiags (ones(nz+1,1)*[-1,1],[-1,0],nz,nz); 


%% BACKWARD - Create cell center distance
dxmb = [dxm;dxm(end)];
dXmb = spdiags(1./dxmb,0,nx,nx);

dymb = [dym;dym(end)];
dYmb = spdiags(1./dymb,0,ny,ny);

dzmb = [dzm;dzm(end)];
dZmb = spdiags(1./dzmb,0,nz,nz);

% Create gradient operators in 3D 
Wxf =  kron( kron( speye(ny) , spdiags(dxmb,0,nx,nx) * d_dx ), speye(nz) );
Wyf = kron( kron( spdiags(dymb,0,ny,ny) * d_dy , speye(nx) ), speye(nz) );
Wzf = kron( kron( speye(ny) , speye(nx) ), spdiags(dzmb,0,nz,nz) * d_dz );

% Gradient operators in 1D
d_dx = spdiags (ones(nx+1,1)*[-1,1],[0,1],nx,nx); d_dx(end,end)=0;
d_dy = spdiags (ones(ny+1,1)*[-1,1],[0,1],ny,ny); d_dy(end,end)=0;
d_dz = spdiags (ones(nz+1,1)*[-1,1],[0,1],nz,nz); d_dz(end,end)=0;

% Create gradient operators in 3D 
Wxb = kron( kron( speye(ny) , spdiags(dxmf,0,nx,nx) * d_dx ), speye(nz) );
Wyb = kron( kron( spdiags(dymf,0,ny,ny) * d_dy , speye(nx) ), speye(nz) );
Wzb = kron( kron( speye(ny) , speye(nx) ), spdiags(dzmf,0,nz,nz) * d_dz );

%% Average forward and backward
Wx = (Wxb + Wxf) / 2;
Wy = (Wyb + Wyf) / 2;
Wz = (Wzb + Wzf) / 2;

%% Square root dimensions for dV
dx = sqrt(dx);
dy = sqrt(dy);
dz = sqrt(dz);

% Create dimensions matrix Wx
dX = kron( kron( speye(ny) ,  spdiags(sqrt(dxmb+dxmf),0,nx,nx) * dXmb * dXmf  ), speye(nz) );
dY = kron( kron( spdiags(dy,0,ny,ny)  , speye(nx) ), speye(nz));
dZ = kron( kron( speye(ny) , speye(nx) ), spdiags(dz,0,nz,nz) );
Vx = dX .* dY .* dZ;

% Create dimensions matrix Wy
dX = kron( kron( speye(ny) , spdiags(dx,0,nx,nx) ), speye(nz) );
dY = kron( kron( spdiags(sqrt(dymb+dymf),0,ny,ny) * dYmb * dYmf   , speye(nx) ), speye(nz));
dZ = kron( kron( speye(ny) , speye(nx) ), spdiags(dz,0,nz,nz) );
Vy = dX .* dY .* dZ;

% Create dimensions matrix Wz
dX = kron( kron( speye(ny) , spdiags(dx,0,nx,nx) ), speye(nz) );
dY = kron( kron( spdiags(dy,0,ny,ny)  , speye(nx) ), speye(nz));
dZ = kron( kron( speye(ny) , speye(nx) ), spdiags(sqrt(dzmb+dzmf),0,nz,nz) * dZmb * dZmf );
Vz = dX .* dY .* dZ;
