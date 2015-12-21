function [ Wx , Wy, Wz, Vx, Vy, Vz ]=get_GRAD_op3D_v3(dx,dy,dz,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

nullcell = reshape(nullcell,nz,ny,nx);

% Create cell center distance
dxm =  dx(1:end-1)/2 + dx(2:end)/2 ;
dXm = spdiags(1./dxm(:),0,nx-1,nx-1);

dym =  dy(1:end-1)/2 + dy(2:end)/2 ;
dYm = spdiags(1./dym(:),0,ny-1,ny-1);

dzm = dz(1:end-1)/2 + dz(2:end)/2 ;%dzm = [dzm(1) dzm];
dZm = spdiags(1./dzm(:),0,nz-1,nz-1);

% Derivative operators in 1D
d_dx = spdiags (ones(nx+1,1)*[-1,1],[0,1],nx-1,nx);
d_dy = spdiags (ones(ny+1,1)*[-1,1],[0,1],ny-1,ny);
d_dz = spdiags (ones(nz+1,1)*[-1,1],[0,1],nz-1,nz);

% Create gradient operators in 3D
Wx = kron( kron( speye(ny) , d_dx ), speye(nz) );
Wy = kron( kron( d_dy , speye(nx) ), speye(nz) );
Wz = kron( kron( speye(ny) , speye(nx) ), d_dz );

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
dz = sqrt(dz(:));

% Kill air cells from the gradient Wx
index = nullcell(:,1:end-1,:) - nullcell(:,2:end,:);
mask = (index~=-1 & index~=1 & nullcell(:,1:end-1,:)~=0); mask = reshape(mask,nz*ny*(nx-1),1);
Wx = spdiags(mask,0,nz*ny*(nx-1),nz*ny*(nx-1)) * Wx;

% Kill air cells from the gradient Wy
index = nullcell(:,:,1:end-1) - nullcell(:,:,2:end);
mask = (index~=-1 & index~=1 & nullcell(:,:,1:end-1)~=0); mask = reshape(mask,nz*(ny-1)*nx,1);
Wy = spdiags(mask,0,nz*(ny-1)*nx,nz*(ny-1)*nx) * Wy;

% Kill air cells from the gradient Wz
index = nullcell(1:end-1,:,:) - nullcell(2:end,:,:);
mask = (index~=-1 & index~=1 & nullcell(1:end-1,:,:)~=0); mask = reshape(mask,(nz-1)*ny*nx,1);
Wz = spdiags(mask,0,(nz-1)*ny*nx,(nz-1)*ny*nx) * Wz;

% Create dimensions matrix Wx
dX = kron( kron( speye(ny) , spdiags(dx(1:end-1),0,nx-1,nx-1) * dXm ), speye(nz) );
dY = kron( kron( spdiags(dy,0,ny,ny)  , speye(nx-1) ), speye(nz));
dZ = kron( kron( speye(ny) , speye(nx-1) ), spdiags(dz,0,nz,nz) );
Vx = dX * dY * dZ ;

% Create dimensions matrix Wy
dX = kron( kron( speye(ny-1) , spdiags(dx,0,nx,nx) ), speye(nz) );
dY = kron( kron( spdiags(dy(1:end-1),0,ny-1,ny-1) * dYm , speye(nx) ), speye(nz));
dZ = kron( kron( speye(ny-1) , speye(nx) ), spdiags(dz,0,nz,nz) );
Vy = dX * dY * dZ ;

% Create dimensions matrix Wz
dX = kron( kron( speye(ny) , spdiags(dx,0,nx,nx) ), speye(nz-1) );
dY = kron( kron( spdiags(dy,0,ny,ny)  , speye(nx) ), speye(nz-1));
dZ = kron( kron( speye(ny) , speye(nx) ), spdiags(dz(1:end-1),0,nz-1,nz-1) * dZm );
Vz = dX * dY * dZ ;

end