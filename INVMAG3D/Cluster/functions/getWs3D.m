function [ Ws , v ] = getWs3D(dx,dy,dz)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

% Square root dimensions for dV
dx = sqrt(dx);
dy = sqrt(dy);
dz = sqrt(dz);

% nullcell = reshape(nullcell,nz,nx,ny);

mcell = nx * ny * nz;

% Create hmid dimensions matrix
dX = kron( kron( ones(ny,1) , dx(:) ), ones(nz,1) );
dY = kron( kron( dy(:)  , ones(nx,1) ), ones(nz,1));
dZ = kron( kron( ones(ny,1) , ones(nx,1) ), dz(:) );

v = dX .* dY .* dZ;

Ws = speye(mcell);

end