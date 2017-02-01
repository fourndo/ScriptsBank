function [ Ws, Gx , Gy, Gz, V, Vx, Vy, Vz ]=get_GRAD_op3D_SQUARE_Kron(dx,dy,dz,nullcell,FLAG)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

% Compute half-cell distance 
dxm =  dx(1:end-1)/2 + dx(2:end)/2 ; dxm = [dxm(:);dxm(end)];
dym =  dy(1:end-1)/2 + dy(2:end)/2 ; dym = [dym(:);dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2 ; dzm = [dzm(:);dzm(end)];

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
dz = sqrt(dz(:));

mcell = nx*ny*nz;

% nullcell = reshape(nullcell,nz,nx,ny);

% Generate Gradient Operators and volume matrices
% mcellx = (nx-1)*ny*nz;
% mcelly = nx*(ny-1)*nz;
% mcellz = nx*ny*(nz-1);

% Gx = spalloc(mcell,mcell,2*mcell);
% Gy = spalloc(mcell,mcell,2*mcell);
% Gz = spalloc(mcell,mcell,2*mcell);
% 
% Vx = speye(mcell);
% Vy = speye(mcell);
% Vz = speye(mcell);



if strcmp(FLAG,'BACK')
    ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[-1,0],n,n);
else
    ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n,n);
end

%% Create gradient operator and volume matrices for x
d_dx = ddx(nx); 

if strcmp(FLAG,'BACK')
    d_dx(1,2) = -1;
else
    d_dx(end,end-1) = 1;
end

Gx = kron(kron(speye(ny),d_dx),speye(nz));
Vx = kron(kron(spdiags(dy,0,ny,ny),spdiags(1./dxm,0,nx,nx)),spdiags(dz,0,nz,nz));

Gx = Gx(nullcell==1,nullcell==1);

% Find row to remove
indx = sum(abs(Gx),2) ~= 2;

Gx(indx,:) = 0;
Vx = Vx(nullcell==1,nullcell==1);

%% Create gradient operator and volume matrices for y
d_dx = ddx(ny); 

if strcmp(FLAG,'BACK')
    d_dx(1,2) = -1;
else
    d_dx(end,end-1) = 1;
end

Gy = kron(kron(d_dx,speye(nx)),speye(nz));
Vy = kron(kron(spdiags(1./dym,0,ny,ny),spdiags(dx,0,nx,nx)),spdiags(dz,0,nz,nz));

Gy = Gy(nullcell==1,nullcell==1);

% Find row to remove
indx = sum(abs(Gy),2) ~= 2;

Gy(indx,:) = 0 ;
Vy = Vy(nullcell==1,nullcell==1);

%% Create gradient operator and volume matrices for z
d_dx = ddx(nz); 

if strcmp(FLAG,'BACK')
    d_dx(1,2) = -1;
else
    d_dx(end,end-1) = 1;
end

Gz = kron(kron(speye(ny),speye(nx)),d_dx);
Vz = kron(kron(spdiags(dy,0,ny,ny),spdiags(dx,0,nx,nx)),spdiags(1./dzm,0,nz,nz));

Gz = Gz(nullcell==1,nullcell==1);

% Find row to remove
indx = sum(abs(Gz),2) ~= 2;

Gz(indx,:) = 0;
Vz = Vz(nullcell==1,nullcell==1);

%% Create hmid dimensions matrix
dX = kron( kron( ones(ny,1) , dx(:) ), ones(nz,1) );
dY = kron( kron( dy(:)  , ones(nx,1) ), ones(nz,1));
dZ = kron( kron( ones(ny,1) , ones(nx,1) ), dz(:) );

v = dX .* dY .* dZ;

Ws = speye(mcell);
V = spdiags( v , 0 , mcell,mcell);

Ws = Ws(nullcell==1,nullcell==1);
V = V(nullcell==1,nullcell==1);
% Gx = X * Gx * X';
% Gy = X * Gy * X';
% Gz = X * Gz * X';
% 
% Vx = X * Vx * X';
% Vy = X * Vy * X';
% Vz = X * Vz * X';


end