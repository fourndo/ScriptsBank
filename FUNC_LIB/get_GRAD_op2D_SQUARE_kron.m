function [ Ws, Gx, Gz, V, Vx, Vz ]=get_GRAD_op2D_SQUARE_FB(dx,dz,nullcell,FLAG)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
nz = length(dz);

% Compute half-cell distance 
dxm = sqrt( dx(1:end-1)/2 + dx(2:end)/2 ) ; dxm = [dxm(:);dxm(end)];
dzm = sqrt( dz(1:end-1)/2 + dz(2:end)/2 ) ; dzm = [dzm(:);dzm(end)];

% Square root dimensions for dV
dx = sqrt(dx(:));
dz = sqrt(dz(:));

mcell = nx*nz;

nullcell = reshape(nullcell,nz,nx);

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

Gx = kron(d_dx,speye(nz));
Vx = kron(spdiags(1./dxm,0,nx,nx),spdiags(dz,0,nz,nz));

Gx = Gx(nullcell==1,nullcell==1);

% Find row to remove
indx = sum(abs(Gx),2) ~= 2;

Gx(indx,:) = 0;
Vx = Vx(nullcell==1,nullcell==1);

%% Create gradient operator and volume matrices for z
d_dx = ddx(nz); 
if strcmp(FLAG,'BACK')
    d_dx(1,2) = -1;
else
    d_dx(end,end-1) = 1;
end

Gz = kron(speye(nx),d_dx);
Vz = kron(spdiags(dx,0,nx,nx),spdiags(1./dzm,0,nz,nz));

Gz = Gz(nullcell==1,nullcell==1);

% Find row to remove
indx = sum(abs(Gz),2) ~= 2;

Gz(indx,:) = 0;
Vz = Vz(nullcell==1,nullcell==1);

%% Create hmid dimensions matrix
dX = kron( dx(:) , ones(nz,1) );
dZ = kron(  ones(nx,1) , dz(:) );

v = dX .* dZ;

Ws = speye(mcell);
V = spdiags( v , 0 , mcell,mcell);

Ws = Ws(nullcell==1,nullcell==1);
V = V(nullcell==1,nullcell==1);
    

end