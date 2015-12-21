function DIV_test(argin)
% Curl operator test
% 
% Computes the residual between the numerical and analytical solution for
% the curl of a function:
% f(x,y,z) = <


% clear all
% close all
% define mesh
count=1;



%% Discretization Iterations
            fprintf('**DIV TEST**\nUniform mesh | | | | |\n')
            fprintf('Cell size\t|residual|\n')
for ii = 1:6
   
n = 2 ^ ii;

%% Create mesh for center and faces
X0 = -2; Y0 = -2; Z0 = -2; %Origin
Xmax = 2; Ymax = 2; Zmax = 2; %Origin

%Create mesh size (cell center)
xf = linspace(X0,Xmax,n+1); xf=xf(:); 
yf = linspace(Y0,Ymax,n+1); yf=yf(:);
zf = linspace(Z0,Zmax,n+1); zf=zf(:);

% Create cell size
dx = abs(xf(2:end) - xf(1:end-1));
dy = abs(yf(2:end) - yf(1:end-1));
dz = abs(zf(2:end) - zf(1:end-1));

% Create hmid vectors
dxn = dx(1:end-1)/2 + dx(2:end)/2; %dxn=[dx(1);dxn;dx(end)];
dyn = dy(1:end-1)/2 + dy(2:end)/2; %dyn=[dy(1);dyn;dy(end)];
dzn = dz(1:end-1)/2 + dz(2:end)/2; %dzn=[dz(1);dzn;dz(end)];

%Create mesh location (center)
xc = X0 + dx(1)/2 + [0;cumsum(dxn)]; xc=xc(:); 
yc = Y0 + dy(1)/2 + [0;cumsum(dyn)]; yc=yc(:);
zc = Z0 + dz(1)/2 + [0;cumsum(dyn)]; zc=zc(:);



% dxn = dx(1:end-1)/2 + dx(2:end)/2;
% dyn = dy(1:end-1)/2 + dy(2:end)/2;
% dzn = dz(1:end-1)/2 + dz(2:end)/2;

nxn = length(dxn); nyn = length(dyn) ; nzn = length(dzn);
nx = length(dx); ny = length(dy) ; nz = length(dz);

% Create diagonal matrices for dimensions
dX = spdiags(1./dx,0,nx,nx);
dY = spdiags(1./dy,0,ny,ny);
dZ = spdiags(1./dz,0,nz,nz);

% Create cell-center dimension matrix
dXn = spdiags(1./[dx(1);dxn;dx(end)],0,nxn+2,nxn+2);
dYn = spdiags(1./[dy(1);dyn;dy(end)],0,nyn+2,nyn+2);
dZn = spdiags(1./[dz(1);dzn;dz(end)],0,nzn+2,nzn+2);

% Create cell-center dimension matrix with first and last dimension
% repeated for BCs.
dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);

% Number of faces in every directions
nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

% Number of edges
nhx = (nx) * (ny+1) * (nz+1);
nhy = (nx+1) * (ny) * (nz+1);
nhz = (nx+1) * (ny+1) * (nz);

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);


%% Create primary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV = [Dx Dy Dz];


 
%% Analytical solution of the form E = A , where A if divergence free

% DIV test

    % Edge variables
    [Xie, Yie, Zie] = ndgrid(xc, yf, zf);
    [Xje, Yje, Zje] = ndgrid(xf, yc, zf);
    [Xke, Yke, Zke] = ndgrid(xf, yf, zc);

    % Face Variables
    [Xif, Yif, Zif] = ndgrid(xf, yc, zc);
    [Xjf, Yjf, Zjf] = ndgrid(xc, yf, zc);
    [Xkf, Ykf, Zkf] = ndgrid(xc, yc, zf);

    % Create cell center variables
    [X,Y,Z] = ndgrid(xc, yc, zc);
    
    % Function for the field
    Eif = @(X,Y,Z) -Z .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    Ejf = @(X,Y,Z) -X .* Z .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    Ekf = @(X,Y,Z) -X .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    
    Aif = Eif(Xif, Yif, Zif);
    Ajf = Ejf(Xjf, Yjf, Zjf);
    Akf = Ekf(Xkf, Ykf, Zkf);
    
    Af = [Aif(:); Ajf(:); Akf(:)];


    % NUMERICAL: CURL A
    
    DIVAf = DIV * Af;
        

    % ANALYTICAL: V x A
    Be = 30*(X .* Z .* Y) .* exp(-5*(X.^2 + Y.^2 + Z.^2));
 
    Be = Be(:);

    % |residual|
    curl3n = norm(DIVAf - Be, 'inf'); 
    
    % Print to screen the result
    fprintf('%3.2e   %3.2e   \n',1/n,curl3n)
    
   

end
    fprintf('**end of test**\n\n')