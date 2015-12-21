function CURL_test(argin)
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
            fprintf('**CURL TEST**\nUniform mesh | | | | |\n')
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


%% Derivative matrices for curl (faces to edge)
% The first and last vectors in the y and z-direction and set to 0
% The last vector in the x-direction is set to 0;

d_dz = dZn * ddx(nz); d_dz=d_dz(:,2:end-1);  d_dz([1,end])=d_dz([1,end])*2;
Dzy = kron( kron( d_dz , speye(ny+1) ), speye(nx) );

d_dy = dYn * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dyz = kron( kron( speye(nz+1) , d_dy ), speye(nx) );

d_dz = dZn * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=d_dz([1,end])*2;
Dzx = kron( kron( d_dz, speye(ny) ) , speye(nx+1) );

d_dx = dXn * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dxz = kron( kron( speye(nz+1) , speye(ny) ), d_dx );

d_dy = dYn * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dyx = kron( kron( speye(nz) , d_dy ), speye(nx+1) );

d_dx = dXn * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dxy = kron( kron( speye(nz) , speye(ny+1) ), d_dx );

Ox = sparse ( nhx , nfx );
Oy = sparse ( nhy , nfy );
Oz = sparse ( nhz , nfz );

CURL = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz];


 
%% Analytical solution of the form E = A , where A if divergence free


%% CURL test

    % Edge variables
    [Xie, Yie, Zie] = ndgrid(xc, yf, zf);
    [Xje, Yje, Zje] = ndgrid(xf, yc, zf);
    [Xke, Yke, Zke] = ndgrid(xf, yf, zc);

    % Face Variables
    [Xif, Yif, Zif] = ndgrid(xf, yc, zc);
    [Xjf, Yjf, Zjf] = ndgrid(xc, yf, zc);
    [Xkf, Ykf, Zkf] = ndgrid(xc, yc, zf);

    % Function for the field
    Eif = @(X,Y,Z) -Z .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    Ejf = @(X,Y,Z) -X .* Z .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    Ekf = @(X,Y,Z) -X .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    
    Aif = Eif(Xif, Yif, Zif);
    Ajf = Ejf(Xjf, Yjf, Zjf);
    Akf = Ekf(Xkf, Ykf, Zkf);
    
    Af = [Aif(:); Ajf(:); Akf(:)];
    
    

%     % Define curl function A
%     Ax = @(x,y,z)(cos(y));
%     Ay = @(x,y,z)(cos(z));
%     Az = @(x,y,z)(cos(x));
% 
    % Need to compute BC for Derichelet
    e = @(n) ones(n, 1);
    bcx = sparse(nx+1, 1); bcx([1,end]) = [-2 2];
    bcy = sparse(ny+1, 1); bcy([1,end]) = [-2 2];
    bcz = sparse(nz+1, 1); bcz([1,end]) = [-2 2];

    % Get the x components of BC on edges
    [Y,X,Z] = meshgrid (yf,xc,zf);
    A_y = Ejf(X,Y,Z);

    BCz = kron(bcz, kron( e(ny+1), e(nx)));

    Hzy = BCz.*A_y(:) * dZ(1);


    % Compute only face values for Dyz
    [Y,X,Z] = meshgrid (yf,xc,zf);
    A_z = Ekf(X,Y,Z);

    BCy = kron(e(nz+1), kron( bcy, e(nx)));
    Hyz = BCy.*A_z(:) * dY(1);

    % Compute the total Hx component, minus sign on the Hyz like the curl.
    Hx = Hzy - Hyz;

    % Get the y components of BC
    [Y,X,Z] = meshgrid (yc,xf,zf);
    A_z = Ekf(X,Y,Z);

    BCx = kron(e(nz+1), kron( e(ny), bcx));
    Hxz = BCx.*A_z(:) * dX(1);



    [Y,X,Z] = meshgrid (yc,xf,zf);
    A_x = Eif(X,Y,Z);

    BCz = kron(bcz, kron( e(ny), e(nx+1)));
    Hzx = BCz.*A_x(:) * dZ(1);

    % Compute the total Hy component, minus sign on the Hzx like the curl.
    Hy = Hxz - Hzx;

    % Get the z components of BC
    % [Y,X,Z] = meshgrid ([yf(1:end-2);yf(end)],xf,zc);
    [Y,X,Z] = meshgrid (yf,xf,zc);
    A_x = Eif(X,Y,Z);

    BCy = kron(e(nz), kron( bcy, e(nx+1)));
    Hyx = BCy.*A_x(:)*dY(1);

    % [Y,X,Z] = meshgrid (yf,[xf(1:end-2);xf(end)],zc);
    [Y,X,Z] = meshgrid (yf,xf,zc);
    A_y = Ejf(X,Y,Z);

    BCx = kron(e(nz), kron( e(ny+1), bcx));
    Hxy = BCx.*A_y(:) * dX(1);

    % Compute the total Hz component, minus sign on the Hxy like the curl.
    Hz = -Hxy + Hyx;


    BC =[Hx;Hy;Hz];

    % NUMERICAL: CURL A
    
    CURLAf = CURL * Af;
    
    Ae = CURLAf + BC;
    

    % ANALYTICAL: V x A
    Bie = 10*(Xie .* Zie.^2 - Xie .* Yie.^2) .* exp(-5*(Xie.^2 + Yie.^2 + Zie.^2));
    
    Bje = 10*(Yje .* Xje.^2 - Yje .* Zje.^2) .* exp(-5*(Xje.^2 + Yje.^2 + Zje.^2));
    
    Bke = 10*(Zke .* Yke.^2 - Zke .* Xke.^2) .* exp(-5*(Xke.^2 + Yke.^2 + Zke.^2));

    Be = [Bie(:); Bje(:); Bke(:)];

    % |residual|
    curl3n = norm(Ae - Be, 'inf'); 
    
    % Print to screen the result
    fprintf('%3.2e   %3.2e   \n',1/n,curl3n)
    
   

end

fprintf('**end of test**\n\n')