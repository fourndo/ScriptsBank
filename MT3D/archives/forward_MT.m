% EOSC 555 - Operators
% 
% Authors
% Operator: Eldad
% Test : Dom
% 
% Date: March 9th, 2013
%
% INTRO:
% Script creating the GRAD, DIV and CURL operators for nodal
% discretization. The code has been tested on a simple vector field
% function:
% E = A + GRAD(phi)
% where A is a divergence-free vector function, 
% A = < -y^3 , x^3 >
% and phi is a curl-free scalar function, 
% phi = 3x + x^2*y - y^4 + z^2
%
% Both have a simple analytical solution and can be used seperately to test
% the GRAD, DIV and CURL operators.
% 
% The model space is 1 x 1 x 1 cube discretized by n x (n+1) x (n+2) cells.
% A solution is computed using the operators and compared to the
% analytical solution. Figure 1 shows that the error on the operators
% is approximatively of order h^2.
%
% The following properties have also been tested: 
% CURLxGRAD(phi) ==0 (Concervative fields are curl-free)
% DIV ( CURL x (E)) == 0 (Divergence of the B-field is 0)

% function [GRAD,DIV,CURL]=ops_primary(dx,dy,dz)


clear all
close all
% define mesh
count=1;
% for n = [8,16,32,64]%,64,128,256,512]
n=3;   
% h = 1/n;
% nodal grid
% tF = linspace(0,1,n+1); tF = tF(:);
% tC = linspace(0,1,n); tC = tC(:);

%% Create mesh for center and faces
X0 = 0; Y0 = 0; Z0 = 0; %Origin

%Create mesh size (faces)
xf = linspace(0,1,n+1); xf=xf(:); 
yf = linspace(0,1,n+1); yf=yf(:);
zf = linspace(0,1,n+1); zf=zf(:);

dx = abs(xf(2:end) - xf(1:end-1));
dy = abs(yf(2:end) - yf(1:end-1));
dz = abs(zf(2:end) - zf(1:end-1));
% dz = [0.02;0.03;0.05;0.08;0.12;0.14;0.15];
%Create mesh (center to center)
% xc = dx(1:end-1)/2 + dx(2:end)/2; xc=xc(:); 
% yc = yf(1:end-1) + dy/2; yc=yc(:);
% zc = zf(1:end-1) + dz/2; zc=zc(:);

dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% dxn = dx(1:end-1)/2 + dx(2:end)/2;
% dyn = dy(1:end-1)/2 + dy(2:end)/2;
% dzn = dz(1:end-1)/2 + dz(2:end)/2;

nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);
nx = length(dx); ny = length(dy) ; nz = length(dz);

% Create diagonal matrices for hmid dimensions
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

% Create cell-center dimension matrix
dX = spdiags(1./dx,0,nx,nx); dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
dY = spdiags(1./dy,0,ny,ny); dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
dZ = spdiags(1./dz,0,nz,nz); dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);

nex = (nx+1) * (ny) * (nz);
ney = (nx) * (ny+1) * (nz);
nez = (nx) * (ny) * (nz+1);

nhx = (nx) * (ny+1) * (nz+1);
nhy = (nx+1) * (ny) * (nz+1);
nhz = (nx+1) * (ny+1) * (nz);

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);

d_dx = ddx(n);

%% Derivative matrices for curl (faces to edge)
% The first and last vectors in the y and z-direction and set to 0
% The last vector in the x-direction is set to 0;

% d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1);  d_dz([1,end])=d_dz([1,end])*2;
% Dzy = kron( kron( d_dz , speye(ny+1) ), speye(nx) );
% 
% d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
% Dyz = kron( kron( speye(nz+1) , d_dy ), speye(nx) );
% 
% d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=d_dz([1,end])*2;
% Dzx = kron( kron( d_dz, speye(ny) ) , speye(nx+1) );
% 
% d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
% Dxz = kron( kron( speye(nz+1) , speye(ny) ), d_dx );
% 
% d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
% Dyx = kron( kron( speye(nz) , d_dy ), speye(nx+1) );
% 
% d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
% Dxy = kron( kron( speye(nz) , speye(ny+1) ), d_dx );
% 
% Ox = sparse ( nhx , nex );
% Oy = sparse ( nhy , ney );
% Oz = sparse ( nhz , nez );
% 
% CURLf = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz];

% clear Dzy Dyz Dxz Dzx Dyx Dxy
% Curl operator from edges to faces

% Dzy = kron( kron( dZ * ddx(nz-1) , speye(ny) ), speye(nx+1) );
% Dyz = kron( kron( speye(nz) , dY * ddx(ny-1) ), speye(nx+1) );
% 
% Dzx = kron( kron( dZ * ddx(nz-1), speye(ny+1) ) , speye(nx) );
% Dxz = kron( kron( speye(nz) , speye(ny+1) ), dX * ddx(nx-1) );
% 
% Dyx = kron( kron( speye(nz+1) , dY * ddx(ny-1) ), speye(nx) );
% Dxy = kron( kron( speye(nz+1) , speye(ny) ), dX * ddx(nx-1) );
% 
% Ox = sparse ( nex , nhx );
% Oy = sparse ( ney , nhy );
% Oz = sparse ( nez , nhz );
% 
% CURLe = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz];
% 
% L = CURLe * CURLf;
% L = CURLf' * CURLf;

%% Attempt at Laplacian operator
% Second partial Derivative in 3 directions

uo = 4*pi*10^-7;

mb = 0.01; % model background [S/m]
mc = 10; % model conductor [S/m]

w = 5^2; % Frequency [hz]


% %Partial derivatives for x-component
d_dx=  dXxl * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; 
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

%Boundary COnditions: 0 derivative on y-plane
DDyx(1:(nx+1):end,:)=0;
DDyx((nx+1):(nx+1):end,:)=0;

d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; 

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
kappa=sqrt(1i*w*uo*mb); 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end) + 1i*kappa/dZm(end,end);

DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );



% %Partial derivatives for y-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

%Boundary COnditions: 0 derivative on x-plane
DDxy(1:(ny+1)*nx:end,:)=0;
DDxy(2:(ny+1)*nx:end,:)=0;
DDxy(3:(ny+1)*nx:end,:)=0;

DDxy((ny+1)*nx-2:(ny+1)*nx:end,:)=0;
DDxy((ny+1)*nx-1:(ny+1)*nx:end,:)=0;
DDxy((ny+1)*nx:(ny+1)*nx:end,:)=0;


d_dy= dYxl * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; 

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
kappa=sqrt(1i*w*uo*mb); 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end) + 1i*kappa/dZm(end,end);

DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );

% %Partial derivatives for z-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz;
DDzz = kron( kron( dd_dz , speye(ny) ), speye(nx) );

Oyx = sparse ( nex , ney );
Ozx = sparse ( nex , nez );

Oxy = sparse ( ney , nex );
Ozy = sparse ( ney , nez );

Oxz = sparse ( nez , nex );
Oyz = sparse ( nez , ney );

L = [DDxx+DDyx+DDzx Oyx Ozx;
     Oxy DDxy+DDyy+DDzy Ozy;
     Oxz Oyz DDxz+DDyz+DDzz];


%% Create center gradient (center to face)

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

GRAD = [Dx; Dy; Dz];

%% Create divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV = [Dx Dy Dz];


%% Averaging operator to faces

av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avx = av(nx); avx([1,end]) = 1;
avy = av(ny); avy([1,end]) = 1;
avz = av(nz); avz([1,end]) = 1;

Avx = kron ( kron ( speye (nz) , speye (ny)), avx);
Avy = kron ( kron ( speye (nz) , avy), speye (nx));
Avz = kron ( kron ( avz , speye (ny)), speye (nx));

AVF = [Avx;Avy;Avz];
 
%% Averaging for edges
% % Averaging two direction perp to the averaging direction
% Avx = kron ( kron ( avz , avy) , speye (nx));
% Avy = kron ( kron ( avz , speye (ny)), avx);
% Avz = kron ( kron ( speye (nz) , avy), avx);
% 
% AVE = [Avx;Avy;Avz];

%% Construct Boundary Condition Matrix for A_vert Dirichlet condition
% Note since we have Az+L = g(x, y, z = +L) = C0 and
% Az-L = g(x, y, z =- L) = C1 we need to institute a boundary condition matrix
% and hit it with the Laplacian.



e = @(n) ones(n, 1);
bcs = [-2000, 2*0]; %Decays to zero at depth, some constant in x and y at the
% surface. The BC's for the z (vertical case) direction are just zeros.
bcx = zeros((nx + 1) * ny * nz, 1); 
bcx(1:(nx + 1) * ny) = bcs(1)*dZm(1)*dZxl(1);
% bcx((end -((nx + 1) * ny) +1 ):end) = bcs(2);

bcy = zeros(nx  * (ny+ 1) * nz, 1);
bcy(1:(nx + 1) * ny) = bcs(1)*dZm(1)*dZxl(1);
% bcy((end-((nx + 1) * ny) +1):end) = bcs(2);

bcz = zeros(nx  * ny * (nz+ 1), 1);% bcz([1,end]) = bcs;

BCx = sparse(bcx);

BCy = sparse(bcy);

BCz = sparse(bcz);

% BCx = kron(e(nz), kron( e(ny), bcx));
% BCy = kron(e(nz), kron( bcy, e(nx)));
% BCz = kron(bcz, kron( e(ny), e(nx)));
% 
% B = [BCx ; BCy; BCz];

% BCx = kron(bcx, kron( e(ny), e(nx)));
% BCy = kron(bcy, kron( e(ny), e(nx)));
% BCz = kron(bcz, kron( e(ny), e(nx)));

% BCx = kron( kron (e(nz),  e(ny)), bcx);
% BCy = kron( kron (e(nz),  bcy), e(nx));
% BCz = kron( kron (bcz,  e(ny)),e(nx));

% BCx = kron( kron (e(nz),  e(ny)), bcx);
% BCy = kron( kron (e(nz),  bcy), e(nx));
% BCz = kron( kron (bcz,  e(ny)),e(nx));
B = [BCx ; BCy; BCz];

A = [L + 1i*w*uo, 1i*uo*GRAD
     DIV      , DIV*GRAD];

% A([1,1]) = A([1,1]) + 1;
 
q = [ B; zeros(size(DIV,1), 1)];

u = A\q;

A = u(1: size(GRAD, 1));
Phi = u( size(GRAD, 1) + 1: end);


E = A + GRAD*Phi;

vecplot(E, dx, dy, dz)
