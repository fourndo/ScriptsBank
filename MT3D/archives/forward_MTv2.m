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
n=4;   
% h = 1/n;
% nodal grid
% tF = linspace(0,1,n+1); tF = tF(:);
% tC = linspace(0,1,n); tC = tC(:);

%% Create mesh for center and faces
X0 = 0; Y0 = 0; Z0 = 0; %Origin

%Create mesh position (faces)
xf = linspace(0,300,n+1); xf=xf(:); 
yf = linspace(0,300,n+1); yf=yf(:);
zf = linspace(0,300,n+1); zf=zf(:);


% Create mesh size
dx = abs(xf(2:end) - xf(1:end-1));
dy = abs(yf(2:end) - yf(1:end-1));
dz = abs(zf(2:end) - zf(1:end-1));
% dz = 25*1.1.^[1:n+10]';
% Create center-center mesh size (hmid)
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% Compute number of faces and cells in 3D
nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);

nx = length(dx); ny = length(dy) ; nz = length(dz);

nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

nface = nfx + nfy + nfz;

mcell = nx * ny *nz;

% Create diagonal matrices for hmid dimensions
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

% Create cell-center dimension matrix
dX = spdiags(1./dx,0,nx,nx); dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
dY = spdiags(1./dy,0,ny,ny); dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
dZ = spdiags(1./dz,0,nz,nz); dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);


% nhx = (nx) * (ny+1) * (nz+1);
% nhy = (nx+1) * (ny) * (nz+1);
% nhz = (nx+1) * (ny+1) * (nz);

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);

%% Model parameters
uo = 4*pi*10^-7;

s = ones(mcell,1) * 1e-2; % Background conductivity [S/m]

w = 1e+3; % Frequency [hz]

skindepth = sqrt( 2/ (s(1) * uo * w));

%% Laplacian operator for primary field
% Second partial Derivative in 3 directions
% %Partial derivatives for x-component
d_dx=  dXxl * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end],:)=0;
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz;

kappa=sqrt(1i*w*uo*s); 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end)/2 + 1i*kappa(1)/dZm(end,end);

DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

% %%Partial derivatives for y-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

d_dy= dYxl * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end],:)=0;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; 

kappa=sqrt(1i*w*uo*s); 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end)/2 + 1i*kappa(1)/dZm(end,end);

DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


% %%Partial derivatives for z-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; dd_dz([1,end],:)=0;
DDzz = kron( kron( dd_dz , speye(ny) ), speye(nx) );

Oyx = sparse ( nfx , nfy );
Ozx = sparse ( nfx , nfz );

Oxy = sparse ( nfy , nfx );
Ozy = sparse ( nfy , nfz );

Oxz = sparse ( nfz , nfx );
Oyz = sparse ( nfz , nfy );

L_p = [DDxx+DDyx+DDzx Oyx Ozx;
     Oxy DDxy+DDyy+DDzy Ozy;%
     Oxz Oyz DDxz+DDyz+DDzz];

%% Laplacian operator for secondary field
% Second partial Derivative in 3 directions
% %Partial derivatives for x-component
d_dx=  dXxl * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end])=0;
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])*3/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; dd_dz([1,end])=dd_dz([1,end])*3/2;
DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

% %%Partial derivatives for y-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])*3/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

d_dy= dYxl * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end])=0;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

d_dz= dZm * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; dd_dz([1,end])=dd_dz([1,end])*3/2; 
DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


% %%Partial derivatives for z-component
d_dx= dXm * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])*3/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYm * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])*3/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= d_dz'*d_dz; dd_dz([1,end])=0;
DDzz = kron( kron( dd_dz , speye(ny) ), speye(nx) );

Oyx = sparse ( nfx , nfy );
Ozx = sparse ( nfx , nfz );

Oxy = sparse ( nfy , nfx );
Ozy = sparse ( nfy , nfz );

Oxz = sparse ( nfz , nfx );
Oyz = sparse ( nfz , nfy );

L_s = [DDxx+DDyx+DDzx Oyx Ozx;
     Oxy DDxy+DDyy+DDzy Ozy;%
     Oxz Oyz DDxz+DDyz+DDzz];
 
%% Create primary gradient (center to face)

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

GRAD_p = [Dx; Dy; Dz];

%% Create primary gradient (center to face)

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=d_dz([1,end])*2;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

GRAD_s = [Dx; Dy; Dz];

%% Create primary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV_p = [Dx Dy Dz];

%% Create primary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1); d_dx([1,end])=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1); d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1); d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV_s = [Dx Dy Dz];

%% Averaging operator to faces

av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avx = av(nx); avx([1,end]) = 1;
avy = av(ny); avy([1,end]) = 1;
avz = av(nz); avz([1,end]) = 1;

Avx = kron ( kron ( speye (nz) , speye (ny)), avx);
Avy = kron ( kron ( speye (nz) , avy), speye (nx));
Avz = kron ( kron ( avz , speye (ny)), speye (nx));

AVF = [Avx;Avy;Avz];

S = @(m) spdiags( (AVF * (m.^-1)).^-1, 0, nface, nface) ;
 

%% Construct Boundary Condition Matrix for A_vert Dirichlet condition
% Note since we have Az+L = g(x, y, z = +L) = C0 and
% Az-L = g(x, y, z =- L) = C1 we need to institute a boundary condition matrix
% and hit it with the Laplacian.

% e = @(n) ones(n, 1);

%Decays to zero at depth, some constant in x and y at the
% surface. The BC's for the z (vertical case) direction are just zeros.
bcs = [-2*1000, 2*0]; 
bcx = zeros((nx + 1) * ny * nz, 1); 
bcx(1:(nx + 1) * ny) = bcs(1)*dZm(1)^2;
% bcx((end -((nx + 1) * ny) +1 ):end) = bcs(2);

bcy = zeros(nx  * (ny+ 1) * nz, 1);
bcy(1:nx * ( ny + 1)) = bcs(1)*dZm(1)^2;
% bcy((end-((nx + 1) * ny) +1):end) = bcs(2);

bcz = zeros(nx  * ny * (nz+ 1), 1);% bcz([1,end]) = bcs;

BCx = sparse(bcx);

BCy = sparse(bcy);

BCz = sparse(bcz);


B = [BCx ; BCy; BCz];

%% Derivative test for primary field
a = randn(nface, 1);
phi = randn(mcell, 1);
u = [a;phi];
m = randn(mcell, 1);
q = randn(mcell+nface, 1);
v = randn(mcell, 1);

sdiag = @(XX) spdiags( XX, 0, length(XX), length(XX));

% GRAD(1) = 2/dx(1);
% Build forward operator
G_p = @(m) [L_p + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD_p
            DIV_p*S(m)      , DIV_p*S(m)*GRAD_p ];

% Build forward operator
G_s = @(m) [L_s + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD_s
            DIV_s*S(m)      , DIV_s*S(m)*GRAD_s ];


% C_p = @(m,u,q) G_p(m)*u - q;
% % Derivative with respect to model parameter
%     dSdm = @(m) S(m).^2 * AVF * sdiag(1./m.^2);
%     % Derivatives
%     dCdm_p = @(a,phi,m) [1i*w*uo*sdiag(a) * dSdm(m) + 1i*w*uo * sdiag(GRAD_p * phi) * dSdm(m)
%                        DIV_p * sdiag(a) * dSdm(m) + DIV_p * sdiag(GRAD_p * phi) * dSdm(m)];
% 
% % Sensitivity Matrix
% J_p = @(m,a,phi) -G_p(m)\dCdm_p(a,phi,m);
% Jacob=J_p(m,a,phi);
% 
% for ii=1:10
% h = 10^(-ii);
% 
% m2=( m + v*h );
% 
% 
% Oh = norm(C_p(m2,u,q) - C_p(m,u,q));
% Oh2 = norm(C_p(m2,u,q) - C_p(m,u,q) - dCdm_p(a,phi,m)*h*v);
% 
% fprintf('%3.2e    %3.2e    %3.2e\n',h,Oh,Oh2)
% end

%% Experiment using simple anomaly and a background conductivity
% Compute primary field
 
q = [ -B; zeros(mcell, 1)];

u = G_p(s)\q;

A = u(1: nface);
Phi = u( nface + 1: end);


E = A + GRAD_p*Phi;

figure;vecplot(E, dx, dy, dz)