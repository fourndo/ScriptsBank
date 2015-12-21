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
% 
% We want to compute 
% J = (iwu)^-1 * curl x curl x E - u * s * E
% With the right boundary matrix


clear all
close all
% define mesh
count=1;



%% Discretization Iterations

for ii = 1%:4
   
n = 2 ^ ii;

%% Create mesh for center and faces
X0 = 0; Y0 = 0; Z0 = 0; %Origin

%Create mesh size (cell center)
xf = linspace(0,1,n+1); xf=xf(:); 
yf = linspace(0,1,n+1); yf=yf(:);
zf = linspace(0,1,n+1); zf=zf(:);

dx = abs(xf(2:end) - xf(1:end-1));
dy = abs(yf(2:end) - yf(1:end-1));
dz = abs(zf(2:end) - zf(1:end-1));


%Create mesh (center to center)
xc = xf(1:end-1) + dx/2; xc=xc(:); 
yc = yf(1:end-1) + dy/2; yc=yc(:);
zc = zf(1:end-1) + dz/2; zc=zc(:);

dxn = dx(1:end-1)/2 + dx(2:end)/2; dxn=[dx(1);dxn;dx(end)];
dyn = dy(1:end-1)/2 + dy(2:end)/2; dyn=[dy(1);dyn;dy(end)];
dzn = dz(1:end-1)/2 + dz(2:end)/2; dzn=[dz(1);dzn;dz(end)];

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
dXn = spdiags(1./dxn,0,nxn,nxn);
dYn = spdiags(1./dyn,0,nyn,nyn);
dZn = spdiags(1./dzn,0,nzn,nzn);

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

d_dx = ddx(n);

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

%% Derivative matrices for curl (faces to edge)
% The first and last vectors in the y and z-direction and set to 0
% The last vector in the x-direction is set to 0;

d_dz = dZn * ddx(nz); d_dz=d_dz(:,2:end-1);  %d_dz([1,end])=d_dz([1,end])*2;
Dzy = kron( kron( d_dz , speye(ny+1) ), speye(nx) );

d_dy = dYn * ddx(ny); d_dy=d_dy(:,2:end-1); %d_dy([1,end])=d_dy([1,end])*2;
Dyz = kron( kron( speye(nz+1) , d_dy ), speye(nx) );

d_dz = dZn * ddx(nz); d_dz=d_dz(:,2:end-1); %d_dz([1,end])=d_dz([1,end])*2;
Dzx = kron( kron( d_dz, speye(ny) ) , speye(nx+1) );

d_dx = dXn * ddx(nx); d_dx=d_dx(:,2:end-1); %d_dx([1,end])=d_dx([1,end])*2;
Dxz = kron( kron( speye(nz+1) , speye(ny) ), d_dx );

d_dy = dYn * ddx(ny); d_dy=d_dy(:,2:end-1); %d_dy([1,end])=d_dy([1,end])*2;
Dyx = kron( kron( speye(nz) , d_dy ), speye(nx+1) );

d_dx = dXn * ddx(nx); d_dx=d_dx(:,2:end-1); %d_dx([1,end])=d_dx([1,end])*2;
Dxy = kron( kron( speye(nz) , speye(ny+1) ), d_dx );

Ox = sparse ( nhx , nfx );
Oy = sparse ( nhy , nfy );
Oz = sparse ( nhz , nfz );

CURLt = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz]';


%% Create center gradient (center to face)

d_dx = dXn * ddx(nx); d_dx=d_dx(:,2:end-1);
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy = dYn * ddx(ny); d_dy=d_dy(:,2:end-1);
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz = dZn * ddx(nz); d_dz=d_dz(:,2:end-1);
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

GRAD = [Dx; Dy; Dz];

%% Create divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1); d_dx(end,end)=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1); d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1); d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV = [Dx Dy Dz];


%% Laplacian operator for primary field
% Second partial Derivative in 3 directions
% %Partial derivatives for x-component
d_dx=  dXxl.^0.5 * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= dXn * d_dx'*d_dx; dd_dx([1,end],:)=0;
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYn.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZn.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz;

kappa=sqrt(1i*c.w*c.uo*m(1)); 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end)/2 + 1i*kappa/dZn(end,end);

DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

% %%Partial derivatives for y-component
d_dx= dXn.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX * d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

d_dy= dYxl.^0.5 * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= dYn * d_dy'*d_dy; dd_dy([1,end],:)=0;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZn.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz; 
 
dd_dz(1) = 3/2*dd_dz(1);
dd_dz(end) = dd_dz(end)/2 + 1i*kappa/dZn(end,end);

DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


% %%Partial derivatives for z-component
d_dx= dXn.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX * d_dx'*d_dx;  dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYn.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl.^0.5 * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= dZn * d_dz'*d_dz; dd_dz([1,end],:)=0;
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
% Averaging two direction perp to the averaging direction
Avx = kron ( kron ( avz , avy) , speye (nx));
Avy = kron ( kron ( avz , speye (ny)), avx);
Avz = kron ( kron ( speye (nz) , avy), avx);

AVE = [Avx;Avy;Avz];


 
%% Analytical solution of the form E = A + Vphi, where A if divergence free
%% Note to Ben: Start copying here. Might have to switch from f to n (e.g.
%% xn == xf)
% Define curl function
Ax = @(x,y,z)(cos(y));
Ay = @(x,y,z)(cos(z));
Az = @(x,y,z)(cos(x));

% Need to compute BC for Derichelet
e = @(n) ones(n, 1);
bcx = sparse(nx+1, 1); bcx([1,end]) = [-2 2];
bcy = sparse(ny+1, 1); bcy([1,end]) = [-2 2];
bcz = sparse(nz+1, 1); bcz([1,end]) = [-2 2];
    
% Get the x components of BC on edges
[Y,X,Z] = meshgrid (yf,xc,zf);
A_y = Ay(X,Y,Z);

BCz = kron(bcz, kron( e(nyn), e(nx)));

Hzy = BCz.*A_y(:) * dZ(1);


% Compute only face values for Dyz
[Y,X,Z] = meshgrid (yf,xc,zf);
A_z = Az(X,Y,Z);

BCy = kron(e(nzn), kron( bcy, e(nx)));
Hyz = BCy.*A_z(:) * dY(1);

% Compute the total Hx component, minus sign on the Hyz like the curl.
Hx = Hzy - Hyz;

% Get the y components of BC
[Y,X,Z] = meshgrid (yc,xf,zf);
A_z = Az(X,Y,Z);

BCx = kron(e(nzn), kron( e(ny), bcx));
Hxz = BCx.*A_z(:) * dX(1);



[Y,X,Z] = meshgrid (yc,xf,zf);
A_x = Ax(X,Y,Z);

BCz = kron(bcz, kron( e(ny), e(nxn)));
Hzx = BCz.*A_x(:) * dZ(1);

% Compute the total Hy component, minus sign on the Hzx like the curl.
Hy = Hxz - Hzx;

% Get the z components of BC
% [Y,X,Z] = meshgrid ([yf(1:end-2);yf(end)],xf,zc);
[Y,X,Z] = meshgrid (yf,xf,zc);
A_x = Ax(X,Y,Z);

BCy = kron(e(nz), kron( bcy, e(nxn)));
Hyx = BCy.*A_x(:)*dY(1);

% [Y,X,Z] = meshgrid (yf,[xf(1:end-2);xf(end)],zc);
[Y,X,Z] = meshgrid (yf,xf,zc);
A_y = Ay(X,Y,Z);

BCx = kron(e(nz), kron( e(nyn), bcx));
Hxy = BCx.*A_y(:) * dX(1);

% Compute the total Hz component, minus sign on the Hxy like the curl.
Hz = -Hxy + Hyx;


BC =[Hx;Hy;Hz];

%% Ben's test

    % Edge variables
    [Xie, Yie, Zie] = ndgrid(xc, yf, zf);
    [Xje, Yje, Zje] = ndgrid(xf, yc, zf);
    [Xke, Yke, Zke] = ndgrid(xf, yf, zc);

    % Face Variables
    [Xif, Yif, Zif] = ndgrid(xf, yc, zc);
    [Xjf, Yjf, Zjf] = ndgrid(xc, yf, zc);
    [Xkf, Ykf, Zkf] = ndgrid(xc, yc, zf);

    % Function for the field
    Aif = cos(Yif);
    Ajf = cos(Zjf);
    Akf = cos(Xkf);
    
    Af = [Aif(:); Ajf(:); Akf(:)];
    CURLAf = CURL * Af;
    
    % NUMERICAL: CURL A
    Ae = CURLAf + BC;

    % ANALYTICAL: V x A
    Bie = -sin(Zie);
    BieX = Bie;
    Bje = -sin(Xje);
    Bke = -sin(Yke);

    Be = [Bie(:); Bje(:); Bke(:)];

    % |residual|
    curl3n = norm(Ae - Be, 'inf'); 
    
    % Print to screen the result
    fprintf('%3.2e   %3.2e   \n',1/n,curl3n)
    
   

%% Eldad's test
%% Constants
% Centre variables

para = 10;

[Xc, Yc, Zc] = ndgrid(xc, yc, zc);

psi = @(x) tanh(para * (x + 1/4)) - tanh(para * (x - 1/4)) + 1/100;
m = psi(Xc(:)) .* psi(Yc(:)) .* psi(Zc(:));

c.uo = 4*pi*10^-7;
c.w = 1e7; % Frequency [hz]
k = 1i*c.w*c.uo;
skindepth = sqrt( 2/ (max(m) * c.uo * c.w));
diffusion = c.w * c.uo * max(m);

Ex = @(X,Y,Z) -Z .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2)) ./ psi(X);
Ey = @(X,Y,Z) -X .* Z .* exp(-5*(X.^2 + Y.^2 + Z.^2)) ./ psi(Y);
Ez = @(X,Y,Z) -X .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2)) ./ psi(Z);




%%
% A*u = b
%
% |L + iwuS  iwuSG |*| A   | = | -iwJ |
% |    DS     DSG  | | phi |   | -DJ  |
%
% where Js = -(iwu)^(-1) C'C*E - S(m)*E
% Harmonic averaging since m is discontinuous on the faces.

%% Boundary condition
% Get the x components of BC on edges
[Y,X,Z] = meshgrid (yf,xc,zf);
E_y = Ey(X,Y,Z);

BCz = kron(bcz, kron( e(nyn), e(nx)));

Hzy = BCz.*E_y(:) * dZ(1);


% Compute only face values for Dyz
[Y,X,Z] = meshgrid (yf,xc,zf);
E_z = Ez(X,Y,Z);

BCy = kron(e(nzn), kron( bcy, e(nx)));
Hyz = BCy.*E_z(:) * dY(1);

% Compute the total Hx component, minus sign on the Hyz like the curl.
Hx = Hzy - Hyz;


% Get the y components of BC
[Y,X,Z] = meshgrid (yc,xf,zf);
E_z = Ez(X,Y,Z);

BCx = kron(e(nzn), kron( e(ny), bcx));
Hxz = BCx.*E_z(:) * dX(1);



[Y,X,Z] = meshgrid (yc,xf,zf);
E_x = Ex(X,Y,Z);

BCz = kron(bcz, kron( e(ny), e(nxn)));
Hzx = BCz.*E_x(:) * dZ(1);

% Compute the total Hy component, minus sign on the Hzx like the curl.
Hy = Hxz - Hzx;


% Get the z components of BC
% [Y,X,Z] = meshgrid ([yf(1:end-2);yf(end)],xf,zc);
[Y,X,Z] = meshgrid (yf,xf,zc);
E_x = Ex(X,Y,Z);

BCy = kron(e(nz), kron( bcy, e(nxn)));
Hyx = BCy.*E_x(:)*dY(1);

% [Y,X,Z] = meshgrid (yf,[xf(1:end-2);xf(end)],zc);
[Y,X,Z] = meshgrid (yf,xf,zc);
E_y = Ey(X,Y,Z);

BCx = kron(e(nz), kron( e(nyn), bcx));
Hxy = BCx.*E_y(:) * dX(1);

% Compute the total Hz component, minus sign on the Hxy like the curl.
Hz = -Hxy + Hyx;


BC =[Hx;Hy;Hz];

%% 
nfaces = size(CURL,2);

S = @(m) spdiags( (AVF * (1./m)).^-1, 0, nfaces, nfaces ) ;

E_x = Ex(Xif, Yif, Zif);
E_y = Ex(Xjf, Yjf, Zjf);
E_z = Ex(Xkf, Ykf, Zkf);

E = [E_x(:); E_y(:); E_z(:)];

% EM PDE Matrix
A = @(m) [L_p + k*S(m), k*S(m)*GRAD
          DIV * S(m)  , DIV*S(m)*GRAD]; 
      
Js = -1/k * CURLt * (CURL * E + BC)+ - S(m) * E;     
% Js = S(m) * E;

b = [ -k*Js; -DIV * Js];

% Forward Operator -> Primary Fields
u = A(m)\b;
% Obtain fields from output
[a_n, phi_n] = vecsplit(u, nfaces, ncells);

%% Pseudo Analytic solution

u = [speye(nfaces) G] \ E;
[a_r, phi_r] = vecsplit(u, nfaces, ncells);

da(ii) = norm(a_n - a_r, 'inf'); %#ok<*SAGROW>
dphi(ii) = norm(phi_r - phi_n, 'inf');
diva(ii) = norm(D * a_r, 'inf');

end 

loglog(1./h, [da', dphi', diva']);
legend('a','phi','DIV a','Location','Best')