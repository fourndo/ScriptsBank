function [L_p,L_s,DIV_p,DIV_s,GRAD,CURL_e,CURL_f,AVF,AVC,AVHx,AVHy,dZm]=get_ops_test(dx,dy,dz,w,s_background)

uo = 4*pi*10^-7; % Magnetic permeability of free-space



% % Create node & face location
% xf = [0;cumsum(dx)]; xf=xf(:); 
% yf = [0;cumsum(dy)]; yf=yf(:);
% zf = [0;cumsum(dz)]; zf=zf(:);

% Create center-center mesh size (hmid)
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% Compute number of face&nodes
nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);

% Compute number of cells
nx=length(dx) ; ny=length(dy) ; nz=length(dz);  

% Number of faces in every directions
nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

% Number of edges
nhx = (nx) * (ny+1) * (nz+1);
nhy = (nx+1) * (ny) * (nz+1);
nhz = (nx+1) * (ny+1) * (nz);

% Create hmid dimensions matrix
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

% Create cell dimension matrix
dX = spdiags(1./dx,0,nx,nx); 
dY = spdiags(1./dy,0,ny,ny); 
dZ = spdiags(1./dz,0,nz,nz); 

% Create over-sized cell dimension matrix with first and last dimension
% repeated for BCs.
dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);


%% Laplacian operator for primary field
% Second partial Derivative in 3 directions
% %Partial derivatives for x-component
d_dx=  dXxl.^0.5 * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= dXm * d_dx'*d_dx; dd_dx([1,end],:)=0;
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz; dd_dz([1,end],:)=0;

% kappa=sqrt(1i*w*uo*s_background); 
% dd_dz(1) = 3/2*dd_dz(1);
% dd_dz(end) = dd_dz(end)/2 + 1i*kappa/dZm(end,end);

DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

% %%Partial derivatives for y-component
d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX * d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

d_dy= dYxl.^0.5 * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= dYm * d_dy'*d_dy; dd_dy([1,end],:)=0;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

%Boundary Conditions: Derichlet on top and Robin at bottom
%Will eventually need to compute kappa for individual cells since
%conductivity will change at every location
d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz; dd_dz([1,end],:)=0;
 
% dd_dz(1) = 3/2*dd_dz(1);
% dd_dz(end) = dd_dz(end)/2 + 1i*kappa/dZm(end,end);

DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


% %%Partial derivatives for z-component
d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX * d_dx'*d_dx;  dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl.^0.5 * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= dZm * d_dz'*d_dz; dd_dz([1,end],:)=0;
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
d_dx=  dXxl.^0.5 * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= dXm * d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])*3/2;
DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])*3/2;
DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz; dd_dz([1,end])=dd_dz([1,end])*3/2;
DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

% %%Partial derivatives for y-component
d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX *d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])*3/2;
DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

d_dy= dYxl.^0.5 * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= dYm * d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])*3/2;
DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ * d_dz'*d_dz; dd_dz([1,end])=dd_dz([1,end])*3/2; 
DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


% %%Partial derivatives for z-component
d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX * d_dx'*d_dx; dd_dx([1,end])=dd_dx([1,end])*3/2;
DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY * d_dy'*d_dy; dd_dy([1,end])=dd_dy([1,end])*3/2;
DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

d_dz= dZxl.^0.5 * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= dZm * d_dz'*d_dz; dd_dz([1,end])=dd_dz([1,end])*3/2;
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

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=d_dz([1,end])*2;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

GRAD = [Dx; Dy; Dz];

%% Create primary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV_p = [Dx Dy Dz];

%% Create secondary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1); d_dx([1,end])=0;
Dx = kron(kron(speye(nz),speye(ny)),d_dx);

d_dy=dY * ddx(ny-1); d_dy([1,end])=0;
Dy = kron(kron(speye(nz),d_dy),speye(nx));

d_dz=dZ * ddx(nz-1); d_dz([1,end])=0;
Dz = kron(kron(d_dz,speye(ny)),speye(nx));

DIV_s = [Dx Dy Dz];

%% Averaging operator cell center to faces

av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avx = av(nx); avx([1,end]) = 1;
avy = av(ny); avy([1,end]) = 1;
avz = av(nz); avz([1,end]) = 1;

Avx = kron ( kron ( speye (nz) , speye (ny)), avx);
Avy = kron ( kron ( speye (nz) , avy), speye (nx));
Avz = kron ( kron ( avz , speye (ny)), speye (nx));

AVF = [Avx;Avy;Avz];


 
%% Averaging operator face to center

av = @(n) spdiags (ones (n+2,1)*[0.5,0.5],[-1,0],n+2,n+1);

avfx = av(nx); avfx = avfx(2:end-1,:);
avfy = av(ny); avfy = avfy(2:end-1,:); 
avfz = av(nz); avfz = avfz(2:end-1,:);

Avfx = kron ( kron ( speye (nz) , speye (ny)), avfx);
Avfy = kron ( kron ( speye (nz) , avfy), speye (nx));
Avfz = kron ( kron ( avfz , speye (ny)), speye (nx));

Ox = sparse ( nx*ny*nz , nfx);
Oy = sparse ( nx*ny*nz , nfy);
Oz = sparse ( nx*ny*nz , nfz);

AVC =   [Avfx Oy Oz
        Ox Avfy Oz
        Ox Oy Avfz];

%% Averaging operator for Hx and Hy on edge

av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avfx = av(nx+1); avfx = avfx(2:end-1,:);
avfy = av(ny+1); avfy = avfy(2:end-1,:); 
% avfz = av(nz+1); avfz = avfz(2:end-1,:);

Avyx = kron ( kron ( speye (nz+1) , avfy), speye (nx));
Avxy = kron ( kron ( speye (nz+1) , speye (ny)), avfx);
% Avfz = kron ( kron ( speye (nz) , avfy), speye (nx+1));

Oxx = sparse ( nfz , nhx );
Oxy = sparse ( nfz , nhy );
Oxz = sparse ( nfz , nhz );
Oyx = sparse ( nfy , nhx );
Oyy = sparse ( nfy , nhy );
Oyz = sparse ( nfy , nhz );
Ozx = sparse ( nfx , nhx );
Ozy = sparse ( nfx , nhy );
Ozz = sparse ( nfx , nhz );

AVHx =   [Avyx Oxy Oxz
           Oyx Oyy Oyz
           Ozx Ozy Ozz];
       
AVHy =   [Oxx Avxy Oxz
           Oyx Oyy Oyz
           Ozx Ozy Ozz];

%% Derivative matrices for curl (edges to faces)

d_dz = dZ * ddx(nz-1); %d_dz=d_dz(:,2:end-1); 
Dzy = kron( kron( d_dz , speye(ny) ), speye(nx+1) );

d_dy = dY * ddx(ny-1); %d_dy=d_dy(:,2:end-1); 
Dyz = kron( kron( speye(nz) , d_dy ), speye(nx+1) );

d_dz = dZ * ddx(nz-1); %d_dz=d_dz(:,2:end-1); 
Dzx = kron( kron( d_dz, speye(ny+1) ) , speye(nx) );

d_dx = dX * ddx(nx-1); %d_dx=d_dx(:,2:end-1); 
Dxz = kron( kron( speye(nz) , speye(ny+1) ), d_dx );

d_dy = dY * ddx(ny-1); %d_dy=d_dy(:,2:end-1);
Dyx = kron( kron( speye(nz+1) , d_dy ), speye(nx) );

d_dx = dX * ddx(nx-1); %d_dx=d_dx(:,2:end-1); 
Dxy = kron( kron( speye(nz+1) , speye(ny) ), d_dx );

Ox = sparse ( nfx , nhx );
Oy = sparse ( nfy , nhy );
Oz = sparse ( nfz , nhz );

CURL_e = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz];

%% Curl operator with Derichlet boundary conditions (faces to edge)
% The first and last vectors in the y and z-direction and set to 0
% The last vector in the x-direction is set to 0;


d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1);  d_dz([1,end])=d_dz([1,end])*2;
Dzy = kron( kron( d_dz , speye(ny+1) ), speye(nx) );

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dyz = kron( kron( speye(nz+1) , d_dy ), speye(nx) );

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=d_dz([1,end])*2;
Dzx = kron( kron( d_dz, speye(ny) ) , speye(nx+1) );

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dxz = kron( kron( speye(nz+1) , speye(ny) ), d_dx );

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=d_dy([1,end])*2;
Dyx = kron( kron( speye(nz) , d_dy ), speye(nx+1) );

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=d_dx([1,end])*2;
Dxy = kron( kron( speye(nz) , speye(ny+1) ), d_dx );

Ox = sparse ( nhx , nfx );
Oy = sparse ( nhy , nfy );
Oz = sparse ( nhz , nfz );

CURL_f = [Ox Dzy -Dyz;-Dzx Oy Dxz;Dyx -Dxy Oz];
end