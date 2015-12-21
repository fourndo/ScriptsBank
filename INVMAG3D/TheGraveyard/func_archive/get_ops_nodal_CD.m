function [AVEN,AVCE,AVCN,GRAD,DIV,DIVbc,Q,R,RI,D,Sz,Sx,Sy,H0,Wr,V] = get_ops_nodal(X0, Y0, Z0, dx, dy, dz, ObsX, ObsY, ObsZ, H, I, D)
% [AVF,GRAD,DIV,B0,Wr,dV] = get_ops(X0, Y0, Z0, dx, dy, dz, H, I, D)
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Use sparse operators
% Inputs:
% Xo, Yo, Zo : Coordinates of the South-West-Top corner of the mesh
% dX, dY, dZ : Vectors of cell sizes for the East, North and Vertical axis
% X,  Y,  Z  : Coordinates of the observation points 
%      H     : Magnitude of inducing field (nT)
%      I     : Inclinaison of H (degrees from horizontal)
%      D     : Declinaison of H (degrees from North)
% nullcell   : 1D vector for active(1) and inactive (0) cells
%
% Reference: 
% Lelievre & Oldenburg 2006 Magnetic forward modelling and inversion for
% high susceptibility. Geophysics J. Int. 166: 76-90

%% Create mesh for testing

% clear all 
% close all
% dx = [1 1 1]*50;
% dy = [1 1 1]*50;
% dz = [1 1 1]*50;
% nx= length(dx);
% ny= length(dy);
% nz= length(dz);
% X0 = 0;
% Y0 = 0;
% Z0 = 0;
% 
% %3D density contrast model
% 
% % figure (1)
% 
% [model] = zeros(nz,nx,ny);
% model(2,2,2)=1.0;
% m=reshape(model,nx*ny*nz,1);
% 
% % Discretized topography
% topo_model = ones(nx,ny)*0;
% 
% % 
% % Create active cell matrix (will later need to discretize topography beforehand)
% nullcell = make_nullcell(dx,dy,dz,X0,Y0,Z0,topo_model);
% 
% mcell = nullcell(end);
% 
% %Create data points
% % Center observations for simplicity
% H = 50000; % Inducing field
% I = 45 ; % Field inclinaison (degree from h-plane)
% D = 0 ;  % Field declinaison (degree from North)
% cellx = cumsum(dx); celly = cumsum(dy); ndx = length(cellx); ndy = length(celly);
% [ObsX, ObsY, ObsZ] = meshgrid(cellx((floor(ndx/2)):(ceil(ndx/2)))-min(dx)/3,celly((floor(ndy/2)):(ceil(ndy/2)))-min(dy)/3,5);
% ObsX = reshape(ObsX,size(ObsX,1) * size(ObsX,2), 1);
% ObsY = reshape(ObsY,length(ObsX),1);
% ObsZ = reshape(ObsZ,length(ObsX),1);
% 
% 
% ndata = length(ObsX);

%% ///Script starts here/// %%
R0 = min( [min(dx) min(dy) min(dz)] ) / 4;

ndata = length(ObsX);

% Inducing field parameters
hy = H * cosd(I)*cosd(D);
hx = H * cosd(I)*sind(D);
hz = H * sind(I);

% Cell size vectors
dx = dx(:);
dy = dy(:);
dz = dz(:);

% Number of cells
nx =length(dx);
ny =length(dy);
nz =length(dz);

mcell = nx * ny *nz;

% Number of cell faces
nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

nfaces = nfx + nfy + nfz;

% Number of cell edges
nex = (nx) * (ny+1) * (nz+1);
ney = (nx+1) * (ny) * (nz+1);
nez = (nx+1) * (ny+1) * (nz);

nedges = nex + ney + nez;
nnodes = (nx+1)*(ny+1)*(nz+1);

% Create center-center mesh size (hmid)
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% Compute number of face&nodes
nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);



%% Create hmid dimensions matrix
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

%% Create cell dimension matrix
dX = spdiags(1./dx,0,nx,nx); 
dY = spdiags(1./dy,0,ny,ny); 
dZ = spdiags(1./dz,0,nz,nz); 


% Cell-center location
x_cc = X0 + cumsum(dx) - dx/2;
y_cc = Y0 + cumsum(dy) - dy/2;
z_cc = Z0 - (cumsum(dz) - dz/2);

% Cell-node location
x_cn = [X0;X0 + cumsum(dx)];
y_cn = [Y0;Y0 + cumsum(dy)];
z_cn = [Z0;Z0 - cumsum(dz)];


%% Select a plane of edges for data (Research only)
b = zeros(1,nz); b(7)=1;
Sz = kron(speye(ny+1),kron(speye(nx+1),b));

b = zeros(1,nz+1); b(7)=1;
Sy = kron(kron(speye(ny),speye(nx+1)),b);

b = zeros(1,nz+1); b(7)=1;
Sx = kron(kron(speye(ny+1),speye(nx)),b);

Sz = [Sz sparse( size(Sz,1) , nex ) sparse( size(Sz,1) , ney )];
Sx = [sparse( size(Sx,1) , nez ) Sx sparse( size(Sx,1) , ney )];
Sy = [sparse( size(Sy,1) , nez ) sparse( size(Sy,1) , nex ) Sy];

%% Sparse matrix functions
e = @ (n) ones(1, n);

Ox = sparse ( nnodes , nex );
Oy = sparse ( nnodes , ney );
Oz = sparse ( nnodes , nez );

O = sparse(ndata,nnodes);

kron3 = @(a,b,c) kron( a , kron( b , c ) );

% %% Averaging operator from face to center
% av = @(n) spdiags (ones (n+2,1)*[0.5,0.5],[-1,0],n+2,n+1);
% 
% avfx = av(nx); avfx = avfx(2:end-1,:);
% avfy = av(ny); avfy = avfy(2:end-1,:); 
% avfz = av(nz); avfz = avfz(2:end-1,:);
% 
% Avfx = kron ( kron ( speye (ny), avfx), speye (nz) );
% Avfy = kron ( kron ( avfy , speye (nx)), speye (nz));
% Avfz = kron ( kron ( speye (ny) , speye (nx)), avfz );
% 
% 
% AVC =   [Avfz Ox Oy
%         Oz Avfx Oy
%         Oz Ox Avfy];

%% Averaging operator center to face, then face to edges
av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avcfx = av(nx); avcfx([1,end]) = 1;
avcfy = av(ny); avcfy([1,end]) = 1;
avcfz = av(nz); avcfz([1,end]) = 1;

Avcfx = kron ( kron ( speye (ny) , avcfx ), speye (nz));
Avcfy = kron ( kron ( avcfy , speye (nx) ), speye (nz));
Avcfz = kron ( kron ( speye (ny) , speye (nx) ), avcfz);

avfex = av(nx); avfex([1,end]) = 1;
avfey = av(ny); avfey([1,end]) = 1;
avfez = av(nz); avfez([1,end]) = 1;

Avfex = kron ( kron ( speye (ny) , avfex), speye (nz+1) );
Avfey = kron ( kron ( avfey , speye (nx+1) ), speye (nz) );
Avfez = kron ( kron ( speye (ny+1) , speye (nx)), avfez );

AVCE = [Avfey*Avcfx;Avfez*Avcfy;Avfex*Avcfz];

%% Averaging operator edge to node

avenx = av(nx); avenx([1,end]) = 1;
aveny = av(ny); aveny([1,end]) = 1;
avenz = av(nz); avenz([1,end]) = 1;

Avenx = kron ( kron ( speye (ny+1) , avenx ), speye (nz+1) );
Aveny = kron ( kron ( aveny , speye (nx+1) ), speye (nz+1) );
Avenz = kron ( kron ( speye (ny+1) , speye (nx+1)), avenz );


AVEN = [Avenz Ox Oy;
        Oz  Avenx Oy;
        Oz Ox Aveny];

%% Averaging operator center to node

AVCN = [Avenz*Avfey*Avcfx + Avenx*Avfez*Avcfy + Aveny*Avfex*Avcfz]/2;

%% Create inducing field matrix
Hx= ones(1,nex) * hx;
Hy= ones(1,ney) * hy;
Hz= ones(1,nez) * hz;

H0 = [Hz';Hx';Hy']; 
% H0 = spdiags(h0,0,nedges,nedges);

% hx= ones(1,mcell) * (p);
% hy= ones(1,mcell) * (q);
% hz= ones(1,mcell) * (r);
% h1 = [hz hx hy];
% 
% Hx = spdiags(hx'*H,0,mcell,mcell);
% Hy = spdiags(hy'*H,0,mcell,mcell);
% Hz = spdiags(hz'*H,0,mcell,mcell);
% 
% H1 = [Hz Hy Hx];
%% Create gradient op (node to edge)
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n,n+1);

d_dx = dX * ddx(nx); %d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
Dx = kron(kron(speye(ny+1),d_dx),speye(nz+1));

d_dy = dY * ddx(ny); %d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx+1)),speye(nz+1));

d_dz = dZ * ddx(nz); %d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
Dz = kron(kron(speye(ny+1),speye(nx+1)),d_dz);

GRAD = [Dz; Dx; Dy];

% %% Create gradient op (face to center)
% ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);
% d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
% Dx = kron(kron(speye(ny),d_dx),speye(nz));
% 
% d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
% Dy = kron(kron(d_dy,speye(nx)),speye(nz));
% 
% d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
% Dz = kron(kron(speye(ny),speye(nx)),d_dz);
% 
% 
% 
% GRAD_ftc = [Dz Ox Ox;
%             Oz Dx Oy;
%             Oz Ox Dy];

%% Create divergence (edge to node)
% First and last column of divergence are set to 1 for BC
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);

d_dx=dXm * ddx(nx); d_dx=d_dx(:,2:end-1);
Dx = kron(kron(speye(ny+1),d_dx),speye(nz+1));

d_dy=dYm * ddx(ny); d_dy=d_dy(:,2:end-1);
Dy = kron(kron(d_dy,speye(nx+1)),speye(nz+1));

d_dz=dZm * ddx(nz); d_dz=d_dz(:,2:end-1);
Dz = kron(kron(speye(ny+1),speye(nx+1)),d_dz);

DIV = [Dz Dx Dy];

%% Create divergence operator for BC (edge to node)
% First and last column of divergence are set to 1 for BC
ddx = @(n) spdiags (zeros (n+1,1)*[1,-1],[0,1],n+1,n+2);

d_dx=ddx(nx); d_dx=d_dx(:,2:end-1); d_dx(1)=1/dxm(1) ; d_dx(end)=-1/dxm(end);
Dx = kron(kron(speye(ny+1),d_dx),speye(nz+1));

d_dy=ddx(ny); d_dy=d_dy(:,2:end-1); d_dy(1)=1/dym(1) ; d_dy(end)=-1/dym(end);
Dy = kron(kron(d_dy,speye(nx+1)),speye(nz+1));

d_dz=ddx(nz); d_dz=d_dz(:,2:end-1); d_dz(1)=1/dzm(1) ; d_dz(end)=-1/dzm(end);
Dz = kron(kron(speye(ny+1),speye(nx+1)),d_dz);

DIVbc = [Dz Dx Dy];
 
%% Create cell dimension matrix

dX = kron3( e(ny) , dx' , e(nz) ); 
dY = kron3( dy' , e(nx), e(nz) ); 
dZ = kron3( e(ny) , e(nx) , dz' );

V = dX.*dY.*dZ;

% V = spdiags(dV',0,mcell,mcell);
% 
% VV = [V O O;
%       O V O;
%       O O V];
  
% % Volume descritization for faces 
% dxf = dx'/2;dyf = dy'/2;dzf = dz'/2;
% dXf = kron3( e(ny) , [dxf(1)*2;(dxf(1:end-1) + dxf(2:end));dxf(end)*2] , e(nz) );
% dYf = kron3( [dyf(1)*2;(dyf(1:end-1) + dyf(2:end));dyf(end)*2] , e(nx) , e(nz) );
% dZf = kron3( e(ny) , e(nx) , [dzf(1)*2;(dzf(1:end-1) + dzf(2:end));dzf(end)*2] );

%% Construct Boundary Condition Matrix (H on nodes)
% Boundary matrix 
qx = zeros(1,nx); qx([1,end]) = hx ; 
Qx = kron3( e(ny+1) , qx, e(nz+1) );

qy = zeros(1,ny); qy([1,end]) = hy ; 
Qy = kron3( qy , e(nx+1), e(nz+1) );

qz = zeros(1,nz); qz([1,end]) = hz ; 
Qz = kron3( e(ny+1) , e(nx+1) , qz );

% Boundary matrix - cell center
Q = [Qz';Qx';Qy'] ;

% Q = spdiags(Q,0,mcell,mcell);
%% Create R: distance obs to center faces 
% Rz = @(ObX,ObY,ObZ)...
%        ( (ObX - kron3( e(ny) , x_cc', e(nz+1) )).^2 +...
%        (ObY - kron3( y_cc' , e(nx), e(nz+1) )).^2 +...
%        (ObZ - kron3( e(ny) , e(nx), z_ce' )).^2 ) .^0.5 ;
% 
% Rx = @(ObX,ObY,ObZ)...
%        ( (ObX - kron3( e(ny) , x_ce', e(nz) )).^2 +...
%        (ObY - kron3( y_cc' , e(nx+1), e(nz) )).^2 +...
%        (ObZ - kron3( e(ny) , e(nx+1), z_cc' )).^2 ) .^0.5 ;
% 
% Ry = @(ObX,ObY,ObZ)...
%        ( (ObX - kron3( e(ny+1) , x_cc', e(nz) )).^2 +...
%        (ObY - kron3( y_ce' , e(nx), e(nz) )).^2 +...
%        (ObZ - kron3( e(ny+1) , e(nx), z_cc' )).^2 ) .^0.5 ;
% 
% Rf = [Rz(ObsX(1),ObsY(1),ObsZ(1))';Ry(ObsX(1),ObsY(1),ObsZ(1))';...
%     Rx(ObsX(1),ObsY(1),ObsZ(1))'] ; 
% 
% Rf = spdiags(1./Rf,0,nfaces,nfaces);




%% Compute depth weighting and harmonic averaging of data.

Wr=zeros(1,mcell);
Rx = spalloc( ndata , nnodes , 3*ndata );
Ry = spalloc( ndata , nnodes , 3*ndata );
Rz = spalloc( ndata , nnodes , 3*ndata );

for ii = 1 : ndata
    
    Rc = ( (ObsX(ii) - kron3( e(ny) , x_cc', e(nz) )).^2 +...
           (ObsY(ii) - kron3( y_cc' , e(nx), e(nz) )).^2 +...
           (ObsZ(ii) - kron3( e(ny) , e(nx), z_cc' )).^2 ) .^0.5 ;
    
    Rn = ( (ObsX(ii) - kron3( e(ny+1) , x_cn', e(nz+1) )).^2 +...
           (ObsY(ii) - kron3( y_cn' , e(nx+1), e(nz+1) )).^2 +...
           (ObsZ(ii) - kron3( e(ny+1) , e(nx+1), z_cn' )).^2 ) .^0.5 ;
   
    % Find closest cell and weighted harmonic average the components
    % Weights of factor 1/R 
    id = min(find(Rn==min(Rn)));
    Rx(ii,id)= 1/Rn(id); Rx(ii,id+(nz+1))= 1/Rn(id+(nz+1)); Rx(ii,id-(nz+1))= 1/Rn(id-(nz+1));
    Ry(ii,id)= 1/Rn(id); Ry(ii,id+(nx+1)*(nz+1))= 1/Rn(id+(nx+1)*(nz+1)); Ry(ii,id-(nx+1)*(nz+1))= 1/Rn(id-(nx+1)*(nz+1));
    Rz(ii,id)= 1/Rn(id); Rz(ii,id+1)= 1/Rn(id+1); Rz(ii,id-1)= 1/Rn(id-1);
        
    %% Build depth weighting matrix
    Wr = Wr + ( V ./ (Rc + R0) .^ 3 ) .^ 2 ;
    
    
end

Wr = Wr(:);

R = [Rx O O;
    O Ry O;
    O O Rz];

D = kron ( [hz hx hy] /50000,speye(ndata)  );
% D = kron ( [1 1 1] ,speye(ndata)  );
RI = spdiags(1./(R*ones(3*nnodes,1)),0,ndata*3,ndata*3);


%% Build forward operator 

% Rc = spdiags(1./Rc,0,mcell,mcell);
% Rc2 = [Rc O O;
%       O Rc O;
%       O O Rc];
  
% E = spdiags([ones(nfz,1)*-1;ones(nfx,1);ones(nfy,1)],0,nfaces,nfaces);
% FwrOp = -ones(1,mcell*3) /(4*pi) * VV * AVC * GRAD_ctf  * Rc * H1 * AVC *  GRAD_ctf;
% FwrOp = -h1 /(4*pi) * VV * Rf *  GRAD_ctf * H1 * AVC *  GRAD_ctf;
% FwrOp = -ones(1,mcell*3) * VV * GRAD_ftc *  H0 * GRAD_ctf * Rc;
% FwrOp = FwrOp';