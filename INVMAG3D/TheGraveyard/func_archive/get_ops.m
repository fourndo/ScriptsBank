function [AVF,AVC,GRAD,DIVp,DIVs,B0,Q,P,PI,D,S] = get_ops(X0, Y0, Z0, dx, dy, dz, ObsX, ObsY, ObsZ, H, I, D)
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

%% ///Script starts here/// %%
R0 = min( [min(dx) min(dy) min(dz)] ) / 4;
% mcell = nX * nY * nZ;
uo = 4 * pi * 10^-7;

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

% Create center-center mesh size (hmid)
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% Compute number of face&nodes
nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);

%% Select a plane of surfaces for data (Research only)
b = zeros(1,nz+1); b(3)=1;
Sz = kron(speye(ny),kron(speye(nx),b));

b = zeros(1,ny+1); b(3)=1;
Sy = kron(b,kron(speye(nx),speye(nz)));

b = zeros(1,nx+1); b(3)=1;
Sx = kron(speye(ny),kron(b,speye(nz)));



S = [Sz zeros(nx*ny,nfx) zeros(nx*ny,nfy)];

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

% Cell-edge location
x_ce = [X0;X0 + cumsum(dx)];
y_ce = [Y0;Y0 + cumsum(dy)];
z_ce = [Z0;Z0 - cumsum(dz)];

%% Sparse matrix functions
e = @ (n) ones(1, n);

Ox = sparse ( mcell , nfx );
Oy = sparse ( mcell , nfy );
Oz = sparse ( mcell , nfz );

O = sparse(ndata,mcell);

kron3 = @(a,b,c) kron( a , kron( b , c ) );

%% Averaging operator from face to center
av = @(n) spdiags (ones (n+2,1)*[0.5,0.5],[-1,0],n+2,n+1);

avfx = av(nx); avfx = avfx(2:end-1,:);
avfy = av(ny); avfy = avfy(2:end-1,:); 
avfz = av(nz); avfz = avfz(2:end-1,:);

Avfx = kron ( kron ( speye (ny), avfx), speye (nz) );
Avfy = kron ( kron ( avfy , speye (nx)), speye (nz));
Avfz = kron ( kron ( speye (ny) , speye (nx)), avfz );


AVC =   [Avfz Ox Oy
        Oz Avfx Oy
        Oz Ox Avfy];

%% Averaging operator cell center to faces (for cell-center scalar)

av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n+1,n);

avx = av(nx); avx([1,end]) = 1;
avy = av(ny); avy([1,end]) = 1;
avz = av(nz); avz([1,end]) = 1;

Avx = kron ( kron ( speye (ny) , avx ), speye (nz));
Avy = kron ( kron ( avy , speye (nx) ), speye (nz));
Avz = kron ( kron ( speye (ny) , speye (nx)), avz );

AVF = [Avz;Avx;Avy];

%% Create inducing field matrix
Hx= ones(1,nfx) * hx;
Hy= ones(1,nfy) * hy;
Hz= ones(1,nfz) * hz;

B0 = [Hz';Hx';Hy']; 
%B0 = spdiags(h0',0,nfaces,nfaces);

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
%% Create gradient op (center to face)
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);

d_dx = dXm * ddx(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
Dx = kron(kron(speye(ny),d_dx),speye(nz));

d_dy = dYm * ddx(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx)),speye(nz));

d_dz = dZm * ddx(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
Dz = kron(kron(speye(ny),speye(nx)),d_dz);

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

%% Create primary divergence (face to center)
% First and last column of divergence are set to 1 for BC
d_dx=dX * ddx(nx-1); d_dx([1,end])=0;
Dx = kron(kron(speye(ny),d_dx),speye(nz));

d_dy=dY * ddx(ny-1); d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx)),speye(nz));

d_dz=dZ * ddx(nz-1); d_dz([1,end])=0;
Dz = kron(kron(speye(ny),speye(nx)),d_dz);

DIVp = [Dz Dx Dy];

%% Create secondary divergence (face to center)
% First and last column of divergence are set to 1 for BC
d_dx=dX * ddx(nx-1);% d_dx([1,end])=0;
Dx = kron(kron(speye(ny),d_dx),speye(nz));

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx)),speye(nz));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(speye(ny),speye(nx)),d_dz);

DIVs = [Dz Dx Dy];


 
%% Create cell dimension matrix

dX = kron3( e(ny) , dx' , e(nz) ); 
dY = kron3( dy' , e(nx), e(nz) ); 
dZ = kron3( e(ny) , e(nx) , dz' );

dV = dX.*dY.*dZ;

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

%% Construct Boundary Condition Matrix
% Boundary matrix 
qx = zeros(1,nx); qx(1) = -hx ./ dx(1)'; qx(end) = hx ./ dx(1)';
Qx = kron3( e(ny) , qx, e(nz) );

qy = zeros(1,ny); qy(1) = -hy ./ dy(1)'; qy(end) = hy ./ dy(1)';
Qy = kron3( qy , e(nx), e(nz) );

qz = zeros(1,nz); qz(1) = -hz ./ dz(1)'; qz(end) = hz ./ dz(1)';
Qz = kron3( e(ny) , e(nx) , qz );

% Boundary matrix - cell center
Q = [Qz'+Qx'+Qy'] ;

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
Px = spalloc( ndata , mcell , 3*ndata );
Py = spalloc( ndata , mcell , 3*ndata );
Pz = spalloc( ndata , mcell , 3*ndata );

for ii = 1 : ndata
    
    Rc = ( (ObsX(ii) - kron3( e(ny) , x_cc', e(nz) )).^2 +...
           (ObsY(ii) - kron3( y_cc' , e(nx), e(nz) )).^2 +...
           (ObsZ(ii) - kron3( e(ny) , e(nx), z_cc' )).^2 ) .^0.5 ;
   
    % Find closest cell and weighted harmonic average the components
    % Weights of factor 1/R 
    id = min(find(Rc==min(Rc)));
    Px(ii,id)= 1/Rc(id); Px(ii,id+nz)= 1/Rc(id+nz); Px(ii,id-nz)= 1/Rc(id-nz);
    Py(ii,id)= 1/Rc(id); Py(ii,id+nx*nz)= 1/Rc(id+nx*nz); Py(ii,id-nx*nz)= 1/Rc(id-nx*nz);
    Pz(ii,id)= 1/Rc(id); Pz(ii,id+1)= 1/Rc(id+1); Pz(ii,id-1)= 1/Rc(id-1);
        
    %% Build depth weighting matrix
    Wr = Wr + ( dV ./ (Rc + R0) .^ 3 ) .^ 2 ;
    
    
end

Wr = Wr(:);

P = [Px O O;
    O Py O;
    O O Pz];

D = kron ( [hz hx hy] /50000,speye(ndata)  );
% D = kron ( [1 1 1] ,speye(ndata)  );
PI = spdiags(1./(P*ones(3*mcell,1)),0,ndata*3,ndata*3);


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