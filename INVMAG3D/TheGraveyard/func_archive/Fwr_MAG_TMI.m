function [FwrOp,Wr,dV,del,H0] = FVD_Fwr_MAG(mcell, X0, Y0, Z0, dx, dy, dz,...
                           ObsX, ObsY, ObsZ, H, I, D,m, nullcell)
% FwdMag(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Use sparse operators for linear approximation
% Inputs:
% Xo, Yo, Zo : Coordinates of the South-West-Top corner of the mesh
% dX, dY, dZ : Vectors of cell sizes for the East, North and Vertical axis
% X,  Y,  Z  : Coordinates of the observation points 
%      H     : Magnitude of inducing field (nT)
%      I     : Inclinaison of H (degrees from horizontal)
%      D     : Declinaison of H (degrees from North)
% nullcell   : 1D vector for active(1) and inactive (0) cells
% 
% Last Update: September 17, 2013

%% ///Script starts here/// %%
e = @ (n) ones(1, n);
kron3 = @(a,b,c) kron( a , kron( b , c ) );


%% Mesh parameters
% Number of cells
nx =length(dx);
ny =length(dy);
nz =length(dz);

mcell = nx * ny *nz;

% Number of faces in each direction
nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

% Create center-center mesh size (hmid)
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% Compute number of face&nodes
nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);

% Create hmid dimensions matrix
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

% Create cell dimension matrix
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

%% Create cell dimension matrix
dimX = kron3( e(ny) , dx' , e(nz) ); 
dimY = kron3( dy' , e(nx), e(nz) ); 
dimZ = kron3( e(ny) , e(nx) , dz' );

dV = dimX.*dimY.*dimZ;

V = spdiags(dV',0,mcell,mcell);

%% Create inducing field matrix
% Decompose field in orthogonal components
hy = cosd(I)*cosd(D);
hx = cosd(I)*sind(D);
hz = sind(I);

H0 = [ones(nfz,1) * hz ; ones(nfx,1) * hx ;ones(nfy,1) * hy] * H; 

%% Averaging operator from face to center
av = @(n) spdiags (ones (n+2,1)*[0.5,0.5],[-1,0],n+2,n+1);

avfx = av(nx); avfx = avfx(2:end-1,:);
avfy = av(ny); avfy = avfy(2:end-1,:); 
avfz = av(nz); avfz = avfz(2:end-1,:);

Avfx = kron ( kron ( speye (ny), avfx), speye (nz) );
Avfy = kron ( kron ( avfy , speye (nx)), speye (nz));
Avfz = kron ( kron ( speye (ny) , speye (nx)), avfz );

AVC =   [V*Avfz*hz V*Avfx*hx V*Avfy*hy];

%% Create gradient op (face to center)
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(ny),d_dx),speye(nz));

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx)),speye(nz));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(speye(ny),speye(nx)),d_dz);
        
%% Create R: distance obs to center faces 
rz = @(ObX,ObY,ObZ)...
       ( (ObX - kron3( e(ny) , x_cc', e(nz+1) )).^2 +...
       (ObY - kron3( y_cc' , e(nx), e(nz+1) )).^2 +...
       (ObZ - kron3( e(ny) , e(nx), z_ce' )).^2 ) .^0.5 ;

rx = @(ObX,ObY,ObZ)...
       ( (ObX - kron3( e(ny) , x_ce', e(nz) )).^2 +...
       (ObY - kron3( y_cc' , e(nx+1), e(nz) )).^2 +...
       (ObZ - kron3( e(ny) , e(nx+1), z_cc' )).^2 ) .^0.5 ;

ry = @(ObX,ObY,ObZ)...
       ( (ObX - kron3( e(ny+1) , x_cc', e(nz) )).^2 +...
       (ObY - kron3( y_ce' , e(nx), e(nz) )).^2 +...
       (ObZ - kron3( e(ny+1) , e(nx), z_cc' )).^2 ) .^0.5 ;

Rz=rz(ObsX(1),ObsY(1),ObsZ(1))';
Rx=rx(ObsX(1),ObsY(1),ObsZ(1))';
Ry=ry(ObsX(1),ObsY(1),ObsZ(1))';

% Radial distances (1/r) for each face components
Rfz = spdiags(1./Rz,0,nfz,nfz);
Rfx = spdiags(1./Rx,0,nfx,nfx);
Rfy = spdiags(1./Ry,0,nfy,nfy);


% R from center cell to obs for distance weighting matrix (Wr)
Rc =    ( (ObsX(1) - kron3( e(ny) , x_cc', e(nz) )).^2 +...
       (ObsY(1) - kron3( y_cc' , e(nx), e(nz) )).^2 +...
       (ObsZ(1) - kron3( e(ny) , e(nx), z_cc' )).^2 ) .^0.5 ;
   
Rc =Rc(:);
%% Operator for DEL (DEL R)
% Z-derivative
dz_dz=  dZm * ddx(nz); dz_dz=dz_dz(:,2:end-1); %dz_dz([1,end])=0;
Dzz = kron( kron( speye(ny) , speye(nx) ), dz_dz );

% X-component
dx_dx=  dXm * ddx(nx); dx_dx=dx_dx(:,2:end-1); %dx_dx([1,end])=0;
Dxx = kron( kron( speye(ny) , dx_dx ), speye(nz) );

% Y-component
dy_dy=  dYm * ddx(ny); dy_dy=dy_dy(:,2:end-1); %dy_dy([1,end])=0;
Dyy = kron( kron( dy_dy  , speye(nx) ), speye(nz) );


del = [Dzz*Dz*Rfz Dzz*Dx*Rfx Dzz*Dy*Rfy;
     Dxx*Dz*Rfz Dxx*Dx*Rfx Dxx*Dy*Rfy;
     Dyy*Dz*Rfz Dyy*Dx*Rfx Dyy*Dy*Rfy];
 


%% Build forward operator 

FwrOp = 1 /(4*pi) * H0' * (AVC  * (del))' ;


%% Build depth weighting matrix
R0 = min( [min(dx) min(dy) min(dz)] ) / 4;
Wr = ( dV ./ ((Rc' + R0) .^ 3) ) .^ 2 ;
Wr = Wr(:);
