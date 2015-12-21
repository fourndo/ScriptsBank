function [FwrOp,Wr,dV] = FwdMag_linear(mcell, X0, Y0, Z0, dx, dy, dz,...
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
R0 = min( [min(dx) min(dy) min(dz)] ) / 4;

dx = dx(:);
dy = dy(:);
dz = dz(:);

%%
nx =length(dx);
ny =length(dy);
nz =length(dz);

mcell = nx * ny *nz;

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

%% Create hmid dimensions matrix
dXm = spdiags(1./dxm,0,nxm,nxm);
dYm = spdiags(1./dym,0,nym,nym);
dZm = spdiags(1./dzm,0,nzm,nzm);

dXmxl = spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
dYmxl = spdiags(1./dym,0,nym,nym);
dZmxl = spdiags(1./dzm,0,nzm,nzm);

% Create cell dimension matrix
dX = spdiags(1./dx,0,nx,nx); 
dY = spdiags(1./dy,0,ny,ny); 
dZ = spdiags(1./dz,0,nz,nz); 
% 
% % Create over-sized cell dimension matrix with first and last dimension
% % repeated for BCs.
% dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
% dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
% dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);

%% Create cell dimension matrix
% dX = spdiags(1./dx,0,nx,nx); 
% dY = spdiags(1./dy,0,ny,ny); 
% dZ = spdiags(1./dz,0,nz,nz); 

e = @ (n) ones(1, n);
kron3 = @(a,b,c) kron( a , kron( b , c ) );

% Cell-center location
x_cc = X0 + cumsum(dx) - dx/2;
y_cc = Y0 + cumsum(dy) - dy/2;
z_cc = Z0 - (cumsum(dz) - dz/2);

% Cell-edge location
x_ce = [X0;X0 + cumsum(dx)];
y_ce = [Y0;Y0 + cumsum(dy)];
z_ce = [Z0;Z0 - cumsum(dz)];
% 
% 
% 
Ox = sparse ( mcell , nfx );
Oy = sparse ( mcell , nfy );
Oz = sparse ( mcell , nfz );
% 
O = sparse(mcell,mcell);

%% Create inducing field matrix
hy = cosd(I)*cosd(D);
hx = cosd(I)*sind(D);
hz = sind(I);
% 
H0 = [ones(nfz,1) * hz ; ones(nfx,1) * hx ;ones(nfy,1) * hy] * H; 
% H0 = spdiags(h0',0,nfaces,nfaces);
% H0 = [H * (r);H * (q);H * (p)];
% 
% hx= ones(1,mcell) * H * (p);
% hy= ones(1,mcell) * H * (q);
% hz= ones(1,mcell) * H * (r);
% h1 = [hz hx hy];
% 
% Hx = spdiags(hx',0,mcell,mcell);
% Hy = spdiags(hy',0,mcell,mcell);
% Hz = spdiags(hz',0,mcell,mcell);
% 
% H1 = [Hz Hx Hy];

%% Averaging operator from face to center
av = @(n) spdiags (ones (n+2,1)*[0.5,0.5],[-1,0],n+2,n+1);

avfx = av(nx); avfx = avfx(2:end-1,:);
avfy = av(ny); avfy = avfy(2:end-1,:); 
avfz = av(nz); avfz = avfz(2:end-1,:);

Avfx = kron ( kron ( speye (ny), avfx), speye (nz) );
Avfy = kron ( kron ( avfy , speye (nx)), speye (nz));
Avfz = kron ( kron ( speye (ny) , speye (nx)), avfz );

% 
% AVC =   [Avfz Ox Oy
%         Oz Avfx Oy
%         Oz Ox Avfy];
AVC =   [Avfz*hz Avfx*hx Avfy*hy];
% 
% %% Averaging operator from center to face
% av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[0,1],n+1,n+2);
% 
% avcx = av(nx); avcx=avcx(:,2:end-1); avcx([1,end])=1;
% Avcx = kron(kron(speye(ny),avcx),speye(nz));
% 
% avcy = av(ny); avcy=avcy(:,2:end-1); avcy([1,end])=1;
% Avcy = kron(kron(avcy,speye(nx)),speye(nz));
% 
% avcz = av(nz); avcz=avcz(:,2:end-1); avcz([1,end])=1;
% Avcz = kron(kron(speye(ny),speye(nx)),avcz);
% 
% AVF = [Avcz; Avcx; Avcy];
%     



%% Create gradient op (center to face)
% ddxm = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);
% 
% d_dx = dXm * ddxm(nx); d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
% Dx = kron(kron(speye(ny),d_dx),speye(nz));
% 
% d_dy = dYm * ddxm(ny); d_dy=d_dy(:,2:end-1); d_dy([1,end])=0;
% Dy = kron(kron(d_dy,speye(nx)),speye(nz));
% 
% d_dz = dZm * ddxm(nz); d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
% Dz = kron(kron(speye(ny),speye(nx)),d_dz);
% 
% GRAD_ctf = [Dz; Dx; Dy];

%% Create gradient op (face to center)
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n+1,n+2);
d_dx=dX * ddx(nx-1);% d_dx(end,end)=0;
Dx = kron(kron(speye(ny),d_dx),speye(nz));

d_dy=dY * ddx(ny-1);% d_dy([1,end])=0;
Dy = kron(kron(d_dy,speye(nx)),speye(nz));

d_dz=dZ * ddx(nz-1);% d_dz([1,end])=0;
Dz = kron(kron(speye(ny),speye(nx)),d_dz);



% GRAD_ftc = [Dz Ox Ox;
%             Oz Dx Oy;
%             Oz Ox Dy];
        
%% Create divergence op (face to center)
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
% DIV = [Dz Dx Dy];
%             Oz Dx Oy;
%             Oz Ox Dy];
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

Rfz = spdiags(1./Rz,0,nfz,nfz);
Rfx = spdiags(1./Rx,0,nfx,nfx);
Rfy = spdiags(1./Ry,0,nfy,nfy);

% Rf = spdiags(1./Rf,0,nfaces,nfaces);

Rc =    ( (ObsX(1) - kron3( e(ny) , x_cc', e(nz) )).^2 +...
       (ObsY(1) - kron3( y_cc' , e(nx), e(nz) )).^2 +...
       (ObsZ(1) - kron3( e(ny) , e(nx), z_cc' )).^2 ) .^-0.5 ;
   
Rc =Rc(:);
%% Operator for DEL (DEL phi)

% Partial derivatives for x-component
% ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n,n+1);
% Averaging operator cell center to faces (for cell-center scalar)
av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[0,1],n,n+1);

% Z-component from center to face. Hits a column vector [Zvec;Xvec;Yvec]
dz_dz=  dZm * ddx(nz); dz_dz=dz_dz(:,2:end-1); dz_dz([1,end])=0;
Dzz = kron( kron( speye(ny) , speye(nx) ), dz_dz );

% dz_dx= dZm * ddx(nz); dz_dx=dz_dx(:,2:end-1); dz_dx([1,end])=0;
% Dzx = kron( kron( speye(ny) , dz_dx ),speye(nz) );
% avzx = av(nx);
% Avzx = kron ( kron ( speye (ny) , avzx), speye (nz+1));

% dz_dy= dZm * ddx(ny); dz_dy = dz_dy(:,2:end-1); dz_dy([1,end])=0;
% Dzy = kron( kron( dz_dy , speye(nx) ),speye(nz) );
% avzy = av(ny);
% Avzy = kron ( kron ( avzy , speye (nx) ), speye (nz+1));

% Average back from cell faces to center
% avzz = av(nz);
% Avzz = kron ( kron ( speye (ny) , speye (nx)), avzz);

% X-component
dx_dx=  dXm * ddx(nx); dx_dx=dx_dx(:,2:end-1); dx_dx([1,end])=0;
Dxx = kron( kron( speye(ny) , dx_dx ), speye(nz) );

% dx_dz= dZm * ddx(nz); dx_dz=dx_dz(:,2:end-1); dx_dz([1,end])=0;
% Dxz = kron( kron( speye(ny) , speye(nx)),dx_dz  );
% avxz = av(nz);
% Avxz = kron ( kron ( speye (ny) , speye (nx+1)), avxz);

% dx_dy= dYm * ddx(ny); dx_dy = dx_dy(:,2:end-1); dx_dy([1,end])=0;
% Dxy = kron( kron( dx_dy , speye(nx) ),speye(nz) );
% avxy = av(ny);
% Avxy = kron ( kron ( avxy , speye (nx+1) ), speye (nz));

% Average back from cell faces to center
% avxx = av(nx);
% Avxx = kron ( kron ( speye (ny) , avxx), speye (nz) );

% Y-component
dy_dy=  dYm * ddx(ny); dy_dy=dy_dy(:,2:end-1); dy_dy([1,end])=0;
Dyy = kron( kron( dy_dy  , speye(nx) ), speye(nz) );

% dy_dz= dZm * ddx(nz); dy_dz=dy_dz(:,2:end-1); dy_dz([1,end])=0;
% Dyz = kron( kron( speye(ny) , speye(nx)),dy_dz  );
% avyz = av(nz);
% Avyz = kron ( kron ( speye (ny+1) , speye (nx)), avyz);

% dy_dx= dXm * ddx(nx); dy_dx = dy_dx(:,2:end-1); dy_dx([1,end])=0;
% Dyx = kron( kron(  speye(ny), dy_dx ),speye(nz) );
% avyx = av(nx);
% Avyx = kron ( kron ( speye (ny+1) , avyx ), speye (nz));

% Average back from cell faces to center
% avyy = av(ny);
% Avyy = kron ( kron ( avyy , speye (nx) ), speye (nz) );



del = [Dzz*Dz*Rfz Dzz*Dx*Rfx Dzz*Dy*Rfy;
     Dxx*Dz*Rfz Dxx*Dx*Rfx Dxx*Dy*Rfy;
     Dyy*Dz*Rfz Dyy*Dx*Rfx Dyy*Dy*Rfy];
 
% del = [Dzz*Dz*Rc Avzz*Avzx*Dzx*Dz*Rc Avzz*Avzy*Dzy*Dz*Rc;
%       Avxx*Avxz*Dxz*Dx*Rc Dxx*Dx*Rc Avxx*Avxy*Dxy*Dx*Rc;
%      Avyy*Avyz*Dyz*Dy*Rc Avyy*Avyx*Dyx*Dy*Rc Dyy*Dy*Rc];
 
%  del = [Dzz Ox Oy;
%      Oz Dxx Oy;
%      Oz Ox Dyy];

%% Create cell dimension matrix

dX = kron3( e(ny) , dx' , e(nz) ); 
dY = kron3( dy' , e(nx), e(nz) ); 
dZ = kron3( e(ny) , e(nx) , dz' );

dV = [dX.*dY.*dZ];

V = spdiags(dV',0,mcell,mcell);

% VV = [V O O;
%       O V O;
%       O O V];
% VV = [dV dV dV];

% Volume descritization for faces 
% dxf = dx'/2;dyf = dy'/2;dzf = dz'/2;
% dXf = kron3( e(ny) , [dxf(1)*2;(dxf(1:end-1) + dxf(2:end));dxf(end)*2] , e(nz) );
% dYf = kron3( [dyf(1)*2;(dyf(1:end-1) + dyf(2:end));dyf(end)*2] , e(nx) , e(nz) );
% dZf = kron3( e(ny) , e(nx) , [dzf(1)*2;(dzf(1:end-1) + dzf(2:end));dzf(end)*2] );


%% TMI Matrix

s = speye(mcell);
TMI = kron([1 1 1],s);

%% Build depth weighting matrix
Wr = ( dV ./ (Rc' + R0) .^ 3 ) .^ 2 ;
Wr = Wr(:);
%% Build forward operator 

% Rc = spdiags(1./Rc,0,mcell,mcell);
% Rc2 = [Rc O O;
%       O Rc O;
%       O O Rc];
  
E = spdiags([ones(nfz,1)*-1;ones(nfx,1);ones(nfy,1)],0,nfaces,nfaces);
% FwrOp = -ones(1,mcell*3) /(4*pi) * VV * AVC * GRAD_ctf  * Rc * H1 * AVC *  GRAD_ctf;
% FwrOp = 1 /(4*pi) * VV * AVC * GRAD_ctf * Rc * DIV * H0 * AVF;
% FwrOp = 1 /(4*pi) * V * (H1 * del) * GRAD_ctf * Rc;
FwrOp = 1 /(4*pi) * V * AVC  * (del) * H0;
% FwrOp = -[ones(1,mcell) zeros(1,mcell) zeros(1,mcell)]  * VV * GRAD_ftc * Rf *  H0 * GRAD_ctf;
FwrOp = FwrOp';