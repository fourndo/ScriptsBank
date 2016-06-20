function [A, G, V ]=get_GRAD_op3D_TENSIL_Kron(dx,dy,dz,nullcell,FLAG)
%Build the derivative weighting matrices on a 17 points stencil
%
%      12-15
%     /13-16
%    / 14-17
%   7-9
%  /x-10
% / 8-11
%1-4
%2-5
%3-6



%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

% Compute half-cell distance 
dxm =  dx(1:end-1)/2 + dx(2:end)/2 ; dxm = [dxm(:);dxm(end)];
dym =  dy(1:end-1)/2 + dy(2:end)/2 ; dym = [dym(:);dym(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2 ; dzm = [dzm(:);dzm(end)];

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
dz = sqrt(dz(:));

mcell = nx*ny*nz;

T=kron(kron(spdiags(ones(ny,1)*[-10,1,10],[-1,0,1],ny,ny),spdiags(ones(nx,1)*[1,2,3],[-1,0,1],nx,nx)),spdiags(ones(nz,1)*[3,5,7],[-1,0,1],nz,nz));

ind = find(T);
val = nonzeros(T);

% Gradient g-yz
g1 = T;
g1( ind( val ~=-140 & val ~= 10 ) ) = 0;
g1( ind( val ==-140 )) = 1;
g1(ind( val ==10 ) ) = -1;

% Remove nullcells
g1 = g1(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g1),2) ~= 2;
g1(indx,:) = 0;

% Gradient g-y
g2 = T;
g2( ind( val ~=-100 & val ~= 10 ) ) = 0;
g2( ind( val ==-100 )) = 1;
g2(ind( val ==10 ) ) = -1;

% Remove nullcells
g2 = g2(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g2),2) ~= 2;
g2(indx,:) = 0;

% Gradient g-y-z
g3 = T;
g3( ind( val ~=-60 & val ~= 10 ) ) = 0;
g3( ind( val ==-60 )) = 1;
g3(ind( val ==10 ) ) = -1;

% Remove nullcells
g3 = g3(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g3),2) ~= 2;
g3(indx,:) = 0;

% Gradient g-y-z
g4 = T;
g4( ind( val ~=-70 & val ~= 10 ) ) = 0;
g4( ind( val ==-70 )) = 1;
g4(ind( val ==10 ) ) = -1;

% Remove nullcells
g4 = g4(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g4),2) ~= 2;
g4(indx,:) = 0;

% Gradient g-y-z
g5 = T;
g5( ind( val ~=-50 & val ~= 10 ) ) = 0;
g5( ind( val ==-50 )) = 1;
g5(ind( val ==10 ) ) = -1;

% Remove nullcells
g5 = g5(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g5),2) ~= 2;
g5(indx,:) = 0;

% Gradient g-y-z
g6 = T;
g6( ind( val ~=-30 & val ~= 10 ) ) = 0;
g6( ind( val ==-30 )) = 1;
g6(ind( val ==10 ) ) = -1;

% Remove nullcells
g6 = g6(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g6),2) ~= 2;
g6(indx,:) = 0;

% Gradient g-y-z
g7 = T;
g7( ind( val ~=14 & val ~= 10 ) ) = 0;
g7( ind( val ==14 )) = 1;
g7(ind( val ==10 ) ) = -1;

% Remove nullcells
g7 = g7(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g7),2) ~= 2;
g7(indx,:) = 0;


% Gradient g-y-z
g8 = T;
g8( ind( val ~=6 & val ~= 10 ) ) = 0;
g8( ind( val ==6 )) = 1;
g8(ind( val ==10 ) ) = -1;

% Remove nullcells
g8 = g8(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g8),2) ~= 2;
g8(indx,:) = 0;

% Gradient g-y-z
g9 = T;
g9( ind( val ~=7 & val ~= 10 ) ) = 0;
g9( ind( val ==7 )) = 1;
g9(ind( val ==10 ) ) = -1;

% Remove nullcells
g9 = g9(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g9),2) ~= 2;
g9(indx,:) = 0;

% Gradient g-y-z
g10 = T;
g10( ind( val ~=5 & val ~= 10 ) ) = 0;
g10( ind( val ==5 )) = 1;
g10(ind( val ==10 ) ) = -1;

% Remove nullcells
g10 = g10(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g10),2) ~= 2;
g10(indx,:) = 0;

% Gradient g-y-z
g11 = T;
g11( ind( val ~=3 & val ~= 10 ) ) = 0;
g11( ind( val ==3 )) = 1;
g11(ind( val ==10 ) ) = -1;

% Remove nullcells
g11 = g11(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g11),2) ~= 2;
g11(indx,:) = 0;

% Gradient g-y-z
g12 = T;
g12( ind( val ~=140 & val ~= 10 ) ) = 0;
g12( ind( val ==140 )) = 1;
g12(ind( val ==10 ) ) = -1;

% Remove nullcells
g12 = g12(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g12),2) ~= 2;
g12(indx,:) = 0;

% Gradient g-y-z
g13 = T;
g13( ind( val ~=100 & val ~= 10 ) ) = 0;
g13( ind( val ==100 )) = 1;
g13(ind( val ==10 ) ) = -1;

% Remove nullcells
g13 = g13(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g13),2) ~= 2;
g13(indx,:) = 0;

% Gradient g-y-z
g14 = T;
g14( ind( val ~=60 & val ~= 10 ) ) = 0;
g14( ind( val ==60 )) = 1;
g14(ind( val ==10 ) ) = -1;

% Remove nullcells
g14 = g14(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g14),2) ~= 2;
g14(indx,:) = 0;

% Gradient g-y-z
g15 = T;
g15( ind( val ~=70 & val ~= 10 ) ) = 0;
g15( ind( val ==70 )) = 1;
g15(ind( val ==10 ) ) = -1;

% Remove nullcells
g15 = g15(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g15),2) ~= 2;
g15(indx,:) = 0;

% Gradient g-y-z
g16 = T;
g16( ind( val ~=50 & val ~= 10 ) ) = 0;
g16( ind( val ==50 )) = 1;
g16(ind( val ==10 ) ) = -1;

% Remove nullcells
g16 = g16(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g16),2) ~= 2;
g16(indx,:) = 0;

% Gradient g-y-z
g17 = T;
g17( ind( val ~=30 & val ~= 10 ) ) = 0;
g17( ind( val ==30 )) = 1;
g17(ind( val ==10 ) ) = -1;

% Remove nullcells
g17 = g17(nullcell==1,nullcell==1);
% Find row to remove
indx = sum(abs(g17),2) ~= 2;
g17(indx,:) = 0;

G{1} = g1;G{2} = g2;G{3} = g3;G{4} = g4;G{5} = g5;G{6} = g6;
G{7} = g7;G{8} = g8;G{9} = g9;G{10} = g10;G{11} = g11;G{12} = g12;
G{13} = g13;G{14} = g14;G{15} = g15;G{16} = g16;G{17} = g17;

%% Create hmid dimensions matrix
dX = kron( kron( ones(ny,1) , dx(:) ), ones(nz,1) );
dY = kron( kron( dy(:)  , ones(nx,1) ), ones(nz,1));
dZ = kron( kron( ones(ny,1) , ones(nx,1) ), dz(:) );

v = dX .* dY .* dZ;

V = spdiags( v , 0 , mcell,mcell);

V = V(nullcell==1,nullcell==1);
% Gx = X * Gx * X';
% Gy = X * Gy * X';
% Gz = X * Gz * X';
% 
% Vx = X * Vx * X';
% Vy = X * Vy * X';
% Vz = X * Vz * X';

%% Define spherical address for each gradients
A = [0 cosd(45) cosd(45);
    0 1 0;
    0 cosd(45) -cosd(45);
    cosd(60) cosd(60) cosd(45);
    cosd(45) cosd(45) 0;
    cosd(60) cosd(60) -cosd(45);
    0 0 1;
    0 0 -1;
    cosd(45) 0 cosd(45);
    1 0 0;
    cosd(45) 0 -cosd(45);
    0 -cosd(45) cosd(45);
    0 -1 0;
    0 -cosd(45) -cosd(45);
    cosd(60) -cosd(60) cosd(45);
    cosd(45) -cosd(45) 0;
    cosd(60) -cosd(60) -cosd(45)];

end