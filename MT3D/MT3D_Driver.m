% EOSC 555 - 3DMT_Driver
% Forward modeling of a conductive body in a uniform half-space model
%
% Author:
% Dominique Fournier
% 
% Co-author:
% Ben Postlethwaite
% Last update: April 20th, 2013
%
% INTRO:
% Program for the forward simulation of 3D MT problem as described by Haber
%  et al. (2000). The algorithm uses the Helmholtz decomposition to express 
% the E-field as the combinaison of a vector potential (A) and scalar 
% potential (phi). See the project report for more details. 
% 
% INPUTS:
% The program takes for input the mesh and conductivity model as prescribed
% by the user. The program allows for both uniform and tensor meshes.
% From the input the propgram creates sparse operator for GRAD, DIV and 
% the Laplacien operators. Default values below solve for a 4x4x4 uniform
% mesh with 3 padding cells (expansion factor=1.2)
%
% OUPUTS:
% The program will procede with a series of test for each operators and
% display the result in the command terminal.
%
% The program also outputs the mesh and cell-center pointset in UBC
% format.
%
% Impedances are computed an displayed to screen.
%
% SUBFUNCTIONS:
%
% get_ops.m : Generates all the operators required for the forward
% modelling with the appropriate boundary conditions for this problem.
% (Note: This function should be modified to be more flexible)
%
% Laplacian_test.m : Create the same Laplacian as used in the code and
% procede with a numerical test for the second partial partial derivative
% on a simple function. 
%
% CURL_test.m : Generate the CURL operator with Derichlet boundary
% conditions and test it on a simple function with known analytical
% solution. 
%
% DIV_test.m : Same test as for the Laplacian but using a diffrent
% function.
%
% GRAD_test.m : Same test as for the DIV_test but using a diffrent
% function.
% 
% Fwr_test.m : Procede with a test on the forward operator as prescribed by
% Haber (2000). Uses a conductivity model and E-field vanishing at the
% boundaries and computes a speudo-analytical solution for A, phi and a 
% source term J. A numerical solution is found using the inverse of the 
% forward operator.
%
% Derivative_test.m : Procede with the standard derivative test.
%
% make_UBC_mesh.m
% make_UBC_model.m : Generate UBC file formats for visualization.

clear all
close all


addpath functions

%% User Input - Build mesh and set model parameters

XYZo = [0 0 0];             % Origin SW top corner (x ,y ,z)
XYZmax = [300 300 300];     % Model limits NE bottom corner (x, y ,z)

nx= 4 ; ny=4 ; nz=4;  % Number of core cells before padding

% Uniform mesh
% dx = ones(nx,1)*(XYZmax(1) - XYZo(1))/nx;
% dy = ones(ny,1)*(XYZmax(2) - XYZo(2))/ny;
% dz = ones(nz,1)*(XYZmax(3) - XYZo(2))/nz;

% Tensor mesh
expan_fact = 1.4; % Padding expansion factor
dx = (XYZmax(1) - XYZo(1))/nx; dx = [dx*expan_fact.^[3:-1:1]';ones(nx,1)*dx;dx*expan_fact.^[1:3]'];
dy = (XYZmax(2) - XYZo(2))/ny; dy = [dy*expan_fact.^[3:-1:1]';ones(ny,1)*dy;dy*expan_fact.^[1:3]'];
dz = (XYZmax(3) - XYZo(2))/nz; dz = [dz*expan_fact;ones(nz,1)*dz;dz*expan_fact.^[1:3]'];

make_UBC_mesh(XYZo,dx,dy,dz); % Create mesh file for UBC viewer

nx=length(dx) ; ny=length(dy) ; nz=length(dz);  % Number of cells

uo = 4*pi*10^-7;            % Magnetic permeability of free-space
s_air = 1e-15;              % Air conductivity  
s_background = 1e-0;        % Background conductivity            
s_anomaly = 1e+3;           % Anomalous conductivity
w = 10^3.0;                 % Frequency [hz]

mcell = nx * ny * nz;
nface = (nx+1) * ny * nz + nx * (ny+1) * nz + nx * ny * (nz+1);

skindepth = sqrt( 2/ (s_background * uo * w));

fprintf('Input background conductivity: %7.3e\n',s_background)
fprintf('Input anomalous conductivity: %7.3e\n\n',s_anomaly)

%% Model parameters
% Build conductivity model for primary field calculations
s_p = ones(nx,ny,nz)*s_background;
s_p(:,:,1:2) = s_air;

s_s = s_p; % Anomaly conductivity [S/m]


% % Build anomalous 3D model (2x2x1 prism at top center)

for ii=ceil(nx/2):ceil(nx/2)+1
    for jj=ceil(ny/2):ceil(ny/2)+1
        for kk=3
            s_s(ii,jj,kk)= s_anomaly;
        end
    end
end

% % Build anomalous 2D model (2x2x1 prism at top center)
% 
% for ii=2:nx/2
%     for jj=2:ny-1
%         for kk=3:nz
%             s_s(ii,jj,kk)= s_anomaly;
%         end
%     end
% end

s_p=reshape(s_p,mcell,1);

s_s=reshape(s_s,mcell,1);

UBC_model = make_UBC_model(s_s,nx,ny,nz);   % Create model in UBC format
save('3DMT_UBC.con','-ascii','UBC_model');

%% Get LAPLACIAN, DIVERGENCE, GRADIENT and AVERAGING operators for the
% primary and secondary field calculations.

[L_px,L_py,L_s,DIV_p,DIV_s,GRAD,CURL_e,CURL_f,AVFsc,AVFvcz,AVC,AVHx,AVHy,dZm]=get_ops(dx,dy,dz,w,s_background);


%% Create functions

% Harmonic averaging
S = @(m) spdiags( (AVFsc * (m.^-1)).^-1, 0,nface , nface) ;

% Build forward operator
G_p = @(m,L) [L + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD
            DIV_p*S(m)      , DIV_p*S(m)*GRAD ];

G_s = @(m) [L_s + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD
            DIV_s*S(m)      , DIV_s*S(m)*GRAD ];
        
C_p = @(m,u,q,L) G_p(m,L)*u - q;

sdiag = @(XX) spdiags( XX, 0, length(XX), length(XX));

dSdm = @(m) S(m).^2 * AVFsc * sdiag(1./m.^2);

% Derivatives
dCdm_p = @(a,phi,m) [1i*w*uo*sdiag(a) * dSdm(m) + 1i*w*uo * sdiag(GRAD * phi) * dSdm(m)
                       DIV_p * sdiag(a) * dSdm(m) + DIV_p * sdiag(GRAD * phi) * dSdm(m)];

% Sensitivity Matrix
J_p = @(m,a,phi) -G_p(m)\dCdm_p(a,phi,m);




%% Perform tests on operators

% Laplacian_test('tensor');
% 
% CURL_test;
% 
% DIV_test;
% 
% GRAD_test;
% 
% Derivative_test(C_p,dCdm_p,L_px,nface,mcell);

Fwr_test;

%% Construct Boundary Condition Matrix for the Laplacian
% Boundary matrix for x-polarization
bce = [-2*1, 2*0]; 
bcx = zeros((nx + 1) * ny * nz, 1); 
bcx(1:(nx + 1) * ny) = bce(1)*dZm(1)^2;

bcy = zeros(nx  * (ny+ 1) * nz, 1);

bcz = zeros(nx  * ny * (nz+ 1), 1);

Bcx = sparse(bcx);

Bcy = sparse(bcy);

Bcz = sparse(bcz);


BCx = [Bcx ; Bcy; Bcz];

% Boundary matrix for y-polarization
bcx = zeros((nx + 1) * ny * nz, 1); 

bcy = zeros(nx  * (ny+ 1) * nz, 1);
bcy(1:nx * ( ny + 1)) = bce(1)*dZm(1)^2;

bcz = zeros(nx  * ny * (nz+ 1), 1);

Bcx = sparse(bcx);

Bcy = sparse(bcy);

Bcz = sparse(bcz);


BCy = [Bcx ; Bcy; Bcz];


q_px = [ -BCx; zeros(mcell, 1)];
q_py = [ -BCy; zeros(mcell, 1)];

%% Forward modeling
% Build forward operator for background conductivity
Gpx = G_p(s_p,L_px);
Gpy = G_p(s_p,L_py);

% Compute vector field A and scalar field Phi: u = [A;phi]
u_px = Gpx\q_px;
u_py = Gpy\q_py;

% Extract potentials for two polarization
Ax = u_px(1: nface);
Phix = u_px( nface + 1: end);

Ay = u_py(1: nface);
Phiy = u_py( nface + 1: end);

% Compute E-field = (A + GRAD*phi)
E_px = Ax + GRAD*Phix;
E_py = Ay + GRAD*Phiy;

% Plot primary field
% figure;vecplot(E_p,AVC, dx, dy, dz)
% title('Primary Field')

S_p=S(s_p);

% Use primary field 
q_sx = @(m) [-1i * w * uo * (S(m)-S_p) * E_px; -DIV_p * (S(m)-S_p) * E_px];

% Compute secondary fields
u_sx = G_s(s_s)\q_sx(s_s);

Ax = u_sx(1: nface);
Phix = u_sx( nface + 1: end);

E_sx = Ax + GRAD*Phix;

% figure;vecplot(E_s,AVC, dx, dy, dz)
% title('Secondary Field')

E_tx = E_sx + E_px;
% figure;vecplot(E_t, dx, dy, dz)
% title('Total Field')

J_tx = S(s_s) * E_tx;
figure;vecplot(J_tx,AVC, dx, dy, dz)
title(['\bfCurrent vectors (x-polarized) for \omega: 10^' num2str(log10(w)) 'Hz'])

% Use primary field y-polarized
q_sy = @(m) [-1i * w * uo * (S(m)-S_p) * E_py; -DIV_p * (S(m)-S_p) * E_py];

% Compute secondary fields
u_sy = G_s(s_s)\q_sy(s_s);

Ay = u_sy(1: nface);
Phiy = u_sy( nface + 1: end);

E_sy = Ay + GRAD*Phiy;

% figure;vecplot(E_s,AVC, dx, dy, dz)
% title('Secondary Field')

E_ty = E_sy + E_py;
% figure;vecplot(E_t, dx, dy, dz)
% title('Total Field')

J_ty = S(s_s) * E_ty;
figure;vecplot(J_ty,AVC, dx, dy, dz)
title(['\bfCurrent vectors (y-polarized) for \omega: 10^' num2str(log10(w)) 'Hz'])


%% Merge primary and secondary fields
u_tx = u_px + u_sx;

u_ty = u_py + u_sy;

E_txy = E_ty + E_tx;
E_sxy = E_sx + E_sy;
E_pxy = E_px + E_py;

J_txy = J_ty + J_tx;
figure;vecplot(J_txy,AVC, dx, dy, dz)
title(['\bfCurrent vectors (total) for \omega: 10^' num2str(log10(w)) 'Hz'])

% Create vector files for Gocad import
make_vector_file(AVC*J_txy,XYZo,dx,dy,dz,['J_total_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_txy,XYZo,dx,dy,dz,['E_total_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_sxy,XYZo,dx,dy,dz,['E_secon_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_px,XYZo,dx,dy,dz,['E_px_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_py,XYZo,dx,dy,dz,['E_py_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_sx,XYZo,dx,dy,dz,['E_sx_w' num2str(w) '.dat'])
% make_vector_file(AVC*E_sy,XYZo,dx,dy,dz,['E_sy_w' num2str(w) '.dat'])

%% Compute Impedances
% Compute H-field by selecting A, curling, averging to face and selecting
% observation points only (cell centers for now)

SA = [speye( nface ), sparse( nface, mcell)]; % Select A

% Create matrix selector for the top cells in center (should be modified to
% take Obs locations)
obs_loc = zeros(nx,1); obs_loc(4:7) = 1;
Obsx = kron([0;0;1;zeros(nz-2,1)],kron(speye (ny),obs_loc));
SObsH = [Obsx' sparse(length(obs_loc),nz*nx*(ny+1)) sparse(length(obs_loc),nz*ny*(nx+1))];
SObsEx = [Obsx' sparse(length(obs_loc),(nz+1)*nx*(ny)) sparse(length(obs_loc),(nz+1)*ny*(nx))];
SObsEy = [sparse(length(obs_loc),(nz+1)*nx*(ny)) Obsx' sparse(length(obs_loc),(nz+1)*ny*(nx))];

Qhx = SObsH * AVHx * CURL_f * SA / (1i*w*uo);
Qhy = SObsH * AVHy * CURL_f * SA / (1i*w*uo);

Qex = SObsEx * AVFvcz * AVC * SA;
Qey = SObsEy * AVFvcz * AVC * SA;

Zxx = @(ux, uy) ( (Qex*ux) .* (Qhy*uy) - (Qex*uy) .* (Qhx*ux) ) ./...
                ( (Qhy*uy) .* (Qhy*ux) - (Qhy*ux) .* (Qhx*uy) );
            
Zxy = @(ux, uy) ( (Qex*uy) .* (Qhx*ux) - (Qex*ux) .* (Qhx*uy) ) ./...
                ( (Qhy*uy) .* (Qhy*ux) - (Qhy*ux) .* (Qhx*uy) ); 
            
Zyx = @(ux, uy) ( (Qey*uy) .* (Qhy*ux) - (Qey*ux) .* (Qhy*uy) ) ./...
                ( (Qhx*uy) .* (Qhy*ux) - (Qhx*ux) .* (Qhy*uy) );  


Zyy = @(ux, uy) ( (Qey*ux) .* (Qhx*uy) - (Qex*uy) .* (Qhx*ux) ) ./...
                ( (Qhx*uy) .* (Qhy*ux) - (Qhx*ux) .* (Qhy*uy) );  

Impxx = Zxx(u_tx,u_ty);
Impxy = Zxy(u_tx,u_ty);
Impyx = Zyx(u_tx,u_ty);
Impyy = Zyy(u_tx,u_ty);

app_res = 1/(w * uo) * (real(Impxy).^2 + imag(Impxy).^2);

fprintf('Computed apparent conductivity from impedance')
app_con = 1./app_res

residual = s_background - app_con ;

fprintf('*** END OF PROGRAM ***\n')