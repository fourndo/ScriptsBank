% EOSC 555 - Operators
% 
% Authors
% Operator: Eldad
% Test : Dom
% 
% Date: March 9th, 2013
%
% INTRO:
% Program for the simulation of 3D MT problem as described by Haber et al.
% 2000. The algorithm uses the Helmholtz decomposition to express the
% E-field as the combinaison of a vector potential (A) and scalar potential
% (phi). 
% The program takes for input the mesh and conductivity model as prescribed
% by the user. The program allows for both uniform and tensor meshes.
% From the input the propgram creates sparse operator for GRAD, DIV and 
% the Laplacien operators. Numerical experiments show that a precision of 
% O(h^2) is achieved on the uniform mesh, but only O(h) on the tensor mesh.
%
% A derivative test is also computed showing O(h^2) precision.

clear all
close all


addpath functions
for ii=2:0.25:3
%% Build mesh and set model parameters
% User Input
XYZo = [0 0 0];             % Origin SW top corner (x ,y ,z)
XYZmax = [300 300 300];     % Model limits NE bottom corner (x, y ,z)

nx= 4 ; ny=4 ; nz=5;  % Number of cells before padding

% Uniform mesh
% dx = ones(nx,1)*(XYZmax(1) - XYZo(1))/nx;
% dy = ones(ny,1)*(XYZmax(2) - XYZo(2))/ny;
% dz = ones(nz,1)*(XYZmax(3) - XYZo(2))/nz;

% Tensor mesh
dx = (XYZmax(1) - XYZo(1))/nx; dx = [dx*1.2.^[3:-1:1]';ones(nx,1)*dx;dx*1.2.^[1:3]'];
dy = (XYZmax(2) - XYZo(2))/ny; dy = [dy*1.2.^[3:-1:1]';ones(ny,1)*dy;dy*1.2.^[1:3]'];
dz = (XYZmax(3) - XYZo(2))/nz; dz = [ones(nz,1)*dz;dz*1.4.^[1:3]'];

make_UBC_mesh(XYZo,dx,dy,dz);

nx=length(dx) ; ny=length(dy) ; nz=length(dz);  % Number of cells


uo = 4*pi*10^-7;            % Magnetic permeability of free-space
s_air = 1e-15;
s_background = 1e-0;        % Background conductivity            
s_anomaly = 1e+3;           % Anomalous conductivity
w = 10^ii;                   % Frequency [hz]

skindepth = sqrt( 2/ (s_background * uo * w));

%% Get LAPLACIAN, DIVERGENCE, GRADIENT and AVERAGING operators for the
% primary and secondary field calculations.

[L_p,L_s,DIV_p,DIV_s,GRAD_p,GRAD_s,CURL_e,CURL_f,AVF,AVC,AVHx,AVHy,dZm]=get_ops(dx,dy,dz,w,s_background);

%% Perform test on Laplacian for different meshes and output result
% 'help Laplacian_test' for more details

% Laplacian_test('uniform')
% Laplacian_test('tensor')
%%
mcell = size(GRAD_p,2) ;
nface = size(DIV_s,2);

% Harmonic averaging
S = @(m) spdiags( (AVF * (m.^-1)).^-1, 0,nface , nface) ;

% Build forward operator
G_p = @(m) [L_p + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD_p
            DIV_p*S(m)      , DIV_p*S(m)*GRAD_p ];

G_s = @(m) [L_s + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD_s
            DIV_s*S(m)      , DIV_s*S(m)*GRAD_s ];
        
C_p = @(m,u,q) G_p(m)*u - q;

sdiag = @(XX) spdiags( XX, 0, length(XX), length(XX));

dSdm = @(m) S(m).^2 * AVF * sdiag(1./m.^2);

% Derivatives


dCdm_p = @(a,phi,m) [1i*w*uo*sdiag(a) * dSdm(m) + 1i*w*uo * sdiag(GRAD_p * phi) * dSdm(m)
                       DIV_p * sdiag(a) * dSdm(m) + DIV_p * sdiag(GRAD_p * phi) * dSdm(m)];

% Sensitivity Matrix
J_p = @(m,a,phi) -G_p(m)\dCdm_p(a,phi,m);

% %% Derivative test for primary field operator
% a = randn(nface, 1);
% phi = randn(mcell, 1);
% u = [a;phi];
% m = randn(mcell, 1);
% q = randn(mcell+nface, 1);
% v = randn(mcell, 1);
% 
% fprintf('Derivative Test')
% fprintf('\t1/h\t\t|f(x)-f(x+h)|\t|f(x+h) - f(x)- dfdx*h|\n')
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
% fprintf('%3.2e    %3.2e \t\t\t %3.2e\n',h,Oh,Oh2)
% end
% 
% fprintf ('\n')


%% Forward modeling of uniform half-space conductivity and unit cube anomaly

%% Model parameters

% Build background model with 2 cells of air

% s_p = kron ([ones(2,1)*s_air;ones(nz-2,1)*s_background],kron(ones(ny,1),ones(nx,1))) ; % Background conductivity [S/m]
s_p = ones(mcell,1)*s_background;

s_s = s_p; % Anomaly conductivity [S/m]

%Build anomalous model
s_s=reshape(s_s,nx,ny,nz);
for ii=ceil(nx/2):ceil(nx/2)+1
    for jj=ceil(ny/2):ceil(ny/2)+1
        for kk=1
            s_s(ii,jj,kk)= s_anomaly;
        end
    end
end

s_s=reshape(s_s,mcell,1);

UBC_model = make_UBC_model(s_s,nx,ny,nz);
save('3DMT_UBC.con','-ascii','UBC_model');




%% Construct Boundary Condition Matrix for Robin boundary conditions
% Need 
bcs = [-2*1, 2*0]; 
bcx = zeros((nx + 1) * ny * nz, 1); 
bcx(1:(nx + 1) * ny) = bcs(1)*dZm(1)^2;


bcy = zeros(nx  * (ny+ 1) * nz, 1);
bcy(1:nx * ( ny + 1)) = bcs(1)*dZm(1)^2;


bcz = zeros(nx  * ny * (nz+ 1), 1);% bcz([1,end]) = bcs;

BCx = sparse(bcx);

BCy = sparse(bcy);

BCz = sparse(bcz);


BC = [BCx ; BCy; BCz];

q_p = [ -BC; zeros(mcell, 1)];

% Build forward operator for background conductivity
Gp = G_p(s_p);

% Compute vector field A and scalar field Phi: u = [A;phi]
u_p = Gp\q_p;
A = u_p(1: nface);
Phi = u_p( nface + 1: end);

% Compute E-field = (A + GRAD*phi)
E_p = A + GRAD_p*Phi;

% Plot primary field
% figure;vecplot(E_p,AVC, dx, dy, dz)
% title('Primary Field')

S_p=S(s_p);

% Use primary field 
q_s = @(m) [-1i * w * uo * (S(m)-S_p) * E_p; -DIV_p * (S(m)-S_p) * E_p];

% Compute secondary fields
u_s = G_s(s_s)\q_s(s_s);

A = u_s(1: nface);
Phi = u_s( nface + 1: end);

E_s = A + GRAD_s*Phi;

% figure;vecplot(E_s,AVC, dx, dy, dz)
% title('Secondary Field')

E_t = E_s + E_p;
% figure;vecplot(E_t, dx, dy, dz)
% title('Total Field')

J_t = S(s_s) * E_t;
figure;vecplot(J_t,AVC, dx, dy, dz)
title(['\bfCurrent vectors for \omega 10^' num2str(log10(w)) 'Hz'])

make_vector_file(AVC*J_t,XYZo,dx,dy,dz)
end
% %% Inverse problem
% Jacob=J_p(s_s,A,Phi);
% % Compute H-field
% H = (CURL*E_t - CURL*(BC/dZm(1)^2))/(-1i * w * uo);
% 
% % Extract Bx and By at surface (cell center)
% Hx = AVHx * H; Hx=Hx(1:nx*ny);
% Hy = AVHy * H; Hy=Hy(1:nx*ny);
% 
% % % Right now the surface Ex and Ey = BC
% % Avc_E_t = AVC * E_t;
% % Ex = Avc_E_t(1:nx*ny);
% % Ey = Avc_E_t((nx*ny*nz+1):nx*ny*nz*2);
% data=[];
% 
% % Need to build Q
% 
% % Compute impedance as data points
% for ii=1:nx*ny
%     
%     Zxy = Bcx(ii)/Hy(ii);
%     Zyx = Bcy(ii)/Hx(ii);
%     Zxx = Bcx(ii)/Hx(ii);
%     Zyy = Bcy(ii)/Hy(ii);
% %     rhoxy(ii) = abs(Zxy)^2/ ( w * uo);
% %     rhoyx(ii) = abs(Zyx)^2/ ( w * uo);
% %     
% %     phixy(ii) = atan (imag(Zxy)/real(Zxy));
% %     phiyx(ii) = atan (imag(Zyx)/real(Zyx));
% end