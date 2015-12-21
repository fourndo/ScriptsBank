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

%% Build mesh and operators
% User Input
XYZo = [0 0 0];             % Origin SW top corner (x ,y ,z)
XYZmax = [300 300 300];     % Model limits NE bottom corner (x, y ,z)

nx= 6 ; ny=7 ; nz=8;  % Number of cells before padding

% Uniform mesh
% dx = ones(nx,1)*(XYZmax(1) - XYZo(1))/nx;
% dy = ones(ny,1)*(XYZmax(2) - XYZo(2))/ny;
% dz = ones(nz,1)*(XYZmax(3) - XYZo(2))/nz;

% Tensor mesh
dx = (XYZmax(1) - XYZo(1))/nx; dx = [dx*1.2.^[3:-1:1]';ones(nx,1)*dx;dx*1.2.^[1:3]'];
dy = (XYZmax(2) - XYZo(2))/ny; dy = [dy*1.2.^[3:-1:1]';ones(ny,1)*dy;dy*1.2.^[1:3]'];
dz = (XYZmax(3) - XYZo(2))/nz; dz = [ones(nz,1)*dz;dz*1.4.^[1:3]'];

nx=length(dx) ; ny=length(dy) ; nz=length(dz);  % Number of cells


uo = 4*pi*10^-7;            % Magnetic permeability of free-space
s_air = 1e-15;
s_background = 1e-0;        % Background conductivity            
s_anomaly = 1e+3;           % Anomalous conductivity
w = 1e-1 * 10^2;                   % Frequency [hz]

skindepth = sqrt( 2/ (s_background * uo * w));

% Get LAPLACIAN, DIVERGENCE, GRADIENT and AVERAGING operators for the
% primary and secondary field calculations.
[L_p,L_s,DIV_p,DIV_s,GRAD,AVF,AVC,dZm]=get_ops(dx,dy,dz,w,s_background);

mcell = size(GRAD,2) ;
nface = size(DIV_s,2);

% Harmonic averaging
S = @(m) spdiags( (AVF * (m.^-1)).^-1, 0,nface , nface) ;

% Build forward operator
G_p = @(m) [L_p + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD
            DIV_p*S(m)      , DIV_p*S(m)*GRAD ];

G_s = @(m) [L_s + 1i*w*uo*S(m), 1i*w*uo*S(m)*GRAD
            DIV_s*S(m)      , DIV_s*S(m)*GRAD ];
        
% %% Derivative test for primary field operator
% a = randn(nface, 1);
% phi = randn(mcell, 1);
% u = [a;phi];
% m = randn(mcell, 1);
% q = randn(mcell+nface, 1);
% v = randn(mcell, 1);
% 
% sdiag = @(XX) spdiags( XX, 0, length(XX), length(XX));
% 
% % GRAD(1) = 2/dx(1);
% % Build forward operator
% 
% C_p = @(m,u,q) G_p(m)*u - q;
% 
%     dSdm = @(m) S(m).^2 * AVF * sdiag(1./m.^2);
%     % Derivatives
%     dCdm_p = @(a,phi,m) [1i*w*uo*sdiag(a) * dSdm(m) + 1i*w*uo * sdiag(GRAD * phi) * dSdm(m)
%                        DIV_p * sdiag(a) * dSdm(m) + DIV_p * sdiag(GRAD * phi) * dSdm(m)];
% 
% % Sensitivity Matrix
% % J_p = @(m,a,phi) -G_p(m)\dCdm_p(a,phi,m);
% % Jacob=J_p(m,a,phi);
% 
% fprintf('\t1/h\t\t|f(x)-f(x+h)|\t|f(x+h) - f(x)- dfdx*h|\n')
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



%% Forward modeling of uniform half-space conductivity and unit cube anomaly

%% Model parameters

% Build background model with 2 cells of air

% s_p = kron ([ones(2,1)*s_air;ones(nz-2,1)*s_background],kron(ones(ny,1),ones(nx,1))) ; % Background conductivity [S/m]
s_p = ones(mcell,1)*s_background;

s_s = s_p; % Anomaly conductivity [S/m]

%Build anomalous model
s_s=reshape(s_s,nx,ny,nz);
for ii=ceil(nx/2)-1:ceil(nx/2)+1
    for jj=ceil(ny/2)-1:ceil(ny/2)+1
        for kk=2
            s_s(ii,jj,kk)= s_anomaly;
        end
    end
end

s_s=reshape(s_s,mcell,1);


%% Construct Boundary Condition Matrix for Robin boundary conditions
% Need 
bcs = [-2*1, 2*0]; 
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

q_p = [ -B; zeros(mcell, 1)];

% Build forward operator for background conductivity
Gp = G_p(s_p);

% Compute vector field A and scalar field Phi: u = [A;phi]
u_p = Gp\q_p;
A = u_p(1: nface);
Phi = u_p( nface + 1: end);

% Compute E-field = (A + GRAD*phi)
E_p = A + GRAD*Phi;

% Plot primary field
figure;vecplot(E_p, dx, dy, dz)
title('Primary Field')

S_p=S(s_p);

% Use primary field 
q_s = @(m) [-1i * w * uo * (S(m)-S_p) * E_p; -DIV_p * (S(m)-S_p) * E_p];

% Compute secondary fields
u_s = G_p(s_s)\q_s(s_s);


A = u_s(1: nface);
Phi = u_s( nface + 1: end);

E_s = A + GRAD*Phi;

% figure;vecplot(E_s, dx, dy, dz)
% title('Secondary Field')

E_t = E_s + E_p;
% figure;vecplot(E_t, dx, dy, dz)
% title('Total Field')

J = S(s_s) * E_t;
figure;vecplot(J, dx, dy, dz)
title(['\bfCurrent vectors for \omega log' num2str(log(w))])


