% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

addpath C:\Users\dominiquef\Dropbox\MASTER\INVMAG3D\data\Block
addpath C:\Users\dominiquef\Dropbox\MASTER\INVMAG3D\functions
addpath C:\Users\dominiquef\Dropbox\Master\Miscellaneous
work_dir = 'C:\Users\dominiquef\Dropbox\MASTER\INVMAG3D\data';
Obsfile = 'C:\Users\dominiquef\Dropbox\MASTER\INVMAG3D\data\Block\UBC_FWR.obs';
u0 = 4 * pi * 10^-7;

% [meshfile]=get_UBCmesh('Mesh_20m.msh');
% [ObsX, ObsY, ObsZ, data, Wd, H, I, D] = read_UBCobs('Obs_block.dat');
% D = 45;
% I = 0;
% dx = meshfile(3,meshfile(3,:)~=0)';
% dy = meshfile(4,meshfile(4,:)~=0)';
% dz = meshfile(5,meshfile(5,:)~=0)';
% 
% nx = meshfile(1,1); %size(X,1);    %number of cell in X
% ny = meshfile(1,2); %size(X,2);    %number of cell in Y
% nz = meshfile(1,3); %size(X,3);    %number of cell in Z
% 
% X0 = meshfile(2,1);
% Y0 = meshfile(2,2);
% Z0 = meshfile(2,3);
% m = load('C:\Users\dominiquef\Dropbox\DOM_Projects\INVMAG3D\data\Block\UBC_refmodal\Block_center.sus');
% m(m==-100)=0;
%  m(m~=0) = (1 + m(m~=0)) * u0;
% m(m==0)=u0;


%% 3D density contrast model

% Size of core
core = 100;
r = 10;

% Pre-allocate space for results of different solvers
% #iter, cpu time, residual, data error
ntest = 5;
spec{1} = zeros(ntest,4);
spec{2} = zeros(ntest,4);
spec{3} = zeros(ntest,4);
spec{4} = zeros(ntest,4);

for ii = 4%1:ntest
nx = 10*ii;
ny = 10*ii;
nz = 10*ii;

fprintf('Iteration: %i  Core mesh: %i\n',ii,nx*ny*nz);

dx = round(100 / nx);
dy = round(100 / ny);
dz = round(100 / nz);

dx = [dx*1.4.^[5:-1:1] ones(1,nx)*dx dx*1.4.^[1:5]];
dy = [dy*1.4.^[5:-1:1] ones(1,ny)*dy dy*1.4.^[1:5]];
dz = [dz*1.4.^[5:-1:1] ones(1,nz)*dz dz*1.4.^[1:5]];

nx= length(dx);
ny= length(dy);
nz= length(dz);

% Center mesh at zero
X0 = -sum(dx) / 2;
Y0 = -sum(dy) / 2;
Z0 =  sum(dz) / 2;

% Position of cell centers
x = X0 + cumsum(dx) - dx/2;
y = Y0 + cumsum(dy) - dy/2;
z = Z0 - cumsum(dz) + dz/2;

[X,Y,Z] = ndgrid(x,y,z);

R = (X.^2 + Y.^2 + Z.^2).^0.5;

[model] = ones(nz,nx,ny) * u0;

model(R<r)=2*u0;

m=reshape(model,nx*ny*nz,1);

make_UBC_mesh(work_dir,X0,Y0,Z0,dx,dy,dz);
save([work_dir '\model1.sus'],'-ascii','m');
% Create data points
% Center observations for simplicity
H = 25000; % Inducing field
I = 45 ; % Field inclinaison (degree from h-plane)
D = 45 ;  % Field declinaison (degree from North)

xloc = [-50:10:50] + 1e-3;
yloc = [-50:10:50] + 1e-3;
[ObsX, ObsY, ObsZ] = meshgrid(xloc,yloc,50 +1e-3);
ObsX = reshape(ObsX,size(ObsX,1) * size(ObsX,2), 1);
ObsY = reshape(ObsY,length(ObsX),1);
ObsZ = reshape(ObsZ,length(ObsX),1);

%% Compute analytical solution for a dipole
mdipole = [cosd(I)*cosd(D);cosd(I)*sind(D);sind(I)];

% Compute distance from center to obs (since anomaly at [0,0,0]
robs = (ObsX.^2 + ObsY.^2 + ObsZ.^2).^0.5;
rxy = (ObsX.^2 + ObsY.^2).^0.5;

% Compute field at obs location
costheta = ObsZ./robs;
sintheta = rxy./robs;
cosphi = ObsX./rxy;
sinphi = ObsY./rxy;

rvec = [ sintheta.*cosphi sintheta.*sinphi costheta ] ;
fwr_x = 4*pi*H^2*u0*(r.^3)./(robs.^3) .* (3 * (rvec * mdipole) .* rvec(:,1) - mdipole(1));
fwr_y = 4*pi*H^2*u0*(r.^3)./(robs.^3) .* (3 * (rvec * mdipole) .* rvec(:,2) - mdipole(2));
fwr_z = 4*pi*H^2*u0*(r.^3)./(robs.^3) .* (3 * (rvec * mdipole) .* rvec(:,3) - mdipole(3));

% FWR_data = [fwr_z;fwr_x;fwr_y];
[arg, arg, arg, FWR_data, arg, H, I, D] = read_MAG3D_obs(Obsfile);

% mcell = nx * ny *nz;

% Number of cell faces
nfx = (nx+1) * (ny) * (nz);
nfy = (nx) * (ny+1) * (nz);
nfz = (nx) * (ny) * (nz+1);

% Number of cell edges
nex = (nx) * (ny+1) * (nz+1);
ney = (nx+1) * (ny) * (nz+1);
nez = (nx+1) * (ny+1) * (nz);

nedges = nex + ney + nez;
nnodes = (nx+1)*(ny+1)*(nz+1);

nfaces = nfx + nfy + nfz;

%% Compute forward operator (G) and distance weighting
% Compute depth weigthing matrix
% mode 0: distance weigthing , 1: depth weighting
% pow: Power of the exponentiel decay (default 2 for grav, 3 for mag)

% Discretized topography
topo_model = ones(nx,ny)*200;
ndata = length(ObsX);
% 
% Create active cell matrix (will later need to discretize topography beforehand)
[nullcell,P] = make_nullcell(dx,dy,dz,X0,Y0,Z0,topo_model);

mcell = nx*nz*ny;

Wr=zeros(mcell,1);
V=zeros(mcell,1);



[AVEN,AVCE,AVCN,AVFC,AVEC,GRAD,DIV,DIVbc,CURL,Q,P,PI,DD,Sz,Sx,Sy,H0,Wr,V] = get_ops_nodal(X0, Y0, Z0, dx, dy, dz, ObsX, ObsY, ObsZ, H, I, D);
M = spdiags(AVCE * m,0,size(AVCE,1),size(AVCE,1));

% Primary field formulation
% phi = (DIV * (M * GRAD)) \ (DIVbc * M * Q);
% d =D * PI * (P * AVEN * (H0 + GRAD * ( phi )));

% % Secondary field formulation
% phi = (DIV * (M * GRAD)) \ (-DIV*(M - u0*speye(nedges))*H0);

% Compute phi using Bi-CGStab
A = DIV * (M * GRAD);
b = -DIV*(M - u0*speye(nedges))*H0;
xo = ones( size(A,2),1 ) * 1e-4;



[phi_BiCG,spec{1}(ii,1:3)] = Bi_CG(A,b,xo);
Hs = GRAD * ( phi_BiCG );
% TMI_BiCG =DD * PI * (P * AVEN * ( Hs ));
TMI_BiCG =DD * (P * AVEN * ( Hs ));

[phi_BiCGStab,spec{2}(ii,1:3)] = BiCG_STAB(A,b,xo);
Hs = GRAD * ( phi_BiCGStab );
% TMI_BiCGStab =DD * PI * (P * AVEN * ( Hs ));
TMI_BiCGStab =DD * (P * AVEN * ( Hs ));

[phi_CGNE,spec{3}(ii,1:3)] = CGNE(A,b,xo);
Hs = GRAD * ( phi_CGNE );
% TMI_CGNE =DD * PI * (P * AVEN * ( Hs ));
TMI_CGNE =DD * (P * AVEN * ( Hs ));

[phi_CGSQR,spec{4}(ii,1:3)] = CGSQR(A,b,xo);
Hs = GRAD * ( phi_CGSQR );
% TMI_CGSQR =DD * PI * (P * AVEN * ( Hs ));
TMI_CGSQR =DD * (P * AVEN * ( Hs ));

% Hs = GRAD * ( phi );
% dTMI =DD * PI * (P * AVEN * ( Hs ));

% fwr_d = PI * (P * AVEN * ( Hs ));
% fwr_dz = fwr_d(1:ndata);
% fwr_dx = fwr_d(ndata+1:2*ndata);
% fwr_dy = fwr_d(2*ndata+1:end);

% figure; 
% vecplot(Hs,AVEC,dx',dy',dz');

% figure; 
% vecplot(CURL*Hs,AVFC,dx',dy',dz');

% figure;imagesc(reshape(fwr_dz,length(unique(ObsX)),length(unique(ObsY))));
% % caxis([-1000 1500])
% title('Hz component');colorbar
% 
% figure;imagesc(reshape(fwr_dx,length(unique(ObsX)),length(unique(ObsY))));
% % caxis([-1000 1500])
% title('Hx component');colorbar
% 
% figure;imagesc(reshape(fwr_dy,length(unique(ObsX)),length(unique(ObsY))));
% % caxis([-1000 1500])
% title('Hy component');colorbar
% 
% figure;imagesc(reshape((fwr_dz.^2+fwr_dx.^2+fwr_dy.^2).^0.5,length(unique(ObsX)),length(unique(ObsY))));
% % caxis([-1000 1500])
% title('Hx+Hy+Hz component');colorbar
% 
% figure;imagesc(reshape(Sz*Hs,nx+1,ny+1));
% 
% title('Hz on edge component');colorbar
% 
% figure;imagesc(reshape(Sx*Hs,nx,ny+1));
% 
% title('Hx on edge component');colorbar
% 
% figure;imagesc(reshape(Sy*Hs,nx+1,ny));
% 
% title('Hy on edge component');colorbar


% x = randn(mcell, 1);
% % Random perturbation f(x+a) where a = h*v
% v = randn(mcell, 1);
% 
% N = 10;
% decrease h by an order of magnitude each iteration
% Take norm of left hand side of eqn 2) and 3), take slope, should act
% as right hand side of 2) and 3).

% for jj = 1 : 10
%     h = 10^(-jj);
% 
%     err2 = norm(G*(x + h*v) - G*x); %#ok<*AGROW>
%     err3 = norm(G*(x + h*v) - G*x - h * v);
% 
%     fprintf('%3.2e    %3.2e    %3.2e\n',h , err2, err3)
%     
% end

TMI_data = reshape(FWR_data,length(unique(ObsX)),length(unique(ObsY)));
figure;
imagesc((TMI_data));
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
title('TMI data - UBC MAG3D')
caxis([min(min(TMI_data)) max(max(TMI_data))]);
colorbar


figure;
imagesc(reshape(TMI_CGNE,length(unique(ObsX)),length(unique(ObsY))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
title('TMI data - CGNE')
caxis([min(min(TMI_data)) max(max(TMI_data))]);
colorbar

figure;
imagesc(reshape(TMI_BiCG,length(unique(ObsX)),length(unique(ObsY))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
title('TMI data - BiCG')
caxis([min(min(TMI_data)) max(max(TMI_data))]);
colorbar

figure;
imagesc(reshape(TMI_CGSQR,length(unique(ObsX)),length(unique(ObsY))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
title('TMI data - CGSQR')
caxis([min(min(TMI_data)) max(max(TMI_data))]);
colorbar

figure;
imagesc(reshape(TMI_BiCGStab,length(unique(ObsX)),length(unique(ObsY))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
title('TMI data - BiCGStab')
caxis([min(min(TMI_data)) max(max(TMI_data))]);
colorbar


% fprintf('Residual between Data and CGNE: %10.7e\n',norm(TMI_CGNE-FWR_data));
% fprintf('Residual between Data and BiCG: %10.7e\n',norm(TMI_BiCG-FWR_data));
% fprintf('Residual between Data and CGSQR: %10.7e\n',norm(TMI_CGSQR-FWR_data));
% fprintf('Residual between Data and BiCGStab: %10.7e\n',norm(TMI_BiCGStab-FWR_data));
spec{1}(ii,4) = norm(TMI_CGNE-FWR_data);
spec{2}(ii,4) = norm(TMI_BiCG-FWR_data);
spec{3}(ii,4) = norm(TMI_CGSQR-FWR_data);
spec{4}(ii,4) = norm(TMI_BiCGStab-FWR_data);

end

% figure
% imagesc(reshape(data,length(unique(ObsX)),length(unique(ObsY))))
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);colorbar
% save G
% save Wr
% load GtG
% load Wr
% Square root for the sum of the squares
% Plus another square root of result because inside the objective function,

break

Wr=Wr.^(1/2);

% Normalize depth weighting with the largest value

Wr = Wr./(V');

Wr = Wr./(max(Wr));

Wr=Wr.^(1/2);

% Wr = ones(mcell,1);
% save('Wr.dat','-ascii','Wr')
IWr = spdiags(1./Wr,0,mcell,mcell);
Wr = spdiags(Wr,0,mcell,mcell);

% Create gradient matrices and corresponding volume vectors
[Wx, Vx] = getWx_3D(mcell,dx,dy,dz,nullcell);
[Wy, Vy] = getWy_3D(mcell,dx,dy,dz,nullcell);
[Wz, Vz] = getWz_3D(mcell,dx,dy,dz,nullcell);
[Ws, V ] = getWs_3D(mcell,dx,dy,dz,nullcell);


% Vectorized model


%% Generate data and noise
%Corrupt with 5% random noise
pct_noise = 0.05;
noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
save noise
Wd = spdiags(1./abs(noise),0,ndata,ndata);
d = d + noise;


%% Compute sensitivity
G = Wd * G * IWr;
GtG = G'*G;


% d = data;
% G = G * IWr;
target = ndata;%sum((Wd*(data - d)).^2);

d = Wd * d;

RHS=(G')*d;

% save ('original','data');
% save('data/data.mat','data');
% save('data/kernel.mat','G');
% save('data/model.mat','m');
% save('data/Wr.mat','Wr');
% save ('kernel2','G2');


%% Inversion
Lx = 4*min(dx);
Ly = 4*min(dy);
Lz = 4*min(dz);


as = 1 / Lx^2;
ax = 1;
ay = 1;
az = 1;

pvec= [0.0 1.0 2.0]; 
qvec= [0.0 1.0 2.0]; 
lvec= [0.5 1.0 1.5];

nl=length(lvec);
nq=length(qvec);
np=length(pvec);

finalmodels = zeros(length(pvec),length(qvec),length(lvec),mcell);
finalphid = zeros(length(pvec),length(qvec),length(lvec));
iter_max = 15;
bounds = [0 1];
for pp=3%1:length(pvec)
    
for qq=3%1:length(qvec)   
    
for ll=2%1:length(lvec)
 
beta = 1e+0;
fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl); 
invmod=ones(mcell,1)*1e-4;
phi_d = 99999;
iter = 1;
    while iter<=iter_max && phi_d(end)>target



%             if ii==1                 %First iteration
% 
%                 modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,Vx,Vy,Vz,alphaX,alphaY,alphaZ,alphaC,2,2,lvec(ll),ii);
% 
%             else              
% 
%                 modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,Vx,Vy,Vz,alphaX,alphaY,alphaZ,alphaC,pvec(pp),qvec(qq),lvec(ll),ii);
% 
%             end
%     modfunc = (alphaX*(Wx*Ws)'*(Wx*Ws) + alphaY*(Wy*Ws)'*(Wy*Ws) + alphaZ*(Wz*Ws)'*(Wz*Ws) + alphaC*(Ws)'*(Ws));

    modfunc =  ax*(Vx * Wx * Wr)' * (Vx * Wx * Wr) +...
               ay*(Vy * Wy * Wr)' * (Vy * Wy * Wr) +...
               az*(Vz * Wz * Wr)' * (Vz * Wz * Wr) +...
               as*(Ws * Wr)' * (Ws * Wr);

    % CG solver
    A=(GtG + beta * ((modfunc)));
    phim = invmod' * A * invmod;
    [invmod]=CG_Solver(invmod,A,RHS,bounds);figure(3);hold off
%     invmod=A\RHS;
    d_pred = G*(invmod);
    phi_d(iter) = sum((d_pred-d).^2);
    
    
        if (phi_d(iter) <= target) && (iter > 2)
           break;
        else
           if phi_d(end) < target*2
              beta = 0.75*beta;
           else
              beta = 0.5*beta;
           end
        end

        fprintf('Iteration %i of %i completed, phi_d: %f\n', iter , iter_max,phi_d(iter) )

    iter = iter+1;
    end
    
        counter =1;
    model_out = zeros(nz,nx,ny);
    for jj = 1:ny
    
        for ii = 1:nx

            for kk = 1:nz

                if nullcell(kk,ii,jj) == 0 
                
                    model_out(kk,ii,jj) = -99999;
                    
                else
                    model_out(kk,ii,jj) = invmod(counter);
                    counter = counter +1;
                end
            end
        end
    end
    model_out = reshape(model_out,nz*nx*ny,1);
    save(['results/lambda' num2str(lvec(ll)) 'q' num2str(qvec(qq))  'p' num2str(pvec(pp)) '.den'],'-ascii','model_out')
   
    
end
                  
end

end

% save('3DMAG_lplq.mat','finalmodels','lvec','pvec','qvec','finalphid');


%%
%Create UBC mesh file
fid1=fopen('data/UBC_mesh.msh', 'w');
fprintf(fid1, '%i %i %i\n', nX, nY, nZ);
fprintf(fid1, '%i %i %i\n', X0, Y0, Z0);

for jj=1:nX
    fprintf(fid1, '%4.2f ', dx(jj));    
end

fprintf(fid1,'\n');

for ii=1:nY
           fprintf(fid1,'%4.2f ', dy(ii));
end

fprintf(fid1, '\n');

for kk=1 : nZ
       fprintf(fid1, '%4.2f ', dz(kk));
end

fclose(fid1);

%Create UBC observation file
count = 1;
fid3 = fopen('data/UBC_obs.obs','w');
fprintf(fid3,'%i\n',ndata);

for ii=1:ndata

    fprintf(fid3,'%4.2f %4.2f %4.2f %8.5f %8.5f\n',ObsX(ii),ObsY(ii),ObsZ(ii),Wd(ii,ii)\d(ii),Wd);
    count = count + 1;

end
fclose(fid3);

figure (2)
imagesc(reshape(Wd\d_pred,length(xloc),length(yloc)))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')