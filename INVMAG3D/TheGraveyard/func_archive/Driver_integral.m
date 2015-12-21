% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath C:\Users\Thomasc\Dropbox\Master\Miscellaneous\;
addpath functions;

% Project folders
work_dir = 'C:\Projects\UBC_Research\Modelling\Inversion\Synthetic\Benshtest';
result_dir = 'C:\Projects\UBC_Research\Modelling\Inversion\Synthetic\Benshtest';

obsfile = 'magpre3d.mag';
meshfile = 'Mesh_4m_core.msh';
model_sus = 'Sphere.sus';
% model_azm = 'Synthetic_Band_azm.sus';
% model_dip = 'Synthetic_Band_dip.sus';
topofile = 'Topo_flat_1585m.topo';

u0 = 4 * pi * 10^-7;



% Create weighting matrix from assigned noise
% Wd = spdiags( wd , 0, ndata, ndata ) ;



% Load mesh file and convert to vectors (UBC format)
[mesh]=get_UBC_mesh([work_dir '\' meshfile]);
nx = mesh(1,1); %size(X,1);    %number of cell in X
ny = mesh(1,2); %size(X,2);    %number of cell in Y
nz = mesh(1,3); %size(X,3);    %number of cell in Z

mcell = nx*ny*nz;

dx = mesh(3,1:nx);
dy = mesh(4,1:ny);
dz = mesh(5,1:nz);

x0 = mesh(2,1);
y0 = mesh(2,2);
z0 = mesh(2,3);

% Load synthetic model
m = load([work_dir '\' model_sus]);

% Create nullcell
[nullcell] = topo_2_model(work_dir,meshfile,topofile);
save([work_dir '\nullcell.dat'],'-ascii','nullcell')
load([work_dir '\nullcell.dat']);


% Load observation file (UBC-MAG format)
[H, I, D, obsx, obsy, obsz, data, arg] = read_MAG3D_TMI([work_dir '\' obsfile]);
% [H, I, D, obsx, obsy, obsz, MAG3D_Bx, MAG3D_By, MAG3D_Bz, wdx, wdy, wdz] = read_MAG3D_3C([work_dir '\' obsfile]);
ndata = length(obsx);

% Create model magnetization vectors

m_azm = ones(mcell,1)*D;
m_dip = ones(mcell,1)*I;
% m_azm = load([work_dir '\' model_azm]);
% m_dip = load([work_dir '\' model_dip]);
M = azmdip_2_xyz(m_azm,m_dip,dx,dy,dz);

% Compute forward model
[ Bx, By, Bz, TMI, magB, obsx, obsy, obsz ] = MAG3D_FWR( m, M, H, D, I, nullcell, obsx, obsy, obsz, x0, y0, z0, dx, dy, dz);


%% Compute forward operator (G) and distance weighting
% Compute depth weigthing matrix
% pow: Power of the exponentiel decay (default 2 for grav, 3 for mag)
%
% Space allocation for cellsize matrix and distance weighting
% Wr=zeros(mcell,1);
% V=zeros(mcell,1);


%% Initiate G matrix - TMI integral equation
% G=zeros(ndata,mcell);
% count = 1;
% tic
% for ii=1:ndata;
%     
%        [G(count,:),wr,V] = Fwr_Mag_Integral_DEVB(mcell, X0, Y0, Z0, dx', dy', dz',...
%            ObsX(ii), ObsY(ii), ObsZ(ii), H, I, D, nullcell);
%         
%         Wr = Wr + wr;
%         
% count = count + 1;
% 
% end
% toc
% 
% 
% save([disk_dir '\FWR_op'],'G');
% save([disk_dir '\Wr'],'Wr');
% save([disk_dir '\V'],'V');
% load([disk_dir '\FWR_op']);
% load([disk_dir '\Wr']);
% load([disk_dir '\V']);
% save(['C:\Users\dominiquef\Desktop' 'FWR_op'],'G');
% d = G*m;

% Plot data in 2D
plot_mag3C(obsx,obsy,Bx,By,Bz,TMI,magB,'Forward')

figure;
scatter(obsx,obsy,30,data,'filled')
% imagesc(reshape(TMI,length(unique(obsx)),length(unique(obsy))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-800 2000]);
colorbar
title('\bfMAG3D Forward data - TMI')
write_MAG3D_3C(work_dir,'MAG3D_3C.obs',H,I,D,obsx,obsy,obsz,Bx,By,Bz);

break
%% Scale the depth weighting
% Square root for the sum of the squares
% Plus another square root of result because inside the objective function,
Wr = Wr.^(1/2);

% Normalize depth weighting with the largest value

Wr = Wr./V;

Wr = Wr./(max(Wr));

Wr=Wr.^(1/2);

% Wr = ones(mcell,1);
% save('Wr.dat','-ascii','Wr')
IWr = spdiags(1./Wr,0,mcell,mcell);
Wr = spdiags(Wr,0,mcell,mcell);

%% Create gradient matrices and corresponding volume vectors
[Wx, Vx] = getWx_3D(mcell,dx,dy,dz,numercell);
[Wy, Vy] = getWy_3D(mcell,dx,dy,dz,numercell);
[Wz, Vz] = getWz_3D(mcell,dx,dy,dz,numercell);
[Ws, V ] = getWs_3D(mcell,dx,dy,dz,numercell);

%% Generate data and noise
%Corrupt with 5% random noise
pct_noise   = 0.05;
noise       = (pct_noise.*max(abs(data))).*randn(ndata,1);
save([func_dir '\noise'],'noise')
load([func_dir '\noise'])
Wd          = spdiags(1./abs(noise),0,ndata,ndata);
d           = d + noise;


%% Compute sensitivity
G   = Wd * G * IWr;
GtG = G'*G;


% d = data;
% G = G * IWr;
target = ndata;%sum((Wd*(data - d)).^2);

d = Wd * d;

RHS=(G')*d;


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
    
for qq=1%1:length(qvec)   
    
for ll=1:length(lvec)
 
beta = 1e+2;
fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl); 
invmod=ones(mcell,1)*1e-4;
phi_d = 99999;
iter = 1;
    while iter<=iter_max && phi_d(end)>target



            if iter==1                 %First iteration

                modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,Vx,Vy,Vz,ax,ay,az,as,2,2,lvec(ll),iter);

            else              

                modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,Vx,Vy,Vz,ax,ay,az,as,pvec(pp),qvec(qq),lvec(ll),iter);

            end
 
    modfunc2 =  ax*(Vx * Wx * Wr)' * (Vx * Wx * Wr) +...
               ay*(Vy * Wy * Wr)' * (Vy * Wy * Wr) +...
               az*(Vz * Wz * Wr)' * (Vz * Wz * Wr) +...
               as*(Ws * Wr)' * (Ws * Wr);

    % CG solver
    A=(GtG + beta * ((modfunc)));
    phim = invmod' * A * invmod;
%     invmod = ConjGrad(invmod,A,RHS,bounds); figure(100);hold off
    invmod=A\RHS;
    d_pred = G*(invmod);
    phi_d(iter) = sum((d_pred-d).^2);
    
    
        if (phi_d(iter) <= target) && (iter > 2)
            fprintf('Iteration %i of %i completed, phi_d: %f\n', iter , iter_max,phi_d(iter) )
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

    model_out = invmod;
    model_out(nullcell==0) = -1;
    save([result_dir '\' 'p' num2str(pvec(pp)) 'q' num2str(qvec(qq)) 'l' num2str(lvec(ll)) '.sus'],'-ascii','model_out')
   
    
end
                  
end

end

% save('3DMAG_lplq.mat','finalmodels','lvec','pvec','qvec','finalphid');

figure
imagesc(reshape(Wd\d,length(unique(obsx)),length(unique(obsy))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);colorbar
title('\bfNoisy data')
colorbar

figure
imagesc(reshape(Wd\d_pred,length(unique(obsx)),length(unique(obsy))))
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);colorbar
title('\bfPredicted data')
colorbar
