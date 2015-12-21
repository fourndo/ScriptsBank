% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23

clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath /tera23/dfournier/Code/MAG3CampB/functions/

% Project folders
work_dir = '/tera_raid/dfournier/TKC/Mag/MAG3C_Inv1';
inpfile   = 'MAG3CvecB.inp'; 

[meshfile,obsfile,susfile,wr_flag,chi_target,alphas,beta,pvec,qvec,lvec] = MAGvec3C_read_inp([work_dir '/' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '/' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

% Load effective susceptibility model
esus = load([work_dir '/' susfile]);
P = kron( speye(3),spdiags(esus,0,mcell,mcell) );

% Initialize dynamic cells
load([work_dir '/nullcell.dat']);
Di = spdiags(nullcell,0,mcell,mcell);

% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '/' obsfile]);

ndata = length(obsx);

if mod(ndata,1)~=0
    
    fprintf('Data does not appear to be 3C. Please revise...\n');
    break
    
end

Wd = spdiags(1./wd,0,3*ndata,3*ndata);

%% Create model magnetization vectors
m_azm = ones(mcell,1)*Dazm;
m_dip = ones(mcell,1)*I;
mv = azmdip_2_xyz(m_azm,m_dip);

M = H *[mv(:,1);mv(:,2);mv(:,3)];

%% Depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
% save([work_dir '\wr.dat'],'-ascii','wr');
% 
% IWr = kron(speye(3),spdiags(1./wr,0,mcell,mcell));
% Wr = kron(speye(3),spdiags(wr,0,mcell,mcell));

%% Compute T matrix
load([work_dir '/Tx']);
load([work_dir '/Ty']);
load([work_dir '/Tz']);

G = [Tx;Ty;Tz];

clear Tx Ty Tz

G = Wd * G * P;

d = Wd * d;

%% Create gradient matrices and corresponding volume vectors
[Wx, Wy, Wz, Vx, Vy, Vz] = get_GRAD_op3D_v4(dx,dy,dz,nullcell);
[Ws, v ] = getWs3D(dx,dy,dz);
Ws = kron(speye(3), Ws);
Wx = kron(speye(3), Wx);
Wy = kron(speye(3), Wy);
Wz = kron(speye(3), Wz);
Vx = kron(speye(3), Vx);
Vy = kron(speye(3), Vy);
Vz = kron(speye(3), Vz);
V = kron( speye(3),spdiags((v),0,mcell,mcell) );
WxtWx = ( Vx * Wx )' * ( Vx * Wx ) ;
WytWy = ( Vy * Wy )' * ( Vy * Wy ) ;
WztWz = ( Vz * Wz )' * ( Vz * Wz ) ;
WstWs = ( V * Ws )' * ( V * Ws ) ;
            

%% Inversion
count=1;
target = chi_target * ndata;

comp_phi = @(m,phi,l) sum( ( G * m - d ).^2 ) +...
    (m)' * l * phi * (m);

for ll= 1:length(lvec)

    for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)   
            
            
            mref        = M;       % Reference model 
            invmod      = M;

            phi_init    = sum((G * invmod - d).^2);   % Initial misfit
            phi_d       = phi_init;
            phi_m       = [];  
            
                        % Message prompt
            head = ['lp' num2str(pvec(pp)) '_lq' num2str(qvec(qq)) '_mu' num2str(lvec(ll))];
            logfile = [work_dir '/Log_lBl_' head '.log'];
            fid = fopen(logfile,'w');
            fprintf(fid,'Starting lp inversion %s\n',head);
            fprintf(fid,'Starting misfit %e\n',phi_init);
            fprintf(fid,'Target misfit %e\n',target);
            fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
            fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
            fprintf(fid,'#cut cells \t # CG Iter\n');
            
            count=1;
           
                    
            while phi_d(end)>target

                
                if count==1                 

                    [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD(invmod,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alphas,2,2,1);
              
                    MOF_start = MOF;
                    
                    if isempty(beta)==1
                        
                        beta = full( sum(sum(G.^2,1)) / sum(diag(MOF,0)) * 1e+4 );
                        
                    end
                        
                    lambda = beta ;
                    phi =  comp_phi(invmod,MOF,lambda);
                    
                                        
                else
                    
                    [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD(invmod,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alphas,pvec(pp),qvec(qq),lvec(ll));

                    mu = ( invmod'* MOF_start * invmod ) / (invmod'*MOF*invmod ) ;

                    lambda = beta(count) * mu;
                    

                end
                
               %% PCG steps
                ncg = 0;
                for npcg = 1:3
                diagA = sum(G.^2,1) + lambda*spdiags(MOF,0)';
%                 PreC     = spdiags(1./diagA(:),0,3*mcell,3*mcell); 
                    
                 A = [ G  ;...
                    sqrt( lambda ) * aVRWs ;...
                    sqrt( lambda ) * aVRWx ;...
                    sqrt( lambda ) * aVRWy ;...
                    sqrt( lambda ) * aVRWz ];

                g = [- (G * invmod - d) ; ...
                    - sqrt( lambda ) * ( aVRWs * (invmod - mref) ) ;...
                    - sqrt( lambda ) * ( aVRWx * invmod ) ;...
                    - sqrt( lambda ) * ( aVRWy * invmod ) ;...
                    - sqrt( lambda ) * ( aVRWz * invmod ) ];


                %% Projected steepest descent
                dm = zeros(3*mcell,1);
                [dm,r,iter] = CGLSQ( dm, A , g);
                ncg = ncg + iter;
                
                %% Step length, line search
                gamma = 2;

                % Initialise phi^k
                phi_temp = 0;   
                while phi_temp > phi(end) || gamma == 2
                    
                    gamma = 0.5 * gamma;
                    
                    m_temp = invmod + gamma * dm;
                    phi_temp = comp_phi(m_temp,MOF,lambda);
                    
                end
                
                % Update model
                invmod = m_temp;
                end
                
                phi(count) = comp_phi(invmod,MOF,lambda);
                phi_d(count) = sum((G*invmod-d).^2);

                % Cool beta
                if phi_d(count) < target*2 && count~=1
                    
                  beta(count+1) = 0.5 * beta(count);

                else

                  beta(count+1) = 0.5 * beta(count);

                end


                fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
                fprintf('Iteration: \t %i  \nBeta: \t %8.5e \n',count,beta(count));
                fprintf(fid,' \t %8.5e ',phi_d(count));
                fprintf('phid:\t %8.5e\n',phi_d(count))
                fprintf(fid,' \t %8.5e ',invmod'*alphas(1)*WstWs*invmod);
                fprintf(fid,' \t %8.5e ',invmod'*alphas(2)*WxtWx*invmod);
                fprintf(fid,' \t %8.5e ',invmod'*alphas(3)*WytWy*invmod);
                fprintf(fid,' \t %8.5e ',invmod'*alphas(4)*WztWz*invmod);
                fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
                fprintf(fid,' \t %8.5e ',phi(count));
                fprintf(fid,' \t\t %i\n',ncg);
                fprintf('Number of CGS iterations: %i\n',ncg);

                
                model_out = P*invmod;
                model_out = reshape(model_out,mcell,3);
                model_out(nullcell==0,:) = -100;
                
                save([work_dir '/Mvec' 'p' num2str(pvec(pp)) 'q' num2str(qvec(qq)) 'l' num2str(lvec(ll)) '.mod'],'-ascii','model_out')
   

            	count = count+1;
            end
            
           
        end
        
    end
    
end

% plot_TMI(obsx,obsy,d,pred,wd,'Observed vs Predicted')

% invmod = IWr*invmod;
% pred = G*invmod;
% plot_mag3C(obsx,obsy,pred,I,D,'Final PRED data')

% plot_mag3C(obsx,obsy,pred3C./[wdx;wdy;wdz],I,D,' Inversion - Predicted')

