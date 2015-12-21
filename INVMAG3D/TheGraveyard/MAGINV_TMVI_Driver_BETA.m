% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23

clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Synthetic\Nut_Cracker\MagInv_mvec';
inpfile   = 'MAGINV_TMVI.inp'; 

% [meshfile,obsfile,susfile,chi_target,alphas,beta,pvec,qvec,lvec,FLAG] = MAG3CM_read_inp([work_dir '\' inpfile]);
[meshfile,obsfile,mstart,mref,magfile,chi_target,alphas,beta,p_bound,s_bound,t_bound,pvec,qvec,lvec,FLAG] = MAGINV_TMVI_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

load([work_dir '\nullcell.dat']);


% Initialize dynamic cells
Di = kron( speye(3), spdiags(nullcell,0,mcell,mcell));

% Initialize selector matrix for the p,s,t components
Sp = kron([1;0;0],ones(mcell,1));
Ss = kron([0;1;0],ones(mcell,1));
St = kron([0;0;1],ones(mcell,1));

lowBp = zeros(3*mcell,1) == 1;
lowBs = zeros(3*mcell,1) == 1;
lowBt = zeros(3*mcell,1) == 1;


% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);

ndata = length(obsx);

if mod(ndata,1)~=0
    
    fprintf('Data does not appear to be 3C. Please revise...\n');
    break
    
end


%% Depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
% save([work_dir '\wr.dat'],'-ascii','wr');

wr = load([work_dir '\wr.dat']);
IWr = kron(speye(3),spdiags(1./wr,0,mcell,mcell));
Wr = kron(speye(3),spdiags(wr,0,mcell,mcell));

%% Compute T matrix
% load([work_dir '\Tx']);
% load([work_dir '\Ty']);
% load([work_dir '\Tz']);
% 
% G = [Tx;Ty;Tz];
% 
% clear Tx Ty Tz T

load([work_dir '\Gp']);
load([work_dir '\Gs']);
load([work_dir '\Gt']);

G = [Gp Gs Gt];
clear Gp Gs Gt

% wd = abs(d)*0.05 + 0.05*std(d);
Wd = spdiags(1./wd,0,ndata,ndata);

% G = Wd * G * Esus * H ;

G = Wd * G * IWr;
d = Wd * d;

%% Create gradient matrices and corresponding volume vectors
[Wx, Wy, Wz, Vx, Vy, Vz] = get_GRAD_op3D_v4(dx,dy,dz,nullcell);
[Ws, V ] = getWs3D(dx,dy,dz);
% Ws = kron(speye(3), Ws);
% Wx = kron(speye(3), Wx);
% Wy = kron(speye(3), Wy);
% Wz = kron(speye(3), Wz);
% Vx = kron(speye(3), Vx);
% Vy = kron(speye(3), Vy);
% Vz = kron(speye(3), Vz);
% V = kron( speye(3), V );
% WxtWx = ( Vx * Wx )' * ( Vx * Wx ) ;
% WytWy = ( Vy * Wy )' * ( Vy * Wy ) ;
% WztWz = ( Vz * Wz )' * ( Vz * Wz ) ;
% WstWs = ( V * Ws )' * ( V * Ws ) ;
            

%% Load reference model

if length(mref) == 1
    
    mref = ones(mcell,1) * mref;

else
    
    mref = load([work_dir '\' mref]);
                
end

mref = [mref;kron([1;1],zeros(mcell,1))];

%% Load starting model model

if length(mstart) == 1
    
    mstart = ones(mcell,1) * mstart;

else
    
    mstart = load([work_dir '\' mstart]);
                
end

mstart = [mstart;kron([1;1],ones(mcell,1)*1e-3)];

%% Inversion
count=1;
target = chi_target * ndata;

comp_phi = @(m,phi,l) sum( ( G * m - d ).^2 ) +...
    (m)' * l * phi * (m);

nullcell = kron(ones(3,1),nullcell); 
                        
mref        =  Wr*mref ;       % Reference model 
invmod      =  Wr*mstart ;

phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phi_init;
phi_m       = [];  
Pac = spdiags(nullcell==1 & lowBp==0 &...
    lowBs==0 & lowBt==0 , 0 , 3*mcell, 3*mcell);

            % Message prompt
head = ['lp' num2str(pvec) '_lq' num2str(qvec) '_mu' num2str(lvec)];
logfile = [work_dir '\Log_lBl_' head '.log'];
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

            [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD_TMVI_v2(invmod,mref,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alphas,[2 2 2],[2 2 2],[1 1 1],FLAG);

            MOF_start = MOF;

            if isempty(beta)==1

                beta = full( sum(sum(G.^2,1)) / sum(diag(MOF,0)) * 1e+4 );

            end

            lambda = beta ;
            phi =  comp_phi(invmod,MOF,lambda);


        elseif count~=1

            [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD_TMVI_v2(invmod,mref,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alphas,pvec,qvec,lvec,FLAG);

            mu = ( invmod'* MOF_start * invmod ) / (invmod'*MOF*invmod ) ;

            lambda = beta(count) * mu;


        end

        diagA = sum(G.^2,1) + lambda*spdiags(MOF,0)';
        PreC  = spdiags(1./diagA(:),0,3*mcell,3*mcell); 


         A = [ G  ;...
            sqrt( lambda ) * aVRWs ;...
            sqrt( lambda ) * aVRWx ;...
            sqrt( lambda ) * aVRWy ;...
            sqrt( lambda ) * aVRWz ];

    switch FLAG

        case 'SMOOTH_MOD'
            g = [- (G * invmod - d) ; ...
        - sqrt( lambda ) * ( aVRWs * (invmod-mref) ) ;...
        - sqrt( lambda ) * ( aVRWx * (invmod) ) ;...
        - sqrt( lambda ) * ( aVRWy * (invmod) ) ;...
        - sqrt( lambda ) * ( aVRWz * (invmod) ) ];

        case 'SMOOTH_MOD_DIF'
            g = [- (G * invmod - d) ; ...
        - sqrt( lambda ) * ( aVRWs * (invmod-mref) ) ;...
        - sqrt( lambda ) * ( aVRWx * (invmod-mref) ) ;...
        - sqrt( lambda ) * ( aVRWy * (invmod-mref) ) ;...
        - sqrt( lambda ) * ( aVRWz * (invmod-mref) ) ];
    end

    %% PCG steps
    ncg = 0;
    for npcg = 1:3
    % Gaussian Newton solver
    % Form Hessian and gradient (only half since solving using
    % CGLSQ)
    A = [ G ;...
        sqrt( lambda ) * aVRWs ;...
        sqrt( lambda ) * aVRWx ;...
        sqrt( lambda ) * aVRWy ;...
        sqrt( lambda ) * aVRWz ];

     switch FLAG

        case 'SMOOTH_MOD'

        g = [- (G *invmod - d) ; ...
            - sqrt( lambda ) * ( aVRWs * (invmod - mref) ) ;...
            - sqrt( lambda ) * ( aVRWx * (invmod) ) ;...
            - sqrt( lambda ) * ( aVRWy * (invmod) ) ;...
            - sqrt( lambda ) * ( aVRWz * (invmod) ) ];

        case 'SMOOTH_MOD_DIF'
                g = [- (G *invmod - d) ; ...
            - sqrt( lambda ) * ( aVRWs * (invmod-mref) ) ;...
            - sqrt( lambda ) * ( aVRWx * (invmod-mref) ) ;...
            - sqrt( lambda ) * ( aVRWy * (invmod-mref) ) ;...
            - sqrt( lambda ) * ( aVRWz * (invmod-mref) ) ];
     end

    dm = zeros(3*mcell,1);
    [dm,r,iter] = PCGLSQ( dm, A , g, PreC, Pac);

    %% Step length, line search                
    ncg = ncg+iter; % Record the number of CG iterations

    % Combine active and inactive cells step if active bounds
    if (sum(lowBp) + sum(lowBs) + sum(lowBt))~=0
        rhs_a = ( Di - Pac ) * (A'*g);

        dm_i = max( abs( dm ) );
        dm_a = max( abs(rhs_a) );                
        dm = dm + rhs_a * dm_i / dm_a /10 ;

    end
    gamma = 2;

    % Initialise phi^k
    phi_temp = 0;   
    while phi_temp > phi(end) || gamma == 2

        gamma = 0.5 * gamma;

        m_temp = invmod + gamma * dm;

        lowBp(nullcell==1 & Sp) = m_temp(nullcell==1 & Sp) <= p_bound(1);
        lowBs(nullcell==1 & Ss) = m_temp(nullcell==1 & Ss) <= s_bound(1);
        lowBt(nullcell==1 & St) = m_temp(nullcell==1 & St) <= t_bound(1);

        % Apply bound on model
        m_temp(lowBp==1) = p_bound(1);
        m_temp(lowBs==1) = s_bound(1);
        m_temp(lowBt==1) = t_bound(1);

        % Update projection matrix
        Pac = spdiags(nullcell==1 & lowBp==0 &...
        lowBs==0 & lowBt==0 , 0 , 3*mcell,3*mcell);

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
%                 fprintf(fid,' \t %8.5e ',invmod'*alphas(1)*WstWs*invmod);
%                 fprintf(fid,' \t %8.5e ',invmod'*alphas(2)*WxtWx*invmod);
%                 fprintf(fid,' \t %8.5e ',invmod'*alphas(3)*WytWy*invmod);
%                 fprintf(fid,' \t %8.5e ',invmod'*alphas(4)*WztWz*invmod);
    fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
    fprintf(fid,' \t %8.5e ',phi(count));
    fprintf(fid,' \t\t %i\n',ncg);
    fprintf('Number of CGS iterations: %i\n',ncg);

   

    count = count+1;
end
%% Ouput result
  model_out = IWr*invmod;
    % Create orthogonal magnetization vectors
    Mp =  model_out(1:mcell);%M + IWr * Esus * invmod;
    Ms =  model_out((1+mcell) : 2*mcell);
    Mt =  model_out((1+2*mcell) : 3*mcell);

    % Create total magnetization vector
%     M = [Mp(1:mcell)+Ms(1:mcell)+Mt(1:mcell),...
%         Mp((1+mcell):2*mcell)+Ms((1+mcell):2*mcell)+Mt((1+mcell):2*mcell),...
%         Mp((1+2*mcell):3*mcell)+Ms((1+2*mcell):3*mcell)+Mt((1+2*mcell):3*mcell)];

M = [-Ms -Mt Mp];

    % Create absolute magnetization
%                 model_out(nullcell==0,:) = -100;

    Mamp = sqrt( Mp.^2 + Ms.^2 + Mt.^2 );
    Mrem = sqrt( Ms.^2 + Mt.^2 );
    Mind = sqrt( Mp.^2);

%                 lMl = sqrt(sum(model_out.^2,2));
%                 lml = [model_out(:,1)./lMl  model_out(:,2)./lMl model_out(:,3)./lMl];

    save([work_dir '\Mvec' 'p' num2str(pvec) 'q' num2str(qvec) 'l' num2str(lvec) '.fld'],'-ascii','M')
    save([work_dir '\M_' num2str(pvec) 'q' num2str(qvec) 'l' num2str(lvec) '.amp'],'-ascii','Mamp')
    save([work_dir '\M_' num2str(pvec) 'q' num2str(qvec) 'l' num2str(lvec) '.ind'],'-ascii','Mind')
    save([work_dir '\M_' num2str(pvec) 'q' num2str(qvec) 'l' num2str(lvec) '.rem'],'-ascii','Mrem')           

% plot_TMI(obsx,obsy,d,pred,wd,'Observed vs Predicted')
%% Plot resulst
% invmod = IWr*invmod;
pred = (G*invmod)  ;
plot_TMI(obsx,obsy,d,pred,wd,'Final PRED data')
write_MAG3D_TMI([work_dir '\TMI_iter_' num2str(count) '.pre'],H,I,Dazm,obsx,obsy,obsz,pred.*wd,wd)
                
% plot_mag3C(obsx,obsy,dtot - pred,I,D,'Final residual data')

% plot_mag3C(obsx,obsy,pred3C./[wdx;wdy;wdz],I,D,' Inversion - Predicted')

