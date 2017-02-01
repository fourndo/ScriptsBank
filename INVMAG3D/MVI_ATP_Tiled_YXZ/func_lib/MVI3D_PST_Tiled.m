function [M,beta_out] = MVI3D_PST_Tiled(work_dir,out_dir,dsep,xn,yn,zn,H, HI, HD, obsx, obsy, obsz, G, d, wd,mstart,mref,esus,chi_target,alphas,beta,bounds,LP,t,eps_FLAG,eps_tresh,FLAG1,FLAG2,max_iter)
% Magnetic Vector Invertion in cartesian coordinates
% Written by: D Fournier 
% Created: 2014/07/23

fprintf('Setting up the inversion...\n');
% Load mesh file and convert to vectors (UBC format)
% [xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = nx*ny*nz;

% pct_cutoff = 50;
if eps_FLAG == 1
    eps_q = eps_tresh(2);
    eps_p = eps_tresh(1);
    
else
    
    eps_q = ones(1,3)*1e-1;
    eps_p = ones(1,3)*1e-1;
end
        
%% Create model magnetization vectors
% normalize effective susceptibility
if length(esus)==1
esus = ones(mcell,1)*esus;
end

% % esus = sqrt(esus / max(esus) + 1e-1);
% esus = esus  / norm(esus);
% esus = kron(ones(3,1), esus );

mref = kron(ones(3,1), mref );
% 
mstart = kron([0;0;1], mstart );

% Need to create I/O for cell base weights
w = ones(4*mcell,1);

% Create bound vector
lowBvec = [ones(mcell,1) * bounds(1,1);ones(mcell,1) * bounds(2,1);ones(mcell,1) * bounds(3,1)];
uppBvec = [ones(mcell,1) * bounds(1,2);ones(mcell,1) * bounds(2,2);ones(mcell,1) * bounds(3,2)];

%% Initialize dynamic cells
% Create selector matrix for active cells
load([work_dir dsep 'nullcell.dat']);
x = spdiags(nullcell,0,mcell,mcell);
x = x(nullcell==1,:);
X = kron(speye(3),x);

mactv = sum(nullcell);

mstart = X * mstart;
mref = X * mref;

% Normalize effective susceptibility weight
esus(esus==-100) = 0;
esus = x*esus;
wj = esus  / max(esus) + 1e-2;

% Wj = spdiags( kron(ones(3,1), wj ) ,0 , 3*mactv , 3*mactv );

lowBvec = X * lowBvec;
uppBvec = X * uppBvec;
t = t(nullcell==1);
%% Load observation file (3C UBC-MAG format)
% [H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);

%% Calculating sensitivity weigthing
tic
wr = zeros(1,3*mactv);
for gg = 1 : ndata
    temp = G{1}(gg,:);
    wr = wr + temp.^2;
end
    
wr = sqrt(wr);
wr = (wr'/max(wr));
Wr = spdiags(wr,0,3*mactv,3*mactv);

Wd = spdiags(1./wd,0,ndata,ndata);
toc

d = Wd * d;
%% Create gradient matrices and corresponding volume vectors
[~, Gx, Gy, Gz, ~, ~, ~, ~] = get_GRAD_op3D_SQUARE_Kron(dx,dy,dz,nullcell,x);
% [Ws, V ] = getWs3D(dx,dy,dz,X);

Ws =  spdiags(x * ( w(1:mcell) ) ./wj ,0,mactv,mactv);
Wx =  spdiags(x * ( w(1+mcell:2*mcell) ) ./wj ,0,mactv,mactv);
Wy =  spdiags(x * ( w(1+2*mcell:3*mcell) ) ./wj ,0,mactv,mactv);
Wz =  spdiags(x * ( w(1+3*mcell:4*mcell) ) ./wj ,0,mactv,mactv);

mactv = sum(nullcell);
%% START INVERSION
    
target = chi_target * ndata;     % Target misifit

% Compute total objective function
objfunc = @(m,phim,b) sum( ((Gvec(G,Wd,m) - d) ).^2 ) + (m)' * b * phim * (m);

invmod      =  mstart ;

% phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = sum(((Gvec(G,Wd,invmod) - d) ).^2);

Pac = speye(3*mactv);

% Message prompt
logfile = [work_dir dsep 'Log_TMVI.log'];
fid = fopen(logfile,'w');
fprintf(fid,'Starting lp inversion\n');
fprintf(fid,'Starting misfit %e\n',phi_d);
fprintf(fid,'Target misfit %e\n',target);
fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
fprintf(fid,'#cut cells \t # CG Iter\n');

count= 0; % Initiate iteration count 
tncg = 0; % Compute total number of CG iterations for the whole inversion
lp_count = 0;

% Switcher defines the different modes of the inversion
% switcher 0: Run the usual l2-l2 inversion, beta decreases by 2
% switcher 1: Run the lp-lq inversion, beta decreases by 0.8
% swircher 2: End on inversion - model update < treshold
% if phid_in < target
%     
%     switcher = 1; 
%     
% else
%     
%     switcher = 0; 
%     
% end

dphim = 100;
switcher = 0;
while switcher ~= 3 && count ~= max_iter

    fprintf(['MVI-Cartesian Formulation\n'])
    count=count+1;
    
    if switcher == 0     
        
            
            delta_p(count,1:3) = 1e-1;%delta_tresh(1)*3;
            delta_q(count,1:3) = 1e-1;%delta_tresh(2)*3;
            
            [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_p(count,:),delta_q(count,:));

        

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum(((Gvec(G,Wd,temp))).^2) / (temp'*MOF*temp) * 1e+3 ;

        end
        
        phi(count) = norm((Gvec(G,Wd,invmod) - d)).^2 +...
                    invmod' * beta(count) * MOF * invmod;
%         tresh = dkdt(2,delta(count));

        
    else

        if lp_count == 1 
            for pst = 1 : 3
                
                if eps_FLAG==0 

                    Z = zeros(3,3);
                    Z(pst,pst) = 1;

                    % Create sub-space matrix and grab the right model
                    Z = kron(Z,speye(mactv));
                    m = Z*(invmod-mref);

                    [pp,qq] = get_eps(m,10,Gx,Gy,Gz);

                    delta_p(pst,1:3) = pp;
                    delta_q(pst,1:3) = qq;

                else

                    delta_p(pst,1:3) = eps_p(pst);
                    delta_q(pst,1:3) = eps_q(pst);

                end
                
            end

            model_out = X'*(invmod);
  
            M = reshape(model_out,mcell,3);
            Mamp = sum(M.^2,2).^0.5;
            Mamp(nullcell==0) = -100;

            pred_TMI = Gvec(G,speye(ndata),invmod);
            
            write_MAG3D_TMI([work_dir dsep 'Tile_MVI_PST_l2l2.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
            save([work_dir dsep 'Tile_MVI_PST_l2l2.fld'],'-ascii','M')
            save([work_dir dsep 'Tile_MVI_PST_l2l2.amp'],'-ascii','Mamp')

            figure;
            subplot(1,3,1)
            temp = invmod(1:mactv) - mref(1:mactv);
            [n, xout] =hist((temp),100); hold off
            bar(xout,n);  
            set(gca,'yscale','log')
            title('L2 amp values')

            subplot(1,3,2)
            [n, xout] =hist((Gx * invmod(1+mactv:2*mactv)),100); hold off
            bar(xout,n);  
            set(gca,'yscale','log')
            title('L2 \Delta \theta values')

            subplot(1,3,3)
            [n, xout] =hist((Gx * invmod(1+2*mactv:3*mactv)),100); hold off
            bar(xout,n);  
            set(gca,'yscale','log')
            title('L2 \Delta \phi values')


        end 
        
        if dphim(end)  < 5%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1
            
            fprintf('\n# # ADJUST BETA # #\n');
            switcher = 2;

        else

            fprintf('\n# # LP-LQ ITER# #\n');

        end
                    
        fprintf('# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p,delta_q);

        
    end
    
    % Save previous model before GN step
    m_in = invmod; 
    
    mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;
        
    %% Pre-conditionner
    diagG = zeros(1,3*mactv);
    for gg = 1 : ndata
        diagG = diagG + (G{1}(gg,:)/wd(gg)).^2;
    end
        
    diagA = diagG + beta(count)*spdiags(mof,0)';
    PreC     = spdiags(1./diagA(:),0,3*mactv,3*mactv);

    %% Gauss-Newton steps
    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
%     fprintf('eps_q: \t %8.5e \t eps_q*: \t %8.5e\n',delta_p(count,1),eps_p)
%     fprintf('eps_p: \t %8.5e \t eps_p*: \t %8.5e\n',delta_q(count,1),eps_q)
    
     
    [invmod, ncg, Pac] = GN_CG_Lin_solver( G, Wd, invmod, mref, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, mof, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );
    
    tncg = tncg + ncg;
    
    %% Save iteration and continue    
    % Measure the update length
    if count==1 
        
        rdm(count) =  1;
        gdm(1) = norm(m_in - invmod);
        
    else
        
        gdm(2) = norm(m_in - invmod);
        rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);
    
        gdm(1) = gdm(2);
        
    end
    
    phi_d(count) = sum(( ( Gvec(G,Wd,invmod) - d ) ).^2);
    phi_m(count) = (invmod)'*(MOF)*(invmod);
    phi(count) = objfunc(invmod,mof,beta(count));
    

    if count ~= 1
        
        dphim(count) = abs(phi_m(count) - phi_m(count-1)) / phi_m(count) *100;
        
    end

    fprintf(' phi_d:\t %8.5e \n',phi_d(count))
    fprintf(' phi_m:\t %8.5e \n',phi_m(count))
    fprintf(' dphi_m:\t %8.5e \n',dphim(count))

    % Get next beta
    [switcher,beta(count+1)] = cool_beta(beta(count),phi_d(count),dphim(count),target,switcher,0.25,2);

    % Right log file
    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf(fid,' \t %8.5e ',sum( (aVRWs*invmod).^2 ) );
    fprintf(fid,' \t %8.5e ',sum( (aVRWx*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWy*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWz*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
    fprintf(fid,' \t %8.5e ',phi(count));
    fprintf(fid,' \t\t %i\n',ncg);

    model_out = X'*(invmod);
  
    M = reshape(model_out,mcell,3);
    Mamp = sum(M.^2,2).^0.5;
    Mamp(nullcell==0) = -100;
      
    pred_TMI = Gvec(G,speye(ndata),invmod);
end
beta_out = beta(count);
save([out_dir dsep 'Tile_MVI_PST_Iter0.fld'],'-ascii','M')
save([out_dir dsep 'Tile_MVI_PST_Iter0.amp'],'-ascii','Mamp')
write_MAG3D_TMI([out_dir dsep 'Tile_MVI_PST_Iter0.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
fclose(fid);
