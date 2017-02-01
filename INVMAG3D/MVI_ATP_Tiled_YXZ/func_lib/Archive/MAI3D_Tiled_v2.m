function [m_esus,pre_ampB,beta_out,phid_out,MOF_out,switcher] = MAI3D_Tiled(work_dir,out_dir,dsep,idx,xn,yn,zn,H,HI,HD, obsx, obsy,obsz, G ,d,wd,mstart,mref,mag_xyz,w,chi_target,phid_in,MOF_in,alphas,beta,bounds,LP,t,delta_FLAG,delta_tresh,switcher,FLAG1,FLAG2,max_iter,beta_tol,ROT)
% Amplitude inversion code
% Dominique Fournier 
% Last update: May 31th, 2015


% Load mesh file and convert to vectors (UBC format)
% [xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);



%% Load observation file (3C UBC-MAG format)
% [H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);
% wd = ones(ndata,1);
Wd   = spdiags(1./wd,0,ndata,ndata);


load([work_dir dsep 'nullcell.dat']);

% Create bound vector
lowBvec = ones(mcell,1) * bounds(1);
uppBvec = ones(mcell,1) * bounds(2);

%% Initialize dynamic cells
% load([work_dir '\nullcell.dat']);

% Create selector matrix for active cells
X = spdiags(nullcell,0,mcell,mcell);
X = X(nullcell==1,:);

mactv = sum(nullcell);

mstart = X * mstart;
mref = X * mref;
mag_xyz = X * mag_xyz;
lowBvec = X* lowBvec;
uppBvec = X * uppBvec;
t = X * t;

%% Create model magnetization vectors
% Azimuth and dip of magnitization

% M = [spdiags(H * mag_xyz(:,1),0,mactv,mactv);spdiags(H * mag_xyz(:,2),0,mactv,mactv);spdiags(H * mag_xyz(:,3),0,mactv,mactv)];


%% Compute depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, 'DISTANCE', 2.75, min(dx)/4);
% wr = save([work_dir '\wr.dat'],'-ascii','wr');
% wr = load([work_dir dsep 'wr.dat']);
% wr = get_wrAmp(obsx, obsy, obsz, d, D, I, xn, yn, zn, nullcell, FLAG, 3, min(dx)/4);
% wr = wr.^(4/3);
% wr = X * wr;
% Wr = spdiags(wr,0,mactv,mactv);



%% Forward operator and apply depth weighting
fprintf('Loading Sensitivity...\n');
% load([work_dir '\Tx']);
% Fx = G{1} * M;
% clear Tx
% 
% % load([work_dir '\Ty']);
% Fy = G{2} * M;
% clear Ty
% 
% % load([work_dir '\Tz']);
% Fz = Tz * M;
% clear Tz


%% Apply data weighting
% Fx = Wd * Fx ;
% Fy = Wd * Fy ;
% Fz = Wd * Fz ;

d = Wd * d;

% ampB = @(m) ( (Wd*(Fx * m)).^2 + (Wd*(Fy * m)).^2 + (Wd*(Fz * m)).^2 ) .^ 0.5;

% bx = Fx *  ones(mactv,1); %bx = bx/norm(bx);
% by = Fy * ones(mactv,1); %by = by/norm(by);
% bz = Fz * ones(mactv,1); %bz = bz/norm(bz);
% 
% bx  = spdiags( bx , 0 , ndata, ndata);
% by  = spdiags( by , 0 , ndata, ndata);
% bz  = spdiags( bz , 0 , ndata, ndata);

% magB   = ampB(invmod);

% lBl   = spdiags( d.^-1 , 0 , ndata, ndata);
% 
% J   = (lBl * [bx by bz]) *  ([Fx;Fy;Fz]) ;
%         j = sum(J.^2,1);
%         j = j.^0.25;
% % 
%         wr = (j' / max(j)).^0.5;
% 
%         save([work_dir dsep 'ApproxJ.dat'],'-ascii','temp');
%% Create gradient matrices and corresponding volume vectors
% [~, Gx, Gy, Gz, V, ~, ~, ~] = get_GRAD_op3D_SQUARE_Kron(dx,dy,dz,nullcell);
% [Ws, V ] = getWs3D(dx,dy,dz,X);

[A, GRAD, V] = get_GRAD_op3D_TENSIL_Kron(dx,dy,dz,nullcell,'FWR');

Rz = @(x)   [cosd(x) -sind(x) 0;
            sind(x) cosd(x) 0;
            0 0 1];

Ry = @(x)   [cosd(x) 0 -sind(x);
            0 1 0;
            sind(x) 0 cosd(x)];

Rx = @(x)   [1 0 0;
            0 cosd(x) -sind(x);
            0 sind(x) cosd(x)];
        
rz = Rz(ROT(1));
%     rz = rz*spdiags(sum(abs(rz),2).^-1,0,3,3);

ry = Ry(ROT(2));
%     ry = ry*spdiags(sum(abs(ry),2).^-1,0,3,3);

rx = Rx(ROT(3));

rot = rx*ry*rz;

% Get index and weight for gradients Gx
[val,ind] = sort(A*rot(1,:)') ;

indx = find(acosd(val)<=30);

Gx = GRAD{ind(indx(1))} * val(indx(1));
denom = val(indx(1));
for ii = 2 : length(indx)
    
    Gx = Gx + GRAD{ind(indx(ii))} * val(indx(ii));
    denom = denom + val(indx(ii));
    
end

Gx = Gx * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gx),2)) ~= 2;
% Gx(indx,:) = 0;

% Get index and weight for gradients Gy
[val,ind] = sort(A*rot(2,:)') ;

indx = find(acosd(val)<=30);

Gy = GRAD{ind(indx(1))} * val(indx(1));
denom = val(indx(1));
for ii = 2 : length(indx)
    
    Gy = Gy + GRAD{ind(indx(ii))} * val(indx(ii));
    denom = denom + val(indx(ii));
    
end

Gy = Gy * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gy),2)) ~= 2;
% Gy(indx,:) = 0;


% Get index and weight for gradients Gz
[val,ind] = sort(A*rot(3,:)') ;

indx = find(acosd(val)<=30);

Gz = GRAD{ind(indx(1))} * val(indx(1));
denom = val(indx(1));
for ii = 2 : length(indx)
    
    Gz = Gz + GRAD{ind(indx(ii))} * val(indx(ii));
    denom = denom + val(indx(ii));
    
end

Gz = Gz * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gz),2)) ~= 2;
% Gz(indx,:) = 0;


Ws =  spdiags(X * ( w(1:mcell))  ,0,mactv,mactv);
Wx =  spdiags(X * ( w(1+mcell:2*mcell) )  ,0,mactv,mactv);
Wy =  spdiags(X * ( w(1+2*mcell:3*mcell))  ,0,mactv,mactv);
Wz =  spdiags(X * ( w(1+3*mcell:4*mcell))  ,0,mactv,mactv);

% v = X * spdiags(V);
alphas = [1 1 1 1];

%% Inversion
target = chi_target * ndata;

comp_phi = @(m,mof,l) sum( ( ampB(G,Wd,m) - d ).^2 ) +...
    (m)' * l * mof * (m);
           
% Initialize inversion
invmod      = mstart;       % Initial model        

% phi_init    = sum((ampB(invmod) - d).^2);   % Initial misfit
phi_d       = phid_in;
phi_m       = 1;
phi         = phid_in + beta*phi_m;

% Message prompt
%             head = ['lp' num2str(pvec(pp)) '_lq' num2str(qvec(qq)) '_mu' num2str(lvec(ll))];
logfile = [work_dir dsep 'MAG_AMP_INV_.log'];
fid = fopen(logfile,'w');
fprintf(fid,'Starting Effective Susceptibility inversion\n');
fprintf(fid,'Starting misfit %e\n',phid_in);
fprintf(fid,'Target misfit %e\n',target);
fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
fprintf(fid,'#cut cells \t # CG Iter\n');


% Initiate active cell
Pac = speye(mactv);

count= 0; % Initiate iteration count 

lp_count = 0;


% Switcher defines the different modes of the inversion
% switcher 0: Run the usual l2-l2 inversion, beta decreases by 2
% switcher 1: Run the lp-lq inversion, beta decreases by 0.8
% swircher 2: End on inversion - model update < treshold
eps_p = 1;
eps_q = 1;

tresh_p = eps_p;
tresh_q = eps_q;

dphim = 100;
while switcher ~= 3 && count~=max_iter 

    fprintf(['Amplitude Inversion: Tile' num2str(idx) '\n'])
    count=count+1;
      


    if  switcher == 0                

        eps_p = 0.1;
        eps_q = 0.1;
        
        delta_p(count) = eps_p;
        delta_q(count) = eps_q;
%                     [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(invmod,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,2,2,1);
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,[2 2 2 2 1],FLAG1,FLAG2,switcher,delta_p(count),delta_q(count));

        if isempty(beta)==1

            
            magB   = ampB(G,Wd,invmod)  ;
    
            bx  = spdiags( Wd * ( G{1} * (invmod) ) , 0 , ndata, ndata);
            by  = spdiags( Wd * ( G{2} * (invmod) ) , 0 , ndata, ndata);
            bz  = spdiags( Wd * ( G{3} * (invmod) ) , 0 , ndata, ndata);

        %     lBl   = ampB(invmod);

            lBl   = spdiags( magB.^-1 , 0 , ndata, ndata);

            B{1} = Wd * lBl * bx;
            B{2} = Wd * lBl * by;
            B{3} = Wd * lBl * bz;

        
            diagJ = zeros(1,mactv);
            bbx = spdiags(B{1});
            bby = spdiags(B{2});
            bbz = spdiags(B{3});

            for gg = 1 : ndata
                diagJ = diagJ +  ( bbx(gg) * G{1}(gg,:) + bby(gg) * G{2}(gg,:) + bbz(gg) * G{3}(gg,:)).^2;
            end

            % Re-weighted pseudo sensitivity weighting
            wr = sqrt(diagJ');
            wr = (wr/max(wr));


    %         j = sum(J.^2,1);
    %         j = (j'.^0.5)./v.^2;
    % 
    %         wr = (j / max(j)).^0.5;
            Wr = spdiags(wr,0,mactv,mactv);

            MOF = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

            temp = randn(mactv,1);
            beta = sum(Gvec(G,B,temp).^2) / (temp'*MOF*temp) * 1e+2 ;
%             beta = 5e+5;
        end
        
        

%         tresh_s = 1;
%         tresh_xyz = 1;
        
    else

        lp_count = lp_count+1;

        if lp_count == 1
            
            [eps_p,eps_q] = get_eps(invmod,10,Gx,Gy,Gz);                
            eps_p = eps_p/5;
            eps_q = eps_q/5;
            
        end

        if delta_p(end)~= eps_p %(dphi_p(end) > 2 || lp_count ==1) && switch_p == 0
                            
%             delta_p(count) = delta_p(count-1)*.5;
            delta_p(count) = eps_p;
            if delta_p(count) < eps_p

                delta_p(count) = eps_p;

            end
        else 

            delta_p(count) = delta_p(count-1);%delta_p(count-1);

        end

        if delta_q(end) ~= eps_q %&& switcher == 1%(dphi_q(end) > 2 || lp_count ==1) && switch_q == 0

%             delta_q(count) = delta_q(count-1)*.5;
            delta_q(count) = eps_q;
            if delta_q(count) < eps_q

                delta_q(count) = eps_q;

            end

        else

%             switch_q = 1;
            delta_q(count) = delta_q(count-1);%delta_q(count-1);

        end
        
        if dphim(end)  < 1 && delta_p(count) == eps_p && delta_q(count) == eps_q%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1
            fprintf('\n# # ADJUST BETA # #\n');
            switcher = 2;
            
        else
            
            fprintf('\n# # LP-LQ ITER# #\n');
            
        end

%         fprintf('\n# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p(count),delta_q(count));

%         tresh = dkdt(LP(:,1),delta(count));

        
    end


    
    %% Gauss Newton Steps
    
    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
    fprintf('eps_p: \t %8.5e \t eps_p*: \t %8.5e\n',delta_p(count),eps_p)
    fprintf('eps_q: \t %8.5e \t eps_q*: \t %8.5e\n',delta_q(count),eps_q)
    
    
    % Save current model
    m_in = invmod;
%     dmdx = sqrt( (Wx * invmod).^2 + (Wy * invmod).^2 + (Wz * invmod).^2 );
    
    tncg = 0;
    ggdm = 1;       % Mesure the relative change in rel|dm|
    ddm = [1 1];    % Mesure of change in model update |dm|
    solves = 1;
      

%         magB   = ampB(invmod)  ;
%         bx  = spdiags( Fx * invmod , 0 , ndata, ndata);
%         by  = spdiags( Fy * invmod , 0 , ndata, ndata);
%         bz  = spdiags( Fz * invmod , 0 , ndata, ndata);
% 
%         magB   = ampB(invmod)  ;
% 
%         lBl   = spdiags( magB.^-1 , 0 , ndata, ndata);
% 
%         J   = lBl * ([bx by bz] *  ([Fx;Fy;Fz] * spdiags(invmod,0,mactv,mactv)) );    


                  
    
    
    while solves < 10 && ggdm > 1e-3 
%         [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p(count),delta_q(count));

%         if solves == 1
%             invmod = pz*invmod;
%             
%         end
        magB   = ampB(G,Wd,invmod)  ;

        bx  = spdiags( Wd * ( G{1} * (invmod) ) , 0 , ndata, ndata);
        by  = spdiags( Wd * ( G{2} * (invmod) ) , 0 , ndata, ndata);
        bz  = spdiags( Wd * ( G{3} * (invmod) ) , 0 , ndata, ndata);

        %     lBl   = ampB(invmod);

        lBl   = spdiags( magB.^-1 , 0 , ndata, ndata);

        B{1} = Wd * lBl * bx;
        B{2} = Wd * lBl * by;
        B{3} = Wd * lBl * bz;
            
%         bb = lBl * [bx by bz];
        
%         J   = B *  ([Fx;Fy;Fz]) ;
        
        diagJ = zeros(1,mactv);
        bbx = spdiags(B{1});
        bby = spdiags(B{2});
        bbz = spdiags(B{3});
        
        for gg = 1 : ndata
            diagJ = diagJ +  ( bbx(gg) * G{1}(gg,:) + bby(gg) * G{2}(gg,:) + bbz(gg) * G{3}(gg,:)).^2;
        end

        % Re-weighted pseudo sensitivity weighting
        wr = sqrt(diagJ');
        wr = (wr/max(wr));

            
%         j = sum(J.^2,1);
%         j = (j'.^0.5)./v.^2;
% 
%         wr = (j / max(j)).^0.5;
        Wr = spdiags(wr,0,mactv,mactv);

        mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;
% 
%         save([work_dir dsep 'ApproxJ.dat'],'-ascii','temp');
            
%         A = [ J  ;...
%         sqrt( beta(count) ) * Wr* aVRWs ;...
%         sqrt( beta(count) ) * Wr* aVRWx ;...
%         sqrt( beta(count) ) * Wr* aVRWy ;...
%         sqrt( beta(count) ) * Wr* aVRWz ];

%         MOF = Wr* MOF*Wr;
        
        diagA   = diagJ' + beta(count)*spdiags( mof ,0);
        PreC    = Pac * spdiags(1./diagA(:),0,mactv,mactv);
        
        
        
        switch FLAG1

            case 'SMOOTH_MOD'
%                 g = [- (magB - d) ; ...
%             - sqrt( beta(count) ) * (  Wr* aVRWs * (invmod-mref) ) ;...
%             - sqrt( beta(count) ) * (  Wr* aVRWx * (invmod) ) ;...
%             - sqrt( beta(count) ) * (  Wr* aVRWy * (invmod) ) ;...
%             - sqrt( beta(count) ) * (  Wr* aVRWz * (invmod) ) ];

                g = - Gtvec(G, B,(magB - d)) - beta(count) * (mof * invmod);           
        
            case 'SMOOTH_MOD_DIF'%
%                 g = [- (magB - d) ; ...
%             - sqrt( beta(count) ) * ( aVRWs * (invmod-mref) ) ;...
%             - sqrt( beta(count) ) * ( aVRWx * (invmod-mref) ) ;...
%             - sqrt( beta(count) ) * ( aVRWy * (invmod-mref) ) ;...
%             - sqrt( beta(count) ) * ( aVRWz * (invmod-mref) ) ];
        end 
        
        phi(count) = comp_phi(invmod, mof,beta(count));
        phi_in = phi(end);
        %% Projected steepest descent
        dm = zeros(mactv,1);
        %[dm,~,ncg] = PCGLSQ( dm, A , g, PreC, Pac);
        [dm,~,ncg] = CG_Lin( dm, G, B, beta(count) * mof , g,speye(mactv), PreC, Pac );
%         fprintf('CG iter %i \n',ncg);
        
        %% Step length, line search                
        tncg = tncg+ncg; % Record the number of CG iterations

        temp = spdiags(Pac);
        
        % Combine active and inactive cells step if active bounds
        if sum(temp)~=mactv
            
            rhs_a = ( speye(mactv) - Pac ) * g;
            dm_i = max( abs( dm ) );
            dm_a = max( abs(rhs_a) );                
            dm = dm + rhs_a * dm_i / dm_a /10 ;

        end
        gamma = 2;
        
        % Reduce step length in order to reduce phid
        phi_out = phi_in;
        m_temp = invmod;
        while (phi_out >= phi_in || gamma == 2) && gamma > 1e-3

%             phi_temp(2) = phi_temp(1);

            gamma = 0.5 * gamma;

            gdm = gamma * dm;

            ddm(2) = norm(gdm);

            m_temp = invmod + gdm;

            lowb = m_temp <= lowBvec;
            uppb = m_temp >= uppBvec;

            % Apply bound on model
            m_temp(lowb==1) = lowBvec(lowb==1);
            m_temp(uppb==1) = uppBvec(uppb==1);

            % Update projection matrix
            Pac = spdiags((lowb==0).*(uppb==0),0,mactv,mactv);
            phi_out = comp_phi(m_temp, mof,beta(count));


        end

        phi_in = phi_out;
        
        if solves == 1

            ggdm = 1;
            ddm(1) = ddm(2);

        else

            ggdm = ddm(2)/ddm(1);

        end


        % Update model
        invmod = m_temp;

        
%         phi_m(count) = (invmod)'*(MOF)*(invmod);
        
        fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,ggdm);

        solves = solves + 1;
         
         
    end
    
%     dmdx = sqrt( (Wx * invmod).^2 + (Wy * invmod).^2 + (Wz * invmod).^2 );
    
%     group_s(count) = sum(abs(m_in) <= tresh_p);
%     group_xyz(count) = sum(abs(dmdx) <= tresh_q);
    
    clear A 
    %% Save results and update beta       
    % Measure the update length
    if count==1 
        
        rdm(count) =  1;
        gdm(1) = norm(m_in - invmod);
        
    else
        
        gdm(2) = norm(m_in - invmod);
        rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);
    
        gdm(1) = gdm(2);
        
    end
    
    phi_d(count) = sum((ampB(G,Wd,invmod)-d).^2);
    phi_m(count) = (invmod)'*(MOF)*(invmod);
    phi(count) = comp_phi(invmod, mof,beta(count));
    
    if count ~= 1
        
        dphim(count) = abs(phi(count) - phi(count-1)) / phi(count) *100;
        
    end

    %Update depth weighting
%     j = sum(J.^2,1);
%     j = j.^0.5;
%     wr = (j' / max(j)).^0.5;
% 
%     Wr = spdiags(wr,0,mactv,mactv);
    
%     fprintf('---------->\n')
    fprintf(' phi_d:\t %8.5e \n',phi_d(count))
    fprintf(' phi_m:\t %8.5e \n',phi_m(count))
    fprintf(' dphi_m:\t %8.5e \n',dphim(count))

%     fprintf('Final Relative dm:\t %8.5e \n', rdm(count));
%     fprintf('<----------\n')
    
%     fprintf('Number of Inactive cells: %i\n',sum(tcells));
%     fprintf('Number of CGS iterations: %i\n\n',ncg);
    
    [switcher,beta(count+1)] = cool_beta(beta(count),phi_d(count),dphim(count),target,switcher,beta_tol,1);

    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf(fid,' \t %8.5e ', alphas(1)*sum(aVRWs*invmod).^2) ;
    fprintf(fid,' \t %8.5e ', alphas(2)*sum(aVRWx*invmod).^2) ;
    fprintf(fid,' \t %8.5e ', alphas(3)*sum(aVRWy*invmod).^2) ;
    fprintf(fid,' \t %8.5e ', alphas(4)*sum(aVRWz*invmod).^2) ;
    fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
    fprintf(fid,' \t %8.5e ',phi(count));
%     fprintf(fid,' \t\t %i ',sum(tcells));
    fprintf(fid,' \t\t %i\n',ncg);
   % Output interation result


    

m_esus = X'* ( invmod);
m_esus(nullcell==0) = -100;

% if lp_count==0
    
%     save([work_dir dsep 'Tile' num2str(idx) '_MAI_l2l2.sus'],'-ascii','m_esus')
    
% else
    
    save([work_dir dsep 'Tile' num2str(idx) '_MAI_lplq.sus'],'-ascii','m_esus')
    
% end
pre_ampB = ampB(G,Wd,invmod).*wd;

write_MAG3D_TMI([work_dir dsep 'Tile' num2str(idx) '_MAI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pre_ampB,wd);


end

if switcher==2
    
    switcher=3;
    
end
% count=count-1;            

fprintf(fid,'End of lp inversion. Total Number of iterations: %i\n',tncg);
fclose(fid);

%%
save([out_dir dsep 'Tile' num2str(idx) '_MAI_esus.sus'],'-ascii','m_esus')

pre_ampB = ampB(G,Wd,invmod).*wd;

write_MAG3D_TMI([out_dir dsep 'Tile' num2str(idx) '_MAI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pre_ampB,wd);

beta_out = beta(end);
phid_out = phi_d(end);
phim_out = phi_m(end);
MOF_out = MOF;
