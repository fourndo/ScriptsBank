function [M,pred_TMI,beta_out,phid_out,MOF_out,switcher] = MVI3D_APD_Tiled(work_dir,out_dir,dsep,idx,xn,yn,zn,H, HI, HD, obsx, obsy, obsz, G, d, wd,mstart,mref,esus,chi_target,alphas,beta,bounds,LP,t,eps_FLAG,eps_tresh,FLAG1,FLAG2,max_iter)
% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23


% Load mesh file and convert to vectors (UBC format)
% [xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = nx*ny*nz;


ndata = length(d);
            
Wd = spdiags(1./wd,0,ndata,ndata);

d = Wd * d;

scl_t = 2*pi;
scl_p = 2*pi;

% for ll = 1
%     
%     for oo = 0
%         
%         for kk = 0 
            
            switcher = 0;
            
            for pst = 1 : 3
                if eps_FLAG == 1
                    eps_q(pst) = eps_tresh(2);
                    eps_p(pst) = eps_tresh(1);

                else

                    eps_q(pst) = ones(1,3)*1e-1;
                    eps_p(pst) = ones(1,3)*1e-1;
                end
            end

            %% Create model magnetization vectors
            % normalize effective susceptibility
            % if length(esus)==1
            % esus = ones(mcell,1)*esus;
            % end

            % % esus = sqrt(esus / max(esus) + 1e-1);
            % esus = esus  / norm(esus);
            % esus = kron(ones(3,1), esus );


            % 
            % mstart = kron(ones(3,1), mstart );

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

%             mstart = X * mstart;
            % mref = X * mref;



            % Normalize effective susceptibility weight
            esus(esus==-100) = 0;
            esus = x*esus.^0;
            wj = esus  / max(esus) + 1e-2;

            % Wj = spdiags(esus,0,3*mactv,3*mactv);

            lowBvec = X * lowBvec;
            uppBvec = X * uppBvec;

            %% Load observation file (3C UBC-MAG format)
            % [H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
            % plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
            % plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');




            %% Load forward operator
            % TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
            fprintf('Loading Sensitivity...\n');
            % load([work_dir '\Tx']);
            % load([work_dir '\Ty']);
            % load([work_dir '\Tz']);

            % Txyz = zeros(ndata,3*mactv);
            % 
            % % Do it in for loop to save memory
            % for ii = 1 : ndata
            %     
            % Txyz(ii,:) = TMI * [Tx(ii,:);Ty(ii,:);Tz(ii,:)] * H;
            % 
            % end

            % clear Tx Ty Tz

            % Case cartesian coordinates
            % [P,S,T] = azmdip_2_pst(Dazm,I,mactv);

            % Create orthogonal forward operators
            % Primary (inducing)
            % G = [Txyz];
            % clear Txyz P S T Tx Ty Tz
            % wd = abs(d)*0.05 + 0.05*std(d);
            

            % G = Wd * G * Esus * H ;

            % G = Wd * G * Wj;


            % Case sperical
            m_uvw = @(a,t,p) [a.*cos(t).*cos(p);...
                a.*cos(t).*sin(p);...
                a.*sin(t)];

            % v = @(a,theta,phi) a*cos(theta)*sin(phi);
            % w = @(a,theta,phi) a*sin(theta);

            sProj = @(a,theta,phi)[spdiags(cos(theta).*cos(phi),0,mactv,mactv) spdiags(-a.*sin(theta).*cos(phi),0,mactv,mactv) spdiags(-a.*cos(theta).*sin(phi),0,mactv,mactv);
                spdiags(cos(theta).*sin(phi),0,mactv,mactv) spdiags(-a.*sin(theta).*sin(phi),0,mactv,mactv) spdiags(a.*cos(theta).*cos(phi),0,mactv,mactv);
                spdiags(sin(theta),0,mactv,mactv) spdiags(a.*cos(theta),0,mactv,mactv) sparse(mactv,mactv)];

            % comp_phi = @(dof,mof,l) sum( ( G*m(a,t,p) - d ).^2 ) +...
            %     (m)' * l * phi * (m);

            

            %% Depth weighting
            aa = mstart(1:mcell);
            tt = mstart(1+mcell:2*mcell);
            pp = mstart(1+2*mcell:3*mcell);

            %         J   =  ;
            % S = sProj(aa,tt,pp);
            %         
            % diagJ = zeros(1,3*mactv);
            % for gg = 1 : ndata
            %     diagJ = diagJ +  ( G{1}(gg,:) * S).^2;
            % end
            %         
            % % Re-weighted pseudo sensitivity weighting for each component
            % wr = sqrt(diagJ');
            % 
            % wr(1:mcell)= (wr(1:mcell)/max(wr(1:mcell)));
            % wr(1+mcell:2*mcell)= (wr(1+mcell:2*mcell)/max(wr(1+mcell:2*mcell))) * (1e-3/pi).^2;
            % wr(1+2*mcell:3*mcell)= (wr(1+2*mcell:3*mcell)/max(wr(1+2*mcell:3*mcell))) * (1e-3/pi).^2;
            % 
            % Wr = spdiags(wr,0,3*mactv,3*mactv);

            % Wr = spdiags(wr,0,3*mactv,3*mactv);

            %% Create gradient matrices and corresponding volume vectors
            [~, Gx, Gy, Gz, ~, ~, ~, ~] = get_GRAD_op3D_SQUARE(dx,dy,dz,nullcell,x);
            % [Ws, V ] = getWs3D(dx,dy,dz,X);

            Ws =  spdiags(x * ( w(1:mcell) ./ wj )  ,0,mactv,mactv);
            Wx =  spdiags(x * ( w(1+mcell:2*mcell) ./ wj)  ,0,mactv,mactv);
            Wy =  spdiags(x * ( w(1+2*mcell:3*mcell) ./ wj)  ,0,mactv,mactv);
            Wz =  spdiags(x * ( w(1+3*mcell:4*mcell) ./ wj)  ,0,mactv,mactv);


            mactv = sum(nullcell);
            %% START INVERSION

            target = chi_target * ndata;     % Target misifit

            % Compute total objective function
            objfunc = @(m,phim,b) sum( ((Gvec(G,Wd,m) - d) ).^2 ) + (m)' * b * phim * (m);

            invmod      =  mstart ;
            aa = invmod(1:mcell);
            tt = invmod(1+mcell:2*mcell);
            pp = invmod(1+2*mcell:3*mcell);



            % phi_init    = sum((G * invmod - d).^2);   % Initial misfit
            phi_d       = sum(((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d) ).^2);
            % phi_m       =    (invmod-mref)'* ( aVRWs'*aVRWs * (invmod-mref) ) +...
            %                   (invmod)'*( aVRWx'*aVRWx * (invmod) ) +...
            %                   (invmod)'*( aVRWy'*aVRWy * (invmod) ) +...
            %                   (invmod)'*( aVRWz'*aVRWz * (invmod) ) ;  
            % phi         = phi_d  + beta*phi_m;

            % Initiate active cell
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

            % dkdt = @(p,ep) (ep).^(1/(2*(2-p)));
            % traffic_s = 1;    % Measures the number of points/iteration passing the lp corner
            % traffic_xyz = 1;    % Measures the number of points/iteration passing the lp corner
            % group_s = 1;  % Measures the number of lower than the lp corner
            dphim = 100;
            while switcher ~= 3 && count ~= max_iter

                fprintf(['MVI: Tile' num2str(idx) '\n'])
                count=count+1;

                if switcher == 0     

                    delta_p(count,1:3) = 1e-1;%delta_tresh(1)*3;
                    delta_q(count,1:3) = 1e-1;%delta_tresh(2)*3;


                    aa = invmod(1:mcell);
                    tt = invmod(1+mcell:2*mcell);
                    pp = invmod(1+2*mcell:3*mcell);

            %         J   =  ;
                    S = sProj(aa,tt,pp);

                    [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_p(count,:),delta_q(count,:));

                    wr = zeros(1,3*mactv);
                    for gg = 1 : ndata
                        wr = wr +  ( (G{1}(gg,:)* S)).^2;
                    end   

            %         % Re-weighted pseudo sensitivity weighting for each component
                    wr = abs(sqrt(wr+1e-8))';
                    wr = wr / max(wr);

                    % Try scaling the transformation matrix
            %         wrs = sum(S.^2,1)'.^0.5;

%                     wr(1:mcell)= (wr(1:mcell)/10.^kk);
                    wr(1+mcell:2*mcell)= (wr(1+mcell:2*mcell)/scl_t) ;
                    wr(1+2*mcell:3*mcell)= (wr(1+2*mcell:3*mcell)/scl_p);

                    Wr = spdiags(wr,0,3*mactv,3*mactv);

                    mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

                    if isempty(beta)==1

                        temp = randn(3*mactv,1);
                        beta = sum(((Gvec(G,Wd,temp))).^2) / (temp'*mof*temp) * 1e+2 ;

                    end


                    phi_MOF = (invmod-mref)'* ( aVRWs'*Wr*aVRWs * (invmod-mref) ) +...
                          (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
                          (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
                          (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;

                    phi(count) = norm((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d)).^2 + beta(count) + phi_MOF;
            %         tresh = dkdt(2,delta(count));


                else

                    lp_count = lp_count+1;

                    if lp_count == 1 && eps_FLAG==0

                        for pst = 1 : 3

                            Z = zeros(1,3);
                            Z(pst) = 1;

                            % Create sub-space matrix and grab the right model
                            Z = kron(Z,speye(mactv));

                            m = Z*invmod;

                            [pp,qq] = get_eps(m,10,Gx,Gy,Gz);

                            eps_p(pst) = pp;
                            eps_q(pst) = qq;

                        end


                    end

                    for pst = 1 : 3

                        if delta_p(end,pst)> eps_p(pst)

                            delta_p(count,pst) = eps_p(pst);

                            if delta_p(count,pst) < eps_p(pst)

                                delta_p(count,pst) = eps_p(pst);

                            end

                        else

                                delta_p(count,:) = eps_p(pst);

                        end

                        if delta_q(end,pst)> eps_q(pst)

                            delta_q(count,pst) = eps_q(pst);

                            if delta_q(count,pst) < eps_q(pst)

                                delta_q(count,pst) = eps_q(pst);

                            end

                        else

                            delta_q(count,pst) = eps_q(pst);

                        end


                    end 

                    if dphim(end)  < 5%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1

                        fprintf('\n# # ADJUST BETA # #\n');
                        switcher = 2;

                    else

                        fprintf('\n# # LP-LQ ITER# #\n');

                    end

                    %% TEMPORARY CHANGE
                    % Compute regularization on total vector
            %         m_tot = sqrt( invmod(1:mactv).^2 + invmod(1+mactv:2*mactv).^2 + invmod(2*mactv+1:3*mactv).^2);

            %         m_tot = kron(ones(3,1),m_tot);

                    fprintf('# # LP-LQ ITER# #\n');
                    [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p(count,:),delta_q(count,:));

            %         tresh = dkdt(LP(:,1),delta(count));


                end

            %     dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );    
            %     group_xyz(count) = sum(abs(dmdx) <= tresh_xyz);


                %% Gauss-Newton steps

                fprintf('\n# # # # # # # # #\n');
                fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
                str_var = ['P','S','T'];
                for pst = 1 : 3
                    fprintf([str_var(pst) 'eps_p: \t %8.5e \t eps_p*: \t %8.5e\n'],delta_p(count,pst),eps_p(pst))
                    fprintf([str_var(pst) 'eps_q: \t %8.5e \t eps_q*: \t %8.5e\n'],delta_q(count,pst),eps_q(pst))
                end
                
                % Save current model
                m_in = invmod;
            %     dmdx = sqrt( (Wx * invmod).^2 + (Wy * invmod).^2 + (Wz * invmod).^2 );
            %     Wr = spdiags(wr,0,3*mactv,3*mactv);

            %     mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

                tncg = 0;
                ggdm = 1;       % Mesure the relative change in rel|dm|
                ddm = [1 1];    % Mesure of change in model update |dm|
                solves = 1;

                phi_in = phi(end);

                while solves < 5 && ggdm > 1e-3

            %         magB   = ampB(invmod)  ;
            % 
            %         bx  = spdiags( Fx * invmod , 0 , ndata, ndata);
            %         by  = spdiags( Fy * invmod , 0 , ndata, ndata);
            %         bz  = spdiags( Fz * invmod , 0 , ndata, ndata);
            % 
            %         %     lBl   = ampB(invmod);
            % 
            %         lBl   = spdiags( magB.^-1 , 0 , ndata, ndata);
            % 
                    aa = invmod(1:mcell);
                    tt = invmod(1+mcell:2*mcell);
                    pp = invmod(1+2*mcell:3*mcell);

            %         J   =  ;
                    S = sProj(aa,tt,pp);
            %         A = [  G * S(aa,tt,pp) ;...
            %         sqrt( beta(count) ) * aVRWs ;...
            %         sqrt( beta(count) ) * aVRWx ;...
            %         sqrt( beta(count) ) * aVRWy ;...
            %         sqrt( beta(count) ) * aVRWz ];

            %         diagA   = sum((G * S(aa,tt,pp)).^2,1) + beta(count)*spdiags(MOF,0)';
            %         PreC    = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);

                    wr = zeros(1,3*mactv);
                    for gg = 1 : ndata
                        wr = wr +  ( (G{1}(gg,:)* S)).^2;
                    end        
                    wr = abs(sqrt(wr+1e-8))';
                    wr = wr / max(wr);

                    % Try scaling the transformation matrix
            %         wrs = sum(S.^2,1)'.^0.5;

%                     wr(1:mcell)= (wr(1:mcell)/10.^kk);
                    wr(1+mcell:2*mcell)= (wr(1+mcell:2*mcell)/scl_t) ;
                    wr(1+2*mcell:3*mcell)= (wr(1+2*mcell:3*mcell)/scl_p);

                    Wr = spdiags(wr,0,3*mactv,3*mactv);

                    mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;


%                     Wr = spdiags(wr,0,3*mactv,3*mactv);


%                     mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

                    diagJ = zeros(1,3*mactv);
                    for gg = 1 : ndata
                        diagJ = diagJ +  ( G{1}(gg,:)* S).^2;
                    end

                    diagA   = diagJ + beta(count)*spdiags( mof ,0)';
                    PreC    = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);

                    switch FLAG1

                        case 'SMOOTH_MOD'
                            g_MOF = ( aVRWs'*Wr*aVRWs * (invmod-mref) ) +...
                                    ( aVRWx'*Wr*aVRWx * (invmod) ) +...
                                    ( aVRWy'*Wr*aVRWy * (invmod) ) +...
                                    ( aVRWz'*Wr*aVRWz * (invmod) );

                            g = - S'*(Gtvec(G,Wd,(Gvec(G,Wd,m_uvw(aa,tt,pp)) - d))) - beta(count)*g_MOF;

                        case 'SMOOTH_MOD_DIF'
            %                 g = [- (magB - d) ; ...
            %             - sqrt( beta(count) ) * ( aVRWs * (invmod-mref) ) ;...
            %             - sqrt( beta(count) ) * ( aVRWx * (invmod-mref) ) ;...
            %             - sqrt( beta(count) ) * ( aVRWy * (invmod-mref) ) ;...
            %             - sqrt( beta(count) ) * ( aVRWz * (invmod-mref) ) ];
                    end 


                    %% Projected steepest descent
                    dm =zeros(3*mcell,1);
                    [dm,~,ncg] = CG_Lin( dm, G, Wd, beta(count) * mof , g, S, PreC, Pac );
            %         fprintf('CG iter %i \n',ncg);

                    %% Step length, line search                
                    tncg = tncg+ncg; % Record the number of CG iterations

                    temp = spdiags(Pac);
                    phi_out = phi_in;
                    m_temp = invmod;
                    % Combine active and inactive cells step if active bounds
                    if sum(temp)~=3*mactv
                        
                        rhs_a = ( speye(3*mactv) - Pac ) * (g);
                        
                        dm_i = max( abs( dm ) );
                        dm_a = max( abs(rhs_a) ); 

                        lowb = m_temp <= lowBvec;
                        uppb = m_temp >= uppBvec;

                        % Remove gradients going the wrong direction
                        rhs_a( lowb & rhs_a <=0) = 0;
                        rhs_a( uppb & rhs_a >=0) = 0;

                        if dm_i < dm_a
                            dm = dm + rhs_a * dm_i / dm_a /5 ;
                        else
                            dm = dm + rhs_a;
                        end

                    end
                    gamma = 2;
 
                    % Reduce step length in order to reduce phid

                    while (phi_out >= phi_in || gamma == 2) && gamma > 1e-4

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
                        Pac = spdiags((lowb==0).*(uppb==0),0,3*mactv,3*mactv);
                        aa = m_temp(1:mcell);
                        tt = m_temp(1+mcell:2*mcell);
                        pp = m_temp(1+2*mcell:3*mcell);

                        phi_MOF = (m_temp-mref)'* ( aVRWs'*Wr*aVRWs * (m_temp-mref) ) +...
                                  (m_temp)'*( aVRWx'*Wr*aVRWx * (m_temp) ) +...
                                  (m_temp)'*( aVRWy'*Wr*aVRWy * (m_temp) ) +...
                                  (m_temp)'*( aVRWz'*Wr*aVRWz * (m_temp) ) ;

                        phi_out = sum((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d).^2) + beta(count) * phi_MOF;
            %             phi_out = comp_phi(aa,tt,pp,MOF,beta(count));


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

                    fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,ggdm);
   
                    solves = solves + 1;


                end
                %% Save iteration and continue

                % Measure traffic around the lp corner
            %     if lp_count >= 1  
            %         temp = sum(abs(invmod) <= tresh_s);
            %         traffic_s(count) = abs(group_s(count) - temp) / group_s(count);
            %         
            %         dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
            %         temp = sum(abs(dmdx) <= tresh_xyz);
            % %         traffic_xyz(count) = abs(group_xyz(count) - temp) / group_xyz(count);
            %         
            % 
            %     end

                % Measure the update length
                if count==1 

                    rdm(count) =  1;
                    gdm(1) = norm(m_in - invmod);

                else

                    gdm(2) = norm(m_in - invmod);
                    rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);

                    gdm(1) = gdm(2);

                end

                aa = invmod(1:mcell);
                tt = invmod(1+mcell:2*mcell);
                pp = invmod(1+2*mcell:3*mcell);

                phi_d(count) = sum(( Gvec(G,Wd,m_uvw(aa,tt,pp)) - d ).^2);
                phi_m(count) =    (invmod-mref)'* ( aVRWs'*aVRWs * (invmod-mref) ) +...
                                  (invmod)'*( aVRWx'*aVRWx * (invmod) ) +...
                                  (invmod)'*( aVRWy'*aVRWy * (invmod) ) +...
                                  (invmod)'*( aVRWz'*aVRWz * (invmod) ) ;

                phi_MOF = (invmod-mref)'* ( aVRWs'*Wr*aVRWs * (invmod-mref) ) +...
                                  (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
                                  (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
                                  (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;

                phi(count) = phi_d(count) + beta(count) * phi_MOF;

                if count ~= 1

                    dphim(count) = abs(phi_m(count) - phi_m(count-1)) / phi_m(count) *100;

                end

                % Get truncated cells
                tcells = spdiags(Pac);


            %     fprintf('---------->\n')
                fprintf(' phi_d:\t %8.5e \n',phi_d(count))
                fprintf(' phi_m:\t %8.5e \n',phi_m(count))
                fprintf(' dphi_m:\t %8.5e \n',dphim(count))
            %     fprintf('Final Relative dm:\t %8.5e \n', rdm(count));
            %     fprintf('<----------\n')
            %     
            %     fprintf('Number of Inactive cells: %i\n',sum(tcells));
            %     fprintf('Number of CGS iterations: %i\n\n',ncg);

                % Get next beta
                [switcher,beta(count+1)] = cool_beta(beta(count),phi_d(count),dphim(count),target,switcher,0.25,2);


                % Right log file
            %     fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
            %     fprintf(fid,' \t %8.5e ',phi_d(count));
            %     fprintf(fid,' \t %8.5e ',sum( (aVRWs*invmod).^2 ) );
            %     fprintf(fid,' \t %8.5e ',sum( (aVRWx*invmod).^2 ));
            %     fprintf(fid,' \t %8.5e ',sum( (aVRWy*invmod).^2 ));
            %     fprintf(fid,' \t %8.5e ',sum( (aVRWz*invmod).^2 ));
            %     fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
            %     fprintf(fid,' \t %8.5e ',phi(count));
            %     fprintf(fid,' \t\t %i ',sum(tcells));
            %     fprintf(fid,' \t\t %i\n',ncg);


                model_out = X'*(m_uvw(aa,tt,pp));
            %      model_out(kron([1;1;1],nullcell)==0) = -100;
                % Create orthogonal magnetization vectors
            %     Mp =  model_out(1:mcell);%M + IWr * Esus * invmod;
            %     Mpp = Mp;
            %     Mpp(Mp<0) = 0;
            %     Mpm = Mp;
            %     Mpm(Mp>0) = 0;
            %     
            %     Ms =  model_out((1+mcell) : 2*mcell);
            %     Mt =  model_out((1+2*mcell) : 3*mcell);

                % Convert back to cartesian for plotting
                % z is flipped because conventionofcode is z positive down
            %     [mp,ms,mt] = azmdip_2_pst(HD,HI,mcell);

                Mxyz = model_out;
                Mx = Mxyz(1 : mcell);
                My = Mxyz(((1+mcell) : 2*mcell));
                Mz = Mxyz(((1+2*mcell) : 3*mcell));
                M = [Mx My Mz];

                % Create absolute magnetization
            %                 model_out(nullcell==0,:) = -100;

                Mamp = sqrt( Mx.^2 + My.^2 + Mz.^2 );
            %     Mrem = sqrt( Mpm.^2 + Ms.^2 + Mt.^2 );
            %     Mind = sqrt( Mpp.^2);


                pred_TMI = Gvec(G,speye(ndata),m_uvw(aa,tt,pp));
                write_MAG3D_TMI([work_dir dsep 'Tile' num2str(idx) '_MVI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
                save([work_dir dsep 'Tile' num2str(idx) '_MVI.fld'],'-ascii','M')
                save([work_dir dsep 'Tile' num2str(idx) '_MVI.amp'],'-ascii','Mamp')
            %     write_MAG3D_TMI([work_dir dsep 'Tile' num2str(idx) '_iter_' num2str(count) '.pre'],H,I,Dazm,obsx,obsy,obsz,(G*invmod).*wd,wd);
            end
            
%             save([out_dir dsep 'K_' num2str(kk) ' O_' num2str(oo) ' L_' num2str(ll) '.fld'],'-ascii','M')
%   
%         end
%     end
% end



% save([out_dir dsep 'Tile' num2str(idx) '_MVI.ind'],'-ascii','Mind')
% save([out_dir dsep 'Tile' num2str(idx) '_MVI.rem'],'-ascii','Mrem') 
write_MAG3D_TMI([out_dir dsep 'Tile' num2str(idx) '_MVI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
% leml = norm(invmod - mtrue,1);
% fprintf(fid,'End of lp inversion. Number of iterations: %i\n',count);
% fprintf(fid,'Final Number of CG iterations: %i\n',tncg);
% fprintf('Final Number of CG iterations: %i\n',tncg);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
% fclose(fid);

beta_out = beta(end);
phid_out = phi_d(end);
MOF_out = MOF;