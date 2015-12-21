% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23

clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Triple_Block_lined';
inpfile   = 'MAGINV_TMVI.inp'; 

% [meshfile,obsfile,susfile,chi_target,alphas,beta,pvec,qvec,lvec,FLAG] = MAG3CM_read_inp([work_dir '\' inpfile]);
[meshfile,obsfile,mstart,mref,esus,chi_target,alphas,beta,bounds,norm_vec,FLAG1,FLAG2] = MAGINV_TMVI_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

pct_cutoff = 75;

eps_s = 2e-3;
eps_x = 1e-4;

delta_s = eps_s;
delta_xyz = eps_x;

%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);
Wd   = spdiags(1./wd,0,ndata,ndata);

%% Create model magnetization vectors
% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir '\' mref]);
    
    if length(mref) == mcell

        mref = kron([1;1;1],mref);

    end
    
else
    
    mref = [ones(mcell,1)*mref;ones(mcell,1)*D*pi/180;ones(mcell,1)*I*pi/180];
    
end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir '\' mstart]);
    
    if length(mstart) == mcell
        
        mstart = kron([1;1;1],mstart);
        
    end
    
else
    
    mstart = [ones(mcell,1)*mstart;ones(mcell,1)*D*pi/180;ones(mcell,1)*I*pi/180];
    
end

% Create or load reference model
if ischar(esus)==1
    
    esus = load([work_dir '\' esus]);
%     esus = esus / max(esus) + 1e-1;
%     esus = mcell / norm(esus) * esus;


else
    
    esus = ones(mcell,1)*esus;
%     esus = esus / max(esus+1e-1)+ 1e-3;
    
end

% Need to create I/O for cell base weights
w = ones(4*mcell,1);

% Create bound vector
lowBvec = [ones(mcell,1) * bounds(1,1);ones(mcell,1) * bounds(2,1);ones(mcell,1) * bounds(3,1)];
uppBvec = [ones(mcell,1) * bounds(1,2);ones(mcell,1) * bounds(2,2);ones(mcell,1) * bounds(3,2)];

%% Initialize dynamic cells
% Create selector matrix for active cells
load([work_dir '\nullcell.dat']);
x = spdiags(nullcell,0,mcell,mcell);
x = x(nullcell==1,:);
X = kron(speye(3),x);

mactv = sum(nullcell);

mstart = X * mstart;
mref = X * mref;

% Normalize effective susceptibility weight
esus(esus==-100) = 0;
esus = x * esus;
esus = esus  / norm(esus) * norm(ones(mactv,1));
esus = [esus;ones(2*mactv,1)];
% esus = kron(ones(3,1), esus / max(esus) );
Wj = spdiags(esus,0,3*mactv,3*mactv);

lowBvec = X * lowBvec;
uppBvec = X * uppBvec;

%% Initialize selector matrix for the p,s,t components and bounds
% Sp = kron([1 0 0],speye(sum(nullcell)));
% Ss = kron([0 1 0],speye(sum(nullcell)));
% St = kron([0 0 1],speye(sum(nullcell)));


%% Generate selection and transition vectors for the different lp zones
if ischar(norm_vec)==1
    
    % NEED TO FIX THIS PART FOR p, s and t COMPONENTS
    % Assume that the input norm_vec is 3*mcell-by-5
    lpmat = load([work_dir '\' norm_vec]);
    [s,LP] = find_zones(lpmat);

    if length(s) == mcell
        
         s = kron([1;1;1],s);

    end
        
    % Smooth out the regions with 8-point averager
    % Power determine the transition length
    A = get_AVG_8pt(dx,dy,dz);
    A = kron(speye(3), A);
    
    trans = A*A*(A*s);
    t = X * trans;
%     t(t>0) = sqrt(t(t>0));
%     LP = X * LP;
    
    s = trans==1;
    s = X * s;
    
else
    
    LP = norm_vec(1,:);
%     LP{2} = norm_vec(2,:);
%     LP{3} = norm_vec(3,:);
    
    s= X * ones(3*mcell,1);
    t= X * ones(3*mcell,1);

%     s{2}= x * ones(mcell,1);
%     t{2}= x * ones(mcell,1);
%     
%     s{3}= x * ones(mcell,1);
%     t{3}= x * ones(mcell,1);
end


%% Depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
% save([work_dir '\wr.dat'],'-ascii','wr');

wr = load([work_dir '\wr.dat']);

% if length(wr) == mcell
%         
%     wr = kron([1;1;1],wr);
%         
% end

% wr = x * wr;
% wr = esus.*wr;

% wr = kron([1;1;1],wr);

%% Load forward operator

load([work_dir '\G']);

% load([work_dir '\Tx']);
% load([work_dir '\Ty']);
% load([work_dir '\Tz']);
% 
% TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
% fprintf('Loading Sensitivity...\n');
% % load([work_dir '\Tx']);
% % load([work_dir '\Ty']);
% % load([work_dir '\Tz']);
% 
% G = zeros(ndata,3*mactv);
% 
% % Do it in for loop to save memory
% for ii = 1 : ndata
%     
% G(ii,:) = TMI * [Tx(ii,:);Ty(ii,:);Tz(ii,:)] * H;
% 
% end
% clear Tx Ty Tz

G = Wd * G * Wj * H;
d = Wd * d;

% Case sperical
m = @(a,t,p) [a.*cos(p).*cos(t);...
    a.*cos(p).*sin(t);...
    a.*sin(p)];

S = @(a,t,p)[spdiags(cos(p).*cos(t),0,mactv,mactv) spdiags(-a.*cos(p).*sin(t),0,mactv,mactv) spdiags(-a.*sin(p).*cos(t),0,mactv,mactv);
             spdiags(cos(p).*sin(t),0,mactv,mactv) spdiags(a.*cos(p).*cos(t),0,mactv,mactv) spdiags(-a.*sin(p).*sin(t),0,mactv,mactv);
             spdiags(sin(p),0,mactv,mactv) sparse(mactv,mactv) spdiags(a.*cos(p),0,mactv,mactv)];


%% Temporary added - Scale the sensitivity
% aa = mstart(1:mactv);
% tt = mstart(1+mactv:2*mactv);
% pp = mstart(1+2*mactv:3*mactv);
% 
% J   = G * S(aa,tt,pp) ;
% wj = (sum(abs(J)));
% wj = sqrt(wj / max(wj));
% % sc_v = 1 / max(sum(abs(J(:,1+mcell:2*mcell))));
% % sc_w = 1 / max(sum(abs(J(:,1+2*mcell:3*mcell))));
% 
% % sc = [ones(mactv,1) * sc_u ; ones(mactv,1) * sc_v ; ones(mactv,1) * sc_w];
% % sc = spdiags(sc_u,0,3*mactv,3*mactv);
% 
% % Wj = Wj * spdiags(wj',0,3*mactv,3*mactv);
% 
% G =G * spdiags(wj',0,3*mactv,3*mactv);

%% Create gradient matrices and corresponding volume vectors
[~, Gx, Gy, Gz, V, Vx, Vy, Vz] = get_GRAD_op3D_SQUARE(dx,dy,dz,nullcell,x);
% [Ws, V ] = getWs3D(dx,dy,dz,X);

Ws =  V * spdiags(x * ( w(1:mcell) .* wr ) ,0,mactv,mactv);
Wx =  Vx * spdiags(x * ( w(1+mcell:2*mcell) .* wr ) ,0,mactv,mactv);
Wy =  Vy * spdiags(x * ( w(1+2*mcell:3*mcell) .* wr ) ,0,mactv,mactv);
Wz =  Vz * spdiags(x * ( w(1+3*mcell:4*mcell) .* wr ) ,0,mactv,mactv);

mactv = sum(nullcell);
%% START INVERSION
    
target = chi_target * ndata;     % Target misifit

% Compute total objective function
objfunc = @(m,phim,b) sum( ( G * m - d ).^2 ) + (m)' * b * phim * (m);

invmod      =  mstart ;

phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phi_init;
% mref_perp = [zeros(3*mcell,1)];
% Initiate active cell
Pac = speye(3*mactv);

% Message prompt
logfile = [work_dir '\Log_TMVI.log'];
fid = fopen(logfile,'w');
fprintf(fid,'Starting lp inversion\n');
fprintf(fid,'Starting misfit %e\n',phi_init);
fprintf(fid,'Target misfit %e\n',target);
fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
fprintf(fid,'#cut cells \t # CG Iter\n');

count= 0; % Initiate iteration count 
tncg = 0; % Compute total number of CG iterations for the whole inversion
% leml = 0;
phi = [];
phi_m = [];
lp_count = 0;
switcher = 0; % Switcher defines the different modes of the inversion
% switcher 0: Run the usual l2-l2 inversion, beta decreases by 2
% switcher 1: Run the lp-lq inversion, beta decreases by 0.8
% swircher 2: End on inversion - model update < treshold

dkdt = @(p,ep) (ep).^(1/(2*(2-p)));
traffic_s = 1;    % Measures the number of points/iteration passing the lp corner
traffic_xyz = 1;    % Measures the number of points/iteration passing the lp corner
group_s = 1;  % Measures the number of lower than the lp corner

while switcher ~= 3 


    count=count+1;
    
    aa = invmod(1:mactv);
    tt = invmod(1+mactv:2*mactv);
    pp = invmod(1+2*mactv:3*mactv);
        
    if switcher == 0     
        
%         delta_s(count) = (prctile(abs(invmod(invmod > 0)),pct_cutoff))*3;
%         delta_xyz(count) = (prctile(abs(invmod(invmod > 0)),pct_cutoff))*3;


        
%                     [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(invmod,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,2,2,1);
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,eps_s,eps_x);

        

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum((G*temp).^2) / (temp'*MOF*temp) * 1e+1 ;

        end
        
        phi_m(count) = ((aVRWs * (invmod-mref) )' * (aVRWs * (invmod-mref) ))+...
            (( (aVRWx) * invmod )' * ( (aVRWx) * ( invmod ) ))+...
            (( (aVRWy) * invmod )' * ( (aVRWy) * ( invmod ) ))+...
            (( (aVRWz) * invmod )' * ( (aVRWz) * ( invmod ) ));
        
        phi(count) = sum((G*m(aa,tt,pp) - d).^2) +...
                    beta(count) * phi_m(count);
                
%         tresh = dkdt(2,delta(count));
        tresh_s = eps_s;
%         tresh_xyz = eps_x;
    else

        lp_count = lp_count+1;
        

        if traffic_s(end)*100 <= 1 && lp_count > 4
            
            switcher = 2;
            
        end
          
        fprintf('\n# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,eps_s,eps_x);
                                                         
%         tresh = dkdt(LP(:,1),delta(count));

        
    end
    
    m_in = invmod;
   
    group_s(count) = sum(abs(m_in) <= tresh_s);


    %% Gauss-Newton steps
    
    fprintf('\n# ## # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
     
    % Save current model
%     m_in = invmod;
%     dmdx = sqrt( (Wx * invmod).^2 + (Wy * invmod).^2 + (Wz * invmod).^2 );
    
    tncg = 0;
    ggdm = 1;       % Mesure the relative change in rel|dm|
    ddm = [1 1];    % Mesure of change in model update |dm|
    solves = 1;
     
    phi_in = phi(end);
    
    while solves < 5 && ggdm > 1e-5
 
        aa = invmod(1:mactv);
        tt = invmod(1+mactv:2*mactv);
        pp = invmod(1+2*mactv:3*mactv);
        
        J   =  G * S(aa,tt,pp);

%         sc = spdiags([ones(mactv,1);ones(mactv,1)/max(aa);ones(mactv,1)/max(aa)],0,3*mactv,3*mactv);
        
        A = [ J ;...
        sqrt( beta(count) ) * aVRWs ;...
        sqrt( beta(count) ) * aVRWx ;...
        sqrt( beta(count) ) * aVRWy ;...
        sqrt( beta(count) ) * aVRWz ];

        diagA   = sum((J).^2,1) + beta(count)*spdiags(MOF,0)';
%         diagA = ones(3*mactv,1);
        PreC    = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);
        
        switch FLAG1

            case 'SMOOTH_MOD'
                g = [- (G*m(aa,tt,pp) - d) ; ...
            - sqrt( beta(count) ) * ( aVRWs * (invmod-mref) ) ;...
            - sqrt( beta(count) ) * ( aVRWx * (invmod) ) ;...
            - sqrt( beta(count) ) * ( aVRWy * (invmod) ) ;...
            - sqrt( beta(count) ) * ( aVRWz * (invmod) ) ];

            case 'SMOOTH_MOD_DIF'
                g = [- (magB - d) ; ...
            - sqrt( beta(count) ) * ( aVRWs * (invmod-mref) ) ;...
            - sqrt( beta(count) ) * ( aVRWx * (invmod-mref) ) ;...
            - sqrt( beta(count) ) * ( aVRWy * (invmod-mref) ) ;...
            - sqrt( beta(count) ) * ( aVRWz * (invmod-mref) ) ];
        end 
        

        %% Projected steepest descent
        dm = zeros(3*mactv,1);
        [dm,~,ncg] = PCGLSQ( dm, A , g, PreC, Pac);
%         fprintf('CG iter %i \n',ncg);
        
        %% Step length, line search                
        tncg = tncg+ncg; % Record the number of CG iterations

        temp = spdiags(Pac);
        
        % Combine active and inactive cells step if active bounds
        if sum(temp)~=3*mactv
            
            rhs_a = ( speye(3*mactv) - Pac ) * (A'*g);
            dm_i = max( abs( dm ) );
            dm_a = max( abs(rhs_a) );                
            dm = dm + rhs_a * dm_i / dm_a /10 ;

        end
        gamma = 2;
        
        % Reduce step length in order to reduce phid
        phi_out = phi_in;
%         m_temp = invmod;
        phi_out = [phi_in phi_in]; 
        while (phi_out(1) > phi_out(2) || gamma == 2) && gamma > 1e-4

            phi_out(2) = phi_out(1);

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
            
            aa = m_temp(1:mactv);
            tt = m_temp(1+mactv:2*mactv);
            pp = m_temp(1+2*mactv:3*mactv);
            
            d_pre = G*m(aa,tt,pp) ;
            
            phi_temp = ((aVRWs * (m_temp-mref) )' * (aVRWs * (m_temp-mref) ))+...
            (( (aVRWx) * m_temp )' * ( (aVRWx) * ( m_temp ) ))+...
            (( (aVRWy) * m_temp )' * ( (aVRWy) * ( m_temp ) ))+...
            (( (aVRWz) * m_temp )' * ( (aVRWz) * ( m_temp ) ));
        
            phi_out(1) = sum((d_pre - d).^2) + beta(count) * phi_temp;
            
%             phi_out = comp_phi(aa,tt,pp,MOF,beta(count));


        end

        phi_in = phi_out(1);
        
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
    if lp_count >= 1  
        temp = sum(abs(invmod) <= tresh_s);
        traffic_s(count) = abs(group_s(count) - temp) / group_s(count);
        
%         dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
%         temp = sum(abs(dmdx) <= tresh_xyz);
%         traffic_xyz(count) = abs(group_xyz(count) - temp) / group_xyz(count);
        

    end
    
    % Measure the update length
    if count==1 
        
        rdm(count) =  1;
        gdm(1) = norm(m_in - invmod);
        
    else
        
        gdm(2) = norm(m_in - invmod);
        rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);
    
        gdm(1) = gdm(2);
        
    end
    
    phi_d(count) = sum(( G*m(aa,tt,pp) - d ).^2);
%     phi_m(count) = (invmod-mref)'*(MOF)*(invmod-mref);
    phi_m(count) = ((aVRWs * (invmod-mref) )' * (aVRWs * (invmod-mref) ))+...
            (( (aVRWx) * invmod )' * ( (aVRWx) * ( invmod ) ))+...
            (( (aVRWy) * invmod )' * ( (aVRWy) * ( invmod ) ))+...
            (( (aVRWz) * invmod )' * ( (aVRWz) * ( invmod ) ));
    
%     phi(count) = objfunc(invmod,MOF,beta(count));
    phi(count) = phi_d(count) + beta(count) * phi_m(count);
    
    % Get truncated cells
    tcells = spdiags(Pac);
    
    
    fprintf('---------->')
    fprintf(' misfit:\t %8.5e ',phi_d(count))
    fprintf('Final Relative dm:\t %8.5e ', rdm(count));
    fprintf('<----------\n')
    
    fprintf('Number of Inactive cells: %i\n',sum(tcells));
    fprintf('Number of CGS iterations: %i\n\n',ncg);

    % Get next beta
    [switcher,beta(count+1)] = cool_beta(beta(count),phi_d(count),rdm(count),target,switcher,0.1,1e-2);


    % Right log file
    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf(fid,' \t %8.5e ',sum( (aVRWs*invmod).^2 ) );
    fprintf(fid,' \t %8.5e ',sum( (aVRWx*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWy*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWz*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
    fprintf(fid,' \t %8.5e ',phi(count));
    fprintf(fid,' \t\t %i ',sum(tcells));
    fprintf(fid,' \t\t %i\n',ncg);
  
    model_out = X'*(Wj*invmod);
    
    aa = model_out(1:mcell);
    tt = model_out(1+mcell:2*mcell);
    pp = model_out(1+2*mcell:3*mcell);
     
    M = m(aa,tt,pp);

    % Dot product with unit vector to get induced component
    M_ind = spdiags(kron(ones(3,1),aa),0,3*mcell,3*mcell) * (m(ones(mcell,1),tt,pp) .*...
        m(ones(mcell,1),ones(mcell,1)*D*pi/180,ones(mcell,1)*I*pi/180));
    
%     M_ind(aa<=0) = 0;
    
    M_rem = M - M_ind;
    
    % Reshape and save
    M = reshape(M,mcell,3);
    M_ind = reshape(M_ind,mcell,3);
    M_rem = reshape(M_rem,mcell,3);
    
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
%     
%     % Convert back to cartesian for plotting
%     % z is flipped because conventionofcode is z positive down
% %     [mp,ms,mt] = azmdip_2_pst(Dazm,I,mcell);
%     
%     Mxyz = [mp ms mt] *  model_out;
%     Mx = Mxyz(1 : mcell);
%     My = Mxyz(((1+mcell) : 2*mcell));
%     Mz = Mxyz(((1+2*mcell) : 3*mcell));
%     M = [Mx My Mz];
%     
%     % Create absolute magnetization
% %                 model_out(nullcell==0,:) = -100;
% 
%     Mamp = sqrt( Mp.^2 + Ms.^2 + Mt.^2 );
%     Mrem = sqrt( Mpm.^2 + Ms.^2 + Mt.^2 );
%     Mind = sqrt( Mpp.^2);
%     
    save([work_dir '\MVI.fld'],'-ascii','M')
    save([work_dir '\MVI.amp'],'-ascii','aa')
    save([work_dir '\MVI.ind'],'-ascii','M_ind')
    save([work_dir '\MVI.rem'],'-ascii','M_rem')           
    write_MAG3D_TMI([work_dir '\MVI.pre'],H,I,Dazm,obsx,obsy,obsz,d_pre.*wd,wd)

end
          
% leml = norm(invmod - mtrue,1);
fprintf(fid,'End of lp inversion. Number of iterations: %i\n',count);
fprintf(fid,'Final Number of CG iterations: %i\n',tncg);
fprintf('Final Number of CG iterations: %i\n',tncg);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
fclose(fid);