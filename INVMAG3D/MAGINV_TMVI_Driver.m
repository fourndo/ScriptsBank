% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23

clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock';
inpfile   = 'MAGINV_TMVI.inp'; 
dsep = '\';

% [meshfile,obsfile,susfile,chi_target,alphas,beta,pvec,qvec,lvec,FLAG] = MAG3CM_read_inp([work_dir '\' inpfile]);
[meshfile,obsfile,topofile,mstart,mref,esus,chi_target,alphas,beta,bounds,norm_vec,FLAG1,FLAG2] = MAGINV_TMVI_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

pct_cutoff = 75;

%% Create model magnetization vectors
% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir '\' mref]);
    
    if length(mref) == mcell

        mref = kron([1;1;1],mref);

    end
    
else
    
    mref = ones(3*mcell,1)*mref;
    
end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir '\' mstart]);
    
    if length(mstart) == mcell
        
        mstart = kron([1;1;1],mstart);
        
    end
    
else
    
    mstart = ones(3*mcell,1)*mstart;
    
end

% Create or load reference model
if ischar(esus)==1
    
    esus = load([work_dir '\' esus]);
%     esus = esus / max(esus) + 1e-1;
%     esus = mcell / norm(esus) * esus;


else
    
    esus = ones(mcell,1)*esus;
    esus = esus / max(esus+1e-1)+ 1e-3;
    
end

% Need to create I/O for cell base weights
w = ones(4*mcell,1);

% Create bound vector
lowBvec = [ones(mcell,1) * bounds(1,1);ones(mcell,1) * bounds(2,1);ones(mcell,1) * bounds(3,1)];
uppBvec = [ones(mcell,1) * bounds(1,2);ones(mcell,1) * bounds(2,2);ones(mcell,1) * bounds(3,2)];

%% Initialize dynamic cells
% Create selector matrix for active cells
if isempty(topofile)==1
    
    nullcell = ones(nx*ny*nz,1);
    
else
    % Load topo
    topo = read_UBC_topo([work_dir dsep topofile]);
    [nullcell,tcell,~] = topocheck(xn,yn,zn,topo+1e-5);
    
end
x = spdiags(nullcell,0,mcell,mcell);
x = x(nullcell==1,:);
X = kron(speye(3),x);

mactv = sum(nullcell);

mstart = X * mstart;
mref = X * mref;

% Normalize effective susceptibility weight
esus(esus==-100) = 0;
esus = x * esus;
esus = x*esus;
wj = esus  / max(esus) + 1e-2;
% esus = kron(ones(3,1), esus / max(esus) );
% Wj = spdiags(esus,0,3*mactv,3*mactv);

lowBvec = X * lowBvec;
uppBvec = X * uppBvec;

%% Initialize selector matrix for the p,s,t components and bounds
Sp = kron([1 0 0],speye(sum(nullcell)));
Ss = kron([0 1 0],speye(sum(nullcell)));
St = kron([0 0 1],speye(sum(nullcell)));


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
    
    s= x * ones(mcell,1);
    t= x * ones(mcell,1);

%     s{2}= x * ones(mcell,1);
%     t{2}= x * ones(mcell,1);
%     
%     s{3}= x * ones(mcell,1);
%     t{3}= x * ones(mcell,1);
end


%% Load observation file (3C UBC-MAG format)
[H, HI, HD, MI, MD, dtype, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir dsep obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);
Wd   = spdiags(1./wd,0,ndata,ndata);



%% Depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
% save([work_dir '\wr.dat'],'-ascii','wr');



% if length(wr) == mcell
%         
%     wr = kron([1;1;1],wr);
%         
% end

% wr = x * wr;
% wr = esus.*wr;

% wr = kron([1;1;1],wr);

%% Load forward operator
[P,S,T] = azmdip_2_pst(HD,HI,mactv);
if dtype == 1
    [G,~,~] = MAGSEN_Func(work_dir,dsep,xn,yn,zn, H, HI, HD, MI, MD,dtype,...
        obsx, obsy, obsz,nullcell, 'SENS', 'G', 3, min(dx)/4);
    
    G = [G * (H * P) G * (H * S) G * (H * T)];
    
else
    
    [Tx,Ty,Tz] = MAGSEN_Func(work_dir,dsep,xn,yn,zn, H, HI, HD, MI, MD,dtype,...
        obsx, obsy, obsz,nullcell, 'SENS', 'TxTyTz', 3, min(dx)/4);
    
    G = [[Tx;Ty;Tz] * (H * P) [Tx;Ty;Tz] * (H * S) [Tx;Ty;Tz] * (H * T)];
    clear Tx Ty Tz

end
% Case cartesian coordinates


% Create orthogonal forward operators
% Primary (inducing)



% G = Wd * G * Esus * H ;

G = Wd * G;
d = Wd * d;


%% Create gradient matrices and corresponding volume vectors
%% Create gradient matrices and corresponding volume vectors
[~, Gx, Gy, Gz, V, Vx, Vy, Vz] = get_GRAD_op3D_SQUARE(dx,dy,dz,nullcell,x);
% [Ws, V ] = getWs3D(dx,dy,dz,X);

%wr = load([work_dir '\wr.dat']);
v = spdiags(V);
j = sum(G(:,1:mcell).^2,1);
j = (j'.^0.5)./v.^2;

wr = (j / max(j)).^0.5;


Ws =  V * spdiags(x * ( w(1:mcell) ).* wr ./wj ,0,mactv,mactv);
Wx =  V * spdiags(x * ( w(1+mcell:2*mcell) ).* wr ./wj ,0,mactv,mactv);
Wy =  V * spdiags(x * ( w(1+2*mcell:3*mcell) ).* wr ./wj ,0,mactv,mactv);
Wz =  V * spdiags(x * ( w(1+3*mcell:4*mcell) ).* wr ./wj ,0,mactv,mactv);

mactv = sum(nullcell);

%% START INVERSION
    
target = chi_target * ndata;     % Target misifit

% Compute total objective function
objfunc = @(m,phim,b) sum( ( G * m - d ).^2 ) + (m)' * b * phim * (m);

invmod      =  mstart ;

phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phi_init;
phi_m       = [];  
mref_perp = [zeros(3*mcell,1)];
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
eps_FLAG = 0;
dphim = 100;
while switcher ~= 3 


    count=count+1;
    
    if switcher == 0     
        
        delta_p(count,1:3) = 1e-1;%delta_tresh(1)*3;
        delta_q(count,1:3) = 1e-1;%delta_tresh(2)*3;
        
%                     [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(invmod,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,2,2,1);
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_p(count,:),delta_q(count,:));

        

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum((G*temp).^2) / (temp'*MOF*temp) * 1e+3 ;

        end
        
        phi(count) = norm(G*invmod - d).^2 +...
                    invmod' * beta(count) * MOF * invmod;
%         tresh = dkdt(2,delta(count));
        tresh_s = 1;
        tresh_xyz = 1;
    if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum((G*temp).^2) / (temp'*MOF*temp) * 1e+5 ;

        end
        
        phi(count) = norm(G*invmod - d).^2 +...
                    invmod' * beta(count) * MOF * invmod;
%         tresh = dkdt(2,delta(count));

        
    else

        lp_count = lp_count+1;
        
        
        if lp_count == 1 && eps_FLAG==0

            for pst = 1 : 3

                Clp = zeros(1,3);
                Clp(pst) = 1;

                % Create sub-space matrix and grab the right model
                Clp = kron(Clp,speye(mactv));

                m = Clp*invmod;

                [pp,qq] = get_eps(m,10,Gx,Gy,Gz);
                
                eps_p(pst) = pp;
                eps_q(pst) = qq;
                
            end
            
            
        end
        
        for pst = 1 : 3
            
            if delta_p(end,pst)> eps_p(pst)

                delta_p(count,pst) = delta_p(count-1,pst)*.5;

                if delta_p(count,pst) < eps_p(pst)

                    delta_p(count,pst) = eps_p(pst);

                end

            else

                    delta_p(count,:) = delta_p(count-1,:);

            end

            if delta_q(end,pst)> eps_q(pst)

                delta_q(count,pst) = delta_q(count-1,pst)*.5;

                if delta_q(count,pst) < eps_q(pst)

                    delta_q(count,pst) = eps_q(pst);

                end

            else

                delta_q(count,pst) = delta_q(count-1,pst);

            end


        end 
        
        if dphim(end)  < 5%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1
            
            fprintf('\n# # ADJUST BETA # #\n');
            switcher = 2;

        else

            fprintf('\n# # LP-LQ ITER# #\n');

        end
          
        fprintf('\n# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p(count,:),delta_q(count,:));

%         tresh = dkdt(LP(:,1),delta(count));

        
    end
    
    m_in = invmod;
    dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
    
    group_s(count) = sum(abs(m_in) <= tresh_s);
    group_xyz(count) = sum(abs(dmdx) <= tresh_xyz);
    
    %% Pre-conditionner
    diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
    PreC     = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);

    %% Gauss-Newton steps

    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
%     fprintf('eps_q: \t %8.5e \t eps_q*: \t %8.5e\n',delta_p(count),eps_p)
%     fprintf('eps_p: \t %8.5e \t eps_p*: \t %8.5e\n',delta_q(count),eps_q)
    
    % Save previous model before GN step
    
    
    [invmod, ncg, Pac] = GN_PCG_solver( G, invmod, mref, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );
    tncg = tncg + ncg;
    
    %% Save iteration and continue

    % Measure traffic around the lp corner
    if lp_count >= 1  
        temp = sum(abs(invmod) <= tresh_s);
        traffic_s(count) = abs(group_s(count) - temp) / group_s(count);
        
        dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
        temp = sum(abs(dmdx) <= tresh_xyz);
        traffic_xyz(count) = abs(group_xyz(count) - temp) / group_xyz(count);
        

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
    
    phi_d(count) = sum(( G*invmod - d ).^2);
    phi_m(count) = (invmod)'*(MOF)*(invmod);
    phi(count) = objfunc(invmod,MOF,beta(count));
    
    if count ~= 1
        
        dphim(count) = abs(phi_m(count) - phi_m(count-1)) / phi_m(count) *100;
        
    end
    
%     fprintf('---------->\n')
    fprintf(' phi_d:\t %8.5e \n',phi_d(count))
    fprintf(' phi_m:\t %8.5e \n',phi_m(count))
    fprintf(' dphi_m:\t %8.5e \n',dphim(count))
    
    
    fprintf('---------->')
    fprintf(' misfit:\t %8.5e ',phi_d(count))
    fprintf('Final Relative dm:\t %8.5e ', rdm(count));
    fprintf('<----------\n')
    
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
    fprintf(fid,' \t\t %i\n',ncg);


     model_out = X'*(invmod);
%      model_out(kron([1;1;1],nullcell)==0) = -100;
    % Create orthogonal magnetization vectors
    Mp =  model_out(1:mcell);%M + IWr * Esus * invmod;
    Mpp = Mp;
    Mpp(Mp<0) = 0;
    Mpm = Mp;
    Mpm(Mp>0) = 0;
    
    Ms =  model_out((1+mcell) : 2*mcell);
    Mt =  model_out((1+2*mcell) : 3*mcell);
    
    % Convert back to cartesian for plotting
    % z is flipped because conventionofcode is z positive down
    
    Mxyz = [P S T] *  model_out;
    Mx = Mxyz(1 : mcell);
    My = Mxyz(((1+mcell) : 2*mcell));
    Mz = Mxyz(((1+2*mcell) : 3*mcell));
    M = [Mx My Mz];
    
    % Create absolute magnetization
%                 model_out(nullcell==0,:) = -100;

    Mamp = sqrt( Mp.^2 + Ms.^2 + Mt.^2 );
    Mrem = sqrt( Mpm.^2 + Ms.^2 + Mt.^2 );
    Mind = sqrt( Mpp.^2);
    
    save([work_dir '\Mvec_MVI_m.fld'],'-ascii','M')
    save([work_dir '\M_MVI_m.amp'],'-ascii','Mamp')
    save([work_dir '\M_MVI_m.ind'],'-ascii','Mind')
    save([work_dir '\M_MVI_m.rem'],'-ascii','Mrem')
    
%     if dtype == 1
%         write_MAG3D_TMI([work_dir '\MVI_d.pre'],H,I,Dazm,obsx,obsy,obsz,(G*invmod).*wd,wd)
%     else
%         write_MAG3D_3C([work_dir '\' obsfile(1:end-4) '_3C.obs'],H,HI,HD,...
%         obsx,obsy,obsz,(G*invmod).*wd,wd);
%     end

end
          
% leml = norm(invmod - mtrue,1);
fprintf(fid,'End of lp inversion. Number of iterations: %i\n',count);
fprintf(fid,'Final Number of CG iterations: %i\n',tncg);
fprintf('Final Number of CG iterations: %i\n',tncg);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
fclose(fid);