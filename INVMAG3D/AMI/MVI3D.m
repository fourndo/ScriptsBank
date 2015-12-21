function [M,pred_TMI,beta_out,phid_out,MOF_out,switcher] = MVI3D(work_dir,meshfile,obsfile,mstart,mref,esus,chi_target,phid_in,MOF_in,alphas,beta,bounds,LP,t,delta_FLAG,delta_tresh,switcher,FLAG1,FLAG2,max_iter)
% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23


% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

% pct_cutoff = 50;

%% Create model vectors
% normalize effective susceptibility
if length(esus)==1
esus = ones(mcell,1)*esus;
end

mref = kron(ones(3,1), mref );
% 
% mstart = kron(ones(3,1), mstart );

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
esus = x*esus;
esus = esus  / norm(esus);
esus = kron(ones(3,1), esus );

% esus = X * esus;





Wj = spdiags(esus,0,3*mactv,3*mactv);

lowBvec = X * lowBvec;
uppBvec = X * uppBvec;



%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);

%% Depth weighting
wr = load([work_dir '\wr.dat']);


%% Load forward operator
TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
fprintf('Loading Sensitivity...\n');
load([work_dir '\Tx']);
load([work_dir '\Ty']);
load([work_dir '\Tz']);

Txyz = zeros(ndata,3*mactv);

% Do it in for loop to save memory
for ii = 1 : ndata
    
Txyz(ii,:) = TMI * [Tx(ii,:);Ty(ii,:);Tz(ii,:)];

end

clear Tx Ty Tz

[P,S,T] = azmdip_2_pst(Dazm,I,mactv);

% Create orthogonal forward operators
% Primary (inducing)
G = [Txyz * (H * P) Txyz * (H * S) Txyz * (H * T)];
clear Txyz P S T
% wd = abs(d)*0.05 + 0.05*std(d);
Wd = spdiags(1./wd,0,ndata,ndata);

% G = Wd * G * Esus * H ;

G = Wd * G * Wj;
d = Wd * d;



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

% phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phid_in;
phi_m       = invmod'*MOF_in*invmod;  
phi         = phid_in + beta*phi_m;

% Initiate active cell
Pac = speye(3*mactv);

% Message prompt
% logfile = [work_dir '\Log_TMVI.log'];
% fid = fopen(logfile,'w');
% fprintf(fid,'Starting lp inversion\n');
% fprintf(fid,'Starting misfit %e\n',phid_in);
% fprintf(fid,'Target misfit %e\n',target);
% fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
% fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
% fprintf(fid,'#cut cells \t # CG Iter\n');

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

dkdt = @(p,ep) (ep).^(1/(2*(2-p)));
traffic_s = 1;    % Measures the number of points/iteration passing the lp corner
traffic_xyz = 1;    % Measures the number of points/iteration passing the lp corner
group_s = 1;  % Measures the number of lower than the lp corner

while switcher ~= 3 && count ~= max_iter


    count=count+1;
    
    if switcher == 0     
        
        if delta_FLAG == 1
            
            delta_s(count,1:3) = delta_tresh(1)*3;
            delta_xyz(count) = delta_tresh(2)*3;
            
        else
            
            for pst = 1 : 3
            
                S = zeros(1,3);
                S(pst) = 1;

                % Create sub-space matrix and grab the right model
                S = kron(S,speye(mcell));

                m = S*invmod;

                delta_s(count,pst) = prctile(abs(m(m > 0)),delta_tresh(1))*3;
            
            end
            
            dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
            delta_xyz(count) = prctile(abs(dmdx(dmdx > 0)),delta_tresh(2))*3;
            
        end
%                     [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(invmod,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,2,2,1);
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_s(count),delta_xyz(count));

        

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum(G*temp.^2) / (temp'*MOF*temp) * 1e+3 ;

        end
        
        phi(count) = norm(G*invmod - d).^2 +...
                    invmod' * beta(count) * MOF * invmod;
%         tresh = dkdt(2,delta(count));
        tresh_s = delta_tresh(1);
        tresh_xyz = delta_tresh(2);
        
    elseif switcher == 1

        lp_count = lp_count+1;
        
%         if lp_count == 1
% %             target = target /10;
%             % Fix treshold for smallness term
%             remnt = invmod(mcell+1:end);
%             tresh_s = 1e-4;
% %             tresh_s = prctile(abs(remnt(remnt~=0)),pct_cutoff);%0.001;
%             delta_s(count) = tresh_s*2;%x_tresh.^2;
%             
%             
%             
%             % Fix treshold for smoothness term
%             %0.001;%
%             delta_xyz(count) = tresh_xyz*2;%x_tresh.^2;
%             
%             
%         end
            
        if traffic_s(end)*100 > 1 && switcher == 1 && lp_count~=1
              
                delta_s(count,:) = delta_s(count-1,:)/2;
                
        else

                delta_s(count,:) = delta_s(count-1,:);

        end
        
        if traffic_xyz(end)*100 > 1 && switcher == 1 && lp_count~=1

            delta_xyz(count) = delta_xyz(count-1)/2;  
            
        else

            delta_xyz(count) = delta_xyz(count-1);

        end
        
        if traffic_xyz(end)*100 <= 1 || traffic_s(end)*100 <= 1
            
            switcher = 2;
            
        end
          
        fprintf('# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_s(count,:),delta_xyz(count));

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

%     fprintf('\n# # # # # # # # #\n');
%     fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
    
    % Save previous model before GN step
    
    
    [invmod, ncg, Pac] = GN_PCG_solver( G, invmod, mref, nullcell, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );
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
    
    % Get truncated cells
    tcells = spdiags(Pac);
    
    
%     fprintf('---------->\n')
    fprintf(' phi_d:\t %8.5e \n',phi_d(count))
    fprintf(' phi_m:\t %8.5e \n',phi_m(count))
%     fprintf('Final Relative dm:\t %8.5e \n', rdm(count));
%     fprintf('<----------\n')
%     
%     fprintf('Number of Inactive cells: %i\n',sum(tcells));
%     fprintf('Number of CGS iterations: %i\n\n',ncg);

    % Get next beta
    [switcher,beta(count+1)] = cool_beta(beta(count),phi_d(count),rdm(count),target,switcher,0.25,1e-2);


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


     model_out = X'*(Wj*invmod);
     model_out(kron([1;1;1],nullcell)==0) = -100;
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
    [mp,ms,mt] = azmdip_2_pst(Dazm,I,mcell);
    
    Mxyz = [mp ms mt] *  model_out;
    Mx = Mxyz(1 : mcell);
    My = Mxyz(((1+mcell) : 2*mcell));
    Mz = Mxyz(((1+2*mcell) : 3*mcell));
    M = [Mx My Mz];
    
    % Create absolute magnetization
%                 model_out(nullcell==0,:) = -100;

    Mamp = sqrt( Mp.^2 + Ms.^2 + Mt.^2 );
    Mrem = sqrt( Mpm.^2 + Ms.^2 + Mt.^2 );
    Mind = sqrt( Mpp.^2);
    
%     save([work_dir '\Mvec_TMVI_iter_' num2str(count) '.fld'],'-ascii','M')
%     save([work_dir '\M_TMVI_iter_' num2str(count) '.amp'],'-ascii','Mamp')
%     save([work_dir '\M_TMVI_iter_' num2str(count) '.ind'],'-ascii','Mind')
%     save([work_dir '\M_TMVI_iter_' num2str(count) '.rem'],'-ascii','Mrem')           
    

end
          
% leml = norm(invmod - mtrue,1);
% fprintf(fid,'End of lp inversion. Number of iterations: %i\n',count);
% fprintf(fid,'Final Number of CG iterations: %i\n',tncg);
% fprintf('Final Number of CG iterations: %i\n',tncg);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
% fclose(fid);
pred_TMI = (G*invmod).*wd;
beta_out = beta(end);
phid_out = phi_d(end);
MOF_out = MOF;