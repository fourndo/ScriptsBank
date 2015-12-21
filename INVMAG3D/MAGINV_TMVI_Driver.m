% Magnetic Vector Invertion
% Written by: D Fournier 
% Last update: 2014/07/23

clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\Tile_AMI\Tile1';
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
esus = esus  / norm(esus);
esus = kron(ones(3,1), esus );
% esus = kron(ones(3,1), esus / max(esus) );
Wj = spdiags(esus,0,3*mactv,3*mactv);

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
    
    s= X * ones(3*mcell,1);
    t= X * ones(3*mcell,1);

%     s{2}= x * ones(mcell,1);
%     t{2}= x * ones(mcell,1);
%     
%     s{3}= x * ones(mcell,1);
%     t{3}= x * ones(mcell,1);
end


%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);
Wd   = spdiags(1./wd,0,ndata,ndata);



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


load([work_dir '\Gp']);
load([work_dir '\Gs']);
load([work_dir '\Gt']);

G = [Gp Gs Gt];
clear Gp Gs Gt

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

% Ws = kron(speye(3), Ws);
% Wx = kron(speye(3), Wx);
% Wy = kron(speye(3), Wy);
% Wz = kron(speye(3), Wz);
% Vx = kron(speye(3), Vx);
% Vy = kron(speye(3), Vy);
% Vz = kron(speye(3), Vz);
% V = kron( speye(3), V );
% Gx = kron(speye(3), Gx);
% Gy = kron(speye(3), Gy);
% Gz = kron(speye(3), Gz);

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

while switcher ~= 3 


    count=count+1;
    
    if switcher == 0     
        
        delta_s(count) = (prctile(abs(invmod(invmod > 0)),pct_cutoff))*3;
        delta_xyz(count) = (prctile(abs(invmod(invmod > 0)),pct_cutoff))*3;
        
%                     [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(invmod,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,2,2,1);
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_s(count),delta_xyz(count));

        

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            beta = sum((G*temp).^2) / (temp'*MOF*temp) * 1e+3 ;

        end
        
        phi(count) = norm(G*invmod - d).^2 +...
                    invmod' * beta(count) * MOF * invmod;
%         tresh = dkdt(2,delta(count));
        tresh_s = 1;
        tresh_xyz = 1;
    else

        lp_count = lp_count+1;
        
        if lp_count == 1
%             target = target /10;
            % Fix treshold for smallness term
            remnt = invmod(mcell+1:end);
            tresh_s = 1e-4;
%             tresh_s = prctile(abs(remnt(remnt~=0)),pct_cutoff);%0.001;
            delta_s(count) = tresh_s*2;%x_tresh.^2;
            
            
            
            % Fix treshold for smoothness term
    dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
            tresh_xyz = prctile(abs(dmdx(dmdx > 0)),pct_cutoff);%0.001;%
            delta_xyz(count) = tresh_xyz*2;%x_tresh.^2;
            
            
        end
            
        if traffic_s(end)*100 > 1 && switcher == 1
            
            
            
            figure(1)
            if lp_count==1
            
                subplot(1,2,1) 
                %%%% TEMPORARY CHANGE
%                 target = target / 4;
            else
                subplot(1,2,2)   
                delta_s(count) = delta_s(count-1)/2;
            end
            xx = sort(abs(invmod));
            mm = xx./(xx.^2 + delta_s(end).^2).^(1-LP(1)/2);
            [n, xout] =hist(abs(invmod),100); hold off
            [h_plot,h1,h2] = plotyy(xout,n,xx,mm/max(mm),'bar','plot');  
            hold(h_plot(1),'on');
            hold(h_plot(2),'on');
            set(h2,'LineWidth',2)
                
            ylim(h_plot(1),[0 mactv]);
            xlim(h_plot(1),[0 max(abs(invmod))])
            xlim(h_plot(2),[0 max(abs(invmod))])

            axis(h_plot(1),'square')
            axis(h_plot(2),'square')
    %                 [n, xout] =hist(h_plot(1),abs(invmod),100);
            set(h1,'barwidth', 1, 'basevalue', 1,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.5);

            plot(h_plot(2),[tresh_s tresh_s],[0 max(n)],'r--','LineWidth',2)

            set(h_plot(1),'yscale','log')

            hold off

            
        else

            delta_s(count) = delta_s(count-1);

        end
        
        if traffic_xyz(end)*100 > 1 && switcher == 1
            
            
            
    dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );
            
            figure(2)
            if lp_count==1
            
                subplot(1,2,1) 
            else
                subplot(1,2,2) 
                delta_xyz(count) = delta_xyz(count-1)/2;  
            end
            xx = sort(abs(dmdx));
            mm = xx./(xx.^2 + delta_xyz(end).^2).^(1-LP(2)/2);
            [n, xout] =hist(abs(dmdx),100); hold off
            [h_plot,h1,h2] = plotyy(xout,n,xx,mm/max(mm),'bar','plot');  
            hold(h_plot(1),'on');
            hold(h_plot(2),'on');
            set(h2,'LineWidth',2)
                
            ylim(h_plot(1),[0 1e+4]);
            xlim(h_plot(1),[0 max(dmdx)])
            xlim(h_plot(2),[0 max(dmdx)])
    
            axis(h_plot(1),'square')
            axis(h_plot(2),'square')
    %                 [n, xout] =hist(h_plot(1),abs(invmod),100);
            set(h1,'barwidth', 1, 'basevalue', 1,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.5);

            plot(h_plot(2),[tresh_xyz tresh_xyz],[0 max(n)],'r--','LineWidth',2)

            set(h_plot(1),'yscale','log')

            hold off

            
        else

            delta_xyz(count) = delta_xyz(count-1);

        end
        
        if traffic_xyz(end)*100 <= 1 || traffic_s(end)*100 <= 1
            
            switcher = 2;
            
        end
          
        fprintf('\n# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_s(count),delta_xyz(count));

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
    
    save([work_dir '\Mvec_TMVI_iter_' num2str(count) '.fld'],'-ascii','M')
    save([work_dir '\M_TMVI_iter_' num2str(count) '.amp'],'-ascii','Mamp')
    save([work_dir '\M_TMVI_iter_' num2str(count) '.ind'],'-ascii','Mind')
    save([work_dir '\M_TMVI_iter_' num2str(count) '.rem'],'-ascii','Mrem')           
    write_MAG3D_TMI([work_dir '\TMVI_iter_' num2str(count) '.pre'],H,I,Dazm,obsx,obsy,obsz,(G*invmod).*wd,wd)

end
          
% leml = norm(invmod - mtrue,1);
fprintf(fid,'End of lp inversion. Number of iterations: %i\n',count);
fprintf(fid,'Final Number of CG iterations: %i\n',tncg);
fprintf('Final Number of CG iterations: %i\n',tncg);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
fclose(fid);