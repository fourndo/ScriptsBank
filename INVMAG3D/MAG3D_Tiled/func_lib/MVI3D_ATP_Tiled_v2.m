function [M,pred_TMI,beta_out,phid_out,switcher] = MVI3D_ATP_Tiled(work_dir,out_dir,dsep,idx,xn,yn,zn,H, HI, HD, obsx, obsy, obsz, G, d, wd,mstart,mref,esus,chi_target,alphas,beta,bounds,LP,t,eps_FLAG,eps_tresh,FLAG1,FLAG2,max_iter,ROT)
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

scl_t = 0.003;
scl_p = 0.003;
% 
% for ll = 1
%     
%     for oo = 1
%         
%         for kk = 0 
            
switcher = 0;

for pst = 1 : 3
    if eps_FLAG == 1
        eps_q(pst) = eps_tresh(2);
        eps_p(pst) = eps_tresh(1);

    else

        eps_q(pst) = 1e-2;
        eps_p(pst) = 1e-2;
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

t = x*t;
%% Load observation file (3C UBC-MAG format)
% [H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');




%% Load forward operator
% TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
fprintf('Loading Sensitivity...\n');

% Case sperical
m_uvw = @(a,t,p) [a.*cos(pi*t).*cos(pi*p);...
    a.*cos(pi*t).*sin(pi*p);...
    a.*sin(pi*t)];

sProj = @(a,theta,phi)[spdiags(cos(pi*theta).*cos(pi*phi),0,mactv,mactv) spdiags(-a.*sin(pi*theta).*cos(pi*phi)*pi,0,mactv,mactv) spdiags(-a.*cos(pi*theta).*sin(pi*phi)*pi,0,mactv,mactv);
    spdiags(cos(pi*theta).*sin(pi*phi),0,mactv,mactv) spdiags(-a.*sin(pi*theta).*sin(pi*phi)*pi,0,mactv,mactv) spdiags(a.*cos(pi*theta).*cos(pi*phi)*pi,0,mactv,mactv);
    spdiags(sin(pi*theta),0,mactv,mactv) spdiags(a.*cos(pi*theta)*pi,0,mactv,mactv) sparse(mactv,mactv)];

%% Create gradient matrices and corresponding volume vectors
% [~, Gx, Gy, Gz, ~, ~, ~, ~] = get_GRAD_op3D_SQUARE(dx,dy,dz,nullcell,x);
% [ ~, Gx , Gy, Gz, ~, ~, ~, ~ ]= get_GRAD_op3D_SQUARE_Kron(dx,dy,dz,nullcell,x);
[A, GRAD, V] = get_GRAD_op3D_TENSIL_Kron(dx,dy,dz,nullcell,'FWR');
% [Ws, V ] = getWs3D(dx,dy,dz,X);

Ws =  spdiags(x * ( w(1:mcell)  )  ,0,mactv,mactv);
Wx =  spdiags(x * ( w(1+mcell:2*mcell) )  ,0,mactv,mactv);
Wy =  spdiags(x * ( w(1+2*mcell:3*mcell) )  ,0,mactv,mactv);
Wz =  spdiags(x * ( w(1+3*mcell:4*mcell) )  ,0,mactv,mactv);

%% Rotate gradient operators

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
% Scale the rows


% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(1,:)'),'descend') ;

Gx = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while round(denom) < (45 )
    
    Gx = Gx + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gx = Gx * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gx),2)) ~= 2;
% Gx(indx,:) = 0;

% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(2,:)'),'descend') ;

Gy = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while round(denom) < (45 )
    
    Gy = Gy + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gy = Gy * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gy),2)) ~= 2;
% Gy(indx,:) = 0;


% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(3,:)'),'descend') ;

Gz = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while round(denom) < (45)
    
    Gz = Gz + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gz = Gz * spdiags(ones(mactv,1)/denom,0,mactv,mactv);



%% START INVERSION
mactv = sum(nullcell);
target = chi_target * ndata;     % Target misifit

% Compute total objective function
objfunc = @(m,phim,b) sum( ((Gvec(G,Wd,m) - d) ).^2 ) + (m)' * b * phim * (m);

invmod      =  X*mstart ;
mref = X*mref;

aa = invmod(1:mactv);
tt = invmod(1+mactv:2*mactv);
pp = invmod(1+2*mactv:3*mactv);



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
alphas = kron(alphas,ones(3,1));
alphas(2:3,1) = 0;
%             alphas(2:3,2:4) = pi;
% dkdt = @(p,ep) (ep).^(1/(2*(2-p)));
% traffic_s = 1;    % Measures the number of points/iteration passing the lp corner
% traffic_xyz = 1;    % Measures the number of points/iteration passing the lp corner
% group_s = 1;  % Measures the number of lower than the lp corner
dphim = 100;
scl_t = 1;
scl_p = 1;
while switcher ~= 3 && count ~= max_iter

    fprintf(['MVI: Tile' num2str(idx) '\n'])
    count=count+1;

    if switcher == 0 %count == 1   

        delta_p(1:3) = 1e-1;%delta_tresh(1)*3;
        delta_q(1:3) = 1e-1;%delta_tresh(2)*3;


        aa = invmod(1:mactv);
        tt = invmod(1+mactv:2*mactv);
        pp = invmod(1+2*mactv:3*mactv);

%         J   =  ;
        S = sProj(aa,tt,pp);

        [aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_noGamma(invmod,mref,1,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,kron([1 1 1],[2 2 2 2 1]),FLAG1,FLAG2,switcher,delta_p,delta_q);

        wr = zeros(1,3*mactv);
        for gg = 1 : ndata
            wr = wr +  ( (G{1}(gg,:)* S)).^2;
        end   

%         % Re-weighted pseudo sensitivity weighting for each component
        wr = abs(sqrt(wr+1e-8))';
        
        wr = wr / max(wr);

        Wr = spdiags(wr,0,3*mactv,3*mactv);

%         mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;
        
%         if switcher == 0 && solves > 1
            
%             if count == 1
%                 
%                 scl_t = 1;
%                 scl_p = 1;
%                 
%             else
%                 phi_a = invmod(1:mactv)'*(mof((1:mactv),(1:mactv))'*mof((1:mactv),(1:mactv))*invmod(1:mactv));
%                 phi_t = invmod(1+mactv:2*mactv)'*(mof((1+mactv:2*mactv),(1+mactv:2*mactv))'*mof((1+mactv:2*mactv),(1+mactv:2*mactv))*invmod(1+mactv:2*mactv));
%                 phi_p = invmod(1+2*mactv:3*mactv)'*(mof((1+2*mactv:3*mactv),(1+2*mactv:3*mactv))'*mof((1+2*mactv:3*mactv),(1+2*mactv:3*mactv))*invmod(1+2*mactv:3*mactv));
                phi_a = (invmod(1:mactv))'*(Ws'*spdiags(wr(1:mactv),0,mactv,mactv)*Ws)*(invmod(1:mactv)) + invmod(1:mactv)'*(Gx'*Wx'*spdiags(wr(1:mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1:mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1:mactv),0,mactv,mactv)*Wz*Gz)*invmod(1:mactv);
                phi_t = (invmod(1+mactv:2*mactv))'*(Ws'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Ws)*(invmod(1+mactv:2*mactv)) + invmod(1+mactv:2*mactv)'*(Gx'*Wx'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+mactv:2*mactv);
                phi_p = (invmod(1+2*mactv:3*mactv))'*(Ws'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Ws)*(invmod(1+2*mactv:3*mactv)) + invmod(1+2*mactv:3*mactv)'*(Gx'*Wx'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+2*mactv:3*mactv);

                scl_t = (phi_a/phi_t)*pi*2;
                scl_p = (phi_a/phi_p)*pi;
                
                scl = max([scl_t,scl_p]);
%             end
        
            fprintf('Scale t: %f\t',scl_t)
            fprintf('Scale p: %f\n',scl_p)
%         else
%             phi_a = 1;
%             phi_t = 1;
%             phi_p = 1;
%         end
        wr(1:mactv)= wr(1:mactv);
        wr(1+mactv:2*mactv)= wr(1+mactv:2*mactv)*scl;%/(max(wr(1:mcell))/max(wr(1+mcell:2*mcell))) ;
        wr(1+2*mactv:3*mactv)= wr(1+2*mactv:3*mactv)*scl/2;%/(max(wr(1:mcell))/max(wr(1+2*mcell:3*mcell)));

        Wr = spdiags(wr,0,3*mactv,3*mactv);

        mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

        if isempty(beta)==1

            temp = randn(3*mactv,1);
            temp = temp/max(temp)*.5;
            beta = sum(((Gvec(G,Wd,temp))).^2) / (temp'*mof*temp) * 1e+3 ;
%             beta = 4e+4;
        end


        phi_MOF = (invmod-mref)'* ( aVRWs'*Wr*aVRWs * (invmod-mref) ) +...
              (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
              (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
              (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;

        phi(count) = norm((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d)).^2 + beta(count) + phi_MOF;
%         tresh = dkdt(2,delta(count));


    else

        if switcher ~= 0
            lp_count = lp_count+1;

            if lp_count == 1 
%                 aa = invmod(1:mactv);
%                 tt = invmod(1+mactv:2*mactv);
%                 pp = invmod(1+2*mactv:3*mactv);
%                 figure;
%                 subplot(2,1,1);hist(aa);
%                 subplot(2,1,2);hist(abs(Gx*tt),100);
%                 xlim([0 0.5]);
%                 title('L2l2 \theta');
%                 for pst = 1 : 3

                if eps_FLAG==0 

                    Z = zeros(1,3);
                    Z(1) = 1;

                    % Create sub-space matrix and grab the right model
                    Z = kron(Z,speye(mactv));
                    m = Z*(invmod-mref);

                    [pp,qq] = get_eps(m,10,Gx,Gy,Gz);

                    delta_p(1) = pp;
                    delta_q(1) = qq;

                else

                    delta_p(1) = eps_p(pst);
                    delta_q(1) = eps_q(pst);

                end
                
                delta_p(2:3) = delta_p(1); % Doesn't matter since not active
                delta_q(2:3) = 1.5e-1; % Fix value since always angles
 
                aa = invmod(1:mactv);
                tt = invmod(1+mactv:2*mactv);
                pp = invmod(1+2*mactv:3*mactv);
                
                model_out = X'*(m_uvw(aa,tt,pp));
    
                Mxyz = model_out;
                Mx = Mxyz(1 : mcell);
                My = Mxyz(((1+mcell) : 2*mcell));
                Mz = Mxyz(((1+2*mcell) : 3*mcell));
                M = [Mx My Mz];
                Mamp = sqrt( Mx.^2 + My.^2 + Mz.^2 );

                pred_TMI = Gvec(G,speye(ndata),m_uvw(aa,tt,pp));
                write_MAG3D_TMI([work_dir dsep 'Tile' num2str(idx) '_l2l2_MVI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
                save([work_dir dsep 'Tile' num2str(idx) '_l2l2_MVI.fld'],'-ascii','M')
                save([work_dir dsep 'Tile' num2str(idx) '_l2l2_MVI.amp'],'-ascii','Mamp')

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


        end

        if dphim(end)  < 1 && lp_count > 1 && phi_d(count-1) < target * (1+0.25) && phi_d(count-1) > target * (1-0.25)%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1

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
        [aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_noGamma(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p,delta_q);



%                     MOF = aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWy'*aVRWy + aVRWz'*aVRWz;


    end

%     dmdx = sqrt( (kron(speye(3),Wx) * invmod).^2 + (kron(speye(3),Wy) * invmod).^2 + (kron(speye(3),Wz) * invmod).^2 );    
%     group_xyz(count) = sum(abs(dmdx) <= tresh_xyz);


    %% Gauss-Newton steps

    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
    str_var = ['A','T','P'];
    for pst = 1 : 3
        fprintf([str_var(pst) 'eps_p: \t %8.5e \t eps_p*: \t %8.5e\n'],delta_p(pst),eps_p(pst))
        fprintf([str_var(pst) 'eps_q: \t %8.5e \t eps_q*: \t %8.5e\n'],delta_q(pst),eps_q(pst))
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

    phi_d(count) = sum(( Gvec(G,Wd,m_uvw(aa,tt,pp)) - d ).^2);

    if switcher ~=0
        temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(pi*2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(pi*temp(1+2*mactv:3*mactv)); 
        
        gamma = phi_m(end) /(...
                    (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                    (( (aVRWx * invmod) )' * Wr * ( aVRWx * invmod ))+...
                    (( (aVRWy * invmod) )' * Wr *( aVRWy * invmod ))+...
                    (( (aVRWz * invmod) )' * Wr *( aVRWz * invmod )));
                % gamma=1;
            aVRWs = sqrt(gamma) * aVRWs; 
            aVRWx = sqrt(gamma) * aVRWx;
            aVRWy = sqrt(gamma) * aVRWy;
            aVRWz = sqrt(gamma) * aVRWz;
    end
    
    temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(pi*2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(pi*temp(1+2*mactv:3*mactv));
                
    phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                      (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
                      (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
                      (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;
                  
    phi_in = phi_d(count) + beta(count) * phi_MOF;

                        % Try scaling the transformation matrix
    solves_max= 5;

    while solves <= solves_max %&& ggdm > 1e-2

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
        aa = invmod(1:mactv);
        tt = invmod(1+mactv:2*mactv);
        pp = invmod(1+2*mactv:3*mactv);

%         J   =  ;
       S = sProj(aa,tt,pp);
        wr = zeros(1,3*mactv);
        for gg = 1 : ndata
            wr = wr +  ( (G{1}(gg,:)* S)).^2;
        end
        
%         if switcher~=0
%             aa_scl = 1;%prctile(unique(aa),99);
%         else
%         aa_scl = prctile(unique(aa),99);
%         end
%         fprintf('Scale at %f\n',aa_scl);
        wr = abs(sqrt(wr+1e-8))';

        
       wr = wr / max(wr);

        Wr = spdiags(wr,0,3*mactv,3*mactv);

%         mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;
        
        if switcher == 0
            
%             if count == 1
%                 
%                 scl_t = 1;
%                 scl_p = 1;
%                 
%             else
%                 phi_a = invmod(1:mactv)'*(mof((1:mactv),(1:mactv))'*mof((1:mactv),(1:mactv))*invmod(1:mactv));
%                 phi_t = invmod(1+mactv:2*mactv)'*(mof((1+mactv:2*mactv),(1+mactv:2*mactv))'*mof((1+mactv:2*mactv),(1+mactv:2*mactv))*invmod(1+mactv:2*mactv));
%                 phi_p = invmod(1+2*mactv:3*mactv)'*(mof((1+2*mactv:3*mactv),(1+2*mactv:3*mactv))'*mof((1+2*mactv:3*mactv),(1+2*mactv:3*mactv))*invmod(1+2*mactv:3*mactv));
                phi_a = (invmod(1:mactv)-mref(1:mactv))'*(Ws'*spdiags(wr(1:mactv),0,mactv,mactv)*Ws)*(invmod(1:mactv)-mref(1:mactv)) + invmod(1:mactv)'*(Gx'*Wx'*spdiags(wr(1:mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1:mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1:mactv),0,mactv,mactv)*Wz*Gz)*invmod(1:mactv);
                phi_tp = invmod(1+mactv:2*mactv)'*(Gx'*Wx'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+mactv:2*mactv)+...
                 invmod(1+2*mactv:3*mactv)'*(Gx'*Wx'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+2*mactv:3*mactv);

                scl = (phi_a/phi_tp)/2;
%                 scl_p = (phi_a/phi_p);
                
%                 scl = max([scl_t,scl_p]);
%             end
        
            fprintf('Scale t: %f\t',scl)
%             fprintf('Scale p: %f\n',scl_p)
            
%         else
%             
% %             scl_t = 1;
%             scl = 1;
            
        end
%         else
%             phi_a = 1;
%             phi_t = 1;
%             phi_p = 1;
%         end
        wr(1:mactv)= wr(1:mactv);
        wr(1+mactv:2*mactv)= wr(1+mactv:2*mactv)*scl;%/(max(wr(1:mcell))/max(wr(1+mcell:2*mcell))) ;
        wr(1+2*mactv:3*mactv)= wr(1+2*mactv:3*mactv)*scl;%/(max(wr(1:mcell))/max(wr(1+2*mcell:3*mcell)));

        Wr = spdiags(wr,0,3*mactv,3*mactv);

        mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

        diagJ = zeros(1,3*mactv);
        for gg = 1 : ndata
            diagJ = diagJ +  ( G{1}(gg,:)* S).^2;
        end

        diagA   = diagJ + beta(count)*spdiags( mof ,0)';
        PreC    = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);

        switch FLAG1

            case 'SMOOTH_MOD'
                
                temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(pi*2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(pi*temp(1+2*mactv:3*mactv));
                
                g_MOF = ( aVRWs'*Wr*aVRWs * temp ) +...
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
        dm =zeros(3*mactv,1);
        [dm,~,ncg] = CG_Lin( dm, G, Wd, beta(count) * mof , g, S, PreC, Pac );
%         fprintf('CG iter %i \n',ncg);

        if max(abs(dm(1+mactv:2*mactv))) > max(abs(bounds(2,:)))
            dm(1+mactv:2*mactv) = dm(1+mactv:2*mactv)/10;
%             dm(1+2*mactv:3*mactv) = dm(1+2*mactv:3*mactv)/10;
        end
        
        if max(abs(dm(1+2*mactv:3*mactv))) > max(abs(bounds(3,:)))
%             dm(1+mactv:2*mactv) = dm(1+mactv:2*mactv)/10;
            dm(1+2*mactv:3*mactv) = dm(1+2*mactv:3*mactv)/10;
        end

        %% Step length, line search                
        tncg = tncg+ncg; % Record the number of CG iterations

        temp = spdiags(Pac);
        phi_out = phi_in;
        phi_temp = phi_in;
        m_temp = invmod;
        % Combine active and inactive cells step if active bounds
        if sum(temp)~=3*mactv

            for dd = 1 : 3
%                 
                k = zeros(3,3);
                k(dd,dd) = 1;
%                 
%     % Create sub-space matrix and grab the right model
            K = kron(k,speye(mactv));
            
            rhs_a = K*( speye(3*mactv) - Pac ) * (g);



%             lowb = K*m_temp <= K*lowBvec & rhs_a < 0;
%             uppb = K*m_temp >= K*uppBvec & rhs_a > 0;
% 
%             % Remove gradients going the wrong direction
%                         rhs_a( lowb) = 0;
%                         rhs_a( uppb) = 0;
                        
            dm_i = max( abs( K*dm ) );
            dm_a = max( abs(rhs_a) ); 
                        if dm_i < dm_a
                            dm = dm + rhs_a * dm_i / dm_a * 1e-8 ;
                        else
                            dm = dm + rhs_a;
                        end
            end
            
        end
        gamma = 2;

        % Reduce step length in order to reduce phid
        count_LS = 0;
        while (phi_out >= phi_in) && count_LS < 5
            count_LS = count_LS+1;
%             phi_temp(2) = phi_temp(1);
            if phi_temp < phi_out
                gamma = 0.5 * gamma;
            else
                gamma = 0.5 * gamma;
            end
            

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

            temp = (m_temp-mref);
%             temp(1+mactv:2*mactv) = sin(pi*2*temp(1+mactv:2*mactv));
%             temp(1+2*mactv:3*mactv) = sin(pi*temp(1+2*mactv:3*mactv));
            
            phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                      (m_temp)'*( aVRWx'*Wr*aVRWx * (m_temp) ) +...
                      (m_temp)'*( aVRWy'*Wr*aVRWy * (m_temp) ) +...
                      (m_temp)'*( aVRWz'*Wr*aVRWz * (m_temp) ) ;

            phi_temp = phi_out;
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
%         phi_m(count) = phi_MOF;
%         fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,ggdm);
        fprintf('GN iter %i: , Step Length: %f, phi_m :\t\t %8.5e\n',solves,ddm(2),phi_MOF);
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

    aa = invmod(1:mactv);
    tt = invmod(1+mactv:2*mactv);
    pp = invmod(1+2*mactv:3*mactv);

    phi_d(count) = sum(( Gvec(G,Wd,m_uvw(aa,tt,pp)) - d ).^2);

    temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(pi*2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(pi*temp(1+2*mactv:3*mactv));
                
    phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                      (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
                      (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
                      (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;
    phi_m(count) = phi_MOF;

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

%%
    aa = invmod(1:mactv);
    tt = invmod(1+mactv:2*mactv);
    pp = invmod(1+2*mactv:3*mactv);
        
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

% save([out_dir dsep 'K_' num2str(kk) ' O_' num2str(oo) ' L_' num2str(ll) '.fld'],'-ascii','M')
%   
%         end
%     end
% end

% figure;
% subplot(2,1,1);hist(aa);
% subplot(2,1,2);hist(abs(Gx*tt),100);
% xlim([0 0.5]);
% title('Final \theta');
figure;
subplot(1,3,1)
temp = invmod(1:mactv) - mref(1:mactv);
[n, xout] =hist((temp),100); hold off
bar(xout,n);  
set(gca,'yscale','log')
title('Final amp values')

subplot(1,3,2)
[n, xout] =hist((Gx * invmod(1+mactv:2*mactv)),100); hold off
bar(xout,n);  
set(gca,'yscale','log')
title('Final \Delta \theta values')

subplot(1,3,3)
[n, xout] =hist((Gx * invmod(1+2*mactv:3*mactv)),100); hold off
bar(xout,n);  
set(gca,'yscale','log')
title('Final \Delta \phi values')
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
% MOF_out = MOF;