function [M] = MVI3D_ATP_Tiled(work_dir,out_dir,dsep,tile_id,xn,yn,zn,H, HI, HD, obsx, obsy, obsz, G, d, wd,mstart,mref,chi_target,alphas,beta,bounds,LP,t,eps_FLAG,eps_tresh,FLAG1,FLAG2,max_iter,ROT)
% Magnetic Vector Inversion in spherical coordinates 
% Written by: D Fournier 
% Last update: 2016/05/01


% Load mesh file and convert to vectors (UBC format)
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = nx*ny*nz;


ndata = length(d);
            
Wd = spdiags(1./wd,0,ndata,ndata);

d = Wd * d;

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

lowBvec = X * lowBvec;
uppBvec = X * uppBvec;

t = x*t;

%% Load forward operator
% TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
fprintf('Loading Sensitivity...\n');

% Case sperical
m_uvw = @(a,t,p) [a.*cos(t).*cos(p);...
    a.*cos(t).*sin(p);...
    a.*sin(t)];

sProj = @(a,theta,phi)[spdiags(cos(theta).*cos(phi),0,mactv,mactv) spdiags(-a.*sin(theta).*cos(phi),0,mactv,mactv) spdiags(-a.*cos(theta).*sin(phi),0,mactv,mactv);
    spdiags(cos(theta).*sin(phi),0,mactv,mactv) spdiags(-a.*sin(theta).*sin(phi),0,mactv,mactv) spdiags(a.*cos(theta).*cos(phi),0,mactv,mactv);
    spdiags(sin(theta),0,mactv,mactv) spdiags(a.*cos(theta),0,mactv,mactv) sparse(mactv,mactv)];

%% Create gradient matrices and corresponding volume vectors
[A, GRAD, ~] = get_GRAD_op3D_TENSIL_Kron(dx,dy,dz,nullcell,'FWR');
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


GGGx = blkdiag(Gx,Gx,Gx);
GGGy = blkdiag(Gy,Gy,Gy);
GGGz = blkdiag(Gz,Gz,Gz);

%% START INVERSION
mactv = sum(nullcell);
target = chi_target * ndata;     % Target misifit
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
% alphas(2:3,:) = 1/pi;
alphas(2:3,1) = 0;
dphim = 100;
    
while switcher ~= 3 && count ~= max_iter

    fprintf(['MVI: Tile' num2str(tile_id) '\n'])
    count=count+1;

    if switcher == 0 %count == 1   
        if eps_FLAG==0  && count ==1

            [pp,qq] = get_eps(aa,10,Gx,Gy,Gz);

            delta_p(1) = 1e-3;
            delta_q(1) = 1e-3;
            
            delta_p(2:3) = 1e-1;%delta_tresh(1)*3;
            delta_q(2) = 5e-2;%delta_tresh(2)*3;
            delta_q(3) = 5e-2;
        elseif count==1

            delta_p(1) = eps_p(1);
            delta_q(1) = eps_q(1);
            
            delta_p(2:3) = 5e-2;%delta_tresh(1)*3;
            delta_q(2:3) = 5e-2;%delta_tresh(2)*3;

        end
        
        


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
        wr = abs(sqrt(wr+1e-10))';
        
        wr = wr / max(wr);

        scl_t = (delta_p(1)/(pi));%(phi_a/(phi_t+phi_p));
        scl_p = (delta_p(2)/pi);%(phi_a/phi_p);

        wr(1:mactv)= wr(1:mactv);
        wr(1+mactv:2*mactv)= wr(1+mactv:2*mactv)*scl_t;%/(max(wr(1:mcell))/max(wr(1+mcell:2*mcell))) ;
        wr(1+2*mactv:3*mactv)= wr(1+2*mactv:3*mactv)*scl_p;%/(max(wr(1:mcell))/max(wr(1+2*mcell:3*mcell)));

        Wr = spdiags(wr,0,3*mactv,3*mactv);
        
        mof = aVRWs'*Wr*aVRWs + aVRWx'*Wr*aVRWx + aVRWy'*Wr*aVRWy + aVRWz'*Wr*aVRWz;

        % Estimate beta if not provided
%         if count==1
%             
%             [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
%             temp = (invmod-mref);
% 
%             g_MOF = ( aVRWs'*Wr*aVRWs * temp ) +...
%                     (aVRWx'*Wr*  gx ) +...
%                     (aVRWy'*Wr*  gy ) +...
%                     (aVRWz'*Wr*  gz );
% 
% %                 temp = g_MOF(1+2*mactv:end);
% %                 temp(abs(temp)>pi) = -sign(temp(abs(temp)>pi)).*(2*pi-abs(temp(abs(temp)>pi)));
% %                 g_MOF(1+2*mactv:end) = temp;
% 
%             g_d = S'*(Gtvec(G,Wd,(Gvec(G,Wd,S*invmod))));
% 
%             ratio = abs((invmod'*g_d)/(invmod'*g_MOF));
% 
%             beta = ratio *0.1 ;
% 
%         end

        [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
        temp = (invmod-mref);
        phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                        (( gx )' * Wr * ( gx ))+...
                        (( gy )' * Wr *( gy ))+...
                        (( gz )' * Wr *( gz ));

        if count == 1
            
            phi_m = phi_MOF;
            
        end
        
        phi(count) = norm((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d)).^2 + beta(count) + phi_MOF;
%         tresh = dkdt(2,delta(count));


    else

        lp_count = lp_count+1;

        if lp_count == 1 

%             if eps_FLAG==0 
% 
%                 Z = zeros(1,3);
%                 Z(1) = 1;
% 
%                 % Create sub-space matrix and grab the right model
%                 Z = kron(Z,speye(mactv));
%                 m = Z*(invmod-mref);
% 
%                 [pp,qq] = get_eps(m,10,Gx,Gy,Gz);
% 
%                 delta_p(1) = pp;
%                 delta_q(1) = qq;
% 
%             else
% 
%                 delta_p(1) = eps_p(1);
%                 delta_q(1) = eps_q(1);
% 
%             end

%             delta_p(2:3) = delta_p(1); % Doesn't matter since not active
%             delta_q(2) = 1e-2; % Fix value since always angles
%             delta_q(3) = 1e-2;
            
            aa = invmod(1:mactv);
            tt = invmod(1+mactv:2*mactv);
            pp = invmod(1+2*mactv:3*mactv);

            model_out = X'*(m_uvw(aa,tt,pp));

            M = reshape(model_out,mcell,3);
            Mamp = sum(M.^2,2).^0.5;
            Mamp(nullcell==0) = -100;

            pred_TMI = Gvec(G,speye(ndata),m_uvw(aa,tt,pp));
            write_MAG3D_TMI([work_dir dsep 'Tile' num2str(tile_id) '_l2l2_MVI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);
            save([work_dir dsep 'Tile' num2str(tile_id) '_l2l2_MVI.fld'],'-ascii','M')
            save([work_dir dsep 'Tile' num2str(tile_id) '_l2l2_MVI.amp'],'-ascii','Mamp')

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

        if dphim(end)  < 1 && lp_count > 1 && phi_d(count-1) < target * (1+0.25) && phi_d(count-1) > target * (1-0.25)%traffic_s(end)*100 <= 1  && traffic_xyz(end)*100 <= 1

            fprintf('\n# # ADJUST BETA # #\n');
            switcher = 2;

        else

            fprintf('\n# # LP-LQ ITER# #\n');

        end

        %% Update the regularization
        fprintf('# # LP-LQ ITER# #\n');
        [aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_noGamma(invmod,mref,phi_m(end),Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alphas,LP,FLAG1,FLAG2,switcher,delta_p,delta_q);

    end

    %% Gauss-Newton steps

    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
    str_var = ['A','T','P'];
%     for pst = 1 : 3
%         fprintf([str_var(pst) 'eps_p: \t %8.5e \t eps_p*: \t %8.5e\n'],delta_p(pst),eps_p(pst))
%         fprintf([str_var(pst) 'eps_q: \t %8.5e \t eps_q*: \t %8.5e\n'],delta_q(pst),eps_q(pst))
%     end

    % Save current model
    m_in = invmod;

    tncg = 0;
    ggdm = 1;       % Mesure the relative change in rel|dm|
    ddm = [1 1];    % Mesure of change in model update |dm|
    solves = 1;

    phi_d(count) = sum(( Gvec(G,Wd,m_uvw(aa,tt,pp)) - d ).^2);

    


                        % Try scaling the transformation matrix
    solves_max= 5;
    if switcher ~=0
        temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(temp(1+2*mactv:3*mactv)); 
        [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
        gamma = phi_m(end) /(...
                    (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                    (( gx )' * Wr * ( gx ))+...
                    (( gy )' * Wr *( gy ))+...
                    (( gz )' * Wr *( gz )));
                % gamma=1;
        aVRWs = sqrt(gamma) * aVRWs; 
        aVRWx = sqrt(gamma) * aVRWx;
        aVRWy = sqrt(gamma) * aVRWy;
        aVRWz = sqrt(gamma) * aVRWz;
    end
        
    while solves <= solves_max %&& ggdm > 1e-2
        
%         avrws = aVRWs;
%         avrwx = aVRWx;
%         avrwy = aVRWy;
%         avrwz = aVRWz;
    
        aa = invmod(1:mactv);
        tt = invmod(1+mactv:2*mactv);
        pp = invmod(1+2*mactv:3*mactv);

        S = sProj(aa,tt,pp);
        wr = zeros(1,3*mactv);
        for gg = 1 : ndata
            wr = wr +  ( (G{1}(gg,:)* S)).^2;
        end
        
        wr = abs(sqrt(wr+1e-10))';
        wr = wr/max(wr);
        
%         wr(1:mactv) = wr(1:mactv)/max(wr(1:mactv));
%         wr(1+mactv:2*mactv)= wr(1+mactv:2*mactv)/max(wr(1+mactv:2*mactv));
%         wr(1+2*mactv:3*mactv)= wr(1+2*mactv:3*mactv)/max(wr(1+2*mactv:3*mactv));

        
%         [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,speye(3*mactv),speye(3*mactv),speye(3*mactv));
        
%         phi_a = (invmod(1:mactv))'*(Ws'*spdiags(wr(1:mactv),0,mactv,mactv)*Ws)*(invmod(1:mactv)) + invmod(1:mactv)'*(Gx'*Wx'*spdiags(wr(1:mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1:mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1:mactv),0,mactv,mactv)*Wz*Gz)*invmod(1:mactv);
% %         phi_t = (invmod(1+mactv:2*mactv))'*(Ws'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Ws)*(invmod(1+mactv:2*mactv)) + invmod(1+mactv:2*mactv)'*(Gx'*Wx'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+mactv:2*mactv);
% %         phi_p = (invmod(1+2*mactv:3*mactv))'*(Ws'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Ws)*invmod(1+2*mactv:3*mactv) + invmod(1+2*mactv:3*mactv)'*(Gx'*Wx'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wx*Gx + Gy'*Wy'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wy*Gy + Gz'*Wz'*spdiags(wr(1+2*mactv:3*mactv),0,mactv,mactv)*Wz*Gz)*invmod(1+2*mactv:3*mactv);
%         phi_t = gx(1+mactv:2*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gx(1+mactv:2*mactv) + gy(1+mactv:2*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gy(1+mactv:2*mactv) + gz(1+mactv:2*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gz(1+mactv:2*mactv);
%         phi_p = gx(1+2*mactv:3*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gx(1+2*mactv:3*mactv) +...
%             gy(1+2*mactv:3*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gy(1+2*mactv:3*mactv) +...
%             gz(1+2*mactv:3*mactv)'*spdiags(wr(1+mactv:2*mactv),0,mactv,mactv)*gz(1+2*mactv:3*mactv);
        %         avrwx = scx*aVRWx;
%         avrwy = scy*aVRWy;
%         avrwz = scz*aVRWz;
        
%         Wr = spdiags(wr,0,3*mactv,3*mactv);
% 
%         g = (aVRWs'*Wr*aVRWs + scx*aVRWx'*Wr*aVRWx + scy*aVRWy'*Wr*aVRWy + scz*aVRWz'*Wr*aVRWz)*invmod;
%         
%         phi_a = max(abs(g(1:mactv)));
%         phi_t = max(abs(g(1+mactv:2*mactv)));
%         phi_p = max(abs(g(1+2*mactv:3*mactv)));

%         phi_a = (invmod(1:mactv))'*g(1:mactv);
%         phi_t = (invmod(1+mactv:3*mactv))'*g(1+mactv:3*mactv);
%         phi_p = (invmod(1+2*mactv:3*mactv))'*g(1+2*mactv:3*mactv);
        
%         scl_t = (phi_a/(phi_t));
%         scl_p = (phi_a/phi_p);
%         sprintf('%e, %e',scl_t,scl_p)
        scl_t = (delta_p(1)/(delta_q(2)));%(phi_a/(phi_t+phi_p));
        scl_p = (delta_p(1)/delta_q(3));%(phi_a/phi_p);

        wr(1:mactv)= wr(1:mactv);
        wr(1+mactv:2*mactv)= wr(1+mactv:2*mactv)*scl_t;%/(max(wr(1:mcell))/max(wr(1+mcell:2*mcell))) ;
        wr(1+2*mactv:3*mactv)= wr(1+2*mactv:3*mactv)*scl_p;%/(max(wr(1:mcell))/max(wr(1+2*mcell:3*mcell)));

        Wr = spdiags(wr,0,3*mactv,3*mactv);


        % Re-adjust the ratio between phid and phim for beta
        temp = (invmod-mref);
%                 temp(1+mactv:2*mactv) = sin(2*temp(1+mactv:2*mactv));
%                 temp(1+2*mactv:3*mactv) = sin(temp(1+2*mactv:3*mactv)); 

%         phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
%                   (invmod)'*( aVRWx'*Wr*aVRWx * (invmod) ) +...
%                   (invmod)'*( aVRWy'*Wr*aVRWy * (invmod) ) +...
%                   (invmod)'*( aVRWz'*Wr*aVRWz * (invmod) ) ;
%                       
%         gamma = phi_m(end) /phi_MOF;
%         
%                 % gamma=1;
%         aVRWs = sqrt(gamma) * avrws; 
%         aVRWx = sqrt(gamma) * avrwx;
%         aVRWy = sqrt(gamma) * avrwy;
%         aVRWz = sqrt(gamma) * avrwz;


        temp = (invmod-mref);
        [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
        phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                        (( gx )' * Wr * ( gx ))+...
                        (( gy )' * Wr *( gy ))+...
                        (( gz )' * Wr *( gz ));

    
        diagJ = zeros(1,3*mactv);
        for gg = 1 : ndata
            diagJ = diagJ +  ( G{1}(gg,:)* S).^2;
        end
        
        mof = aVRWs'*Wr*aVRWs +...
            (aVRWx * GGGx)'*Wr*(aVRWx * GGGx) +...
            (aVRWy * GGGy)'*Wr*(aVRWy * GGGy) +...
            (aVRWz * GGGz)'*Wr*(aVRWz * GGGz);
        
        


        switch FLAG1

            case 'SMOOTH_MOD'
                
                
%                 temp(1+mactv:end) = abs(invmod(1+mactv:end));
%                 temp(1+2*mactv:end) = abs(invmod(1+2*mactv:end));
% -sign(aa(abs(aa)>pi)).*(2*pi-abs(aa(abs(aa)>pi)));

                [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
                
                temp = (invmod-mref);
                g_MOF = ( aVRWs'*Wr*aVRWs * temp ) +...
                        ((aVRWx * GGGx)'*Wr*  gx ) +...
                        ((aVRWy * GGGy)'*Wr*  gy ) +...
                        ((aVRWz * GGGz)'*Wr*  gz );
                
%                 temp = g_MOF(1+2*mactv:end);
%                 temp(abs(temp)>pi) = -sign(temp(abs(temp)>pi)).*(2*pi-abs(temp(abs(temp)>pi)));
%                 g_MOF(1+2*mactv:end) = temp;

                g_d = S'*(Gtvec(G,Wd,(Gvec(G,Wd,m_uvw(aa,tt,pp)) - d)));
                
%                 if solves==1
%                 scale = ratio/abs((invmod'*g_d)/(invmod'*g_MOF));
%                 
%                 beta(count) = scale*beta(count);
%                 end
                g = -g_d  - beta(count)*g_MOF;

            case 'SMOOTH_MOD_DIF'
                
                fprintf('SMOOTH_MOD_DIF NOT IMPLEMENTED...SORRY!\n')
    %                 g = [- (magB - d) ; ...
    %             - sqrt( beta(count) ) * ( aVRWs * (invmod-mref) ) ;...
    %             - sqrt( beta(count) ) * ( aVRWx * (invmod-mref) ) ;...
    %             - sqrt( beta(count) ) * ( aVRWy * (invmod-mref) ) ;...
    %             - sqrt( beta(count) ) * ( aVRWz * (invmod-mref) ) ];
        end 


        phi_in = phi_d(count) + beta(count) * phi_MOF;
                
        diagA   = diagJ + beta(count)*spdiags( mof ,0)';
        PreC    = Pac * spdiags(1./diagA(:),0,3*mactv,3*mactv);
        %% Projected steepest descent
        dm =zeros(3*mactv,1);
        [dm,~,ncg] = CG_Lin( dm, G, Wd, beta(count) * mof , g, S, PreC, Pac );

        %% TESTING 
        if max(abs(dm(1+mactv:2*mactv))) > max(abs(bounds(2,:)))
            dm = dm/max(abs(dm(1+mactv:2*mactv)))*pi/2;
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
        while (phi_out >= phi_in) && count_LS < 10
            count_LS = count_LS+1;
%             phi_temp(2) = phi_temp(1);
            if phi_temp < phi_out
                gamma = 0.5 * gamma;
            else
                gamma = 0.5 * gamma;
            end
            

            gdm = gamma * dm;
            
%             gdm(1+mactv:end) = rem(gdm(1+mactv:end),pi);
%             
            ddm(2) = norm(gdm);
            
            m_temp = invmod + gdm;
%             
%             aa = m_temp(1:mactv);
%             tt = m_temp(1+mactv:2*mactv);
%             pp = m_temp(1+2*mactv:3*mactv);
%             m_xyz = reshape(m_uvw(aa,tt,pp),mcell,3);
%             
%             aa = sum(m_xyz.^2,2).^0.5;
%     
%             tt = zeros(mactv,1);
%             pp = zeros(mactv,1);
% 
%             tt(aa>0) = asin(m_xyz(aa>0,3)./(aa(aa>0)));
% 
% 
%             dydx = m_xyz(aa>0,2)./m_xyz(aa>0,1);
%             pp(aa>0) = atan(dydx);
%             pp(m_xyz(:,2)>0 & m_xyz(:,1)<0) = pp(m_xyz(:,2)>0 & m_xyz(:,1)<0) + pi;
%             pp(m_xyz(:,2)<0 & m_xyz(:,1)<0) = pp(m_xyz(:,2)<0 & m_xyz(:,1)<0) - pi;
%     
%             m_temp = [aa;tt;pp];
%             t_ind = kron([0;1;0],ones(mactv,1));
%             p_ind = kron([0;0;1],ones(mactv,1));
%             
            lowb = m_temp <= lowBvec;
            uppb = m_temp >= uppBvec;
% 
%             pp_low = m_temp(lowb==1 & p_ind);           
%             pp_low = mod(pp_low,-2*pi);
%             pp_low(pp_low<-pi) = 2*pi+pp_low(pp_low<-pi);
%             
% %             tt_low_ind = lowb==1 & t_ind==1;
% %             tt_low = m_temp(tt_low_ind);
% %             tt_low = -pi/2-rem(tt_low,pi/2);
% %             
% %             tt_low_ind = reshape(tt_low_ind,mactv,3);
% %             tt_low_ind = tt_low_ind(:,[1 3 2]);
% %             tt_low_ind = tt_low_ind(:);
%             
%             
%             pp_high = m_temp(uppb==1 & p_ind);
%             pp_high = mod(pp_high,2*pi);
%             pp_high(pp_high>pi) = -2*pi+pp_high(pp_high>pi);
%             
% %             tt_high_ind = uppb==1 & t_ind==1;
% %             tt_high = m_temp(tt_high_ind);
% %             tt_high = pi/2-rem(tt_high,pi/2);
% %             
% %             tt_high_ind = reshape(tt_high_ind,mactv,3);
% %             tt_high_ind = tt_high_ind(:,[1 3 2]);
% %             tt_high_ind = tt_high_ind(:);
% %             
% 
%             
%             % Apply bound on model
%             m_temp(lowb==1) = lowBvec(lowb==1);
%             m_temp(lowb==1 & p_ind) = pp_low;
% %             m_temp(lowb==1 & t_ind) = tt_low;
%             
%             m_temp(uppb==1) = uppBvec(uppb==1);
%             m_temp(uppb==1 & p_ind) = pp_high;
%             m_temp(uppb==1 & t_ind) = tt_high;
            
%             pp_flip = m_temp(tt_high_ind | tt_low_ind)+pi;
%             pp_flip(pp_flip>pi) = pp_flip(pp_flip>pi)-2*pi;
%             m_temp(tt_high_ind | tt_low_ind) = pp_flip;
            
            lowb(1+mactv:end) = 0;
            uppb(1+mactv:end) = 0;
            
            m_temp(lowb==1) = lowBvec(lowb==1);
            % Update projection matrix
            Pac = spdiags((lowb==0).*(uppb==0),0,3*mactv,3*mactv);

            aa = m_temp(1:mactv);
            tt = m_temp(1+mactv:2*mactv);
            pp = m_temp(1+2*mactv:3*mactv);

            temp = (m_temp-mref);
            
            [gx, gy, gz] = projectAngle(m_temp,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
                    
            phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                        (( gx )' * Wr * ( gx ))+...
                        (( gy )' * Wr *( gy ))+...
                        (( gz )' * Wr *( gz ));
            phi_temp = phi_out;
            phi_out = sum((Gvec(G,Wd,m_uvw(aa,tt,pp)) - d).^2) + beta(count) * phi_MOF;

        end

        phi_in = phi_out;

        if solves == 1

            ggdm = 1;
            ddm(1) = ddm(2);

        else

            ggdm = ddm(2)/ddm(1);

        end


        % Project model back to spherical
        aa = m_temp(1:mactv);
        tt = m_temp(1+mactv:2*mactv);
        pp = m_temp(1+2*mactv:3*mactv);
        m_xyz = reshape(m_uvw(aa,tt,pp),mcell,3);

        aa = sum(m_xyz.^2,2).^0.5;

        tt = zeros(mactv,1);
        pp = zeros(mactv,1);

        tt(aa>0) = asin(m_xyz(aa>0,3)./(aa(aa>0)));
        pp(aa>0) = atan2(m_xyz(aa>0,2),m_xyz(aa>0,1));
        
        invmod = [aa;tt;pp];
%         phi_m(count) = phi_MOF;
%         fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,ggdm);
        fprintf('GN iter %i: , Step Length: %f, phi_m :\t\t %8.5e\n',solves,ddm(2),phi_MOF);
        fprintf('Number of line search %i\n',count_LS);
        solves = solves + 1;


    end
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

    aa = invmod(1:mactv);
    tt = invmod(1+mactv:2*mactv);
    pp = invmod(1+2*mactv:3*mactv);

    phi_d(count) = sum(( Gvec(G,Wd,m_uvw(aa,tt,pp)) - d ).^2);

    temp = (invmod-mref);
                
    [gx, gy, gz] = projectAngle(invmod,GGGx,GGGy,GGGz,aVRWx,aVRWy,aVRWz);
                    
    phi_MOF = (temp)'* ( aVRWs'*Wr*aVRWs * (temp) ) +...
                (( gx )' * Wr * ( gx ))+...
                (( gy )' * Wr *( gy ))+...
                (( gz )' * Wr *( gz ));
              
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


    % Write log file
    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf(fid,' \t %8.5e ',sum( (aVRWs*invmod).^2 ) );
    fprintf(fid,' \t %8.5e ',sum( (aVRWx*GGGx*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWy*GGGy*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',sum( (aVRWz*GGGz*invmod).^2 ));
    fprintf(fid,' \t %8.5e ',phi_MOF);
    fprintf(fid,' \t %8.5e ',phi(count));
    fprintf(fid,' \t\t %i ',sum(tcells));
    fprintf(fid,' \t\t %i\n',ncg);

%%
    aa = invmod(1:mactv);
    tt = invmod(1+mactv:2*mactv);
    pp = invmod(1+2*mactv:3*mactv);
        
    model_out = X'*(m_uvw(aa,tt,pp));
    M = reshape(model_out,mcell,3);
    Mamp = sum(M.^2,2).^0.5;
    Mamp(nullcell==0) = -100;
    save([work_dir dsep 'Tile'  '_MVIatp.fld'],'-ascii','M')
    save([work_dir dsep 'Tile'  '_MVIphi.mod'],'-ascii','pp')
    save([work_dir dsep 'Tile'  '_MVItheta.mod'],'-ascii','tt')
    save([work_dir dsep 'Tile'  '_MVIatp.amp'],'-ascii','Mamp')

    pred_TMI = Gvec(G,speye(ndata),m_uvw(aa,tt,pp));
    write_MAG3D_TMI([work_dir dsep 'Tile' '_MVIatp.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);

%     write_MAG3D_TMI([work_dir dsep 'Tile' num2str(idx) '_iter_' num2str(count) '.pre'],H,I,Dazm,obsx,obsy,obsz,(G*invmod).*wd,wd);
end

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
write_MAG3D_TMI([out_dir dsep 'Tile' num2str(tile_id) '_MVI.pre'],H,HI,HD,HI,HD,obsx,obsy,obsz,pred_TMI,wd);

function [gx, gy, gz] = projectAngle(m,gradx,grady,gradz,wx,wy,wz)
    mactv = length(m)/3;
    gx = gradx *  m;
                
    temp = gx(1+2*mactv:end);
    temp(abs(temp)>=pi) = -sign(temp(abs(temp)>=pi)).*(2*pi-abs(temp(abs(temp)>=pi)));
    gx(1+2*mactv:end) = temp;
    gx = wx * gx;


    gy = grady *  m;
    temp = gy(1+2*mactv:end);
    temp(abs(temp)>=pi) = -sign(temp(abs(temp)>=pi)).*(2*pi-abs(temp(abs(temp)>=pi)));
    gy(1+2*mactv:end) = temp;
    gy = wy * gy;

    gz = gradz *  m;
    temp = gz(1+2*mactv:end);
    temp(abs(temp)>=pi) = -sign(temp(abs(temp)>=pi)).*(2*pi-abs(temp(abs(temp)>=pi)));
    gz(1+2*mactv:end) = temp;
    gz = wz * gz;