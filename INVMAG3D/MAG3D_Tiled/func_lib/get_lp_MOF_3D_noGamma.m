function [aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_noGamma(model,model_ref,phim,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alpha,LP,FLAG1,FLAG2,FLAG3,delta_s,delta_xyz)
% Function get_lp_WRW(Ws,Wx,Wy,Wz,p,q,r)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). The lp-norm is approximated through a lag diffusivity approach
% where:
% 
% V3: ALL SCALED ON PHI_l2 (2014-08-02) 
%
% gradphi = Ws' * Rs * Ws +  Wx' * Rx * Wx + ...
% R = diag(1./ ( f(m).^(2-p) + epsilon )
% f(m) is abs(grad(m)) or abs(m)
%
% INPUT:
% m  :   Model at the kth iteration
% Ws :   Weighting matrix for the smallest component
% Wx :   Weighting matrix for the gradient x-component
% Wy :   Weighting matrix for the gradient y-component
% Wz :   Weighting matrix for the gradient z-component
% delta = 1e-5;
mcell = size(Ws,1);

% Wr = spdiags(wr,0,mcell,mcell);

%% MODIFICATION FOR TMVI
% If LP is size(3*5), then the process is repeated for the three model
% spaces

nspace = size(LP,2)/5;



for jj = 1 : nspace
%     
    % Temporary test - scale the phi and theta component for MVI
%     if jj ==2
%         
%         scs = 1e-1;
%         scx = 1e-1;
%         
%     elseif jj == 3
%         
%         scs = 1e-1;
%         scx = 1e-1;
%         
%     else
        
        scs = 1;
        scx = 1;
        
%     end
    
    d = zeros(1,nspace);
    d(jj) = 1;
    
    % Create sub-space matrix and grab the right model
    D = kron(d,speye(mcell));
    
    m = D * model;
    mref = D * model_ref;
    
%     LP = LP(:,(1:5)+(jj-1)*5);
    
    pvec    = zeros(mcell,1);
    qxvec   = zeros(mcell,1);
    qyvec   = zeros(mcell,1);
    qzvec   = zeros(mcell,1);
    rvec    = zeros(mcell,1);
    
    lp = LP(:,(1:5)+(5*(jj-1)));
    % Generate lp vectors
    for ii = 1 : size(LP,1)

        pvec    =  pvec + t(:,ii)*lp(ii,1);
        qxvec   =  qxvec + t(:,ii)*lp(ii,2);
        qyvec   =  qyvec + t(:,ii)*lp(ii,3);
        qzvec   =  qzvec + t(:,ii)*lp(ii,4);
        rvec    =  rvec + t(:,ii)*lp(ii,5);
        
    end
    
    switch FLAG1
        case 'SMOOTH_MOD'
            % Compute gradients then compute magnitude         
            dmx = Gx * m;
            dmy = Gy * m;
            dmz = Gz * m;

        case 'SMOOTH_MOD_DIF'


            dmx = Gx * (m - mref);
            dmy = Gy * (m - mref);
            dmz = Gz * (m - mref);

        otherwise
            fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
            fprintf('Please revise input\n')
            return
    end


    switch FLAG2
        case 'dmdx'
            % Compute gradients then compute magnitude 
            GRADx =  dmx;
            GRADy =  dmy;
            GRADz =  dmz;

        case 'GRADm'

            GRADx = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
            GRADy = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
            GRADz = sqrt(dmx.^2 + dmy.^2 + dmz.^2);

        otherwise
            fprintf('FLAG must be GRADm | dmdx\n')
            fprintf('Please revise input\n')
            return
    end
      

%         if jj == 2
%             
%             temp = sin(2*pi*(m - mref));
%             
%         elseif jj == 3
%             
%             temp = sin(pi*(m - mref));
%             
%         else
            
            temp = (m - mref);
            scl = max(m);
            
%         end
        rs = sqrt( 1./ ( abs(temp) .^2 + delta_s(jj).^2 ).^ (1 - pvec/2)  );
        rx = sqrt( 1./ ( abs(GRADx) .^2 + delta_xyz(jj).^2 ).^ (1 - qxvec/2)  );
        ry = sqrt( 1./ ( abs(GRADy) .^2 + delta_xyz(jj).^2 ).^ (1 - qyvec/2)  );
        rz = sqrt( 1./ ( abs(GRADz) .^2 + delta_xyz(jj).^2 ).^ (1 - qzvec/2)  );
        
%         Rs = spdiags( rs ,0,length(rs),length(rs));
% 
%         Rx = spdiags( rx ,0,size(Wx,1),size(Wx,1));
% 
%         Ry = spdiags( ry ,0,size(Wy,1),size(Wy,1)); 
%         Rz = spdiags( rz ,0,size(Wz,1),size(Wz,1));                  

        eta_s = delta_s(jj).^(1 - pvec/2);
        eta_x = delta_xyz(jj).^(1 - qxvec/2);
        eta_y = delta_xyz(jj).^(1 - qyvec/2);
        eta_z = delta_xyz(jj).^(1 - qzvec/2);
                 
        tRs = spdiags( sqrt(eta_s) .* rs,0,mcell,mcell);
        tRx = spdiags( sqrt(eta_x) .* rx,0,mcell,mcell);
        tRy = spdiags( sqrt(eta_y) .* ry,0,mcell,mcell);
        tRz = spdiags( sqrt(eta_z) .* rz,0,mcell,mcell);
    
        avrws =  spdiags(sqrt( rvec * alpha(jj,1)  ),0,mcell,mcell) * Ws * tRs;
        avrwx =  spdiags(sqrt( (2 - rvec)  * alpha(jj,2)  ),0,mcell,mcell) * Wx * tRx * Gx  ;
        avrwy =  spdiags(sqrt( (2 - rvec)  * alpha(jj,3)  ),0,mcell,mcell) * Wy * tRy * Gy  ;
        avrwz =  spdiags(sqrt( (2 - rvec)  * alpha(jj,4)  ),0,mcell,mcell) * Wz * tRz * Gz  ;                    
                       

%     Form the final model objective function
    if jj == 1
        
        aVRWs = kron(d,avrws);
        aVRWx = kron(d,avrwx);
        aVRWy = kron(d,avrwy);
        aVRWz = kron(d,avrwz);
        
    else
        
        aVRWs = [aVRWs;kron(d,avrws)];
        aVRWx = [aVRWx;kron(d,avrwx)];
        aVRWy = [aVRWy;kron(d,avrwy)];
        aVRWz = [aVRWz;kron(d,avrwz)];
        
    end

end



