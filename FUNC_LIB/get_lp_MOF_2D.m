function [MOF,aVRWs,aVRWx,aVRWz,gamma] = get_lp_MOF_2D(m,mref,grad_MOF,s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,LP,FLAG1,FLAG2,FLAG3,delta_p,delta_q)
% Function get_lp_MOF_2D(m,mref,phim,s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,LP,FLAG1,FLAG2,FLAG3,delta)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). 
% 
% V3: Pass on cell-based p,q and r 
%
%
% INPUT:
% m  :   Model at the kth iteration
% Ws :   Weighting matrix for the smallest component
% Wx :   Weighting matrix for the gradient x-component
% Wy :   Weighting matrix for the gradient y-component
% Wz :   Weighting matrix for the gradient z-component


mcell = size(V,1);

Wr      = spdiags(wr,0,mcell,mcell);
pvec    = zeros(mcell,1);
qxvec   = zeros(mcell,1);
qzvec   = zeros(mcell,1);
rvec    = zeros(mcell,1);

% Generate lp vectors
for jj = 1 : size(LP,1)
    
    pvec    =  pvec + t(:,jj)*LP(jj,1);
    qxvec   =  qxvec + t(:,jj)*LP(jj,2);
    qzvec   =  qzvec + t(:,jj)*LP(jj,3);
    rvec    =  rvec + t(:,jj)*LP(jj,4);
end


switch FLAG2
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude 
       
        dmx = Wx * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        % Compute gradients then compute magnitude 
        
        dmx = Wx * (m - mref);
        dmz = Wz * (m - mref);
        
    otherwise
        fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
        fprintf('Please revise input\n')
        return
end

% Augmente system if dealing with 3-component kernel
switch FLAG1
    case 'GRADm'
        % Compute gradients then compute magnitude 
        GRADx = dmx;
        GRADz = dmz;
        
    case 'lGRADml'
        
        % Compute gradients then compute magnitude 
        GRADx = sqrt(dmx.^2 + dmz.^2);
        GRADz = sqrt(dmx .^2 + dmz.^2);
        
    otherwise
        fprintf('FLAG must be GRADm | lGRADml\n')
        fprintf('Please revise input\n')
        return
end

% Pre-allocate memory
tRs = zeros(mcell,1);
tRx = zeros(mcell,1);
tRz = zeros(mcell,1);

        
if FLAG3~=0
      

    
%     for ii = 1 : size(s,2)
        
%         tt =t(:,ii)./sum(t,2);
%         T = spdiags((tt),0,mcell,mcell);
%         tt(tt>0) = sqrt(tt(tt>0));

        rs = sqrt(1./ ( abs((m - mref)) .^2 + delta_p.^2 ).^(1 - pvec/2));
        rx = sqrt(1./ ( abs(GRADx) .^2 + (delta_q).^2 ).^(1 - qxvec/2));
        rz = sqrt(1./ ( abs(GRADz) .^2 + (delta_q).^2 ).^(1 - qzvec/2));

        Rs = spdiags( (rs)  ,0,length(rs),length(rs));

        Rx = spdiags( (rx) ,0,size(Wx,1),size(Wx,1));

        Rz = spdiags( (rz) ,0,size(Wz,1),size(Wz,1));    

        
%         S = spdiags(s(:,ii),0,mcell,mcell);

                        
        avrws =  Rs * Ws ;
        avrwx =  Rx * Wx  ;
        avrwz =  Rz * Wz  ;   

%         eta_s = 1 / max(abs( (avrws)'* (avrws) * (m - mref)) ) ;
%         eta_s = 1  / max(rs );
            eta_s = delta_p.^(1 - pvec/2);
%         eta_s = 1;%mcell / (rs'*rs );
        % Scale the gradient p-norm on the smallest q-norm
        switch FLAG2
        case 'SMOOTH_MOD'

            % Compute gradients then compute magnitude 
%             eta_x = 1 / max(abs( (avrwx)'* (avrwx) * (m)) ) ;
%             eta_z = 1 / max(abs( (avrwz)'* (avrwz) * (m)) ) ;
% 
%             eta_x = 1  / max(rx );
%             eta_z = 1  / max(rz );

            eta_x = delta_q.^(1 - qxvec/2);
            eta_z = delta_q.^(1 - qzvec/2);

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs( (avrwx)'* (avrwx) * (m - mref)) ) ;
            grad_z = 1 / max(abs( (avrwz)'* (avrwz) * (m - mref)) ) ;


        end

        tRs = tRs + sqrt(eta_s) .* rs;
        tRx = tRx + sqrt(eta_x) .* rx;
        tRz = tRz + sqrt(eta_z) .* rz;

%     end
    
%     eta_s = mcell / (tRs'*tRs );
%     eta_x = mcell / (tRx'*tRx );
%     eta_z = mcell / (tRz'*tRz );
%     
    tRs = spdiags( tRs,0,mcell,mcell);
    tRx = spdiags( tRx,0,mcell,mcell);
    tRz = spdiags( tRz,0,mcell,mcell);
    
    aVRWs = spdiags(sqrt( rvec * alpha(1) ),0,mcell,mcell) * Wr * tRs * V * Ws;
    aVRWx = spdiags(sqrt( (2 - rvec) * alpha(2) ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
    aVRWz = spdiags(sqrt( (2 - rvec) * alpha(3) ),0,mcell,mcell) * Wr * tRz * Vz * Wz;
    
%% ALTERNATIVE FORMULATION
%     if size(LP,1) ~= 1;

%         aVRWs = spdiags(sqrt( rvec ),0,mcell,mcell) * aVRWs;
%         aVRWx = spdiags(sqrt( ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx;
%         aVRWz = spdiags(sqrt( ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz;

%     % Otherwise just a scaler
%     else
% 
%         aVRWs = sqrt( LP(1,4) ) * aVRWs;
%         aVRWx = sqrt(( 2.0 - LP(1,4) ) ) * aVRWx;
%         aVRWz = sqrt(( 2.0 - LP(1,4) ) ) * aVRWz;
% 
%     end
%% Scale final phi_m     
    gamma = grad_MOF /...
            ( m'*( aVRWs' * aVRWs + aVRWx' * aVRWx + aVRWz'*aVRWz) * m ) ;
         
aVRWs = sqrt(gamma) * aVRWs; 
aVRWx = sqrt(gamma) * aVRWx;
aVRWz = sqrt(gamma) * aVRWz;

else

    % Form the matrices by including the weights and volume
    aVRWs =  sqrt( LP(1,4) * alpha(1) ) * Wr * V * Ws ;
    aVRWx =  sqrt( ( 2.0 - LP(1,4) ) * alpha(2) ) * Wr * Vx * Wx  ;
    aVRWz =  sqrt( ( 2.0 - LP(1,4) ) * alpha(3) ) * Wr * Vz * Wz  ;

end


% Form the final model objective function
MOF = (aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWz'*aVRWz);


