function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D(m,mref,phim,Ws,Wx,Wy,Wz,Gx,Gy,Gz,s,t,alpha,LP,FLAG1,FLAG2,FLAG3,delta_s,delta_xyz)
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
mcell = length(m);

% Wr = spdiags(wr,0,mcell,mcell);

%% MODIFICATION FOR TMVI
% If LP is size(3*4), then the process is repeated for the three model
% spaces


    
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

% % Generate lp vectors
% pvec    = zeros(mcell,1);
% qxvec   = zeros(mcell,1);
% qyvec   = zeros(mcell,1);
% qzvec   = zeros(mcell,1);
% rvec    = zeros(mcell,1);
% 
% for jj = 1 : size(LP,1)
%     
%     pvec    =  pvec + t(:,jj)*LP(jj,1);
%     qxvec   =  qxvec + t(:,jj)*LP(jj,2);
%     qyvec   =  qyvec + t(:,jj)*LP(jj,3);
%     qzvec   =  qzvec + t(:,jj)*LP(jj,4);
%     rvec    =  rvec + t(:,jj)*LP(jj,5);
%     
% end


tRs = spalloc(mcell,mcell,mcell);
tRx = spalloc(mcell,mcell,mcell);
tRy = spalloc(mcell,mcell,mcell);
tRz = spalloc(mcell,mcell,mcell);




if FLAG3~=0
    
    rvec = zeros(mcell,1);

    for ii = 1 : size(s,2)
        
        tt =t(:,ii);        
        
        tt(tt>0) = sqrt(tt(tt>0));
        
        T = spdiags((tt),0,mcell,mcell);
%         S = spdiags(s(:,ii),0,mcell,mcell);
        
        rs = 1./ ( abs((m - mref)) .^2 + delta_s.^2 ).^(1 - LP(ii,1)/2);
        rx = 1./ ( abs(GRADx) .^2 + delta_xyz.^2 ).^(1 - LP(ii,2)/2);
        ry = 1./ ( abs(GRADy) .^2 + delta_xyz.^2 ).^(1 - LP(ii,3)/2);
        rz = 1./ ( abs(GRADz) .^2 + delta_xyz.^2 ).^(1 - LP(ii,4)/2);

        rvec = rvec + (T * ones(mcell,1)) * LP(ii,5);
        
        Rs = spdiags( (rs) .^ 0.5 ,0,length(rs),length(rs));

        Rx = spdiags( (rx) .^ 0.5,0,size(Wx,1),size(Wx,1));

        Ry = spdiags( (ry) .^ 0.5,0,size(Wy,1),size(Wy,1)); 
        Rz = spdiags( (rz) .^ 0.5,0,size(Wz,1),size(Wz,1));                  

        avrws =  T * Rs ;
        avrwx =  T * Rx * Gx  ;
        avrwy =  T * Ry * Gy  ;
        avrwz =  T * Rz * Gz  ;
       
        
        grad_s = 1 / max(abs((((avrws))'* ((avrws) * ((m - mref)) ) ) ) );
        
         switch FLAG1
        case 'SMOOTH_MOD'
            
            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs( ( avrwx'* ( (avrwx *  m) ) ) ) );
            grad_y = 1 / max(abs( ( avrwy'* ( (avrwy *  m) ) ) ) );
            grad_z = 1 / max(abs( ( avrwz'* ( (avrwz *  m) ) ) ) );

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs(  (avrwx'* ( (avrwx * ( m - mref ) ) ) ) ));
            grad_y = 1 / max(abs(  (avrwy'* ( (avrwy * ( m - mref ) ) ) ) ));
            grad_z = 1 / max(abs(  (avrwz'* ( (avrwz * ( m - mref ) ) ) ) ));

         end
        
        tRs = tRs + sqrt( grad_s ) * T * Rs ;
        tRx = tRx + sqrt( grad_x ) * T * Rx ;
        tRy = tRy + sqrt( grad_y ) * T * Ry ;
        tRz = tRz + sqrt( grad_z ) * T * Rz ;
        
        
             
    end
    
    aVRWs =  spdiags(sqrt( rvec * alpha(1) ),0,mcell,mcell) * Ws * tRs ;
    aVRWx =  spdiags(sqrt( (2 - rvec) * alpha(2) ),0,mcell,mcell) * Wx * tRx *  Gx  ;
    aVRWy =  spdiags(sqrt( (2 - rvec) * alpha(3) ),0,mcell,mcell) * Wy * tRy *  Gy  ;
    aVRWz =  spdiags(sqrt( (2 - rvec) * alpha(4) ),0,mcell,mcell) * Wz * tRz *  Gz  ;                    
                        
    gamma = phim /(...
            ((aVRWs * (m-mref) )' * (aVRWs * (m-mref) ))+...
            (( (aVRWx) * m )' * ( (aVRWx) * ( m ) ))+...
            (( (aVRWy) * m )' * ( (aVRWy) * ( m ) ))+...
            (( (aVRWz) * m )' * ( (aVRWz) * ( m ) )));
                 
%     switch FLAG1
%         case 'SMOOTH_MOD'
%             
%         sc_phigx = gamma;
%         sc_phigy = gamma;
%         sc_phigz = gamma;
%         
%         case 'SMOOTH_MOD_DIF'
% 
%         sc_phigx = gamma;
%         sc_phigy = gamma;
%         sc_phigz = gamma;
% 
%     end

    aVRWs = sqrt(gamma) * aVRWs; 
    aVRWx = sqrt(gamma) * aVRWx;
    aVRWy = sqrt(gamma) * aVRWy;
    aVRWz = sqrt(gamma) * aVRWz;
    
else
       
    % Form the matrices by including the weights and volume
    aVRWs =  sqrt( alpha(1) ) * Ws ;
    aVRWx =  sqrt( alpha(2) ) * Wx * Gx;
    aVRWy =  sqrt( alpha(3) ) * Wy * Gy;
    aVRWz =  sqrt( alpha(4) ) * Wz * Gz;
    
end

% If more than one zone, average the scaling between lp,lq on transition
% if size(LP,1) ~= 1;
        
%     aVRWs = spdiags(sqrt( gamma * rvec ),0,mcell,mcell) * aVRWs;
%     aVRWx = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx;
%     aVRWy = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy;
%     aVRWz = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz;
    
% Otherwise just a scaler
% else
%     
%     aVRWs = sqrt( gamma * LP(1,5) ) * aVRWs;
%     aVRWx = sqrt( sc_phigx * ( 2.0 - LP(1,5) ) ) * aVRWx;
%     aVRWy = sqrt( sc_phigy * ( 2.0 - LP(1,5) ) ) * aVRWy;
%     aVRWz = sqrt( sc_phigz * ( 2.0 - LP(1,5) ) ) * aVRWz;
%     
% end





        MOF = aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWy'*aVRWy + aVRWz'*aVRWz;



end

