function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_TMVI(m_3C,mref_3C,phim,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,s,t,wr,alpha,LP,FLAG1,FLAG2,FLAG3,delta)
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

mcell = length(m_3C)/3;

%% P+ COMPONENT
m = m_3C(1:mcell);
mref = mref_3C(1:mcell);

Wr = spdiags(wr(1:mcell),0,mcell,mcell);
% IWR = spdiags(1./wr,0,mcell,mcell);

switch FLAG1
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude         
        dmx = Wx * m;
        dmy = Wy * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        
        dmx = Wx * (m - mref);
        dmy = Wy * (m - mref);
        dmz = Wz * (m - mref);
        
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

sc_phis = 1;
sc_phigx = 1;
sc_phigy = 1;
sc_phigz = 1;

aVRWs_p = spalloc(mcell,mcell,mcell);
aVRWx_p = spalloc(mcell,mcell,mcell);
aVRWy_p = spalloc(mcell,mcell,mcell);
aVRWz_p = spalloc(mcell,mcell,mcell);


if FLAG3~=0
    
    rvec = zeros(mcell,1);
    
    for ii = 1 : size(s{1},2)
        
        tt =t{1}(:,ii);
        
        rvec = rvec+t{1}(:,ii)*LP{1}(ii,5);
        
        tt(tt>0) = sqrt(tt(tt>0));
        
        T = spdiags((tt),0,mcell,mcell);
        S = spdiags(s{1}(:,ii),0,mcell,mcell);
        
        rs = 1 ./ ( abs(((m-mref))) .^( 2.0 - LP{1}(ii,1) ) + delta );
        rx = 1 ./ ( abs(GRADx) .^( 2.0 - LP{1}(ii,2) ) + delta ) ;
        ry = 1 ./ ( abs(GRADy) .^( 2.0 - LP{1}(ii,3) ) + delta ) ;
        rz = 1 ./ ( abs(GRADz) .^( 2.0 - LP{1}(ii,4) ) + delta ) ;
        
        
        Rs = spdiags( (rs) .^0.5,0,mcell,mcell);

        Rx = spdiags( (rx) .^0.5,0,mcell,mcell);

        Ry = spdiags( (ry) .^0.5,0,mcell,mcell);

        Rz = spdiags( (rz) .^0.5,0,mcell,mcell);                    

        avrws =  S * Rs * Ws ;
        avrwx =  S * Rx * Wx  ;
        avrwy =  S * Ry * Wy  ;
        avrwz =  S * Rz * Wz  ;
       
        
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
        
        aVRWs_p = aVRWs_p + sqrt( grad_s ) * T * Rs;
        aVRWx_p = aVRWx_p + sqrt( grad_x ) * T * Rx;
        aVRWy_p = aVRWy_p + sqrt( grad_y ) * T * Ry;
        aVRWz_p = aVRWz_p + sqrt( grad_z ) * T * Rz;
        
        
             
    end
    
    aVRWs_p =  sqrt( alpha(1) ) *aVRWs_p*Wr *   V  * Ws ;
    aVRWx_p =  sqrt( alpha(2) ) *aVRWx_p*Wr *   Vx  * Wx  ;
    aVRWy_p =  sqrt( alpha(3) ) *aVRWy_p*Wr *   Vy  * Wy  ;
    aVRWz_p =  sqrt( alpha(4) ) *aVRWz_p*Wr *   Vz  * Wz  ;                    
                        
    sc_phis = phim(1) /(...
            ((aVRWs_p * (m-mref) )' * (aVRWs_p * (m-mref) ))+...
            (( (aVRWx_p) * m )' * ( (aVRWx_p) * ( m ) ))+...
            (( (aVRWy_p) * m )' * ( (aVRWy_p) * ( m ) ))+...
            (( (aVRWz_p) * m )' * ( (aVRWz_p) * ( m ) )));
                 
    switch FLAG1
        case 'SMOOTH_MOD'
            
        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;
        
        case 'SMOOTH_MOD_DIF'

        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;

    end


else
       
    % Form the matrices by including the weights and volume
    aVRWs_p =  sqrt( alpha(1) ) *Wr *   V * Ws ;
    aVRWx_p =  sqrt( alpha(2) ) *Wr *   Vx * Wx  ;
    aVRWy_p =  sqrt( alpha(3) ) *Wr *   Vy * Wy  ;
    aVRWz_p =  sqrt( alpha(4) ) *Wr *   Vz * Wz  ;
    
end

% If more than one zone, average the scaling between lp,lq on transition
if size(LP{1},1) ~= 1;
        
    aVRWs_p = spdiags(sqrt( sc_phis * rvec ),0,mcell,mcell) * aVRWs_p;
    aVRWx_p = spdiags(sqrt( sc_phigx * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx_p;
    aVRWy_p = spdiags(sqrt( sc_phigy * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy_p;
    aVRWz_p = spdiags(sqrt( sc_phigz * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz_p;
    
% Otherwise just a scaler
else
    
    aVRWs_p = sqrt( sc_phis * LP{1}(1,5) ) * aVRWs_p;
    aVRWx_p = sqrt( sc_phigx * ( 2.0 - LP{1}(1,5) ) ) * aVRWx_p;
    aVRWy_p = sqrt( sc_phigy * ( 2.0 - LP{1}(1,5) ) ) * aVRWy_p;
    aVRWz_p = sqrt( sc_phigz * ( 2.0 - LP{1}(1,5) ) ) * aVRWz_p;
    
end



% Form the final model objective function
MOF_pp = aVRWs_p'*aVRWs_p + aVRWx_p'*aVRWx_p + aVRWy_p'*aVRWy_p + aVRWz_p'*aVRWz_p;


%% P- COMPONENT
% m_p = zeros(mcell,1);
% m_p(m<0) = m(m<0);

% m = sqrt( m_p.^2 +  m_3C((1+mcell):2*mcell).^2 +  m_3C((1+2*mcell):3*mcell).^2 );

% m = m_3C((1+mcell):2*mcell);
% mref = mref_3C((1+mcell):2*mcell);
% 
% % IWR = spdiags(1./wr,0,mcell,mcell);
% 
% switch FLAG1
%     case 'SMOOTH_MOD'
%         % Compute gradients then compute magnitude         
%         dmx = Wx * m_pst;
%         dmy = Wy * m_pst;
%         dmz = Wz * m_pst;
%         
%     case 'SMOOTH_MOD_DIF'
%         
%         
%         dmx = Wx * (m - mref);
%         dmy = Wy * (m - mref);
%         dmz = Wz * (m - mref);
%         
%     otherwise
%         fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
%         fprintf('Please revise input\n')
%         return
% end
% 
% 
% switch FLAG2
%     case 'dmdx'
%         % Compute gradients then compute magnitude 
%         GRADx =  dmx;
%         GRADy =  dmy;
%         GRADz =  dmz;
% 
%     case 'GRADm'
% 
%         GRADx = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
%         GRADy = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
%         GRADz = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
% 
%     otherwise
%         fprintf('FLAG must be GRADm | dmdx\n')
%         fprintf('Please revise input\n')
%         return
% end
% 
% sc_phis = 1;
% sc_phigx = 1;
% sc_phigy = 1;
% sc_phigz = 1;
% 
% aVRWs_pm = spalloc(mcell,mcell,mcell);
% aVRWx_pm = spalloc(mcell,mcell,mcell);
% aVRWy_pm = spalloc(mcell,mcell,mcell);
% aVRWz_pm = spalloc(mcell,mcell,mcell);
% 
% 
% if FLAG3~=0
%     
%     rvec = zeros(mcell,1);
%     
%     for ii = 1 : size(s{2},2)
%         
%         tt =t{2}(:,ii);
%         
%         rvec = rvec+t{2}(:,ii)*LP{2}(ii,5);
%         
%         tt(tt>0) = sqrt(tt(tt>0));
%         
%         T = spdiags((tt),0,mcell,mcell);
%         S = spdiags(s{2}(:,ii),0,mcell,mcell);
%         
%         rs = 1 ./ ( abs((m_pst)) .^( 2.0 - LP{2}(ii,1) ) + delta );
%         rx = 1 ./ ( abs(GRADx) .^( 2.0 - LP{2}(ii,2) ) + delta ) ;
%         ry = 1 ./ ( abs(GRADy) .^( 2.0 - LP{2}(ii,3) ) + delta ) ;
%         rz = 1 ./ ( abs(GRADz) .^( 2.0 - LP{2}(ii,4) ) + delta ) ;
%         
%         
%         Rs = spdiags( (rs) .^0.5,0,mcell,mcell);
% 
%         Rx = spdiags( (rx) .^0.5,0,mcell,mcell);
% 
%         Ry = spdiags( (ry) .^0.5,0,mcell,mcell);
% 
%         Rz = spdiags( (rz) .^0.5,0,mcell,mcell);                    
% 
%         avrws =  S * Rs * Ws ;
%         avrwx =  S * Rx * Wx  ;
%         avrwy =  S * Ry * Wy  ;
%         avrwz =  S * Rz * Wz  ;
%        
%         
%         grad_s = 1 / max(abs((((avrws))'* ((avrws) * ((m - mref)) ) ) ) );
%         
%          switch FLAG1
%         case 'SMOOTH_MOD'
%             
%             % Compute gradients then compute magnitude 
%             grad_x = 1 / max(abs( ( avrwx'* ( (avrwx *  m) ) ) ) );
%             grad_y = 1 / max(abs( ( avrwy'* ( (avrwy *  m) ) ) ) );
%             grad_z = 1 / max(abs( ( avrwz'* ( (avrwz *  m) ) ) ) );
% 
%         case 'SMOOTH_MOD_DIF'
% 
%             % Compute gradients then compute magnitude 
%             grad_x = 1 / max(abs(  (avrwx'* ( (avrwx * ( m - m_ref ) ) ) ) ));
%             grad_y = 1 / max(abs(  (avrwy'* ( (avrwy * ( m - m_ref ) ) ) ) ));
%             grad_z = 1 / max(abs(  (avrwz'* ( (avrwz * ( m - m_ref ) ) ) ) ));
% 
%          end
%         
%         aVRWs_pm = aVRWs_pm + sqrt( grad_s ) * T * Rs;
%         aVRWx_pm = aVRWx_pm + sqrt( grad_x ) * T * Rx;
%         aVRWy_pm = aVRWy_pm + sqrt( grad_y ) * T * Ry;
%         aVRWz_pm = aVRWz_pm + sqrt( grad_z ) * T * Rz;
%         
%         
%              
%     end
%     
%     aVRWs_pm =  sqrt( alpha(1) ) *aVRWs_pm * Wr * V  * Ws ;
%     aVRWx_pm =  sqrt( alpha(2) ) *aVRWx_pm * Wr * Vx * Wx  ;
%     aVRWy_pm =  sqrt( alpha(3) ) *aVRWy_pm * Wr * Vy * Wy  ;
%     aVRWz_pm =  sqrt( alpha(4) ) *aVRWz_pm * Wr * Vz * Wz  ;                       
%                         
%     sc_phis = phim(2) /(...
%             ((aVRWs_pm * m )' * (aVRWs_pm * m ))+...
%             (( (aVRWx_pm) * m )' * ( (aVRWx_pm) * ( m ) ))+...
%             (( (aVRWy_pm) * m )' * ( (aVRWy_pm) * ( m ) ))+...
%             (( (aVRWz_pm) * m )' * ( (aVRWz_pm) * ( m ) )));
%                  
%     switch FLAG1
%         case 'SMOOTH_MOD'
%             
%         sc_phigx = sc_phis;
%         sc_phigy = sc_phis;
%         sc_phigz = sc_phis;
% 
% 
%         case 'SMOOTH_MOD_DIF'
% 
%         sc_phigx = sc_phis;
%         sc_phigy = sc_phis;
%         sc_phigz = sc_phis;
% 
%     end
% 
% 
% else
%     
%     
%     % Form the matrices by including the weights and volume
%     aVRWs_pm=  sqrt( alpha(1) ) * Wr * V * Ws ;
%     aVRWx_pm =  sqrt( alpha(2) ) * Wr * Vx * Wx  ;
%     aVRWy_pm =  sqrt( alpha(3) ) * Wr * Vy * Wy  ;
%     aVRWz_pm =  sqrt( alpha(4) ) * Wr * Vz * Wz  ;
% 
% %     WstRsWs = ( aVRWs )' * ( aVRWs ) ;
% %     WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
% %     WztRzWz = ( aVRWz )' * ( aVRWz ) ;
%     
% end
% 
% % If more than one zone, average the scaling between lp,lq on transition
% if size(LP{2},1) ~= 1;
%         
%     aVRWs_pm = spdiags(sqrt( sc_phis * rvec ),0,mcell,mcell) * aVRWs_pm;
%     aVRWx_pm = spdiags(sqrt( sc_phigx * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx_pm;
%     aVRWy_pm = spdiags(sqrt( sc_phigy * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy_pm;
%     aVRWz_pm = spdiags(sqrt( sc_phigz * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz_pm;
%     
% % Otherwise just a scaler
% else
%     
%     aVRWs_pm = sqrt( sc_phis * LP{2}(1,5) ) * aVRWs_pm;
%     aVRWx_pm = sqrt( sc_phigx * ( 2.0 - LP{2}(1,5) ) ) * aVRWx_pm;
%     aVRWy_pm = sqrt( sc_phigy * ( 2.0 - LP{2}(1,5) ) ) * aVRWy_pm;
%     aVRWz_pm = sqrt( sc_phigz * ( 2.0 - LP{2}(1,5) ) ) * aVRWz_pm;
%     
% end
% 
% 
% 
% % Form the final model objective function
% MOF_pm = aVRWs_pm'*aVRWs_pm + aVRWx_pm'*aVRWx_pm + aVRWy_pm'*aVRWy_pm + aVRWz_pm'*aVRWz_pm;

%% S-COMPONENT

m = m_3C((1+mcell):2*mcell);
mref = mref_3C((1+mcell):2*mcell);

switch FLAG1
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude         
        dmx = Wx * m;
        dmy = Wy * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        
        dmx = Wx * (m - mref);
        dmy = Wy * (m - mref);
        dmz = Wz * (m - mref);
        
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

sc_phis = 1;
sc_phigx = 1;
sc_phigy = 1;
sc_phigz = 1;

aVRWs_s = spalloc(mcell,mcell,mcell);
aVRWx_s = spalloc(mcell,mcell,mcell);
aVRWy_s = spalloc(mcell,mcell,mcell);
aVRWz_s = spalloc(mcell,mcell,mcell);


if FLAG3~=0
    
    rvec = zeros(mcell,1);
    
    for ii = 1 : size(s{2},2)
        
        tt =t{2}(:,ii);
        
        rvec = rvec+t{2}(:,ii)*LP{2}(ii,5);
        
        tt(tt>0) = sqrt(tt(tt>0));
        
        T = spdiags((tt),0,mcell,mcell);
        S = spdiags(s{2}(:,ii),0,mcell,mcell);
        
        rs = 1 ./ ( abs(((m-mref))) .^( 2.0 - LP{2}(ii,1) ) + delta );
        rx = 1 ./ ( abs(GRADx) .^( 2.0 - LP{2}(ii,2) ) + delta ) ;
        ry = 1 ./ ( abs(GRADy) .^( 2.0 - LP{2}(ii,3) ) + delta ) ;
        rz = 1 ./ ( abs(GRADz) .^( 2.0 - LP{2}(ii,4) ) + delta ) ;
        
        
        Rs = spdiags( (rs) .^0.5,0,mcell,mcell);

        Rx = spdiags( (rx) .^0.5,0,mcell,mcell);

        Ry = spdiags( (ry) .^0.5,0,mcell,mcell);

        Rz = spdiags( (rz) .^0.5,0,mcell,mcell);                    

        avrws =  S * Rs * Ws ;
        avrwx =  S * Rx * Wx  ;
        avrwy =  S * Ry * Wy  ;
        avrwz =  S * Rz * Wz  ;
       
        
        grad_s = 1 / max(abs((((avrws))'* ((avrws) * ((m - mref)) ) ) ) );
        
         switch FLAG1
        case 'SMOOTH_MOD'
            
            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs( ( avrwx'* ( (avrwx *  m) ) ) ) );
            grad_y = 1 / max(abs( ( avrwy'* ( (avrwy *  m) ) ) ) );
            grad_z = 1 / max(abs( ( avrwz'* ( (avrwz *  m) ) ) ) );

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs(  (avrwx'* ( (avrwx * ( m - m_ref ) ) ) ) ));
            grad_y = 1 / max(abs(  (avrwy'* ( (avrwy * ( m - m_ref ) ) ) ) ));
            grad_z = 1 / max(abs(  (avrwz'* ( (avrwz * ( m - m_ref ) ) ) ) ));

         end
        
        aVRWs_s = aVRWs_s + sqrt( grad_s ) * T * Rs;
        aVRWx_s = aVRWx_s + sqrt( grad_x ) * T * Rx;
        aVRWy_s = aVRWy_s + sqrt( grad_y ) * T * Ry;
        aVRWz_s = aVRWz_s + sqrt( grad_z ) * T * Rz;
        
        
             
    end
    
    aVRWs_s =  sqrt( alpha(1) ) *aVRWs_s * Wr * V  * Ws ;
    aVRWx_s =  sqrt( alpha(2) ) *aVRWx_s * Wr * Vx * Wx  ;
    aVRWy_s =  sqrt( alpha(3) ) *aVRWy_s * Wr * Vy * Wy  ;
    aVRWz_s =  sqrt( alpha(4) ) *aVRWz_s * Wr * Vz * Wz  ;                       
                        
    sc_phis = phim(2) /(...
            ((aVRWs_s * (m-mref) )' * (aVRWs_s * (m-mref) ))+...
            (( (aVRWx_s) * m )' * ( (aVRWx_s) * ( m ) ))+...
            (( (aVRWy_s) * m )' * ( (aVRWy_s) * ( m ) ))+...
            (( (aVRWz_s) * m )' * ( (aVRWz_s) * ( m ) )));
                 
    switch FLAG1
        case 'SMOOTH_MOD'
            
        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;


        case 'SMOOTH_MOD_DIF'

        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;

    end


else
    
    
    % Form the matrices by including the weights and volume
    aVRWs_s =  sqrt( alpha(1) ) * Wr * V * Ws ;
    aVRWx_s =  sqrt( alpha(2) ) * Wr * Vx * Wx  ;
    aVRWy_s =  sqrt( alpha(3) ) * Wr * Vy * Wy  ;
    aVRWz_s =  sqrt( alpha(4) ) * Wr * Vz * Wz  ;

%     WstRsWs = ( aVRWs )' * ( aVRWs ) ;
%     WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
%     WztRzWz = ( aVRWz )' * ( aVRWz ) ;
    
end

% If more than one zone, average the scaling between lp,lq on transition
if size(LP{2},1) ~= 1;
        
    aVRWs_s = spdiags(sqrt( sc_phis * rvec ),0,mcell,mcell) * aVRWs_s;
    aVRWx_s = spdiags(sqrt( sc_phigx * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx_s;
    aVRWy_s = spdiags(sqrt( sc_phigy * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy_s;
    aVRWz_s = spdiags(sqrt( sc_phigz * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz_s;
    
% Otherwise just a scaler
else
    
    aVRWs_s = sqrt( sc_phis * LP{2}(1,5) ) * aVRWs_s;
    aVRWx_s = sqrt( sc_phigx * ( 2.0 - LP{2}(1,5) ) ) * aVRWx_s;
    aVRWy_s = sqrt( sc_phigy * ( 2.0 - LP{2}(1,5) ) ) * aVRWy_s;
    aVRWz_s = sqrt( sc_phigz * ( 2.0 - LP{2}(1,5) ) ) * aVRWz_s;
    
end



% Form the final model objective function
MOF_s = aVRWs_s'*aVRWs_s + aVRWx_s'*aVRWx_s + aVRWy_s'*aVRWy_s + aVRWz_s'*aVRWz_s;

%% T-COMPONENT

m = m_3C((1+2*mcell):3*mcell);
mref = mref_3C((1+2*mcell):3*mcell);

switch FLAG1
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude         
        dmx = Wx * m;
        dmy = Wy * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        
        dmx = Wx * (m - mref);
        dmy = Wy * (m - mref);
        dmz = Wz * (m - mref);
        
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

sc_phis = 1;
sc_phigx = 1;
sc_phigy = 1;
sc_phigz = 1;

aVRWs_t = spalloc(mcell,mcell,mcell);
aVRWx_t = spalloc(mcell,mcell,mcell);
aVRWy_t = spalloc(mcell,mcell,mcell);
aVRWz_t = spalloc(mcell,mcell,mcell);


if FLAG3~=0
    
    rvec = zeros(mcell,1);
    
    for ii = 1 : size(s{2},2)
        
        tt =t{2}(:,ii);
        
        rvec = rvec+t{2}(:,ii)*LP{2}(ii,5);
        
        tt(tt>0) = sqrt(tt(tt>0));
        
        T = spdiags((tt),0,mcell,mcell);
        S = spdiags(s{2}(:,ii),0,mcell,mcell);
        
        rs = 1 ./ ( abs(((m-mref))) .^( 2.0 - LP{2}(ii,1) ) + delta );
        rx = 1 ./ ( abs(GRADx) .^( 2.0 - LP{2}(ii,2) ) + delta ) ;
        ry = 1 ./ ( abs(GRADy) .^( 2.0 - LP{2}(ii,3) ) + delta ) ;
        rz = 1 ./ ( abs(GRADz) .^( 2.0 - LP{2}(ii,4) ) + delta ) ;
        
        
        Rs = spdiags( (rs) .^0.5,0,mcell,mcell);

        Rx = spdiags( (rx) .^0.5,0,mcell,mcell);

        Ry = spdiags( (ry) .^0.5,0,mcell,mcell);

        Rz = spdiags( (rz) .^0.5,0,mcell,mcell);                    

        avrws =  S * Rs * Ws ;
        avrwx =  S * Rx * Wx  ;
        avrwy =  S * Ry * Wy  ;
        avrwz =  S * Rz * Wz  ;
       
        
        grad_s = 1 / max(abs((((avrws))'* ((avrws) * ((m - mref)) ) ) ) );
        
         switch FLAG1
        case 'SMOOTH_MOD'
            
            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs( ( avrwx'* ( (avrwx *  m) ) ) ) );
            grad_y = 1 / max(abs( ( avrwy'* ( (avrwy *  m) ) ) ) );
            grad_z = 1 / max(abs( ( avrwz'* ( (avrwz *  m) ) ) ) );

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            grad_x = 1 / max(abs(  (avrwx'* ( (avrwx * ( m - m_ref ) ) ) ) ));
            grad_y = 1 / max(abs(  (avrwy'* ( (avrwy * ( m - m_ref ) ) ) ) ));
            grad_z = 1 / max(abs(  (avrwz'* ( (avrwz * ( m - m_ref ) ) ) ) ));

         end
        
        aVRWs_t = aVRWs_t + sqrt( grad_s ) * T * Rs;
        aVRWx_t = aVRWx_t + sqrt( grad_x ) * T * Rx;
        aVRWy_t = aVRWy_t + sqrt( grad_y ) * T * Ry;
        aVRWz_t = aVRWz_t + sqrt( grad_z ) * T * Rz;
        
        
             
    end
    
    aVRWs_t =  sqrt( alpha(1) ) *aVRWs_t * Wr * V  * Ws ;
    aVRWx_t =  sqrt( alpha(2) ) *aVRWx_t * Wr * Vx * Wx  ;
    aVRWy_t =  sqrt( alpha(3) ) *aVRWy_t * Wr * Vy * Wy  ;
    aVRWz_t =  sqrt( alpha(4) ) *aVRWz_t * Wr * Vz * Wz  ;                       
                        
    sc_phis = phim(3) /(...
            ((aVRWs_t * (m-mref) )' * (aVRWs_t * (m-mref) ))+...
            (( (aVRWx_t) * m )' * ( (aVRWx_t) * ( m ) ))+...
            (( (aVRWy_t) * m )' * ( (aVRWy_t) * ( m ) ))+...
            (( (aVRWz_t) * m )' * ( (aVRWz_t) * ( m ) )));
                 
    switch FLAG1
        case 'SMOOTH_MOD'
            
        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;


        case 'SMOOTH_MOD_DIF'

        sc_phigx = sc_phis;
        sc_phigy = sc_phis;
        sc_phigz = sc_phis;

    end


else
    
    
    % Form the matrices by including the weights and volume
    aVRWs_t =  sqrt( alpha(1) ) * Wr * V * Ws ;
    aVRWx_t =  sqrt( alpha(2) ) * Wr * Vx * Wx  ;
    aVRWy_t =  sqrt( alpha(3) ) * Wr * Vy * Wy  ;
    aVRWz_t =  sqrt( alpha(4) ) * Wr * Vz * Wz  ;

%     WstRsWs = ( aVRWs )' * ( aVRWs ) ;
%     WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
%     WztRzWz = ( aVRWz )' * ( aVRWz ) ;
    
end

% If more than one zone, average the scaling between lp,lq on transition
if size(LP{2},1) ~= 1;
        
    aVRWs_t = spdiags(sqrt( sc_phis * rvec ),0,mcell,mcell) * aVRWs_t;
    aVRWx_t = spdiags(sqrt( sc_phigx * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx_t;
    aVRWy_t = spdiags(sqrt( sc_phigy * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy_t;
    aVRWz_t = spdiags(sqrt( sc_phigz * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz_t;
    
% Otherwise just a scaler
else
    
    aVRWs_t = sqrt( sc_phis * LP{2}(1,5) ) * aVRWs_t;
    aVRWx_t = sqrt( sc_phigx * ( 2.0 - LP{2}(1,5) ) ) * aVRWx_t;
    aVRWy_t = sqrt( sc_phigy * ( 2.0 - LP{2}(1,5) ) ) * aVRWy_t;
    aVRWz_t = sqrt( sc_phigz * ( 2.0 - LP{2}(1,5) ) ) * aVRWz_t;
    
end



% Form the final model objective function
MOF_t = aVRWs_t'*aVRWs_t + aVRWx_t'*aVRWx_t + aVRWy_t'*aVRWy_t + aVRWz_t'*aVRWz_t;



%% Build the complete model objective function
MOF = blkdiag( MOF_pp , MOF_s , MOF_t ); 

if FLAG3~=0
    % Scale phim
    sc_phim = sum(phim) / (m_3C'*MOF*m_3C);
    
else
    
    sc_phim = 1;
    
end

MOF = sc_phim * MOF;
aVRWs = sqrt(sc_phim) * blkdiag( aVRWs_p , aVRWs_s , aVRWs_t );
aVRWx = sqrt(sc_phim) * blkdiag( aVRWx_p , aVRWx_s , aVRWx_t );
aVRWy = sqrt(sc_phim) * blkdiag( aVRWy_p , aVRWy_s , aVRWy_t );
aVRWz = sqrt(sc_phim) * blkdiag( aVRWz_p , aVRWz_s , aVRWz_t );
