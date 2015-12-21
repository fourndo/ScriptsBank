function [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D_v3(m,mref,phim,s,T,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,p,qx,qz,r,FLAG1,FLAG2,FLAG3,delta)
% Function get_lp_WRW(Ws,Wx,Wy,Wz,p,q,r)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). The lp-norm is approximated through a lag diffusivity approach
% where:
% 
% V3: Pass on cell-based p,q and r 
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


mcell = nx*nz;

Wr = spdiags(wr,0,mcell,mcell);

% wrx = reshape(wr,nz,nx); wrx = wrx(:,2:end); wrx = wrx(:);
% wrz = reshape(wr,nz,nx); wrz = wrz(2:end,:); wrz = wrz(:);

% V = spdiags((v),0,mcell,mcell);
% Wrx = spdiags(wrx,0,length(wrx),length(wrx));
% Wrz = spdiags(wrz,0,length(wrz),length(wrz));

% rs = 1./ ( abs(Ws * Wr *  (m - mref)) .^( 2.0 - p ) + delta );
rs = 1./ ( abs(Ws * (m - mref)) .^( 2.0 - p ) + delta );
Rs = spdiags( rs.^ 0.5 ,0,length(rs),length(rs));

switch FLAG2
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude 
%         dmx = Wx * Wr * m;
%         dmy = Wy * Wr * m;
%         dmz = Wz * Wr * m;
        
        dmx = Wx * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        % Compute gradients then compute magnitude 
%         dmx = Wx * Wr * (m - mref);
%         dmy = Wy * Wr * (m - mref);
%         dmz = Wz * Wr * (m - mref);
        
        dmx = Wx * (m - mref);
        dmz = Wz * (m - mref);
        
    otherwise
        fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
        fprintf('Please revise input\n')
        return
end

% Create averaging operator center to face
% acfx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx-1,nx);
% Acfx = kron( acfx , speye(nz) );
% 
% 
% acfz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz-1,nz);
% Acfz = kron(  speye(nx) , acfz );
% 
% % Create averaging operator face to center
% afcz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz,nz+1);afcz = afcz(:,2:end-1);
% Afcz = kron( speye(nx) , afcz );
% 
% afcx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx,nx+1);afcx = afcx(:,2:end-1);
% Afcx = kron( afcx , speye(nz) );

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

rx = 1./ ( abs(GRADx) .^( 2.0 - qx ) + delta ) ;
%                     rx = rx / max(rx) + 1e-6;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

rz = 1./ ( abs(GRADz) .^( 2.0 - qz ) + delta ) ;
%                     rz = rz / max(rz) + 1e-6;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                    

aVRWs =  Wr * V * Rs * Ws ;
aVRWx =  Wr * Vx * Rx * Wx  ;
aVRWz =  Wr * Vz * Rz * Wz  ;



% phi = ( alpha(1) * (Wr * V * Ws  )' * (Wr * V * Ws  )  +...
%         alpha(2) * (Wr * Vx * Wx  )' * (Wr * Vx * Wx) +...
%         alpha(3) * (Wr * Vz * Wz  )' * (Wr * Vz * Wz) );



% Scale gradient
% dphi2_dm = ( (V*Ws)' * V*Ws ) * (m - mref);
% dphip_dm = (aVRWs' * aVRWs) * (m - mref);
% 
% scl_dphim = max(dphi2_dm)/max(dphip_dm);
% 
% [mm,I] = sort(m);    
% 
% aVRWs = scl_dphim * V * Rs * Ws + spdiags(sign(dphi2_dm)*min(abs(dphi2_dm)),0,mcell,mcell);    
 

% figure;plot(mm,scl_dphim*dphip_dm(I) + sign(dphi2_dm)*min(abs(dphi2_dm)));

% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs )' * ( aVRWs ) ;
WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
WztRzWz = ( aVRWz )' * ( aVRWz ) ;

scale_s = 1;
scale_g = 1;

scale_gs = zeros(mcell,1) ;
scale_gx  = zeros(mcell,1) ;
scale_gz  = zeros(mcell,1) ;

if FLAG3~=0
        
    for ii = 1 : size(s,2)

        S= spdiags(s{ii},0,mcell,mcell);
        
        max_gs = 1 / max(abs( S * (aVRWs'* (aVRWs * (S * (m - mref)) ) ) ) );
        scale_gs = scale_gs + T{ii} * max_gs;

        % Scale the gradient p-norm on the smallest q-norm
        switch FLAG2
        case 'SMOOTH_MOD'
            % Compute gradients then compute magnitude 
            
%             dmdx = max(abs(Z * alpha(2) * WxtRxWx * Z * m ));
%                 
%             dmdz = max(abs(Z * alpha(3) * WztRzWz * Z * m ));
            
%             scale_gx = ax / max(abs( ax * WxtRxWx * invmod ));
            max_gx = 1 / max(abs(  (S *aVRWx'* (aVRWx * (S * m) ) ) ) );
            max_gz = 1 / max(abs(  (S *aVRWz'* (aVRWz * (S * m) ) ) ) );
            
            scale_gx = scale_gx + T{ii} * max_gx;
            scale_gz = scale_gz + T{ii} * max_gz;
%             scale_x(z{ii}==1) = 1 / dmdx;
%             scale_z(z{ii}==1) = 1 / dmdz;
%             ( m' * Z * (alpha(2) * WxtRxWx +...
%                     alpha(3) * WztRzWz) * Z * m );

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            max_gx = 1 / max(abs(  (S *aVRWx'* (aVRWx * (S * ( m - m_ref ) ) ) ) ));
            max_gz = 1 / max(abs(  (S *aVRWz'* (aVRWz * (S * ( m - m_ref ) ) ) ) ));
            
            scale_gx = scale_gx + T{ii} * max_gx;
            scale_gz = scale_gz + T{ii} * max_gz;


        end
        
    end
    
    aVRWs = spdiags(sqrt( alpha(1) * scale_gs ),0,mcell,mcell) * aVRWs;
    aVRWx = spdiags(sqrt( alpha(2) * scale_gx ),0,mcell,mcell) * aVRWx;
    aVRWz = spdiags(sqrt( alpha(3) * scale_gz ),0,mcell,mcell) * aVRWz;
    
    scale_s = 0.5 * phim /...
            ( m' * ( aVRWs'*aVRWs) * m );
                 
    switch FLAG2
        case 'SMOOTH_MOD'
            
            scale_g = 0.5 * phim /...
            ( m' * ( aVRWx'*aVRWx +...
                     aVRWz' * aVRWz) * m );

        case 'SMOOTH_MOD_DIF'

            % Compute gradients then compute magnitude 
            scale_g = 0.5 * phim /...
            ( (m - mref)' * ( aVRWx'*aVRWx +...
                     aVRWx' * aVRWx) * (m - mref) );

    end


    
end

% Final scales
% scale(1) =  scale_s * alpha(1);
% scale(2) =  scale_g * alpha(2);
% scale(3) =  scale_g * alpha(3);
  
aVRWs = spdiags(sqrt( scale_s .* r),0,mcell,mcell) * aVRWs;
aVRWx = spdiags(sqrt( scale_g .* ( 2.0 - r ) ) ,0,mcell,mcell) * aVRWx;
aVRWz = spdiags(sqrt( scale_g .* ( 2.0 - r ) ) ,0,mcell,mcell) * aVRWz;
    
% end    
set(figure(10), 'Position', [50 200 750 750])
figure(10);imagesc(reshape(WstRsWs*m,nz,nx));colorbar
set(figure(11), 'Position', [50 200 750 750])
figure(11);imagesc(reshape(WxtRxWx*m,nz,nx));colorbar

set(figure(12), 'Position', [50 200 750 750])
figure(12);imagesc(reshape(aVRWs'*aVRWs*m,nz,nx));colorbar
set(figure(13), 'Position', [50 200 750 750])
figure(13);imagesc(reshape((aVRWx'*aVRWx +aVRWz'*aVRWz)*m,nz,nx));colorbar

% Form the final model objective function
MOF = aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWz'*aVRWz;

% Scale MOF
% if FLAG3~=0
%     
%     mu = phim / ((m- mref )' * MOF * (m- mref));
%     
% else
%     
%     mu =1 ;
%     
% end
% MOF = mu * MOF;
% aVRWs = sqrt(mu) * aVRWs; 
% aVRWx = sqrt(mu) * aVRWx; 
% aVRWz = sqrt(mu) * aVRWz; 
