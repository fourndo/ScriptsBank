function [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D_v2(m,mref,phim,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,p,q,r,FLAG1,FLAG2,delta)
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


mcell = nx*nz;

Wr = spdiags(wr,0,mcell,mcell);

wrx = reshape(wr,nz,nx); wrx = wrx(:,2:end); wrx = wrx(:);
wrz = reshape(wr,nz,nx); wrz = wrz(2:end,:); wrz = wrz(:);

% V = spdiags((v),0,mcell,mcell);
Wrx = spdiags(wrx,0,length(wrx),length(wrx));
Wrz = spdiags(wrz,0,length(wrz),length(wrz));

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
acfx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx-1,nx);
Acfx = kron( acfx , speye(nz) );


acfz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz-1,nz);
Acfz = kron(  speye(nx) , acfz );

% Create averaging operator face to center
afcz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz,nz+1);afcz = afcz(:,2:end-1);
Afcz = kron( speye(nx) , afcz );

afcx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx,nx+1);afcx = afcx(:,2:end-1);
Afcx = kron( afcx , speye(nz) );

% Augmente system if dealing with 3-component kernel
switch FLAG1
    case 'GRADm'
        % Compute gradients then compute magnitude 
        GRADx = dmx;
        GRADz = dmz;
        
    case 'lGRADml'
        
        % Compute gradients then compute magnitude 
        GRADx = sqrt(dmx.^2 + (Acfx * Afcz * dmz).^2);
        GRADz = sqrt((Acfz * Afcx * dmx ).^2 + dmz.^2);
        
    otherwise
        fprintf('FLAG must be GRADm | lGRADml\n')
        fprintf('Please revise input\n')
        return
end


rx = 1./ ( abs(GRADx) .^( 2.0 - q ) + delta ) ;
%                     rx = rx / max(rx) + 1e-6;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

rz = 1./ ( abs(GRADz) .^( 2.0 - q ) + delta ) ;
%                     rz = rz / max(rz) + 1e-6;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                    

aVRWs = V * Rs * Ws ;
aVRWx = Vx * Rx * Wx  ;
aVRWz = Vz * Rz * Wz  ;

% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs )' * ( aVRWs ) ;
WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
WztRzWz = ( aVRWz )' * ( aVRWz ) ;

phi = ( m - mref )' *...
      ( alpha(1) * ((Wr * V * Ws  )' * (Wr * V * Ws  ) * (m - mref) )  +...
        alpha(2) * ((Wrx * Vx * Wx  )' * (Wrx * Vx * Wx  ) * (m - mref) ) +...
        alpha(3) * ((Wrz * Vz * Wz  )' * (Wrz * Vz * Wz  ) * (m - mref) ) );


if p==2 && q==2
    
    scale(1) = 1 ;
    scale_g  = 1 ;
    
else
    
    phis = ((m - mref)' * WstRsWs * (m - mref) );

    % Scale the smallest q-norm on the smallest norm
    scale(1) = 0.5 * (phim) / (alpha(1) * phis);

    % Scale the gradient p-norm on the smallest q-norm
    switch FLAG2
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude 
        scale_g = 0.5 * phim /...
        ( m' * (alpha(2) * WxtRxWx +...
                alpha(3) * WztRzWz) * m );
        
    case 'SMOOTH_MOD_DIF'
        
        % Compute gradients then compute magnitude 
        scale_g = 0.5 * phim /...
        ( (m - mref)' * (alpha(2) * WxtRxWx +...
                        alpha(3) * WztRzWz) * (m - mref) );
        

    end
    
end

% Final scales
scale(1) =  r * scale(1) * alpha(1);
scale(2) = (2.0 - r) * scale_g * alpha(2);
scale(3) = (2.0 - r) * scale_g * alpha(3);

aVRWs = sqrt(scale(1)) * Wr * aVRWs;
aVRWx = sqrt(scale(2)) * Wrx * aVRWx;
aVRWz = sqrt(scale(3)) * Wrz * aVRWz;

% Form the final model objective function
MOF =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WztRzWz ) ;

% Scale MOF
mu = phi / ((m- mref )' * MOF * (m- mref));

MOF = mu * MOF;