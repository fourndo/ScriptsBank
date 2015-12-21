function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD_v2(m,mref,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,p,q,s,FLAG)
% Function get_lp_WRW(Ws,Wx,Wy,Wz,p,q,r)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). The lp-norm is approximated through a lag diffusivity approach
% where:
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
delta = 1e-10;
mcell = length(m);

rs = 1./ ( abs(Ws * (m-mref)) .^( 2.0 - p ) + delta );
%                     rs = rs / max(rs) + 1e-6;
Rs = spdiags( rs.^ 0.5 ,0,mcell,mcell);

% Compute gradients then compute magnitude 
switch FLAG
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude 
        dmx = Wx * m;
        dmy = Wy * m;
        dmz = Wz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        % Compute gradients then compute magnitude 
        dmx = Wx * (m - mref);
        dmy = Wy * (m - mref);
        dmz = Wz * (m - mref);
        
    otherwise
        fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
        fprintf('Please revise input\n')
        return
end

% Create averaging operator center to face
% acfx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx-1,nx);
% Acfx = kron( kron( speye(ny) , acfx ), speye(nz) );
% 
% acfy = spdiags (ones(ny+1,1)*[0.5,0.5],[0,1],ny-1,ny);
% Acfy = kron( kron( acfy , speye(nx) ), speye(nz) );
% 
% acfz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz-1,nz);
% Acfz = kron( kron( speye(ny) , speye(nx) ), acfz );
% 
% % Create averaging operator face to center
% afcz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz,nz+1);afcz = afcz(:,2:end-1);
% Afcz = kron( kron( speye(ny) , speye(nx) ), afcz );
% 
% afcy = spdiags (ones(ny+1,1)*[0.5,0.5],[0,1],ny,ny+1);afcy = afcy(:,2:end-1);
% Afcy = kron( kron( afcy , speye(nx) ), speye(nz) );
% 
% afcx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx,nx+1);afcx = afcx(:,2:end-1);
% Afcx = kron( kron( speye(ny) , afcx ), speye(nz) );

% GRADx = sqrt(dmx.^2 + (Acfx * Afcy * dmy).^2 + (Acfx * Afcz * dmz).^2);
% GRADy = sqrt((Acfy * Afcx * dmx ).^2 +  dmy.^2 + (Acfy * Afcz * dmz).^2);
% GRADz = sqrt((Acfz * Afcx * dmx ).^2 + (Acfz * Afcy * dmy).^2 + dmz.^2);

GRAD = sqrt(dmx.^2 + (dmy).^2 + (dmz).^2);

% rx = 1./ ( abs(GRADx) .^( 2.0 - q ) + delta ) ;
% %                     rx = rx / max(rx) + 1e-6;
% Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));
% 
% ry = 1./ ( abs(GRADy) .^( 2.0 - q ) + delta ) ;
% %                     ry = ry / max(ry) + 1e-6;
% Ry = spdiags( ry .^0.5,0,size(Wy,1),size(Wy,1));
% 
% rz = 1./ ( abs(GRADz) .^( 2.0 - q ) + delta ) ;
% %                     rz = rz / max(rz) + 1e-6;
% Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                    

r = 1./ ( abs(GRAD) .^( 2.0 - q ) + delta ) ;
%                     rx = rx / max(rx) + 1e-6;
R = spdiags( r .^0.5,0,mcell,mcell);

aVRWs = V * Rs * Ws;
aVRWx = Vx * R * Wx;
aVRWy = Vy * R * Wy;
aVRWz = Vz * R * Wz;
% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs )' * ( aVRWs ) ;
WxtRxWx = ( aVRWx )' * ( aVRWx ) ;
WytRyWy = ( aVRWy )' * ( aVRWy ) ;
WztRzWz = ( aVRWz )' * ( aVRWz ) ;

if p==2 && q==2
    
    scale(1) = 1 ;
    scale_g  = 1 ;
    
else
    
    phis = (m' * WstRsWs * m );

    % Scale the smallest q-norm on the smallest norm
    scale(1) = (m' * (V*Ws)' * (V * Ws) * m ) / (phis);

    % Scale the gradient p-norm on the smallest q-norm
    switch FLAG
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude 
        scale_g = scale(1) * alpha(1) * phis /...
        ( m' * (alpha(2) * WxtRxWx +...
        alpha(3) * WytRyWy +...
        alpha(4) * WztRzWz) * m );
        
    case 'SMOOTH_MOD_DIF'
        
        % Compute gradients then compute magnitude 
        scale_g = scale(1) * alpha(1) * phis /...
        ( (m - mref)' * (alpha(2) * WxtRxWx +...
        alpha(3) * WytRyWy +...
        alpha(4) * WztRzWz) * (m - mref) );
        

    end
    
end

% Final scales
scale(1) = s * scale(1) * alpha(1);
scale(2) = (2.0 - s) * scale_g * alpha(2);
scale(3) = (2.0 - s) * scale_g * alpha(3);
scale(4) = (2.0 - s) * scale_g * alpha(4);

aVRWs = sqrt(scale(1)) * aVRWs;
aVRWx = sqrt(scale(2)) * aVRWx;
aVRWy = sqrt(scale(3)) * aVRWy;
aVRWz = sqrt(scale(4)) * aVRWz;

% Form the final model objective function
MOF =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WytRyWy +...
scale(4) * WztRzWz ) ;