function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF(m,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,p,q,r)
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

rs = 1./ ( abs(Ws * m) .^( 2.0 - p ) + delta );
%                     rs = rs / max(rs) + 1e-6;
Rs = spdiags( rs.^ 0.5 ,0,mcell,mcell);

rx = 1./ ( abs(Wx * m) .^( 2.0 - q ) + delta ) ;
%                     rx = rx / max(rx) + 1e-6;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

ry = 1./ ( abs(Wy * m) .^( 2.0 - q ) + delta ) ;
%                     ry = ry / max(ry) + 1e-6;
Ry = spdiags( ry .^0.5,0,size(Wy,1),size(Wy,1));

rz = 1./ ( abs(Wz * m) .^( 2.0 - q ) + delta ) ;
%                     rz = rz / max(rz) + 1e-6;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                    

aVRWs = V * Rs * Ws;
aVRWx = Vx * Rx * Wx;
aVRWy = Vy * Ry * Wy;
aVRWz = Vz * Rz * Wz;
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
    scale_g = scale(1) * alpha(1) * phis /...
        ( m' * (alpha(2) * WxtRxWx +...
        alpha(3) * WytRyWy +...
        alpha(4) * WztRzWz) * m );
end

% Final scales
scale(1) = r * scale(1) * alpha(1);
scale(2) = (2.0 - r) * scale_g * alpha(2);
scale(3) = (2.0 - r) * scale_g * alpha(3);
scale(4) = (2.0 - r) * scale_g * alpha(4);

aVRWs = sqrt(scale(1)) * aVRWs;
aVRWx = sqrt(scale(2)) * aVRWx;
aVRWy = sqrt(scale(3)) * aVRWy;
aVRWz = sqrt(scale(4)) * aVRWz;

% Form the final model objective function
MOF =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WytRyWy +...
scale(4) * WztRzWz ) ;