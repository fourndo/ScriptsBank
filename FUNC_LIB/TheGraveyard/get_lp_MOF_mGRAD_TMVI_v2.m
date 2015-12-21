function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_mGRAD_TMVI_v2(m_pst,mref_pst,nx,ny,nz,V,Ws,Vx,Wx,Vy,Wy,Vz,Wz,alpha,p,q,l,FLAG)
% Function get_lp_WRW(Ws,Wx,Wy,Wz,p,q,r)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). The lp-norm is approximated through a lag diffusivity approach
% where:
% 
% V3: ALL SCALED ON PHI_l2 (2014-08-02) 
% TMVI: The s,t orthogonal direction are computed together.
%
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
mcell = nx*ny*nz;

%% Create averaging operator center to face for gradient
acfx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx-1,nx);
Acfx = kron( kron( speye(ny) , acfx ), speye(nz) );

acfy = spdiags (ones(ny+1,1)*[0.5,0.5],[0,1],ny-1,ny);
Acfy = kron( kron( acfy , speye(nx) ), speye(nz) );

acfz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz-1,nz);
Acfz = kron( kron( speye(ny) , speye(nx) ), acfz );

% Create averaging operator face to center
afcz = spdiags (ones(nz+1,1)*[0.5,0.5],[0,1],nz,nz+1);afcz = afcz(:,2:end-1);
Afcz = kron( kron( speye(ny) , speye(nx) ), afcz );

afcy = spdiags (ones(ny+1,1)*[0.5,0.5],[0,1],ny,ny+1);afcy = afcy(:,2:end-1);
Afcy = kron( kron( afcy , speye(nx) ), speye(nz) );

afcx = spdiags (ones(nx+1,1)*[0.5,0.5],[0,1],nx,nx+1);afcx = afcx(:,2:end-1);
Afcx = kron( kron( speye(ny) , afcx ), speye(nz) );

% % Augmente system if dealing with 3-component kernel
% if size(Ws,1) == 3*mcell
%     
%     Acfx = kron(speye(3), Acfx);
%     Acfy = kron(speye(3), Acfy);
%     Acfz = kron(speye(3), Acfz);
%     
%     Afcx = kron(speye(3), Afcx);
%     Afcy = kron(speye(3), Afcy);
%     Afcz = kron(speye(3), Afcz);
%   
% end

%% PRIMARY FIELD COMPONENT
m = m_pst(1:mcell);
mref = mref_pst(1:mcell);

% SMALLNESS TERM
Wsm = Ws * (m - mref);

rs = 1./ ( abs(Wsm) .^( 2.0 - p(1) ) + delta );

Rs = spdiags( rs.^0.5 , 0, mcell,mcell);

% GRADIENT TERM
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

GRADx = sqrt(dmx.^2 + (Acfx * Afcy * dmy).^2 + (Acfx * Afcz * dmz).^2);
GRADy = sqrt((Acfy * Afcx * dmx ).^2 +  dmy.^2 + (Acfy * Afcz * dmz).^2);
GRADz = sqrt((Acfz * Afcx * dmx ).^2 + (Acfz * Afcy * dmy).^2 + dmz.^2);


rx = 1./ ( abs(GRADx) .^( 2 - q(1) ) + delta ) ;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

ry = 1./ ( abs(GRADy) .^( 2 - q(1) ) + delta ) ;
Ry = spdiags( ry .^0.5,0,size(Wy,1),size(Wy,1));

rz = 1./ ( abs(GRADz) .^( 2 - q(1) ) + delta ) ;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));


               

aVRWs_p = V * Rs * Ws;
aVRWx_p = Vx * Rx * Wx;
aVRWy_p = Vy * Ry * Wy;
aVRWz_p = Vz * Rz * Wz;

% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs_p )' * ( aVRWs_p ) ;
WxtRxWx = ( aVRWx_p )' * ( aVRWx_p ) ;
WytRyWy = ( aVRWy_p )' * ( aVRWy_p ) ;
WztRzWz = ( aVRWz_p )' * ( aVRWz_p ) ;

    
phis = alpha(1) * ((m - mref)' * WstRsWs * (m - mref) );

phi = (m - mref)' * ( alpha(1) * ((V * Ws)' * (V * Ws) * (m - mref) )  +...
                    alpha(2) * ((Vx * Wx)' * (Vx * Wx) * (m - mref) ) +...
                    alpha(3) * ((Vy * Wy)' * (Vy * Wy) * (m - mref) ) +...
                    alpha(4) * ((Vz * Wz)' * (Vz * Wz) * (m - mref) ) );

% Scale the smallest q-norm on the smallest norm
scale(1) = phi / (phis);

% Scale the gradient p-norm on the smallest q-norm
switch FLAG
case 'SMOOTH_MOD'
    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( m' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * m );

case 'SMOOTH_MOD_DIF'

    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( (m - mref)' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * (m - mref) );


end
    
% end

% Final scales
scale(1) = l(1) * scale(1) * alpha(1);
scale(2) = (2.0 - l(1)) * scale_g * alpha(2);
scale(3) = (2.0 - l(1)) * scale_g * alpha(3);
scale(4) = (2.0 - l(1)) * scale_g * alpha(4);

aVRWs_p = sqrt(scale(1)) * aVRWs_p;
aVRWx_p = sqrt(scale(2)) * aVRWx_p;
aVRWy_p = sqrt(scale(3)) * aVRWy_p;
aVRWz_p = sqrt(scale(4)) * aVRWz_p;

MOF_p =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WytRyWy +...
scale(4) * WztRzWz ) ;

%% Secondary component
m = m_pst( (1+mcell) : 2*mcell );
mref = mref_pst( (1+mcell) : 2*mcell );

% SMALLNESS TERM
Wsm = Ws * (m - mref);

rs = 1./ ( abs(Wsm) .^( 2.0 - p(2) ) + delta );

Rs = spdiags( rs.^0.5 , 0, mcell, mcell );

% GRADIENT TERM
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

GRADx = sqrt(dmx.^2 + (Acfx * Afcy * dmy).^2 + (Acfx * Afcz * dmz).^2);
GRADy = sqrt((Acfy * Afcx * dmx ).^2 +  dmy.^2 + (Acfy * Afcz * dmz).^2);
GRADz = sqrt((Acfz * Afcx * dmx ).^2 + (Acfz * Afcy * dmy).^2 + dmz.^2);


rx = 1./ ( abs(GRADx) .^( 2.0 - q(2) ) + delta ) ;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

ry = 1./ ( abs(GRADy) .^( 2.0 - q(2) ) + delta ) ;
Ry = spdiags( ry .^0.5,0,size(Wy,1),size(Wy,1));

rz = 1./ ( abs(GRADz) .^( 2.0 - q(2)) + delta ) ;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                

aVRWs_s = V * Rs * Ws;
aVRWx_s = Vx * Rx * Wx;
aVRWy_s = Vy * Ry * Wy;
aVRWz_s = Vz * Rz * Wz;

% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs_s )' * ( aVRWs_s ) ;
WxtRxWx = ( aVRWx_s )' * ( aVRWx_s ) ;
WytRyWy = ( aVRWy_s )' * ( aVRWy_s ) ;
WztRzWz = ( aVRWz_s )' * ( aVRWz_s ) ;

    
phis = alpha(1) * ((m - mref)' * WstRsWs * (m - mref) );

phi = (m - mref)' * ( alpha(1) * ((V * Ws)' * (V * Ws) * (m - mref) )  +...
                    alpha(2) * ((Vx * Wx)' * (Vx * Wx) * (m - mref) ) +...
                    alpha(3) * ((Vy * Wy)' * (Vy * Wy) * (m - mref) ) +...
                    alpha(4) * ((Vz * Wz)' * (Vz * Wz) * (m - mref) ) );

% Scale the smallest q-norm on the smallest norm
scale(1) = phi / (phis);

% Scale the gradient p-norm on the smallest q-norm
switch FLAG
case 'SMOOTH_MOD'
    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( m' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * m );

case 'SMOOTH_MOD_DIF'

    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( (m - mref)' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * (m - mref) );


end
    
% end

% Final scales
scale(1) = l(2) * scale(1) * alpha(1);
scale(2) = (2.0 - l(2)) * scale_g * alpha(2);
scale(3) = (2.0 - l(2)) * scale_g * alpha(3);
scale(4) = (2.0 - l(2)) * scale_g * alpha(4);

aVRWs_s = sqrt(scale(1)) * aVRWs_s;
aVRWx_s = sqrt(scale(2)) * aVRWx_s;
aVRWy_s = sqrt(scale(3)) * aVRWy_s;
aVRWz_s = sqrt(scale(4)) * aVRWz_s;

MOF_s =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WytRyWy +...
scale(4) * WztRzWz ) ;

%% Tertiary component
m = m_pst( (1+2*mcell) : 3*mcell );
mref = mref_pst( (1+2*mcell) : 3*mcell );

% SMALLNESS TERM
Wsm = Ws * (m - mref);

rs = 1./ ( abs(Wsm) .^( 2.0 - p(3) ) + delta );

Rs = spdiags( rs.^0.5 , 0, mcell, mcell );

% GRADIENT TERM
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

GRADx = sqrt(dmx.^2 + (Acfx * Afcy * dmy).^2 + (Acfx * Afcz * dmz).^2);
GRADy = sqrt((Acfy * Afcx * dmx ).^2 +  dmy.^2 + (Acfy * Afcz * dmz).^2);
GRADz = sqrt((Acfz * Afcx * dmx ).^2 + (Acfz * Afcy * dmy).^2 + dmz.^2);


rx = 1./ ( abs(GRADx) .^( 2.0 - q(3) ) + delta ) ;
Rx = spdiags( rx .^0.5,0,size(Wx,1),size(Wx,1));

ry = 1./ ( abs(GRADy) .^( 2.0 - q(3) ) + delta ) ;
Ry = spdiags( ry .^0.5,0,size(Wy,1),size(Wy,1));

rz = 1./ ( abs(GRADz) .^( 2.0 - q(3) ) + delta ) ;
Rz = spdiags( rz .^0.5,0,size(Wz,1),size(Wz,1));                

aVRWs_t = V * Rs * Ws;
aVRWx_t = Vx * Rx * Wx;
aVRWy_t = Vy * Ry * Wy;
aVRWz_t = Vz * Rz * Wz;

% Form the matrices by including the weights and volume
WstRsWs = ( aVRWs_t )' * ( aVRWs_t ) ;
WxtRxWx = ( aVRWx_t )' * ( aVRWx_t ) ;
WytRyWy = ( aVRWy_t )' * ( aVRWy_t ) ;
WztRzWz = ( aVRWz_t )' * ( aVRWz_t ) ;

    
phis = alpha(1) * ((m - mref)' * WstRsWs * (m - mref) );

phi = (m - mref)' * ( alpha(1) * ((V * Ws)' * (V * Ws) * (m - mref) )  +...
                    alpha(2) * ((Vx * Wx)' * (Vx * Wx) * (m - mref) ) +...
                    alpha(3) * ((Vy * Wy)' * (Vy * Wy) * (m - mref) ) +...
                    alpha(4) * ((Vz * Wz)' * (Vz * Wz) * (m - mref) ) );

% Scale the smallest q-norm on the smallest norm
scale(1) = phi / (phis);

% Scale the gradient p-norm on the smallest q-norm
switch FLAG
case 'SMOOTH_MOD'
    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( m' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * m );

case 'SMOOTH_MOD_DIF'

    % Compute gradients then compute magnitude 
    scale_g = phi /...
    ( (m - mref)' * (alpha(2) * WxtRxWx +...
    alpha(3) * WytRyWy +...
    alpha(4) * WztRzWz) * (m - mref) );


end
    
% end

% Final scales
scale(1) = l(3) * scale(1) * alpha(1);
scale(2) = (2.0 - l(3)) * scale_g * alpha(2);
scale(3) = (2.0 - l(3)) * scale_g * alpha(3);
scale(4) = (2.0 - l(3)) * scale_g * alpha(4);

aVRWs_t = sqrt(scale(1)) * aVRWs_t;
aVRWx_t = sqrt(scale(2)) * aVRWx_t;
aVRWy_t = sqrt(scale(3)) * aVRWy_t;
aVRWz_t = sqrt(scale(4)) * aVRWz_t;

MOF_t =  ( scale(1) * WstRsWs +...
scale(2) * WxtRxWx + scale(3) * WytRyWy +...
scale(4) * WztRzWz ) ;

%% Build the complete model objective function
MOF = blkdiag( MOF_p , MOF_s , MOF_t ); 

aVRWs = blkdiag( aVRWs_p , aVRWs_s , aVRWs_t );
aVRWx = blkdiag( aVRWx_p , aVRWx_s , aVRWx_t );
aVRWy = blkdiag( aVRWy_p , aVRWy_s , aVRWy_t );
aVRWz = blkdiag( aVRWz_p , aVRWz_s , aVRWz_t );
