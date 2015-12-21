function modfunc = phim(x,Wx,Wz,Ws,ax,az,ac,p,q,l,var)
% function modfunc = phim(x,p,Wx,Wz,ax,az,ac,Ws,q)
% ---
% Generalized model objective function for the 2D problem given by:
%                Phim = ||x||^p + lambda*||m||^q
%
% Inputs:
%             model: x
%         exponents: p,q
% derivative matrix: Wx, Wz
%            alphas: ax,az
% ---
mcell = length(x);
%% Model derivatives:
Wxm = (Wx) * x;
Wzm = (Wz) * x;

% Small term in compact function
delta=1e-8;

%% Calculate p-norm:
R = 1./( abs(Wxm).^(2-p) + delta );
Rx = spdiags(R,0,size(Wx,1),size(Wx,1));
R = 1./( abs(Wzm).^(2-p) + delta );
Rz = spdiags(R,0,size(Wz,1),size(Wz,1));

phiXZ = ax*(Wx*Ws)'*Rx*(Wx*Ws) + az*(Wz*Ws)'*Rz*(Wz*Ws);

scaleXZ = (x' * ( ax*(Wx*Ws)'*(Wx*Ws) + az*(Wz*Ws)'*(Wz*Ws) ) * x )/ ( x' * phiXZ * x );

%% Add second portion ("q-norm"):
Rc = 1./( abs(x).^(2-q) + delta );
WctWc = Ws'*spdiags(Rc,0,mcell,mcell)*Ws;

% Calculate lambda to trade-off from p-norm:
scaleC = (x'*(scaleXZ * phiXZ)*x)/(x'* (ac * WctWc) * x);

% Final model objective function:
modfunc = var * scaleXZ * phiXZ + l * scaleC * ac * WctWc;


end
