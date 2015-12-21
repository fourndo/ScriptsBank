function modfunc = phim(x,p,Wx,ax,lambda,Ws,q)
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

%% Model derivatives:
Wxm = (Wx) * x;
% Wzm = (Wz) * x;

% Small term in compact function
delta=1e-11;

%% Calculate p-norm:
R = 1./( abs(Wxm).^(2-p) + delta );
Rx = diag(R);
% R = 1./( abs(Wzm).^(2-p) + delta );
% Rz = diag(R);

modfunc = ax*Wx'*Rx*Wx;% + az*Wz'*Rz*Wz;

if (nargin < 7)
   return;
end
%% Add second portion ("q-norm"):
Rc = 1./( abs(x).^(2-q) + delta );
WctWc = Ws'*diag(Rc)*Ws;
% Calculate lambda to trade-off from p-norm:
scale = (x'*modfunc*x)/(x'*WctWc*x);

% Final model objective function:
modfunc = modfunc + scale*lambda*WctWc;

return
end
