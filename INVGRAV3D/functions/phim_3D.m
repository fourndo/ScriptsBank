function modfunc = phim_3D(x,Wx,Wy,Wz,Ws,Wr,ax,ay,az,as,p,q,lambda,iter)
% function modfunc = phim(mode,x,p,Wx,Wz,lambda,q)
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
% mcell = length(V);
% IV = spdiags(1./sqrt(V'),0,mcell,mcell);

%% Model derivatives:
Wxm = Wx * x;
Wzm = Wz * x;
Wym = Wy * x;
% Small term in compact function
delta = 1e-8;
epsilon = 1e-8;
%% Calculate p-norm:
R = 1./( ( abs(Wxm) ) .^(2-p)  + delta );
Rx = spdiags(R,0,length(R),length(R));

R = 1./( ( abs(Wzm) ) .^(2-p)  + delta  );
Rz = spdiags(R,0,length(R),length(R));

R = 1./( ( abs(Wym) ) .^(2-p) + delta  );
Ry = spdiags(R,0,length(R),length(R));

phiXYZ = ax*(Wx*Ws)'*Rx*(Wx*Ws) + ay*(Wy*Ws)'*Ry*(Wy*Ws) + az*(Wz*Ws)'*Rz*(Wz*Ws);

if iter==1
    scaleXYZ = 1;
else
    
scaleXYZ = (x'*(ax*(Wx*Ws)'*(Wx*Ws) + ay*(Wy*Ws)'*(Wy*Ws) + az*(Wz*Ws)'*(Wz*Ws))*x)/(x'*phiXYZ*x);

end

% if ~(nargin == 8)
%    return;
% end
%% Add second portion ("q-norm"):

R = 1./( ( abs(x) ).^(2-q) + epsilon  );
Rc=spdiags(R,0,length(R),length(R)); 

WctWc = as * Ws' * Rc * Ws;

scaleC = (x'*( scaleXYZ * phiXYZ )*x)/(x' * WctWc * x);
% NOTE: Lambda will should be an input, it controles the relative weight
% between WctWc and (WxtWx+WztWz) after it has been scaled.

% Final model objective function:
modfunc = scaleXYZ * (ax*(Wx*Ws*Wr)'*Rx*(Wx*Ws*Wr) +...
    ay*(Wy*Ws*Wr)'*Ry*(Wy*Ws*Wr) +...
    az*(Wz*Ws*Wr)'*Rz*(Wz*Ws*Wr) ) +...
    scaleC * lambda * as * (Ws * Wr)' * Rc * (Ws * Wr);
% modfunc = scaleXYZ * phiXYZ + scaleC*lambda*WctWc;

end
