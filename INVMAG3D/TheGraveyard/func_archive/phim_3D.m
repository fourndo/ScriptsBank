function modfunc = phim_3D(x,Wx,Wy,Wz,Ws,Wr,Vx,Vy,Vz,ax,ay,az,as,p,q,lambda,iter)
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

phiXYZ = ax*(Vx * Wx)'*Rx*(Vx * Wx) + ay*(Vy * Wy)'*Ry*(Vy * Wy) + az*(Vz * Wz)'*Rz*(Vz * Wz);

if iter==1
    scaleXYZ = 1;
else
    
scaleXYZ = (x'*(ax*(Vx * Wx)'*(Vx * Wx) + ay*(Vy * Wy)'*(Vy * Wy) + az*(Vz * Wz)'*(Vz * Wz))*x)/(x'*phiXYZ*x);

end

% if ~(nargin == 8)
%    return;
% end
%% Add second portion ("q-norm"):

R = 1./( ( abs(x) ).^(2-q) + epsilon  );
Rc=spdiags(R,0,length(R),length(R)); 

WctWc = as * Ws' * Rc * Ws;

if iter==1
    scaleC = 1;
else
    scaleC = (x'*( scaleXYZ * phiXYZ )*x)/(x' * WctWc * x);
end
% NOTE: Lambda will should be an input, it controles the relative weight
% between WctWc and (WxtWx+WztWz) after it has been scaled.

% Final model objective function:
modfunc = (2-lambda) * ( scaleXYZ * (ax*(Vx * Wx * Wr)'*Rx*(Vx * Wx * Wr) +...
    ay*(Vy * Wy * Wr)' * Ry * (Vy * Wy * Wr) +...
    az*(Vz * Wz * Wr)'*Rz*(Vz * Wz * Wr) ) ) +...
    scaleC * lambda * (Ws * Wr)' * Rc * (Ws * Wr);
% modfunc = scaleXYZ * phiXYZ + scaleC*lambda*WctWc;

end
