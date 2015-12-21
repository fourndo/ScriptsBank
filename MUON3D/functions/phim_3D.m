function modfunc = phim_3D(x,p,Ws,Wx,Wy,Wz,as,ax,ay,az,lambda,q,iter,mnull)
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

%% Model derivatives:
Wxm = (Wx) * x;
Wzm = (Wz) * x;
Wym = (Wy) * x;
% Small term in compact function
delta=1e-7;

%% Calculate p-norm:
R = 1./( abs(Wxm).^(2-p) + delta );
Rx = spdiags(R,0,length(R),length(R));

R = 1./( abs(Wzm).^(2-p) + delta );
Rz = spdiags(R,0,length(R),length(R));

R = 1./( abs(Wym).^(2-p) + delta );
Ry = spdiags(R,0,length(R),length(R));

phiXYZ = ax*Wx'*Rx*Wx + ay*Wy'*Ry*Wy + az*Wz'*Rz*Wz;

if iter==1
    scaleXYZ = 1;
else
    
scaleXYZ = (x'*(ax*Wx'*Wx + ay*Wy'*Wy + az*Wz'*Wz)*x)/(x'*phiXYZ*x);
end

% if ~(nargin == 8)
%    return;
% end
%% Add second portion ("q-norm"):
R = 1./( abs(x).^(2-q) + delta );

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
% modfunc = ax*Wx'*Rx*Wx + az*Wz'*Rz*Wz + scaleC*ac*WctWc;
modfunc = abs(2-lambda) * scaleXYZ * phiXYZ + scaleC*as*lambda*WctWc;

end
