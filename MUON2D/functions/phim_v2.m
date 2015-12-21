function modfunc = phim_v2(x,p,Wx,Wz,ax,az,ac,q)
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

% Small term in compact function
delta=1e-11;

%% Calculate p-norm:
R = 1./( abs(Wxm).^(2-p) + delta );
Rx = diag(R);
R = 1./( abs(Wzm).^(2-p) + delta );
Rz = diag(R);

phiXYZ = ax*Wx'*Rx*Wx + az*Wz'*Rz*Wz;

% Scale the minimum gradient support on the l2-gradient to keep the
% norm of the objective function roughly the same
scaleXYZ = (x'*(ax*Wx'*Wx + az*Wz'*Wz)*x)/(x'*phiXYZ*x);


%% Add second portion ("q-norm"):

Rc = 1./( abs(x).^(2-q) + delta );
WctWc=diag(Rc); 

scaleC = (x'*(scaleXYZ * phiXYZ)*x)/(x'*WctWc*x);

% Final model objective function:
% modfunc = ax*Wx'*Rx*Wx + az*Wz'*Rz*Wz + scaleC*ac*WctWc;
modfunc = scaleXYZ * phiXYZ + scaleC*ac*WctWc;
return
end
