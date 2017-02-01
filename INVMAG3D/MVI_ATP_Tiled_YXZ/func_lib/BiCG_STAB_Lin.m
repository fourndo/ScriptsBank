function [x,spec] = BiCG_STAB_Lin(x,G, W,MOF,RHS, Proj, PreC, Patv)
% Compute x = A\b using BiConjugateGradient Stabilized method.
% Algorithm from H.A. Van Der Vorst (1992)

% Initial residual
r= Patv*RHS - (Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*x))) + Patv*MOF*Patv*x);    
ro = PreC * r;

rho_in = 1; 
alpha = 1;
omega = 1;
v = 0;
p = 0;

count = 1;
tic
while norm(r) > 1e-4 && count < 100
    
    rho_out = ro'*r;
    
    beta = ( rho_out / rho_in ) * ( alpha / omega );
    
    p = r + beta * ( p - omega * v );
    
    v = Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*p))) + Patv*MOF*Patv*p;
    
    alpha = rho_out / ( ro' * v );
    
    s = r - alpha * v;
    
    t = Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*s))) + Patv*MOF*Patv*s;
    
    h = PreC * t;
    omega = ( h' * s ) / ( h' * t );
    
    x = x + alpha * p + omega * s;
    
    r = s - omega*t;
    
    rho_in = rho_out;
    
%     figure(3)
%     semilogy (count,norm(r),'*');
%     hold on
    
    lrl(count) = norm(r);
    count = count + 1;
    
end

spec = [count, norm(r) , toc];

fprintf('Bi-CGStab steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
% figure;
% plot (lrl,':*');title('Residual curve Bi-CGStab');ylim([0 1e-2])

    