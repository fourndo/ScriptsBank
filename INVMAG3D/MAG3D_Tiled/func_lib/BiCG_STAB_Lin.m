function [x,r,count] = BiCG_STAB_Lin(x,G, W,MOF,RHS, Proj, PreC, Patv)
% Compute x = A\b using BiConjugateGradient Stabilized method.
% Algorithm from H.A. Van Der Vorst (1992)

% Initial residual
r= Patv*RHS - (Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*x))) + Patv*MOF*Patv*x);    
ro = r;

rho0 = 1; 
a = 1;
w0 = 1;
v_0 = 0;
p_0 = 0;

count = 1;
tic
while norm(r) > 1e-4 && count < 20
    
    rho_i = ro'*r;
    
    beta = ( rho_i / rho0 ) * ( a / w0 );
    
    p_i = r + beta * ( p_0 - w0 * v_0 );
    
    y = PreC * p_i;
    
    v_0 = Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*y))) + Patv*MOF*Patv*y;
    
    a = rho_i / ( ro' * v_0 );
    
    h = x + a * y;
    
%     if norm(a*y) < 1e-4
%         
%         return
%         
%     end
    
    s = r - a * v_0;
    
    z = PreC * s;
    
    t = Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*z))) + Patv*MOF*Patv*z;
    
    w0 = ( (PreC*t)' * (PreC*s) ) / ( (PreC*t)' * (PreC*t) );
    
    x = h + w0 * z;
    
    r = s - w0*t;
    
    rho0 = rho_i;
    
%     figure(3)
%     semilogy (count,norm(r),'*');
%     hold on
    
    lrl(count) = norm(r);
    count = count + 1;
    
end

% spec = [count, norm(r) , toc];

% fprintf('Bi-CGStab steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
% figure;
% plot (lrl,':*');title('Residual curve Bi-CGStab');ylim([0 1e-2])

    