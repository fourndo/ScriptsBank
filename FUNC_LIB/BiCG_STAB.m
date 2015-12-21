function [x,spec] = BiCG_STAB(A,b,x)
% Compute x = A\b using BiConjugateGradient Stabilized method.
% Algorithm from H.A. Van Der Vorst (1992)

% Initial residual
r = b - A*x;
ro = r;

rho_in = 1; 
alpha = 1;
omega = 1;
v = 0;
p = 0;

count = 1;
tic
while norm(r) > 1e-4 && count < 200
    
    rho_out = ro'*r;
    
    beta = ( rho_out / rho_in ) * ( alpha / omega );
    
    p = r + beta * ( p - omega * v );
    
    v = A * p;
    
    alpha = rho_out / ( ro' * v );
    
    s = r - alpha * v;
    
    t = A*s;
    
    omega = ( t' * s ) / ( t' * t );
    
    x = x + alpha * p + omega * s;
    
    r = s - omega*t;
    
    rho_in = rho_out;
    
    figure(3)
    semilogy (count,norm(r),'*');
    hold on
    
    lrl(count) = norm(r);
    count = count + 1;
    
end

spec = [count, norm(r) , toc];

fprintf('Bi-CGStab steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
% figure;
% plot (lrl,':*');title('Residual curve Bi-CGStab');ylim([0 1e-2])

    