function [x,spec] = Bi_CG(A,b,x)
% Compute x = A\b using BiConjugateGradient method.
% Algorithm from H.A. Van Der Vorst (1992)

% Initial residual
r = b - A*x;
r2 = b - A'*x;

p = r;
p2 = r2;

count = 1;
tic
while norm(r) > 1e-6 && count < 2000
    
    alpha = ( r2'*r ) / ( p2' * A * p );
    
    x = x + alpha * p;
      
    rout = r - alpha*A*p;
    r2out = r2 - alpha*A'*p2;
    
    beta = ( r2out' * rout ) / ( r2' * r );
    
    p = rout + beta * p;
    p2 = r2out + beta * p2;
    
    
    r = rout;
    r2 = r2out;
    
    lrl(count) = norm(r);
    count = count + 1;
    
end

spec = [count, norm(r) , toc];

fprintf('Bi-CG steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
figure;
plot (lrl,':*');title('Residual curve Bi-CG');ylim([0 1e-2])
hold on