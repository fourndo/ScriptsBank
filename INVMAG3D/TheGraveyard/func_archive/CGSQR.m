function [x,spec] = CGSQR(A,b,x)
% Compute x = A\b using Conjugate Gradient Square method.
% Algorithm from H.A. Van Der Vorst (1992)

% Initial residual
r = b - A*x;
r0 = r;

p = 0;
q = 0;
rhoin = 1;

count = 1;
tic
while norm(r) > 1e-6  && count < 2000 %&& norm(r) < 10
    
    rho = r0' * r;
    
    if rho==0
        
        fprintf('CG-S Stopped before convergence at iteration %i\n',count);
        break
        
    end
    
    if count==1
        u = r;
        p = u;
    else
        beta = rho/rhoin;

        u = r + beta*q;

        p = u + beta*(q + beta*p);
    end
    
    v = A*p;
    
    alpha = rho / (r0'*v);    
    
    q = u - alpha*v;
    
    w = u + q;
        
    rout = r - alpha*(A*w);
    
    % Re-started if iteration is going rogue
    if norm(rout) > norm(r)
        
        r0 = r;
        rho = r0' * r;
        u = r;
        p = u;
        v = A*p;
        alpha = rho / (r0'*v);    
    
        q = u - alpha*v;

        w = u + q;

        r = r - alpha*(A*w);
        fprintf('Restarted CGS\n')
        lrl(count) = norm(r);
        count = count +1;
    else
        
        r = rout;
        
    end
    
    x = x + alpha * w;
    rhoin = rho;
    
    
    
    lrl(count) = norm(r);
    
    count = count +1;
    
    
end

spec = [count, norm(r) , toc];
fprintf('CG_SQR steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
figure;
plot (lrl,':*');title('Residual curve CGSQR');ylim([0 1e-2])
hold on