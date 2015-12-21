function [x0,spec]=CGNE(A,RHS,x0)


r=((A*(x0))-RHS);   
z = A'*r;
p=z;

tic;
count = 1;
while norm(r) > 1e-5 %&& count < 10000
    
    
    Ap=(A*p);
    alpha=norm(z)^2/norm(Ap)^2;
   
    
     x0=x0-alpha*(p);
     
    r=r-alpha*Ap;
    znew = A'*r;
%     rnew=r'*r;
    beta = norm(znew)^2/norm(z)^2;
    p = znew + beta * p;
    
    z = znew;


%     invert(:)=x0(:);
    
    
    lrl(count) = norm(r);
    count = count +1;
    
end

spec = [count, norm(r) , toc];

fprintf('CG-NE steps: %i, Residual: %6.3e, in %8.5f sec\n', count, norm(r) , toc);
figure;
plot (lrl,':*');title('Residual curve CGNE');ylim([0 1e-2])
hold on
    
    