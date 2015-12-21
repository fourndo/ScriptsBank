function [x0]=ConjGrad(x0,A,RHS,bounds)


r=(A*(x0)-RHS);    
p=-r;

% invert=zeros(20,60);
rold=r'*r;
rnorm=norm(r);
dnorm=norm(RHS);
count=1;
error=rnorm/dnorm;
misfit(count)=norm(r)^2;

% Treshhold vectors


%     figure(5)
%     plot (count,misfit(count),'*');
%     hold on
while rnorm>=1e-4
    
    Ap=A*p;
    alpha=rold./(p'*Ap);
   
    
     x0=x0+alpha.*p;
     
     x0(x0<0) = bounds(1);
     p(x0<0) = 0;
     x0(x0>1) = bounds(2);
     p(x0>1) = 0;
     
    r=r+alpha*Ap;       %+(0*drxz').^2);
    rnew=r'*r;
    p=-r+rnew/rold.*p;
    rold=rnew;


%     invert(:)=x0(:);
    
    rnorm=norm(r);
    error=rnorm/dnorm;
    count=count+1;
    misfit(count)=norm(r)^2;
    
    figure(100)
    semilogy(count,misfit(count),'*');
    hold on
    
    
    
end

