function [x0]=conjgrad_v3(x0,A,RHS)


r=(A*(x0)-RHS);    
p=-r;

% invert=zeros(20,60);
rold=r'*r;
rnorm=norm(r);
dnorm=norm(RHS);
count=1;
error=rnorm/dnorm;
misfit(count)=norm(r)^2;


%     figure(5)
%     plot (count,misfit(count),'*');
%     hold on
while count<30%error>=1e-6
    
    Ap=A*p;
    alpha=rold./(p'*Ap);
   
    
     x0=x0+alpha.*p;

    r=r+alpha*Ap;       %+(0*drxz').^2);
    rnew=r'*r;
    p=-r+rnew/rold.*p;
    rold=rnew;


%     invert(:)=x0(:);
    
    rnorm=norm(r);
    error=rnorm/dnorm;
    count=count+1;
    misfit(count)=norm(r)^2;
    
    figure(3)
    plot (count,misfit(count),'*');
    hold on
    
    
    
end

