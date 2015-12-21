function [x0,count]=conjgrad_v2(x0,A,RHS)


r=(A*(x0)-RHS);    
p=-r;

% invert=zeros(20,60);
rold=r'*r;
rnorm=sqrt(sum(r.^2));
dnorm=sqrt(sum(RHS.^2));
count=1;
error=rnorm/dnorm;
% misfit(count)=norm(r)^2;


%     figure(5)
%     plot (count,misfit(count),'*');
%     hold on
while (error(count) >= 10e-4) && (count < 100) 
    
    Ap=A*p;
    alpha=rold/(p'*Ap);
   
    
     x0=x0+alpha*p;

    r=r+alpha*Ap;       %+(0*drxz').^2);
    rnew=r'*r;
    p=-r+(rnew/rold)*p;
    rold=rnew;

   
%     invert(:)=x0(:);
    
    rnorm=sqrt(sum(r.^2));
    
    count=count+1;
%     misfit(count)=norm(r)^2;
    error(count)=rnorm/dnorm;

   

  
end
sprintf('CG iterations:%d',count)
%     figure(5)
%     plot (misfit,'-+');
%     hold on