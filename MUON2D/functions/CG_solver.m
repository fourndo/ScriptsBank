function [x0,count]=CG_solver(x0,A,RHS)


r=(A*(x0)-RHS);    
p=-r;

rold=r'*r;
rnorm=sqrt(sum(r.^2));
dnorm=sqrt(sum(RHS.^2));
count=1;
error=rnorm/dnorm;

while (error(count) >= 10e-4) && (count < 10000) 
    
    Ap=A*p;
    alpha=rold/(p'*Ap);
   
    
     x0=x0+alpha*p;

    r=r+alpha*Ap;       %+(0*drxz').^2);
    rnew=r'*r;
    p=-r+(rnew/rold)*p;
    rold=rnew;

    
    rnorm=sqrt(sum(r.^2));
    
    count=count+1;

    error(count)=rnorm/dnorm;

   

  
end
sprintf('CG iterations:%d',count)
