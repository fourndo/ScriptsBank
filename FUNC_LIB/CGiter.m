function [x,normr,count]=CGiter(x,A,RHS)

r=(A*(x)-RHS);    
p=-r;

rold=r'*r;
count=0;

normr = norm(r);
while count < length(x) && normr > 1e-3
    
    count=count+1;
    Ap=A*p;
    alpha=rold./(p'*Ap);
   
    x=x+alpha.*p;
     
    r=r+alpha*Ap;      
    rnew=r'*r;
    p=-r+rnew/rold.*p;
    rold=rnew;

   
    
    normr = norm(r);
    
%     figure(3)
%     semilogy (count,normr,'*');
%     hold on
    
    
    
end
