function [x,r,count]=ConjGrad(x,A,RHS)

r=(A*(x)-RHS);    
p=-r;

rold=r'*r;
count=1;

while norm(r) > 1e-4 && (count < 500)
    
    Ap=A*p;
    alpha=rold./(p'*Ap);
   
    x=x+alpha.*p;
    
%     x(x<0) = 1e-6;
%     x(x>1) = 1;
    
    r=r+alpha*Ap;      
    rnew=r'*r;
    p=-r+rnew/rold.*p;
    rold=rnew;

    count=count+1;
    
%     figure(100)
%     semilogy (count,norm(r),'*');
%     hold on
    
    
    
end

% hold off