function [x,r,count]=CGLSQ(x,A,RHS)

r=(A*(x)-RHS);    
p=-(A'*r);

s0=p;

sold=s0'*s0;
count=1;

while (count < 400)
    
    q = A*p;
    
    alpha=sold./(q'*q);
   
    x = x + alpha.*p;
%     
%     x(x<0) = 1e-6;
%     x(x>1) = 1;
    
    r = r + alpha * q;  
    
    s = A' * r;
    
    snew = s' * s;
    
    if (snew) / (norm(s0)) < 1e-4
        
        return
        
    end
    
    p = -s + ( snew / sold ) * p;
    
    sold=snew;

    count=count+1;
    
%     figure(100)
%     semilogy (count,(snew) / (norm(s0)),'*');
%     hold on
    
    
    
end

% hold off