function [x,r,count]=PCGLSQ(x,A,RHS,PC, Pa)

A = A*Pa;

r =  RHS - ( A * x );    
p = PC*( A' * r );

s0 = p;

sold= ( A' * r )' * s0;
count=1;

while (count < 200) %&& sold > 1e-6
    
    q = A*p;
    
    alpha = sold / ( q' * q );
   
    x = x + alpha .* p;
    
    r = r - alpha * q;  
    s = A' * r;
    h = PC * s;
    
%     s = (A' * r);
    
    snew = s' * h;
    
    if (snew) / (norm(s0)) < 1e-4
        
        return
        
    end
    
    p = h + ( snew / sold ) * p;
    
    sold=snew;

    count=count+1;
    
%     figure(100)
%     semilogy (count,(snew) / (norm(s0)),'*');
%     hold on
%     
       
end

% hold off