function [x,r,count]=CGLSQ(x,A,RHS)

r = ( RHS - A * (x) );    
p = ( A' * r );

s0 = p;

sold= (A'*r)' * s0;
count=0;

while  count < length(x)
    
    count=count+1;
    
    q = A*p;
    
    alpha = sold / ( q' * q );
   
    x = x + alpha .* p;
    
    r = r - alpha * q;  
    s = A' * r;
    h =  s;
    
%     s = (A' * r);
    
    snew = s' * h;
    
    if (snew) / (norm(s0)) < 1e-4
        
        return
        
    end
    
    p = h + ( snew / sold ) * p;
    
    sold=snew;

    
% % Plot change in residual if needed    
%     figure(100)
%     semilogy (count,(snew) / (norm(s0)),'*');
%     hold on
%     
       
end

% hold off