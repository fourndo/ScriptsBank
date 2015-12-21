function [Wx,Vx]=getWx1D(nx,dx)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
Wx = spalloc(nx,nx,nx*2);
vx = zeros(nx,1);
count=1;
for ii=1:nx
 
    if ii < nx
        Wx(count,count+1) = -1; 
        Wx(count,count) = 1;

        vx(count) = ( ( dx(ii) + dx(ii+1) ) /2 )^-0.5;
        count=count+1;
        
    else
        
        Wx(count,count-1) = 1; 
        Wx(count,count) = -1;

        vx(count) = ( ( dx(ii) + dx(ii-1) ) /2 )^-0.5;
        count=count+1;
        
    end
        
    
end

Vx = spdiags(vx,0,nx,nx);