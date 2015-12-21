function [Wx,Vx]=getWx1D(nx,dx)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
Wx = spalloc(nx-1,nx,nx*2);
vx = zeros(nx-1,1);
count=1;
for ii=1:nx-1
 
    Wx(count,count+1) = -sqrt(1/dx(ii)); 
    Wx(count,count) = sqrt(1/dx(ii));
   
    vx(count) = ( ( dx(ii) + dx(ii+1) ) /2 )^-0.5;
    count=count+1;
    
end

Vx = spdiags(vx,0,nx-1,nx-1);