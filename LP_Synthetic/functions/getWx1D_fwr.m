function [Wx]=getWx1D_fwr(mcell,nx,dx)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
Wx = sparse(mcell-1,mcell);
count=1;
for ii=1:nx-1
  
    Wx(count,count+1) = sqrt(1/dx(ii+1)); 
    Wx(count,count) = -sqrt(1/dx(ii));

    count=count+1;
    
end

