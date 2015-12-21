function [Wx]=getWx1D_cntr(mcell,nX,dX)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
Wx = zeros(mcell,mcell);
count=1;
for ii=1:nX
    if ii == nX   
    Wx(count,count-1) = -sqrt(1/dX(ii-1)/2); 
%     Wx(count,count) = -sqrt(1/dX(ii));
    elseif ii == 1
    Wx(count,count+1) = sqrt(1/dX(ii+1)/2);
    else
    Wx(count,count+1) = -sqrt(1/dX(ii-1)/2);
    Wx(count,count-1) = sqrt(1/dX(ii+1)/2);
    end
    count=count+1;
    
end
