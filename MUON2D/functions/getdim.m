function [ws,wx,wz]=getdim(mcell,dX,dZ)
% Compute dimension matrices (experimental)

nX=length(dX);
nZ=length(dZ);

ws=zeros(mcell,1);
wx=zeros(mcell,1);
wz=zeros(mcell,1);

count=0;
for ii=1:nX
    
    for jj=1:nZ
        
        count=count+1;
        
        ws(count)= dX(ii) * dZ(jj);
        wx(count)= dZ(jj) / dX(ii);
        wz(count)= dX(ii) / dZ(jj);
        
    end
    
end

ws = spdiags(ws,0,mcell,mcell);
wx = spdiags(wx,0,mcell,mcell);
wz = spdiags(wz,0,mcell,mcell);