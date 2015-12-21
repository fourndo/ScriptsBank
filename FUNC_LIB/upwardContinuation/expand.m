function [aout,lx1,lx2,ly1,ly2,x,y]=expand(ain,imean,oldx,oldy,dx,dy)
%
%  function [aout,lx1,lx2,ly1,ly2]=expand(ain,imean,oldx,oldy,dx,dy)
%
%  padd a 2D map to the nearest power of 2
%
%  imean: flag indicating the value to padd the grid to:
%  imean=0: padd to zero
%            =1: padd to mean value of the edge
%

[ny,nx]=size(ain);
dc=(sum(ain(1,:))+sum(ain(ny,:))+sum(ain(2:ny-1,1))+...
       sum(ain(2:ny-1,nx)))/(2*ny+2*nx-4);
   
if imean==1,
   val=dc;
else
    val=0.0;
end

%
% The nearest power of 2 
%
n2x = pow2(ceil(log2(nx)));
n2y = pow2(ceil(log2(ny)));
if (n2x-nx)/nx<(0.1), n2x=n2x*2; end;
if (n2y-ny)/ny<(0.1), n2y=n2y*2; end;


%
% Padding with linear tapering
%
lx1=fix((n2x-nx)/2);
lx2=n2x-nx-lx1;
ly1=fix((n2y-ny)/2);
ly2=n2y-ny-ly1;

aout=zeros(n2y,n2x);
aout(ly1+1:ly1+ny,lx1+1:lx1+nx)=ain;

for ii=lx1+1:lx1+nx
   for jj=-ly1:-1
      wt=(jj+ly1)/ly1;
      wk=aout(ly1+1,ii)+(aout(ly1+1,ii)-aout(ly1+1-jj,ii))*wt;
      aout(ly1+1+jj,ii)=wk+(val-wk)*(1-wt);
   end;
   for jj=1:ly2,
      wt=(ly2-jj)/ly2;
      wk=aout(n2y-ly2,ii)+(aout(n2y-ly2,ii)-aout(ly1+ny-jj,ii))*wt;
      aout(ly1+ny+jj,ii)=wk+(val-wk)*(1-wt);
   end;
end;

for jj=1:n2y,
   for ii=-lx1:-1,
      wt=(ii+lx1)/lx1;
      wk=aout(jj,lx1+1)+(aout(jj,lx1+1)-aout(jj,lx1+1-ii))*wt;
      aout(jj,lx1+1+ii)=wk+(val-wk)*(1-wt);
   end;
   for ii=1:lx2,
      wt=(lx2-ii)/lx2;
      wk=aout(jj,n2x-lx2)+(aout(jj,n2x-lx2)-aout(jj,lx1+nx-ii))*wt;
      aout(jj,lx1+nx+ii)=wk+(val-wk)*(1-wt);
   end;
end;

x = zeros(1,lx1+lx2+nx);
y = zeros(1,ly1+ly2+ny);
kk=1;
for ii = lx1:-1:1
    x(ii)=oldx(1)-(dx*kk);
    kk=kk+1;
end
kk = 1;
for ii = (nx+lx1):length(x)
    x(ii)=oldx(nx)+(dx*kk);
    kk=kk+1;
end

kk=1;
for ii = (lx1+1):(lx1+nx)
    x(ii)=oldx(kk);
    kk=kk+1;
end

kk = 1;
for ii = ly1:-1:1
    y(ii)=oldy(1)-(dy*kk);
    kk=kk+1;
end
kk = 1;
for ii = (ny+ly1):length(y)
    y(ii)=oldy(ny)+(dy*kk);
    kk=kk+1;
end

kk = 1;
for ii = (ly1+1):(ly1+ny)
    y(ii)=oldy(kk);
    kk=kk+1;
end

%
% end of function expand
%