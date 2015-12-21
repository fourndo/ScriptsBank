function aout=unfold2(ain)
%
% function aout=unfold2(ain)
%
% unfold a 2D Fourier transform to shift
% the DC component to the "center"
%
[nx1,ny1]=size(ain);
kx=nx1/2;
ky=ny1/2;
aout(1:kx-1,1:ky-1)=ain(kx+2:nx1,ky+2:ny1);
aout(kx:nx1,1:ky-1)=ain(1:kx+1,ky+2:ny1);
aout(1:kx-1,ky:ny1)=ain(kx+2:nx1,1:ky+1);
aout(kx:nx1,ky:ny1)=ain(1:kx+1,1:ky+1);
%
%end function unfold2
end