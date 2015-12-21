function aout=fold2(ain)

% function aout=fold2(ain)
%
% fold a 2D array to shift DC component
% to the lower-left corner
%
[nx1,ny1]=size(ain);
kx=nx1/2;
ky=ny1/2;
aout(1:kx+1,1:ky+1)=ain(kx:nx1,ky:ny1);
aout(kx+2:nx1,1:ky+1)=ain(1:kx-1,ky:ny1);
aout(1:kx+1,ky+2:ny1)=ain(kx:nx1,1:ky-1);
aout(kx+2:nx1,ky+2:ny1)=ain(1:kx-1,1:ky-1);
end