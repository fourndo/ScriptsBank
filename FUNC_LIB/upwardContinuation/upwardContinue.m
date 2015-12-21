% function newdat = upwardContinue(delz,oldx,oldy,delx,dely,olddat)
% oldx = 0:10:200, then delx = 10;
% Author: Kris Davis
% Date: February 1, 2005
% Modified: 10, April, 2008 (K. Davis)
%
%
% This program sets up the framework for inputting gridded data and apply
% the Fourier transform. It also calculates the power spectrum in 2D.  This
% program calls fold2D.m which is another program designed to fold a
% function in two dimensions. 

% Reading in the data set using the function read2D()
% clear all;
% close all;
function newdat = upwardContinue(delz,oldx,oldy,delx,dely,olddat)

% [nx,ny,delx,dely,oldx,oldy,olddat,inpf] = readgxf;

% Continued height
delz = abs(delz);
[nx,ny] = size(olddat);
[dat,lx1,lx2,ly1,ly2] = expand(olddat,1,oldx,oldy,delx,dely);
nx = nx + lx1 + lx2;
ny = ny + ly1 + ly2;

%Nyquist 
niqx = nx/2;
niqy = ny/2;

% Calculating frequencies for FFT
dox    = 2*pi/(delx*nx);
doy    = 2*pi/(dely*ny);
omegax = dox*(-niqx+1:niqx);
omegay = doy*(-niqy+1:niqy);

%Meshing for the power spectrum
[Omegax, Omegay] = meshgrid(omegax, omegay);

%Radial wavenumber
Omegar = sqrt(Omegax.^2 + Omegay.^2);
Omegar(niqx,niqy) = 1e-10;

% FFT data
fftDat = unfold2(fft2(fold2(dat)));

% Upward continue
ucfftDat = fftDat.*exp(-Omegar*delz);

%Take Inverse FFT
ucDat = unfold2(real(ifft2(fold2(ucfftDat))));

% New data is not the expanded version
newdat = ucDat((lx1+1:(nx-lx2)),(ly1+1:(ny-ly2)));

end
