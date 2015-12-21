close all;
clear all;

%
% Author: Kris Davis
% Date: February 22, 2005
%
%   This program computes the gravity gradient tensors of a sphere in the
%   Fourier domain. The dimensions and properties of the sphere can be
%   determined by the user or the user may use the default parmeters.  The
%   total-field magnetic anomaly is also computed as well as the real and
%   imaginary parts of the FFT of the magnetic data.  The field is computed
%   in the Fourier domain as with the gravity and is inversed transformed
%   to the space domain
%

% default = input('Would you like to use the default parameters (Y or N)?  ', 's');
%% Set up parameters
delz = 50; % Upward continuation
rho = 1; %kg/m^3
r = 200; %meters, radius
%Calculate Mass of Object
m = (4/3)*pi*(r^3)*rho;
%Source Position
a = 0; %x
b = 0; %y
%Z is positive down
d = 250; %depth to source
z = 0; %observation height
xmin = -775; %meters
xmax = 800;
ymin = -775;
ymax = 800;
inc = 25; %increment

g = 6.63e-11; %g constant

%Set up the parameters
x = xmin:inc:xmax;
y = ymin:inc:ymax;
nx = (xmax-xmin)/inc + 1;
ny = (ymax-ymin)/inc + 1;

%Nyquist
niqx = nx/2;
niqy = ny/2;

%Distance between points in x and y directions
dx = x(2) - x(1);
dy = y(2) - y(1);
n = 1;

% Calculating frequencies for FFT
for ii = -niqx + 1:niqx
    omegax(ii + niqx) = ((2*pi)*(ii))/(nx*dx);
end
for jj = -niqy + 1:niqy
    omegay(jj + niqy) = ((2*pi)*(jj))/(ny*dy);
end

dwx = omegax(2) - omegax(1);
dwy = omegay(2) - omegay(1);

%Creating Mesh of frequencies for later use
[OMEGAY,OMEGAX] = meshgrid(omegay,omegax);

% Calcuating h to make it easier for Tensor calculation
h = -(d-z);
h1 = -(d-(z-delz));

% Calculating Wr
Wr=sqrt(OMEGAX.^2+OMEGAY.^2);
Wr(niqx,niqy)=1.0E-10;

%% Vertical gravity anomaly at both heights
Gz = g*m*10^12*exp(h*Wr);
GzUC = g*m*10^12*exp(h1*Wr);

% Folding data
newGz = fold2(Gz);
newGzUC = fold2(GzUC);

%Perform inverse FFT with a scaling factor
inGz = real(ifft2(newGz))*((nx*ny*dwx*dwy)/(2*pi));
inGzUC = real(ifft2(newGzUC))*((nx*ny*dwx*dwy)/(2*pi));

%Unfold the inverse FFT data
Gz = unfold2(inGz);
GzUC = unfold2(inGzUC);

%% Test upward continuation operator for gravity
newdat = upwardContinue(delz,x,y,dx,dy,Gz);
UCdiff = abs(GzUC - newdat);

figure('name','vertical gravity');
subplot(221)
imagesc(x,y,Gz);
title('original');
colorbar;
axis equal;
axis tight;

subplot(223)
imagesc(x,y,GzUC);
title('Fourier modelled');
colorbar;
axis equal;
axis tight;

subplot(222)
imagesc(x,y,newdat);
title('UC');
colorbar;
axis equal;
axis tight;

subplot(224)
imagesc(x,y,UCdiff);
title('UC diff');
colorbar;
axis equal;
axis tight;

%% Magnetics 
incl = 65.0;
decl = 25;

I = -(incl*pi)/180;
D = -(decl*pi)/180;

% Bo components
ux = cos(I)*cos(D);
uy = cos(I)*sin(D);
uz = sin(I);

% Magnetization vector components 
m = (4/3)*pi*(r^3)*2.0;
mx = m*ux;
my = m*uy;
mz = m*uz;

%Vectors km, kb and the total mag field
km = sqrt(-1)*OMEGAX*mx + sqrt(-1)*OMEGAY*my + mz*Wr;
kb = sqrt(-1)*OMEGAX*ux + sqrt(-1)*OMEGAY*uy + uz*Wr;
magT = (1e-7)*sqrt((mx^2)+(my^2)+(mz^2))*km.*kb.*(exp(h*Wr))./Wr;
magT1 = (1e-7)*sqrt((mx^2)+(my^2)+(mz^2))*km.*kb.*(exp(h1*Wr))./Wr;

mag = unfold2(real(ifft2(fold2(magT))*((nx*ny*dwx*dwy)/(2*pi))));
magUC = unfold2(real(ifft2(fold2(magT1))*((nx*ny*dwx*dwy)/(2*pi))));

%% Magnetic upward continuation
newdat = upwardContinue(delz,x,y,dx,dy,mag);
UCdiff = abs(magUC - newdat);


figure('name','magnetic dipole');
subplot(221)
imagesc(y,x,mag);
title('original');
colorbar;
axis equal;
axis tight;

subplot(223)
imagesc(y,x,magUC);
title('Fourier modelled');
colorbar;
axis equal;
axis tight;

subplot(222)
imagesc(y,x,newdat);
title('UC');
colorbar;
axis equal;
axis tight;

subplot(224)
imagesc(y,x,UCdiff);
title('UC diff');
colorbar;
axis equal;
axis tight;