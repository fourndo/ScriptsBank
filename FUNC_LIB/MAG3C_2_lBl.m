% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
% close all
clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath FUNC_LIB;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Amplitude\I45D45\MAG3D_FULL';
obsfile   = 'magfor3d.mag'; 
outfile   = 'EQS_lBl.obs';

% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(obsx);

if mod(ndata,1)~=0
    
    fprintf('Data does not appear to be 3C. Please revise...\n');
    break
    
end

datax = d(1:ndata) ; wdx = wd(1:ndata);
datay = d( (ndata+1) : (2*ndata)) ; wdy = wd( (ndata+1) : (2*ndata) );
dataz = d( (2*ndata+1) : (3*ndata)) ; wdz = wd( (2*ndata+1) : (3*ndata) );

ampd = sqrt( datax.^2 + datay.^2 + dataz.^2 );
wd = sqrt( wdx.^2 + wdy.^2 + wdz.^2 );

write_MAG3D_TMI([work_dir '\' outfile],H,I,Dazm,obsx,obsy,obsz,ampd,wd)
