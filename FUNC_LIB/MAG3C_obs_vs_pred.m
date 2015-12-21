% Load in MAG3D_3C observation files and write out |B|
% Dominique Fournier 2014/11/02
% close all
clear all
close all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath FUNC_LIB;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Equi_source\Checker_board\ES_56x64_ax0';
obsfile   = 'Obs_checkers_56x64_3C.obs';
prefile   = 'magfor3d.mag';
% outfile   = 'EQS_lBl.obs';

% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d_obs, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
plot_mag3C(obsx,obsy,d_obs,I,D,'Observed 3C-data')

% Load predicted file (3C UBC-MAG format)
[~, ~, ~, ~, ~, ~, ~, d_pre, ~] = read_MAG3D_obs([work_dir '\' prefile]);
plot_mag3C(obsx,obsy,d_pre,I,D,'Predicted 3C-data')

% Plot residual
plot_mag3C(obsx,obsy,(d_obs - d_pre)./wd,I,D,'Residual (Normalized)')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

% ndata = length(obsx);
% 
% if mod(ndata,1)~=0
%     
%     fprintf('Data does not appear to be 3C. Please revise...\n');
%     break
%     
% end
% 
% datax = d_obs(1:ndata) ; wdx = wd(1:ndata);
% datay = d_obs( (ndata+1) : (2*ndata)) ; wdy = wd( (ndata+1) : (2*ndata) );
% dataz = d_obs( (2*ndata+1) : (3*ndata)) ; wdz = wd( (2*ndata+1) : (3*ndata) );
% 
% ampd = sqrt( datax.^2 + datay.^2 + dataz.^2 );
% wd = sqrt( wdx.^2 + wdy.^2 + wdz.^2 );
% 
% write_MAG3D_TMI([work_dir '\' outfile],H,I,Dazm,obsx,obsy,obsz,ampd,wd)
