% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

% addpath ..\FUNC_LIB\;
addpath /tera23/dfournier/Code/MAG3CampB/functions/

% Project folders
work_dir = '/tera_raid/dfournier/TKC/Mag/MAG3C_Inv1';
inpfile = 'MAG3Csen.inp';

[meshfile,obsfile,magfile,topofile, wr_flag] = MAG3Csen_read_inp([work_dir '/' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '/' meshfile]);

% Write logfile
fid = fopen([work_dir '/MAG3Csen.log'],'w');
fprintf(fid,'MAG3Csen\n');
fprintf(fid,'Generates sparse matrices for magnetostatic forward modeling: Tx, Ty, Tz\n');
fprintf(fid,'Topographic model: nullcell.dat\n');
fprintf(fid,'DISTANCE | DEPTH weighting: wr.dat\n\n');
fprintf(fid,'Written by: Dominique Fournier\n');
fprintf(fid,'Last update: July 14th, 2014\n\n');
fprintf(fid,'INPUT FILES: \n');
fprintf(fid,'Mesh: \t\t\t %s \n',meshfile);
fprintf(fid,'Obsfile: \t\t\t %s \n',obsfile);
fprintf(fid,'Topography: \t\t\t %s \n',magfile);
fprintf(fid,'Magnetization model: \t\t\t %s \n',topofile);
fprintf(fid,'Weighting: \t\t\t %s \n',wr_flag);
fclose(fid);
%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

if isempty(topofile)==1
    
    nxn = length(xn);
    nyn = length(yn);
    topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    
else
    % Load topo
    topo = read_UBC_topo([work_dir '/' topofile]);
end

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
save([work_dir '/nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
cellID = find(nullcell==1);

% Get nodal discretization for octree levels
celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);

% Load susceptibility model
mcell = length(nullcell);

% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '/' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')
ndata = length(obsx);

% Create depth weighting
wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
save([work_dir '/wr.dat'],'-ascii','wr');

% Load magnetization model
if isempty(magfile) == 1
    
    mag_azmdip = [ones(mcell,1)*Dazm ones(mcell,1)*I];
    
else
    
    mag_azmdip = load([work_dir '\' magfile]);
end

%% Create model magnetization vectors
% Azimuth and dip of magnitization
mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
M = [spdiags(H * mag_xyz(:,1),0,mcell,mcell);spdiags(H * mag_xyz(:,2),0,mcell,mcell);spdiags(H * mag_xyz(:,3),0,mcell,mcell)];

% Compute sensitivities
% Pre-allocate
Tx = zeros(ndata,3*mcell);
Ty = zeros(ndata,3*mcell);
Tz = zeros(ndata,3*mcell);
progress = -1;
tic 
fid = fopen([work_dir '/MAG3Csen.log'],'a');
for ii = 1:ndata
   
    % compute kernel for active cells
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),mcell,celln,cellID);

    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf(fid,'Computed %i pct of data in %8.5f sec\n',d_iter*5,toc);
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end

fprintf(fid,'Sensitivity calculation completed in: %f min\n',toc/60);
fclose(fid);

save([work_dir '/Tx'],'Tx','-v7.3');
save([work_dir '/Ty'],'Ty','-v7.3');
save([work_dir '/Tz'],'Tz','-v7.3');


