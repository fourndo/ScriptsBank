% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\Tile_AMI\Tile1';
inpfile = 'MAG3Csen.inp';

[meshfile,obsfile,topofile, wr_flag, sen_flag] = MAG3Csen_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

% Write logfile
fid = fopen([work_dir '\MAG3Csen.log'],'w');
fprintf(fid,'MAG3Csen\n');
fprintf(fid,'Generates sparse matrices for magnetostatic forward modeling: Tx, Ty, Tz\n');
fprintf(fid,'Topographic model: nullcell.dat\n');
fprintf(fid,'DISTANCE | DEPTH weighting: wr.dat\n\n');
fprintf(fid,'Written by: Dominique Fournier\n');
fprintf(fid,'Last update: July 14th, 2014\n\n');
fprintf(fid,'INPUT FILES: \n');
fprintf(fid,'Mesh: \t\t\t %s \n',meshfile);
fprintf(fid,'Obsfile: \t\t\t %s \n',obsfile);
fprintf(fid,'Topography: \t\t\t %s \n',topofile);
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
    topo = read_UBC_topo([work_dir '\' topofile]);
end

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
cellID = find(nullcell==1);

% Get nodal discretization for octree levels
celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);

% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')
ndata = length(obsx);

% Create depth weighting
wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, wr_flag);
save([work_dir '\wr.dat'],'-ascii','wr');

% Compute sensitivities
% Pre-allocate

%% Create projection matrix for TMI
TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];

mcell = sum(nullcell);


T = zeros(ndata,3*mcell);

progress = -1;
tic 
fid = fopen([work_dir '\MAG3Csen.log'],'a');
for ii = 1:ndata
   
    % compute kernel for active cells
    [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),celln);

    T(ii,:) = TMI * [tx;ty;tz];
    
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf(fid,'Computed %i pct of data in %8.5f sec\n',d_iter*5,toc);
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end

fprintf(fid,'Sensitivity calculation completed in: %f min\n',toc/60);
fclose(fid);

%% Define p,s,t components

% Get magnetization matrix in p,s,t format
[p,s,t] = azmdip_2_pst(Dazm,I,mcell);

% Create orthogonal forward operators
% Primary (inducing)
Gp = T * (H * p);
save([work_dir '\Gp'],'Gp','-v7.3');
clear Gp;

Gs = T * (H * s);
save([work_dir '\Gs'],'Gs','-v7.3');
clear Gs;

Gt = T * (H * t);
save([work_dir '\Gt'],'Gt','-v7.3');
clear Gt;

