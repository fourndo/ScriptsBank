% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\MAG';
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
fprintf('Begin Calculation for TopoCheck\n')
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
save([work_dir '\nullcell.dat'],'-ascii','nullcell');
load([work_dir '\nullcell.dat']);

% Get index of active cells
cellID = find(nullcell==1);

% Get nodal discretization for octree levels
celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);

% Load susceptibility model
mcell = size(celln,1);

%% Load observation loc  
[H, BI, BD, MI, MD, obsx, obsy, obsz, d, uncert] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')

% [XYZd_out] = Filter_xy([obsx obsy obsz d uncert],25);

ndata = length(obsx);

%% Create TMI projection operator
D = mod(450-BD,360);
Ptmi = [(cosd(BI) * cosd(D)) (cosd(BI) * sind(D)) sind(BI)];


%% Create depth weighting
fprintf(['Begin Calculation for' wr_flag 'weighting\n'])
wr = get_wr(obsx, obsy, obsz, D, BI, xn, yn, zn, nullcell, wr_flag, 3 , 0);
wr = wr(:);
save([work_dir '\wr.dat'],'-ascii','wr');

%% Compute sensitivities
% Pre-allocate

switch sen_flag
    case 'TxTyTz'
    
    Tx = zeros(ndata,3*mcell);
    Ty = zeros(ndata,3*mcell);
    Tz = zeros(ndata,3*mcell);
    
    otherwise
    
    G = zeros(ndata,3*mcell);
    
end

progress = -1;
tic 
fid = fopen([work_dir '\MAG3Csen.log'],'a');
for ii = 1:ndata
   
    % compute kernel for active cells
    [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),celln);

    switch sen_flag
        
        case 'TxTyTz'

        Tx(ii,:) = tx;
        Ty(ii,:) = ty;
        Tz(ii,:) = tz;


        
        otherwise

        G(ii,:) = Ptmi * [tx;ty;tz];
    
%                 if ii == 665
%             
%             fprintf('Shit!\n')
%             
%                 end
        
    end 
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf(fid,'Computed %i pct of data in %8.5f sec\n',d_iter*5,toc);
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end

fprintf(fid,'Sensitivity calculation completed in: %f min\n',toc/60);
fclose(fid);

switch sen_flag
        
    case 'TxTyTz'

    save([work_dir '\Tx'],'Tx','-v7.3');
    save([work_dir '\Ty'],'Ty','-v7.3');
    save([work_dir '\Tz'],'Tz','-v7.3');
    
    otherwise

    save([work_dir '\G'],'G','-v7.3');

end 



