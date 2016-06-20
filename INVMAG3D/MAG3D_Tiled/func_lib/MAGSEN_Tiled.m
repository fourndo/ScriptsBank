function [G] = MAGSEN_Tiled(work_dir,dsep,xn,yn,zn, H, BI, BD, MI, MD, obsx, obsy, obsz,nullcell, wr_flag, sen_flag, R, R0, M)
% Function MAGSEN_Tiles(work_dir,meshfile,obsfile,topofile, wr_flag, sen_flag)
% Generate sensitivity matrix for MAG3D style inversion
% Dominique Fournier 2013/01/23
%
% INPUT VARIABLES
% work_dir : Working directory where input files are and sensitivity saved
% meshfile : Mesh file
% obsfile  : Observation file in UBC format
% topofile : Either UBC format file or []
% wr_flaf  : Type of sensitivity weighting, either DEPTH or DISTANCE
% sen_flag : Type of sensitivity matrix, either TxTyTz or G

% Load mesh file and convert to vectors (UBC format)
% [xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

% Write logfile
fid = fopen([work_dir dsep 'MAG3Csen.log'],'w');
fprintf(fid,'MAG3Csen\n');
fprintf(fid,'Generates sparse matrices for magnetostatic forward modeling: Tx, Ty, Tz\n');
fprintf(fid,'Topographic model: nullcell.dat\n');
fprintf(fid,'DISTANCE | DEPTH weighting: wr.dat\n\n');
fprintf(fid,'Written by: Dominique Fournier\n');
fprintf(fid,'Last update: July 14th, 2014\n\n');
fprintf(fid,'INPUT FILES: \n');
fprintf(fid,'Weighting: \t\t\t %s \n',wr_flag);
fclose(fid);
%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

nx = length(xn)-1;
ny = length(yn)-1;
nz = length(zn)-1;

% Create nullcell
% fprintf('Begin Calculation for TopoCheck\n')
% [nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% save([work_dir dsep 'nullcell.dat'],'-ascii','nullcell');
% load([work_dir dsep 'nullcell.dat']);

% Get index of active cells
cellID = find(nullcell==1);

% Get nodal discretization for octree levels
% celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);
% Create node pair for each cell
znzn = kron(kron(ones(ny,1),ones(nx,1)),[zn(1:end-1)' zn(2:end)']);
xnxn = kron(kron(ones(ny,1),[xn(1:end-1)' xn(2:end)']),ones(nz,1));
ynyn = kron(kron([yn(1:end-1)' yn(2:end)'],ones(nx,1)),ones(nz,1));

celln = [znzn(:,1) xnxn(:,1) ynyn(:,1) znzn(:,2) xnxn(:,2) ynyn(:,2)];
celln = celln(nullcell==1,:);

% Load susceptibility model
mcell = size(celln,1);

% Load observation loc

% [H, I, Dazm, D, obsx, obsy, obsz, d, uncert] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')

% [XYZd_out] = Filter_xy([obsx obsy obsz d uncert],25);

ndata = length(obsx);

%% Create TMI projection operator
Ptmi = [(cosd(BI) * cosd(BD)) (cosd(BI) * sind(BD)) sind(BI)];


%% Create depth weighting
if strcmp(wr_flag,'DISTANCE')==1 || strcmp(wr_flag,'DEPTH')==1
    fprintf(['Begin Calculation for ' wr_flag ' weighting\n'])
    wr = get_wr(obsx, obsy, obsz, MD, MI, xn, yn, zn, nullcell, wr_flag, R, R0);
    wr = wr(:);
    save([work_dir dsep 'wr.dat'],'-ascii','wr');
end

%% Compute sensitivities
% Pre-allocate

switch sen_flag
    case 'Guvw'
    
        G{1} = zeros(ndata,3*mcell);
        
    case 'Gpst'
    
        G{1} = zeros(ndata,3*mcell);
        % Case cartesian coordinates
        [P,S,T] = azmdip_2_pst(BD,BI,mcell);
        
    case 'GxGyGz'
        
        G{1} = zeros(ndata,mcell);
        G{2} = zeros(ndata,mcell);
        G{3} = zeros(ndata,mcell);
    
    otherwise
    
        G{1} = zeros(ndata,mcell);
%     Ty = spalloc(ndata,3*mcell,0);
%     Tz = spalloc(ndata,3*mcell,0);
    
end

progress = -1;
tic 
fid = fopen([work_dir dsep 'MAG3Csen.log'],'a');
fprintf(['Begin Sensitivity CALC\n'])
for ii = 1:ndata
   
    % compute kernel for active cells
    [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),celln);

    switch sen_flag
        case 'Guvw'
            G{1}(ii,:) = Ptmi * [tx;ty;tz] * H;
                
        case 'Gpst'

            G{1}(ii,:) = [Ptmi * [tx;ty;tz] * (H * P) Ptmi * [tx;ty;tz] * (H * S) Ptmi * [tx;ty;tz] * (H * T)];

        case 'GxGyGz'
            
            G{1}(ii,:) = tx * M;
            G{2}(ii,:) = ty * M;
            G{3}(ii,:) = tz * M;

        
        otherwise

        G{1}(ii,:) = Ptmi * [tx;ty;tz] * M;
    

    end 
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf(fid,'Computed %i pct of data in %8.5f sec\n',d_iter*5,toc);
        fprintf('Computed %i pct of data in %8.5f sec\r',d_iter*5,toc)
        progress = d_iter;

    end
            
end

fprintf(fid,'Sensitivity calculation completed in: %f min\n',toc/60);
fclose(fid);

% switch sen_flag
%         
%     case 'TxTyTz'
% 
%     save([work_dir dsep 'Tx'],'Tx','-v7.3');
%     save([work_dir dsep 'Ty'],'Ty','-v7.3');
%     save([work_dir dsep 'Tz'],'Tz','-v7.3');
%     
%     otherwise
% 
%     save([work_dir dsep 'G'],'G','-v7.3');
% 
% end 



