% Mira Geoscience
% Voxmerge.mat
%
% INPUT: 
% Models (UBC format): The order of import is not important.
% Mesh file (UBC format): One mesh file per model
% Final mesh (UBC format): Mesh that covers at least all the tiles.
% 
%
% OUTPUT:
% Files the same size in UBC format
% Author: D.Fournier
% Revision: 7
% Last update: January, 2016
%
% Change since version 5:
% - Script meshes that don't have to be co-planar
% - New interpolation scheme using radial distance (much faster but same result)
% - Implement a topocheck to cut air cells as a final step.
% - Uses TriScatteredInterp instead of scatteredInterpolant
%
% Change since version 6:
% -Allow to import rotated grids. Requires rotation parameters

close all
clear all

% root_dir = pwd;
addpath .\functions

%% INPUTS PARAMETERS
% Load files
work_dir ='C:\Users\DominiqueFournier\ownCloud\Research\Kevitsa\Modeling\MAG\Aiborne\100m\Tiles_l2';

outFile = 'SimPEG_Merged.sus';
% cd(work_dir);

% Maximum search radius (m), recommended x2 smallest cell size
rangemax = 1000;

% Small constant for smoothing factor (0:max)
delta = 20;

% Number of neighbours to interpolate from
n_interp = 9;

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=0;

dsep = '\';

%## Input files ##
% Specify final mesh file 
VM_meshfile = [work_dir dsep '..\Mesh_global_100m.msh'];

% Flag for topography: 'no_topo' | 'topofile' | 'nullfile'
% flag1 = 'no_topo';
flag1 = 'topofile'; topofile = '..\..\..\Kevitsa.topo';
% flag1 = 'nullfile'; topofile = '..\Forward\nullcell.dat';

% Flag to remove padding cells from tiles
flag2 = 'rem_pad'; % Either 'rem_pad'  | 'default'
% If flag is 'rem_pad' -> Specify padding to remove in [E,W,S,N]
pad = [1 1 1 1];
% pad{2} = [16 16 16 16];
% pad{3} = [16 16 16 16];

% Flag for interpolation space
flag3 = 'linear'; %Either 'linear' | 'log'



mType = 'fld';
files = dir([work_dir dsep '*.' mType]);


for ii = 1:size(files,1)
    
    id = regexp(files(ii).name,'\d+','match');
    
    if ~isempty(id)
        mod_file{ii} = [work_dir dsep files(ii).name];


        meshfile{ii} = [work_dir dsep 'MAG_Tile' num2str(id{1}) '.msh'];
        T{ii} = [0 0 0 0 0];
    end
end

% mod_file{1} = [work_dir dsep 'maginv3d_002.sus'];
% meshfile{1} = [work_dir dsep '..\Mesh_A2_ROT_v2.msh'];

% mod_file{1} = [work_dir dsep '..\Forward\Mesh_20m_Dens\gzinv3d_001.den'];
% mod_file{2} = [work_dir dsep 'tile2_chg_new.chg'];
% mod_file{3} = [work_dir dsep 'Model13200.chg'];
% mod_file{4} = [work_dir dsep 'Model13700.chg'];
% mod_file{5} = [work_dir dsep 'Model14200.chg'];
% mod_file{6} = [work_dir dsep 'Model14700.chg'];
% mod_file{7} = [work_dir dsep 'Model15200.chg'];
% mod_file{8} = [work_dir dsep 'Model15700.chg'];
% mod_file{9} = [work_dir dsep 'Model15950.chg'];
% mod_file{10} = [work_dir dsep 'Model16200.chg'];

% meshfile{1} = [work_dir dsep 'Mesh_20m.msh'];
% meshfile{2} = [work_dir dsep 'mesh_tile2_new.msh'];
% meshfile{3} = [work_dir dsep 'Mesh13200.msh'];
% meshfile{4} = [work_dir dsep 'Mesh13700.msh'];
% meshfile{5} = [work_dir dsep 'Mesh14200.msh'];
% meshfile{6} = [work_dir dsep 'Mesh14700.msh'];
% meshfile{7} = [work_dir dsep 'Mesh15200.msh'];
% meshfile{8} = [work_dir dsep 'Mesh15700.msh'];
% meshfile{9} = [work_dir dsep 'Mesh15950.msh'];
% meshfile{10} = [work_dir dsep 'Mesh16200.msh'];

% meshfile{1}=[work_dir dsep 'FLIN_ZTEM_100m.msh'];
% meshfile{2}=[work_dir dsep 'rot14_Y2_20m_grid.mesh'];
% meshfile{3}=[work_dir dsep 'rot14_Y3_20m_grid.msh'];


% Define mesh transformation parameters [x0 y0 dx dy theta]
% T{1} = [657000 5973000 0 0 -45];
% T{2} = [0 0 0 0 0];
% T{3} = [313700 6074300 0 0 0];
% T{4} = [313700 6074300 0 0 0];
% T{5} = [313700 6074300 0 0 0];
% T{6} = [313700 6074300 0 0 0];
% T{7} = [313700 6074300 0 0 0];
% T{8} = [313700 6074300 0 0 0];
% T{9} = [313700 6074300 0 0 0];
% T{10} = [313700 6074300 0 0 0];

% or if empty, use cell size
dx = 5;
dy = 5;
dz = 5;


%% \\\\\\\\\\\\\\\\\ %%
% SCRIPT STARTS HERE  %
%%\\\\\\\\\\\\\\\\\\ %%

% Load mesh files for all the tiles
ntiles = size(meshfile,2);
for ii = 1 : ntiles
    
    mtemp = importdata(mod_file{ii});
    mtemp(mtemp==-1) = ndv;
    if strcmp(mType,'fld')
        model{ii}{1} = mtemp(:,1);
        model{ii}{2} = mtemp(:,2);
        model{ii}{3} = mtemp(:,3);
    else
        model{ii}{1} = mtemp;
    end
        
    
    mesh{ii} = get_UBC_mesh(meshfile{ii});
    
    X0(ii) = mesh{ii}(2,1);
    Y0(ii) = mesh{ii}(2,2);
    Z0(ii) = mesh{ii}(2,3);

    Xmax(ii) = mesh{ii}(2,1)+sum(mesh{ii}(3,:));
    Ymax(ii) = mesh{ii}(2,2)+sum(mesh{ii}(4,:));
    Zmax(ii) = mesh{ii}(2,3)-sum(mesh{ii}(5,:));
    
    dxmin(ii) = min(mesh{ii}(3,mesh{ii}(3,:)~=0));
    dymin(ii) = min(mesh{ii}(4,mesh{ii}(4,:)~=0));
    dzmin(ii) = min(mesh{ii}(5,mesh{ii}(5,:)~=0));
    
end


%% Create final voxel merged mesh (VM)
% If mesh file is provided
if isempty(VM_meshfile)==0
    
    VM_mesh = get_UBC_mesh(VM_meshfile);
    
else
% If no meshfile provided then the final mesh is created from extent.
% Find the extent of all the tiles
VM_mesh(2,1) = min(X0);
VM_mesh(2,2) = min(Y0);
VM_mesh(2,3) = max(Z0);

XO_VM = min(X0);
YO_VM = min(Y0);
ZO_VM = max(Z0);

Xmax_VM = max(Xmax);
Ymax_VM = max(Ymax);
Zmax_VM = min(Zmax);

VM_mesh(1,1) = ceil( (Xmax_VM - XO_VM) / dx );
VM_mesh(1,2) = ceil( (Ymax_VM - YO_VM) / dy );
VM_mesh(1,3) = ceil( (ZO_VM - Zmax_VM) / dz );%

VM_mesh(3,1:VM_mesh(1,1)) = ones(1,VM_mesh(1,1)) * dx;
VM_mesh(4,1:VM_mesh(1,2)) = ones(1,VM_mesh(1,2)) * dy; 
VM_mesh(5,1:VM_mesh(1,3)) = ones(1,VM_mesh(1,3)) * dz;

write_UBC_mesh(work_dir,'VOXMERGE.msh',XO_VM,YO_VM,ZO_VM,VM_mesh(3,1:VM_mesh(1,1)),VM_mesh(4,1:VM_mesh(1,2)),VM_mesh(5,1:VM_mesh(1,3)))
end

if strcmp(mType,'fld')
    VM_model{1} = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
    VM_model{2} = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
    VM_model{3} = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
else
    VM_model{1} = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
end
[VMy,VMx]=meshgrid((VM_mesh(2,2)+cumsum(VM_mesh(4,1:VM_mesh(1,2)))-VM_mesh(4,1:VM_mesh(1,2))/2),...
                            (VM_mesh(2,1)+cumsum(VM_mesh(3,1:VM_mesh(1,1)))-VM_mesh(3,1:VM_mesh(1,1))/2));

VMz = VM_mesh(2,3)-cumsum(VM_mesh(5,1:VM_mesh(1,3))) + (VM_mesh(5,1:VM_mesh(1,3)))/2;

%% Assign tile number to interpolate for each column of final mesh
tile_ID = zeros(ntiles , VM_mesh(1,1) , VM_mesh(1,2));
for jj = 1 : VM_mesh(1,2)
    
    for ii = 1 : VM_mesh(1,1)
    
        for kk = 1 : ntiles
            
%             if  VMx(ii,jj) < Xmax(kk) && VMx(ii,jj) > X0(kk) && VMy(ii,jj) < Ymax(kk) && VMy(ii,jj) > Y0(kk)
               
                tile_ID(kk,ii,jj) = 1;
                
%             end
        end
        
    end
    
end
%% Propagate data upward
% This steps deals with conflicting topo between adjacente tiles.
% The finale mesh must be edited using TOPOCHECK to substrat the real topo

for ww = 1 : ntiles
    
    for cc = 1 : length(model{ww})
        model{ww}{cc}=reshape(model{ww}{cc},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));

        switch flag2

            case 'rem_pad'
                % Temporary line to remove padding cells
                model{ww}{cc}(:,[1:pad(1) (end-pad(2)+1):end],:) = ndv;
                model{ww}{cc}(:,:,[1:pad(3) (end-pad(4)+1):end]) = ndv;

        end

        model{ww}{cc}=model{ww}{cc}-ndv;

        for i=1:mesh{ww}(1,1)

            for j=1:mesh{ww}(1,2)

                columnxy=model{ww}{cc}(:,i,j);
                ground = (mesh{ww}(1,3)-nnz(columnxy));

                if ground ~= mesh{ww}(1,3)
                    columnxy(1:ground)=columnxy(ground+1)*ones(ground,1);
                end

                model{ww}{cc}(:,i,j)=columnxy;


            end
        end
        model{ww}{cc}=model{ww}{cc}+ndv;
        
    end
end


%% Create 2D dissolving weight tiles
% 1- Increasing values towards the center of each tile. Levels are defined
%    by the number of cells away from the edge.
% 2- Normalized cosine weights from levels.
% 3- Normalized 2D cell size.

% Create disolving weights from egde:0 towards center:max
for ww = 1 : ntiles
    
    level = 1;
    levmax = 0;
    % Take the top slice of the model
    nullit = reshape(model{ww}{1}(1,:,:),mesh{ww}(1,1),mesh{ww}(1,2));
    nullit = nullit~=ndv;
    
    weights = zeros(mesh{ww}(1,1),mesh{ww}(1,2));

    % Create cell center array for each tile
    x0=mesh{ww}(2,1);
    y0=mesh{ww}(2,2);
    z0=mesh{ww}(2,3);
    nx = mesh{ww}(1,1);
    ny = mesh{ww}(1,2);
    nz = mesh{ww}(1,3);
    dx=mesh{ww}(3,1:nx);
    dy=mesh{ww}(4,1:ny);
    
    
    dz=mesh{ww}(5,1:nz);
    [dY{ww},dX{ww}] = meshgrid(dy,dx);
    [yy,xx]=meshgrid((y0+cumsum(dy)-dy/2),...
                            (x0+cumsum(dx)-dx/2));
    
    % ## Version 7 ##
    % Rotate the coordinate axis
    rot = [cosd(T{ww}(5)) -sind(T{ww}(5));
           sind(T{ww}(5)) cosd(T{ww}(5))];
    
    % Translate to centroid and rotate
    xy = rot*[xx(:)' -  T{ww}(1);yy(:)' - T{ww}(2)];
    
    % Translate tile back to global + shift
    xy(1,:) = xy(1,:) + T{ww}(1) + T{ww}(3);
    xy(2,:) = xy(2,:) + T{ww}(2) + T{ww}(4);
    
    YY{ww} = reshape(xy(2,:),size(dY{ww}));
    XX{ww} = reshape(xy(1,:),size(dX{ww}));
    
    ZZ{ww} = (z0-cumsum(dz)+dz/2);
    % First level along edges
    for ii=1:mesh{ww}(1,1)
        
        for jj=1:mesh{ww}(1,2)
            
           if  nullit(ii,jj)==1
               
               level = min([ii jj (mesh{ww}(1,1)-ii+1) (mesh{ww}(1,2)-jj+1)]);
               weights(ii,jj) = level;
               
           end
           
           if level > levmax
               
               levmax = level;
               
           end
            
        end
        
    end
    
    
    % Disolve from no-data-value cells
    for kk = 1 : levmax;
        
        level = 1;
        Wmax(ww) = 0;
    
        for ii=2:nx-1

            for jj=2:ny-1

               if  nullit(ii,jj)==1

                   level = min(min(weights(ii-1:ii+1,jj-1:jj+1)))+1;
                   weights(ii,jj) = level;

               end

               if level > Wmax(ww)

                   Wmax(ww) = level;

               end

            end

        end
        
    end
    
    
    weights(nullit) = ( weights(nullit) -1 )/...
                                     (Wmax(ww)-1) + 1e-6;
    % Cosine tapper 
    weights = (-0.5*cos(-weights*pi)+0.5)  ;  

%     % Reverse weights if needed
%     if ww == 1
%         
%         weights(nullit) =  (max(max(weights(nullit))) - weights(nullit));
%         
%     end
    
    % Cellsize weights
    weights(nullit) = weights(nullit) ./ (dX{ww}(nullit).*dY{ww}(nullit)) ;  

    % Compute weights with cosine tapper and volume scaled
    % Replicate weights for entire column (Z)
    W{ww} = repmat(weights,[1 1 nz]);
    W{ww} = permute(W{ww},[3 1 2]);
    
    % Display 2D weight pattern
%     figure;imagesc(rot90(reshape(weights,nx,ny)))
%     temp = ('\bfWeights for tile: 1');
% %     set(gca,'XDir','Reverse')
% %     set(gca,'YDir','Reverse')
%     axis equal
%     title(temp);
%     colorbar

    % Find overlapping cells vertically between final mesh and tiles
    % Deals with tiles that are offseted or with different cellsizes in Z
    Z_interp{ww} = zeros( VM_mesh(1,3) ,  mesh{ww}(1,3) );
    
    for jj = 1 : VM_mesh(1,3)
        
        % Find cells within current final cell
        dzz = abs(VMz(jj) - ZZ{ww});
        
        index = dzz < rangemax;
        
        if sum(index) > 0
            
            Z_interp{ww}(jj,index) = (dzz(index) + 1e-2).^-1 / sum((dzz(index) + 1e-2).^-1) ;
            
        end
        
    end
    
    Z_interp{ww} = sparse(Z_interp{ww});
    
end

%% Convert model to log space if required by the flag

if strcmp(flag3,'log')
    for ii = 1 : length(model)
        
        for cc = 1 : length(model{ii})

            model{ii}{cc} = log10(model{ii}{cc});
        end

    end
end

%% Beggining the merging computations
% Iterate over all cells of final mesh and find the closest cells to 
% average. Computes an harmonic weigthed average with inverse distance
% weight. Assumes that all mesh have the same cells vertically.  
dx = min(VM_mesh(3,1:VM_mesh(1,1)));
dy = min(VM_mesh(4,1:VM_mesh(1,2)));
mcell = VM_mesh(1,1)*VM_mesh(1,2);

count = 1;
progress = -1;
tic
for ii = 1 : VM_mesh(1,1)
        
        for jj = 1 : VM_mesh(1,2)
            
            numer = 0;
            denom = 0;
            
            % Cycle through all the tiles to find overlapping cells
            tile_in = find(tile_ID(:,ii,jj) ==1);
            
            for cc = 1 : length(VM_model)
                for kk = 1 : length(tile_in)

                    % Compute distance from current cell to all other cells in
                    % 2D                
                    R = sqrt((VMx(ii,jj) - XX{ tile_in(kk) }).^2 +...
                            (VMy(ii,jj) - YY{ tile_in(kk) }).^2);

                    [r,index] = sort(R(:));

    %                 r = r + 1e-2;


                    % Compute weighted average using weights and distance from
                    % cell center
                    ll = r(1:n_interp) < rangemax;

                    [i,j] = ind2sub([size(R,1) size(R,2)],index(ll));

                    invr = spdiags(1./(r(ll)+delta),0,sum(ll),sum(ll));
                    ww   = reshape(W{ tile_in(kk) }(:,index(ll)),size(W{ tile_in(kk) },1),sum(ll));
                    wght = (invr * ww');

                    % Sum if multiple culumn with same R distance
                    mtemp = reshape(model{ tile_in(kk) }{cc}(:,index(ll)),size(model{ tile_in(kk) }{cc},1),sum(ll));
                    numer = numer +  Z_interp{ tile_in(kk) } *sum( (mtemp .* wght') , 2);
                    denom = denom +  Z_interp{ tile_in(kk) } *sum( wght' , 2);


                end

                if sum(denom)~=0

                   VM_model{cc}(:,ii,jj) =  (numer ./ denom);

                end
            end
            
            % Monitor progress and print to screen
             d_iter = floor(count/mcell*100);
            if  d_iter > progress

                fprintf('Merged %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end
            count = count +1;
            
        end
        
end

%% Inverse log model if flag3==log
if strcmp(flag3,'log')

    for cc = 1 : length(VM_model)
        VM_model{cc}(VM_model~=ndv) = 10.^VM_model{cc}(VM_model{cc}~=ndv);
    end

end

%% Compute topocheck on final mesh to remove aircells
switch flag1
    
    case 'topofile'
        
        % Load topography file (UBC format)
        topo = read_UBC_topo([work_dir '\' topofile]);      
        
        % Get nodal discretization of final mesh
        xn = [VM_mesh(2,1) VM_mesh(2,1) + cumsum( VM_mesh(3,1:VM_mesh(1,1)) ) ];
        yn = [VM_mesh(2,2) VM_mesh(2,2) + cumsum( VM_mesh(4,1:VM_mesh(1,2)) ) ];
        zn = [VM_mesh(2,3) VM_mesh(2,3) - cumsum( VM_mesh(5,1:VM_mesh(1,3)) ) ];
        
        [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
                
        [nullcell,topocell,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
        
        save([work_dir '\nullcell.dat'],'-ascii','nullcell')
        fprintf('Nullcell has been saved to: \n %s \n for future use\n',work_dir);
        
        % Assign air cells
        for cc = 1 : length(VM_model)
            VM_model{cc} = VM_model{cc}(:);
            if strcmp(mType,'fld')
                VM_model{cc}(nullcell==0) = 0;
                VM_model{cc}(VM_model{cc}==ndv) = 0;
            else
                VM_model{cc}(nullcell==0) = ndv; 
        
            end
        end    
    case 'nullfile'

        nullcell = load([work_dir '\' topofile]);
        % Assign air cells
        for cc = 1 : length(VM_model)
            VM_model{cc}(nullcell==0) = ndv; 
        end
        
    otherwise
        
           
end
%% Write back to file in UBC format   
if strcmp(mType,'fld')
    VM_model{1}(isnan(VM_model{1})) = ndv;
    VM_model{2}(isnan(VM_model{2})) = ndv;
    VM_model{3}(isnan(VM_model{3})) = ndv;
    VM_out = [VM_model{1}';VM_model{2}';VM_model{3}']';
    VM_amp = sum([VM_model{1}';VM_model{2}';VM_model{3}']'.^2, 2).^0.5;
    save([work_dir '\AMP_' outFile],'-ascii','VM_amp');
else
    
    
    VM_out=VM_model{1};
    VM_out(isnan(VM_out)) = ndv;
end
save([work_dir '\' outFile],'-ascii','VM_out');
fprintf('Program VOXMERGE v5.1 succesfully ended\n')
fprintf('Merged model saved: %s\\%s\n',work_dir,outFile);



    