% TileMerge
%
% INPUT: 
% work_dir: Working directory
% OPTIONAL
% meshfile : Final mesh
% 
%
% OUTPUT:
% Final merged model
% Author: D.Fournier
%
% Last update: June 9th, 2015
%

clear all
close all


%% INPUTS PARAMETERS
% Load files
work_dir ='C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\MAG\Composite\ALL_Tiles';

% Specify final mesh file 
out_meshfile = '..\Mesh_30m_padded.msh';

% Or leave empty and specify discretization
% out_meshfile = [];
% or if empty, use cell size
% dx = 5;
% dy = 5;
% dz = 5;

% cd(work_dir);

% Maximum search radius (m), recommended x2 smallest cell size
rangemax = 100;

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-100;

dsep = '\';

% Flag for topography: 'no_topo' | 'topofile' | 'nullfile'
% flag1 = 'no_topo';
flag1 = 'topofile'; topofile = '..\..\DEM_30m.topo';
% flag1 = 'nullfile'; nullfile = '..\nullcell.dat';

% Flag to remove padding cells from tiles
flag2 = 'rem_pad'; %'default'; % Either 'rem_pad'  | 'default'
xpad = 5;
ypad = 5;

%% Generate list of input models

mod_list=ls(work_dir);
count = 0;
% Input files
for ii = 3:size(mod_list,1)
    
    if ~isempty(regexp(mod_list(ii,:),'sus','match'))
        count = count + 1;
        
        in_file = strtrim(mod_list(ii,:));
        
        mod_file{count}=[work_dir dsep in_file];
        
        meshfile{count}=[work_dir dsep in_file(1:end-4) '.msh'];
        
    end
    
end
        

% model{3}=importdata([work_dir '\' 'Tile3.den']);
% meshfile{3}=[work_dir '\' 'Tile3.msh'];








%% \\\\\\\\\\\\\\\\\ %%
% SCRIPT STARTS HERE  %
%%\\\\\\\\\\\\\\\\\\ %%

% Load mesh files for all the tiles
ntiles = size(meshfile,2);
for ii = 1 : ntiles
    
    model{ii}=importdata(mod_file{ii});
    
    model{ii}(model{ii}==-1) = ndv;
    
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


%% Create final voxel merged mesh (out)
% If mesh file is provided
if isempty(out_meshfile)==0
    
    out_mesh = get_UBC_mesh([work_dir dsep out_meshfile]);
    
else
% If no meshfile provided then the final mesh is created from extent.
% Find the extent of all the tiles
out_mesh(2,1) = min(X0);
out_mesh(2,2) = min(Y0);
out_mesh(2,3) = max(Z0);

XO_out = min(X0);
YO_out = min(Y0);
ZO_out = max(Z0);

Xmax_out = max(Xmax);
Ymax_out = max(Ymax);
Zmax_out = min(Zmax);

out_mesh(1,1) = ceil( (Xmax_out - XO_out) / dx );
out_mesh(1,2) = ceil( (Ymax_out - YO_out) / dy );
out_mesh(1,3) = ceil( (ZO_out - Zmax_out) / dz );

out_mesh(3,1:out_mesh(1,1)) = ones(1,out_mesh(1,1)) * dx;
out_mesh(4,1:out_mesh(1,2)) = ones(1,out_mesh(1,2)) * dy; 
out_mesh(5,1:out_mesh(1,3)) = ones(1,out_mesh(1,3)) * dz;

write_UBC_mesh(work_dir,'VOXMERGE.msh',XO_out,YO_out,ZO_out,out_mesh(3,1:out_mesh(1,1)),out_mesh(4,1:out_mesh(1,2)),out_mesh(5,1:out_mesh(1,3)))
end

out_model = ones(out_mesh(1,3),out_mesh(1,1),out_mesh(1,2))*ndv;

[outy,outx]=meshgrid((out_mesh(2,2)+cumsum(out_mesh(4,1:out_mesh(1,2)))-out_mesh(4,1:out_mesh(1,2))/2),...
                            (out_mesh(2,1)+cumsum(out_mesh(3,1:out_mesh(1,1)))-out_mesh(3,1:out_mesh(1,1))/2));

outz = out_mesh(2,3)-cumsum(out_mesh(5,1:out_mesh(1,3))) + (out_mesh(5,1:out_mesh(1,3)))/2;

%% Assign tile number to interpolate for each column of final mesh
tile_ID = zeros(ntiles , out_mesh(1,1) , out_mesh(1,2));
for jj = 1 : out_mesh(1,2)
    
    for ii = 1 : out_mesh(1,1)
    
        for kk = 1 : ntiles
            
            if  outx(ii,jj) < Xmax(kk) && outx(ii,jj) > X0(kk) && outy(ii,jj) < Ymax(kk) && outy(ii,jj) > Y0(kk)
               
                tile_ID(kk,ii,jj) = 1;
                
            end
        end
        
    end
    
end
%% Propagate data upward
% This steps deals with conflicting topo between adjacente tiles.
% The finale mesh must be edited using TOPOCHECK to substrat the real topo

for ww = 1 : ntiles
    model{ww}=reshape(model{ww},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));
    
    switch flag2
        
        case 'rem_pad'
            % Temporary line to remove padding cells
            model{ww}(:,[1:xpad (end-xpad+1):end],:) = ndv;
            model{ww}(:,:,[1:ypad (end-ypad+1):end]) = ndv;
    
    end
    
    model{ww}=model{ww}-ndv;

    for i=1:mesh{ww}(1,1)
        
        for j=1:mesh{ww}(1,2)
            
            columnxy=model{ww}(:,i,j);
            ground = (mesh{ww}(1,3)-nnz(columnxy));
            
            if ground ~= mesh{ww}(1,3)
                columnxy(1:ground)=columnxy(ground+1)*ones(ground,1);
            end
            
            model{ww}(:,i,j)=columnxy;


        end
    end
    model{ww}=model{ww}+ndv;
end


%% Create 2D dissolving weight tiles
% 1- Increasing values towards the center of each tile. Levels are defined
%    by the number of cells away from the edge.
% 2- Normalized cosine weights from levels.
% 3- Normalized 2D cell size.

% Create disolving weights from egde:0 towards center:max
for ww = 1 : ntiles
    
    level = 1;
    leoutax = 0;
    % Take the top slice of the model
    nullit = reshape(model{ww}(1,:,:),mesh{ww}(1,1),mesh{ww}(1,2));
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
    [YY{ww},XX{ww}]=meshgrid((y0+cumsum(dy)-dy/2),...
                            (x0+cumsum(dx)-dx/2));
    
     ZZ{ww} = (z0-cumsum(dz)+dz/2);
    % First level along edges
    for ii=1:mesh{ww}(1,1)
        
        for jj=1:mesh{ww}(1,2)
            
           if  nullit(ii,jj)==1
               
               level = min([ii jj (mesh{ww}(1,1)-ii+1) (mesh{ww}(1,2)-jj+1)]);
               weights(ii,jj) = level;
               
           end
           
           if level > leoutax
               
               leoutax = level;
               
           end
            
        end
        
    end
    
    
    % Disolve from no-data-value cells
    for kk = 1 : leoutax;
        
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
%     set(gca,'XDir','Reverse')
%     set(gca,'YDir','Reverse')
%     axis equal
%     title(temp);
%     colorbar

    % Find overlapping cells vertically between final mesh and tiles
    % Deals with tiles that are offseted or with different cellsizes in Z
    Z_interp{ww} = zeros( out_mesh(1,3) ,  mesh{ww}(1,3) );
    
    for jj = 1 : out_mesh(1,3)
        
        % Find cells within current final cell
        dzz = abs(outz(jj) - ZZ{ww});
        
        index = dzz < out_mesh(5,jj);
        
        if sum(index) > 0
            
            Z_interp{ww}(jj,index) = (dzz(index) + 1e-2).^-1 / sum((dzz(index) + 1e-2).^-1) ;
            
        end
        
    end
    
    Z_interp{ww} = sparse(Z_interp{ww});
    
end

 
%% Beggining the merging computations
% Iterate over all cells of final mesh and find the closest cells to 
% average. Computes an harmonic weigthed average with inverse distance
% weight. Assumes that all mesh have the same cells vertically.  
dx = min(out_mesh(3,1:out_mesh(1,1)));
dy = min(out_mesh(4,1:out_mesh(1,2)));
mcell = out_mesh(1,1)*out_mesh(1,2);

count = 1;
progress = -1;
tic
for ii = 1 : out_mesh(1,1)
        
        for jj = 1 : out_mesh(1,2)
            
            numer = 0;
            denom = 0;
            
            % Cycle through all the tiles to find overlapping cells
            tile_in = find(tile_ID(:,ii,jj) ==1);
            
            for kk = 1 : length(tile_in)
                
                % Compute distance from current cell to all other cells in
                % 2D                
                R = sqrt((outx(ii,jj) - XX{ tile_in(kk) }).^2 +...
                        (outy(ii,jj) - YY{ tile_in(kk) }).^2);
                
                [r,index] = min(R(:));
                
%                 r = r + 1e-2;
                [i,j] = ind2sub([size(R,1) size(R,2)],index);
                
                % Compute weighted average using weights and distance from
                % cell center
                if r < rangemax
                    
                    invr = min([1/r 1/(dx/2)]);
                    wght = (1/invr * W{ tile_in(kk) }(:,i,j)')';
                    
                    % Sum if multiple culumn with same R distance
                    numer = numer + sum( Z_interp{ tile_in(kk) } * (model{ tile_in(kk) }(:,i,j) .* wght) , 2);
                    denom = denom + sum( Z_interp{ tile_in(kk) } * wght , 2);

                    
                end
                
            end
            
            if sum(denom)~=0
                
               out_model(:,ii,jj) =  (numer ./ denom);
               
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

%% Compute topocheck on final mesh to remove aircells
switch flag1
    
    case 'topofile'
        
        % Load topography file (UBC format)
        topo = read_UBC_topo([work_dir '\' topofile]);      
        
        % Get nodal discretization of final mesh
        xn = [out_mesh(2,1) out_mesh(2,1) + cumsum( out_mesh(3,1:out_mesh(1,1)) ) ];
        yn = [out_mesh(2,2) out_mesh(2,2) + cumsum( out_mesh(4,1:out_mesh(1,2)) ) ];
        zn = [out_mesh(2,3) out_mesh(2,3) - cumsum( out_mesh(5,1:out_mesh(1,3)) ) ];
        
        [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
                
        [nullcell,topocell,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
        
        save([work_dir '\nullcell.dat'],'-ascii','nullcell')
        fprintf('Nullcell has been saved to: \n %s \n for future use\n',work_dir);
        
        % Assign air cells
        out_model(nullcell==0) = ndv;    
                
    case 'nullfile'

        nullcell = load([work_dir '\' nullfile]);
        out_model(nullcell==0) = ndv;
        
    otherwise
        
           
end
%% Write back to file in UBC format   
out_model=reshape(out_model,out_mesh(1,3)*out_mesh(1,2)*out_mesh(1,1),1);
out_model(isnan(out_model)) = ndv;
filename='VOXMERGE_model.dat';
save([work_dir '\' filename],'-ascii','out_model');
fprintf('Program VOXMERGE v5.1 succesfully ended\n')
fprintf('Merged model saved: %s\\%s\n',work_dir,filename);



    