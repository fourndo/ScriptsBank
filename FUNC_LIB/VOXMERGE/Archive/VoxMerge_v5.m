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
% Revision: 4.1
% Last update: June 19th, 2014
%
% Change since version 3:
% - Script merges 2 co-planar meshes (same number of vertical cells)
% - New interpolation scheme using radial distance (much faster but same result)
% - Implement a topocheck to cut air cells as a final step.
% - Uses TriScatteredInterp instead of scatteredInterpolant

clear all
close all

root_dir = pwd;
addpath(root_dir)


%% INPUTS
% Load files
work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\VOXMERGE\Test';

cd(work_dir);

% Maximum search radius (m)
rangemax = 100;

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-1;

mod_file{1} = 'Tile1.den';
mod_file{2} = 'Tile2.den';
mod_file{3} = 'Tile3.den';

meshfile{1}=[work_dir '\' 'Tile1.msh'];
meshfile{2}=[work_dir '\' 'Tile2.msh'];
meshfile{3}=[work_dir '\' 'Tile3.msh'];

% model{3}=importdata([work_dir '\' 'Tile3.den']);
% meshfile{3}=[work_dir '\' 'Tile3.msh'];

% Specify final mesh file 
VM_meshfile = [];
% or cell size
dx = 10;
dy = 10;
dz = 10;

% Topofile
topofile = [];


%% \\\\\\\\\\\\\\\\\ %%
% SCRIPT STARTS HERE  %
%%\\\\\\\\\\\\\\\\\\ %%

% Load mesh files for all the tiles
ntiles = size(meshfile,2);
for ii = 1 : ntiles
    
    model{ii}=importdata([work_dir '\' mod_file{ii}]);
    
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
VM_mesh(2,3) = min(Z0);

XO_VM = min(X0);
YO_VM = min(Y0);
ZO_VM = max(Z0);

Xmax_VM = max(Xmax);
Ymax_VM = max(Ymax);
Zmax_VM = min(Zmax);

VM_mesh(1,1) = ceil( (Xmax_VM - XO_VM) / dx );
VM_mesh(1,2) = ceil( (Ymax_VM - YO_VM) / dy );
VM_mesh(1,3) = ceil( (ZO_VM - Zmax_VM) / dz );

VM_mesh(3,1:VM_mesh(1,1)) = ones(1,VM_mesh(1,1)) * dx;
VM_mesh(4,1:VM_mesh(1,2)) = ones(1,VM_mesh(1,2)) * dy; 
VM_mesh(5,1:VM_mesh(1,3)) = ones(1,VM_mesh(1,3)) * dz;

write_UBC_mesh(work_dir,XO_VM,YO_VM,ZO_VM,VM_mesh(3,1:VM_mesh(1,1)),VM_mesh(4,1:VM_mesh(1,2)),VM_mesh(5,1:VM_mesh(1,3)))
end

VM_model = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;

[VMy,VMx]=meshgrid((VM_mesh(2,2)+cumsum(VM_mesh(4,1:VM_mesh(1,2)))-VM_mesh(4,1:VM_mesh(1,2))/2),...
                            (VM_mesh(2,1)+cumsum(VM_mesh(3,1:VM_mesh(1,1)))-VM_mesh(3,1:VM_mesh(1,1))/2));

VMz = (VM_mesh(2,3)-cumsum(VM_mesh(5,1:VM_mesh(1,3)))+VM_mesh(5,1:VM_mesh(1,3))/2);                        
%% Propagate data upward
% This steps deals with conflicting topo between adjacente tiles.
% The finale mesh must be edited using TOPOCHECK to substrat the real topo

for ww = 1 : ntiles
    model{ww}=reshape(model{ww},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));
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
    levmax = 0;
    % Take the top slice of the model
    nullit = reshape(model{ww}(1,:,:),mesh{ww}(1,1),mesh{ww}(1,2));
    nullit = nullit~=ndv;
    
    weights{ww} = zeros(mesh{ww}(1,1),mesh{ww}(1,2));

    % Create cell center array for each tile
    x0=mesh{ww}(2,1);
    y0=mesh{ww}(2,2);
    nx = mesh{ww}(1,1);
    ny = mesh{ww}(1,2);
    nz = mesh{ww}(1,3);
    dx=mesh{ww}(3,1:nx);
    dy=mesh{ww}(4,1:ny); 
    [dY{ww},dX{ww}] = meshgrid(dy,dx);
    [YY{ww},XX{ww}]=meshgrid((y0+cumsum(dy)-dy/2),...
                            (x0+cumsum(dx)-dx/2));
    
    % First level along edges
    for ii=1:mesh{ww}(1,1)
        
        for jj=1:mesh{ww}(1,2)
            
           if  nullit(ii,jj)==1
               
               level = min([ii jj (mesh{ww}(1,1)-ii+1) (mesh{ww}(1,2)-jj+1)]);
               weights{ww}(ii,jj) = level;
               
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

                   level = min(min(weights{ww}(ii-1:ii+1,jj-1:jj+1)))+1;
                   weights{ww}(ii,jj) = level;

               end

               if level > Wmax(ww)

                   Wmax(ww) = level;

               end

            end

        end
        
    end
    
    % Display 2D weight pattern
    figure;imagesc(rot90(reshape(weights{ww},nx,ny)))
    temp = (['\bfWeights for tile: ' mod_file{ww}]);
%     set(gca,'XDir','Reverse')
%     set(gca,'YDir','Reverse')
    axis equal
    title(temp);
    colorbar

end

% Compute weights with cosine tapper and volume scaled
for ww = 1 : ntiles
    
    nx = mesh{ww}(1,1);
    ny = mesh{ww}(1,2);
    nz = mesh{ww}(1,3);
    
    active = (weights{ww}~=0);

    % Normalize the weights
    weights{ww}(active) = ( weights{ww}(active) -1 )/...
                                     (Wmax(ww)-1) + 1e-6;
    % Cosine tapper 
    weights{ww} = (-0.5*cos(-weights{ww}*pi)+0.5)  ;  

    % Cellsize weights
    weights{ww}(active) = weights{ww}(active) ./ (dX{ww}(active).*dY{ww}(active)) ;

    % Replicate weights for entire column (Z)
    W{ww} = repmat(weights{ww},[1 1 nz]);
    W{ww} = permute(W{ww},[3 1 2]);
   
end
   
clear weights
 
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

            if ii==114 && jj==217
                
                fprintf('Check\n');
                
            end
            
            % Cycle through all the tiles to find overlapping cells
            for kk = 1 : ntiles
                
                % Compute distance from current cell to all other cells in
                % 2D                
                R = sqrt((VMx(ii,jj) - XX{kk}).^2 +...
                        (VMy(ii,jj) - YY{kk}).^2);
                
                [r,index] = min(R(:));
                r = r + 1e-2;
                [i,j] = ind2sub([size(R,1) size(R,2)],index);
                
                % Get index for Z overlap
                
                
                % Compute weighted average using weights and distance from
                % cell center
                if r < rangemax
                    wght = (1/r * W{kk}(:,i,j)')';
                    
                    % Sum if multiple culumn with same R distance
                    numer = numer + sum( model{kk}(:,i,j) .*...
                                            wght , 2);
                    denom = denom + sum(  wght , 2);
                    
                end
                
            end
            
            if denom~=0
                
               VM_model(:,ii,jj) = numer ./ denom;
               
            end
            
            % Monitor progress and print to screen
             d_iter = floor(count/mcell*100);
            if  d_iter > progress

                fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end
            count = count +1;
            
        end
        
end

%% Compute topocheck on final mesh to remove aircells
% Get nodal discretization of final mesh
% xn = [VM_mesh(2,1) VM_mesh(2,1) + cumsum( VM_mesh(3,1:VM_mesh(1,1)) ) ];
% yn = [VM_mesh(2,2) VM_mesh(2,2) + cumsum( VM_mesh(4,1:VM_mesh(1,2)) ) ];
% zn = [VM_mesh(2,3) VM_mesh(2,3) - cumsum( VM_mesh(5,1:VM_mesh(1,3)) ) ];
% 
% [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
% 
% if isempty(topofile)==1
%     
%     nxn = length(xn);
%     nyn = length(yn);
%     topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
%     
% else
%     % Load topo
%     topo = read_UBC_topo(topofile);
% end
% 
% [nullcell,topocell,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% 
% % Assign air cells
% VM_model(nullcell==0) = ndv;    

%% Write back to file in UBC format   
VM_model=reshape(VM_model,VM_mesh(1,3)*VM_mesh(1,2)*VM_mesh(1,1),1);
filename='VOXMERGE_model.dat';
save([work_dir '\' filename],'-ascii','VM_model');
fprintf('Program VOXMERGE v5.1 succesfully ended\n')
fprintf('Merged model saved: %s\\%s\n',work_dir,filename);



    