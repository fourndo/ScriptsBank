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
% Revision: XX
% Last update: August 23th, 2013
%
%**** CODE NOT FINISHED ****
clear all
close all
root_dir = pwd;
addpath(root_dir)
addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB'

%% INPUTS
% Load files
work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\Desktop\Thomas_merge';

cd(work_dir);

model{1}=importdata([work_dir '\' 'VPmg_Peppler_Hobdad_Den_4Targeting_prop.den']);
meshfile{1}=[work_dir '\' 'VPmg_Peppler_Hobdad_Den_4Targeting_compact.msh'];

model{2}=importdata([work_dir '\' 'VPmg_Lamelee_Den_4Targeting_prop.den']);
meshfile{2}=[work_dir '\' 'VPmg_Lamelee_Den_4Targeting_compact.msh'];

% model{3}=importdata([work_dir '\' 'Inv_IP_Tile3.chg']);
% meshfile{3}=[work_dir '\' 'ROT40_Titan_Tile3_nopad.msh'];

% Specify final mesh file 
VM_meshfile = [work_dir '\LAM_And_PEP_50m_VOI_Merged_compact.msh'];
% or cell size
dx = 50;
dy = 50;

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-99999;

% Maximum search radius (m)
rangemax = 100;

%% Load mesh files for all the tiles
ntiles = size(meshfile,2);
for ii = 1 : ntiles
    
    mesh{ii} = get_UBC_mesh(meshfile{ii});
    X0{ii} = mesh{ii}(2,1);
    Y0{ii} = mesh{ii}(2,2);
    Z0{ii} = mesh{ii}(2,3);

    Xmax{ii} = mesh{ii}(2,1)+sum(mesh{ii}(3,:));
    Ymax{ii} = mesh{ii}(2,2)+sum(mesh{ii}(4,:));
    Zmax{ii} = mesh{ii}(2,3)-sum(mesh{ii}(5,:));
    
    dxmin{ii} = min(mesh{ii}(3,mesh{ii}(3,:)~=0));
    dymin{ii} = min(mesh{ii}(4,mesh{ii}(4,:)~=0));
    dzmin{ii} = min(mesh{ii}(5,mesh{ii}(5,:)~=0));
    
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

Xmax_VM = max(Xmax);
Ymax_VM = max(Xmax);
Zmax_VM = max(Xmax);

VM_mesh(1,1) = ceil( (Xmax_VM - XO_VM) / min(dxmin) );
VM_mesh(1,2) = ceil( (Ymax_VM - YO_VM) / min(dymin) );
VM_mesh(1,3) = mesh{1}(1,3);

VM_mesh(3,1:VM_mesh(1,1)) = ones(1,VM_mesh(1,1)) * dx;
VM_mesh(4,1:VM_mesh(1,2)) = ones(1,VM_mesh(1,2)) * dy; 
VM_mesh(5,1:VM_mesh(1,3)) = mesh{1}(5,:); 

end

VM_model = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
% [VMdy,VMdx]=meshgrid(VM_mesh(4,1:VM_mesh(1,2)),VM_mesh(3,1:VM_mesh(1,1)));

[VMy,VMx]=meshgrid((VM_mesh(2,2)+cumsum(VM_mesh(4,1:VM_mesh(1,2)))-VM_mesh(4,1:VM_mesh(1,2))/2),...
                            (VM_mesh(2,1)+cumsum(VM_mesh(3,1:VM_mesh(1,1)))-VM_mesh(3,1:VM_mesh(1,1))/2));

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
Wmax = 0;
for ww = 1 : ntiles
    
    % Take the top slice of the model
    weights{ww} = reshape(model{ww}(1,:,:),mesh{ww}(1,1),mesh{ww}(1,2));

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
    
    for ii=1:mesh{ww}(1,1)
        
        for jj=1:mesh{ww}(1,2)
            
           if  weights{ww}(ii,jj)~=ndv
               
               level = min([ii jj (mesh{ww}(1,1)-ii+1) (mesh{ww}(1,2)-jj+1)]);
               weights{ww}(ii,jj) = level;
               
           else
               
               weights{ww}(ii,jj) = 0;
               
           end
           
           if level > Wmax
               
               Wmax = level;
               
           end
            
        end
        
    end
    

end


for ww = 1 : ntiles
       % Normalize the weights
   weights{ww}(weights{ww}~=0) = ( weights{ww}(weights{ww}~=0) -1 )/...
                                     (Wmax-1) + 1e-6;
   % Cosine tapper 
   weights{ww} = (-0.5*cos(-weights{ww}*pi)+0.5)  ;  
   
   % Cellsize weights
   weights{ww}(weights{ww}~=0) = weights{ww} ./ (dX{ww}.*dY{ww}) ;
       
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

            
            % Cycle through all the tiles to find overlapping cells
            for kk = 1 : ntiles
                
                logic = zeros(VM_mesh(1,1),VM_mesh(1,2));
                range = -1;
            
                while sum(sum(logic))==0 && (range * dx <= rangemax || range * dy <= rangemax)
                    
                    range = range+1;
                    logic = ( XX{kk} >= ( VMx(ii,jj) - range * dx) ) .*...
                        ( XX{kk} <= ( VMx(ii,jj) + range * dx) ) .*...
                        ( YY{kk} >= ( VMy(ii,jj) - range * dy) ) .*...
                        ( YY{kk} <= ( VMy(ii,jj) + range * dy) );
                    
                end
                
                % Compute harmonic averaging of overlapping cells from tile
                % kk. Scaled by the inverse distance to the centroid of
                % merging cells.
                if sum(sum(logic))~=0
                    r = sqrt( (XX{kk}(logic==1) -  VMx(ii,jj)).^2 +...
                        (YY{kk}(logic==1) -  VMy(ii,jj)).^2 + 1e-1 );
                    
                    r = spdiags(1./r,0,length(r),length(r));
                    
                    wght = (r * W{kk}(:,logic==1)')';
                    
                    numer = numer + sum( model{kk}(:,logic==1) .*...
                                            wght , 2);
                    denom = denom + sum(  wght , 2);
                    
                end
                
            end
            
            if numer~=0
                
               VM_model(:,ii,jj) = numer ./ denom;
               
            end

             d_iter = floor(count/mcell*100);
            if  d_iter > progress

                fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end
            count = count +1;
            
        end
        
end

            
% topo = load('topo_model.txt');

%% Write back to file in UBC format
% for ww=(size(model,2))
    
    
    VM_model=reshape(VM_model,VM_mesh(1,3)*VM_mesh(1,2)*VM_mesh(1,1),1);
%     VM_model(topo==0)=ndv;
%     out=model{ww};
    filename='Merged_tiles_VM_v2.dat';
    save([work_dir '\' filename],'-ascii','VM_model')
% end


    