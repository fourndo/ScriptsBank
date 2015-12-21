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
% Revision: 3
% Last update: December 29th, 2013
%
% CODE IN DEVELOPMENT

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

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-99999;

% Maximum search radius (x/2*cell)
rangemax = 1;

%% Load mesh files for all the tiles
ntiles = size(meshfile,2);
mesh{ntiles,1}= [];
X0 = zeros(ntiles,1);
Y0 = zeros(ntiles,1);
Z0 = zeros(ntiles,1);
Xmax = zeros(ntiles,1);
Ymax = zeros(ntiles,1);
Zmax = zeros(ntiles,1);
dxmin = zeros(ntiles,1);
dymin = zeros(ntiles,1);
dzmin = zeros(ntiles,1);

for ii = 1 : ntiles
    
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

% Make sure that tiles are all coplanar, otherwise switch to full 3D merge.
argin = 'coplanar';
for ii = 1 : ntiles-1
    
    if Z0(ii)~=Z0(ii+1) || Zmax(ii)~=Zmax(ii+1)
        
        argin = 'staggered';
        break
        
    end
    
end

switch argin 
    
    case 'staggered'
    
        fprintf('Program has detected that tiles are %s\n',argin);
        fprintf('Merging will be done in 3D... may take several minutes');
    
    otherwise
    
        fprintf('Program has detected that tiles are %s\n',argin);
        fprintf('Merging coefficients in 2D only');
        
end

%% Create final voxel merged mesh (VM)
% If mesh file is provided
if isempty(VM_meshfile)==0
    
    VM_mesh = get_UBC_mesh(VM_meshfile);
    
else

% User input if mesh not provided
dx = input('Final mesh file not specified. Enter cell size for x (m)\n');
dy = input('Final mesh file not specified. Enter cell size for y (m)\n');
dz = input('Final mesh file not specified. Enter cell size for z (m)\n');

% If no meshfile provided then the final mesh is created from extent.
% Find the extent of all the tiles
VM_mesh(2,1) = min(X0);
VM_mesh(2,2) = min(Y0);
VM_mesh(2,3) = max(Z0);

Xmax_VM = max(Xmax);
Ymax_VM = max(Ymax);
Zmax_VM = min(Zmax);

VM_mesh(1,1) = ceil( (Xmax_VM - VM_mesh(2,1)) / min(dxmin) );
VM_mesh(1,2) = ceil( (Ymax_VM - VM_mesh(2,2)) / min(dymin) );
VM_mesh(1,3) = ceil( (VM_mesh(2,3) - Zmax_VM) / min(dzmin) );

VM_mesh(3,1:VM_mesh(1,1)) = ones(1,VM_mesh(1,1)) * dx;
VM_mesh(4,1:VM_mesh(1,2)) = ones(1,VM_mesh(1,2)) * dy; 
VM_mesh(5,1:VM_mesh(1,3)) = ones(1,VM_mesh(1,3)) * dz; 

end

VM_model = ones(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2))*ndv;
% [VMdy,VMdx]=meshgrid(VM_mesh(4,1:VM_mesh(1,2)),VM_mesh(3,1:VM_mesh(1,1)));

% Create cell center coordinate matrices for final mesh
switch argin
    
    case 'coplanar'
        [VMy,VMx]=meshgrid((VM_mesh(2,2)+cumsum(VM_mesh(4,1:VM_mesh(1,2)))-VM_mesh(4,1:VM_mesh(1,2))/2),...
                            (VM_mesh(2,1)+cumsum(VM_mesh(3,1:VM_mesh(1,1)))-VM_mesh(3,1:VM_mesh(1,1))/2));

    otherwise
        [VMx,VMz,VMy]=meshgrid((VM_mesh(2,1)+cumsum(VM_mesh(3,1:VM_mesh(1,1)))-VM_mesh(3,1:VM_mesh(1,1))/2),...
            (VM_mesh(2,3)-cumsum(VM_mesh(5,1:VM_mesh(1,3)))+VM_mesh(5,1:VM_mesh(1,3))/2),...
            (VM_mesh(2,2)+cumsum(VM_mesh(4,1:VM_mesh(1,2)))-VM_mesh(4,1:VM_mesh(1,2))/2));

end
    

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
% 3- Normalized 3D cell size.
% 4- Vertical weigths on horizontal slices (1 = [top,bottom] , max =
% center)

weights{ntiles,1}= [];
XX{ntiles,1}= [];
YY{ntiles,1}= [];
ZZ{ntiles,1}= [];

for ww = 1 : ntiles
    
    Wmax{ww} = 0;
    % Take the top slice of the model
    weights{ww} = reshape(model{ww},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));

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
    
    switch argin
    
    case 'coplanar'
        [YY{ww},XX{ww}]=meshgrid((y0+cumsum(dy)-dy/2),...
                            (x0+cumsum(dx)-dx/2));
    otherwise
        [XX{ww},ZZ{ww},YY{ww}] = meshgrid((x0+cumsum(dx)-dx/2),...
            (z0-cumsum(dz)+dz/2),(y0+cumsum(dy)-dy/2));
    end
    
    for ii=1:nx
        
        for jj=1:ny
            
           if  weights{ww}(1,ii,jj)~=ndv
               
               level = min([ii jj (nx-ii+1) (ny-jj+1)]);
               weights{ww}(:,ii,jj) = level;
               
           else
               
               weights{ww}(:,ii,jj) = 0;
               
           end
           
            
        end
        
    end
        
    
    if Wmax{ww} < max(max(max(weights{ww})))
        
        Wmax{ww} = max(max(max(weights{ww})));
        
    end
    
end


for ww = 1 : ntiles
    
   
   
    % Normalize the weights
    weights{ww}(weights{ww}~=0) = ( weights{ww}(weights{ww}~=0) -1 )/...
                                     (Wmax{ww}-1) + 1e-6;
    % Cosine tapper 
    weights{ww} = (-0.5*cos(-weights{ww}*pi)+0.5)  ;  

%     % Cellsize weights
%     nx = mesh{ww}(1,1);
%     ny = mesh{ww}(1,2);
%     nz = mesh{ww}(1,3);
%     [dX,dZ,dY] = meshgrid(mesh{ww}(3,1:nx),mesh{ww}(5,1:nz),mesh{ww}(4,1:ny));
% 
%     weights{ww} = weights{ww} ./ (dX.*dY.*dZ) ;

    switch argin
    
    case 'staggered'
        
        % Apply vertical weight
        zlevels = [1:floor(nz/2) ceil(nz/2):-1:1];
        vweight = repmat(zlevels',[1 nx ny]);
        weights{ww} = vweight.*weights{ww};
        
    end
   
end
    
 
%% Beggining the merging computations
% Iterate over all cells of final mesh and find the closest cells to 
% average. Computes an harmonic weigthed average with inverse distance
% weight. 
% OPTION 1: Assumes that all mesh have the same cells vertically. (Faster)
% OPTION 2: Tiles are not co-planar. Algorithm cycles through all cells 
% vertically. (Much slower)

switch argin
    
    case 'coplanar'

                        
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

                        while sum(sum(sum(logic)))==0 && (range < rangemax)

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
                                (YY{kk}(logic==1) -  VMy(ii,jj)).^2 + 1e-2 );

                            r = spdiags(1./(r+1),0,length(r),length(r));

                            wght = (r * weights{kk}(:,logic==1)')';

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
        
    case 'staggered';
                    
        dx = min(VM_mesh(3,1:VM_mesh(1,1)));
        dy = min(VM_mesh(4,1:VM_mesh(1,2)));
        dz = min(VM_mesh(5,1:VM_mesh(1,3)));
        
        for ii = 1 : VM_mesh(1,1)

                for jj = 1 : VM_mesh(1,2)
                    
                    for kk = 1 : VM_mesh(1,3)

                        numer = 0;
                        denom = 0;


                        % Cycle through all the tiles to find overlapping cells
                        for tt = 1 : ntiles

                            logic = zeros(VM_mesh(1,3),VM_mesh(1,1),VM_mesh(1,2));
                            range = -1;

                            while sum(sum(sum(logic)))==0 && (range < rangemax)

                                range = range+1;
                                logic = ( XX{tt} >= ( VMx(kk,ii,jj) - range * dx) ) .*...
                                    ( XX{tt} <= ( VMx(kk,ii,jj) + range * dx) ) .*...
                                    ( YY{tt} >= ( VMy(kk,ii,jj) - range * dy) ) .*...
                                    ( YY{tt} <= ( VMy(kk,ii,jj) + range * dy) ) .*...
                                    ( ZZ{tt} >= ( VMz(kk,ii,jj) - range * dz) ) .*...
                                    ( ZZ{tt} <= ( VMz(kk,ii,jj) + range * dz) );
                            end

                            % Compute harmonic averaging of overlapping cells from tile
                            % kk. Scaled by the inverse distance to the centroid of
                            % merging cells.
                            if sum(sum(sum(logic)))~=0
                                r = sqrt( (XX{tt}(logic==1) -  VMx(kk,ii,jj)).^2 +...
                                    (YY{tt}(logic==1) -  VMy(kk,ii,jj)).^2 +...
                                    (ZZ{tt}(logic==1) -  VMz(kk,ii,jj)).^2 + 1e-2 );

                                r = spdiags(1./(r+1),0,length(r),length(r));

                                wght = (r * weights{ww}(logic==1));

                                numer = numer + sum(sum( model{tt}(logic==1) .*...
                                                        wght , 2));
                                denom = denom + sum(sum(  wght , 2));


                            end

                        end

                        if numer~=0

                           VM_model(kk,ii,jj) = numer ./ denom;

                        end

                    end
                
                end

        end
        
    otherwise

    fprintf('Input argument not recognize\n')
    fprintf('Should be "coplanar" or "staggered" tiles\n')
    break

end
            
% topo = load('topo_model.txt');

%% Write back to file in UBC format  
VM_model=reshape(VM_model,VM_mesh(1,3)*VM_mesh(1,2)*VM_mesh(1,1),1);
%     VM_model(topo==0)=ndv;
%     out=model{ww};
filename='Merged_tiles_out.dat';
save([work_dir '\' filename],'-ascii','VM_model')

make_UBC_mesh(work_dir,VM_mesh);
    
fprintf('Merge Completed!!\n');
fprintf('You can find merged model here:\n%s\n',[work_dir '\' filename]);


    