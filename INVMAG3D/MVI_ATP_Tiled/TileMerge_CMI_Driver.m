% function TileMerge_AMI(work_dir,out_meshfile,tilefile,rangemax,ndv,dsep,flag1,arg1,flag2,arg2)
% function TileMerge(work_dir,out_meshfile,rangemax,ndv,dsep,flag1,arg1,flag2,arg2)
% Merge inverted models from seperate tiles onto a large mesh.
% Tiles can be overlapping in any orientation horizontally, but must be
% overlapping vertically. Use a dissolving weights towards the outside of
% the tiles to compute a weighted average.
%
% INPUT: 
% work_dir: Working directory
% out_meshfile : Final mesh to be interpolated onto
% rangemax : Maximum interpolation distance
% ndv : no-data-value 
% dsep : directory seperator (Windows '\') (Linux '/')
% flag1 : either 'no_topo' || 'topofile' || 'nullcell'
% arg1 : either [] || 'topofile.dat' || 'nullcell.dat'
% flag2 : either 'rem_pad' || 'default'
% arg2 : either [padx pady] || []
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

addpath '.\func_lib'
%% INPUTS
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\CMI_Tiles\ALL_Tiles';
out_meshfile = '..\..\Mesh_50m.msh';
tilefile = '..\..\Tiles_50m_Tight.dat';
rangemax = 200;
ndv = -100;
dsep = '\';
flag1 = 'topofile';
arg1='..\..\CDED_76c_d_RandSUB.dat';
flag2 = '';
arg2 = '';

%% \\\\\\\\\\\\\\\\\ %%
% SCRIPT STARTS HERE  %
%%\\\\\\\\\\\\\\\\\\ %%
%% Generate list of input models from directory

T = load([work_dir dsep tilefile]);
ntiles = size(T,1);

filein = ls(work_dir);
filenum = [];
for ii = 1 : size(filein,1)
    
    temp1 = regexp(filein(ii,:),'\d*','match');
    temp2 = regexp(filein(ii,:),'[.]','split');
    if ~isempty(temp1) && strcmp(strtrim(temp2{2}),'fld')
        
        filenum = [filenum str2num(temp1{:})];
        
    end
end

tiles = unique(filenum);
ntiles = length(tiles);

% Input files
for ii = 1:ntiles
    
        
        in_file = ['Tile' num2str(tiles(ii))];
        
        fld_file{ii}=[work_dir dsep in_file '_MVI.fld'];
        rem_file{ii}=[work_dir dsep in_file '_MVI.rem'];
        ind_file{ii}=[work_dir dsep in_file '_MVI.ind'];
        
        msh_file{ii}=[work_dir dsep in_file '.msh'];
        

    
end

% Load mesh files for all the tiles
% model has dimension [ntiles-by-5]
% where second dimension is [rem ind mx my mz]
for ii = 1 : ntiles
    
    fprintf('Loading Tile %s\n',fld_file{ii});
    fld_model = importdata(fld_file{ii});
    model{ii,1} = importdata(rem_file{ii});
    model{ii,2} = importdata(ind_file{ii});
    
    model{ii,3} = fld_model(:,1);
    model{ii,4} = fld_model(:,2);
    model{ii,5} = fld_model(:,3);

%     model{ii}(model{ii}==-1) = ndv;
    
    mesh{ii} = get_UBC_mesh(msh_file{ii});
    
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

write_UBC_mesh(work_dir,'VOXMERGE.msh',dsep,XO_out,YO_out,ZO_out,out_mesh(3,1:out_mesh(1,1)),out_mesh(4,1:out_mesh(1,2)),out_mesh(5,1:out_mesh(1,3)))
end

for mm = 1 : 5
    
    out_model{mm} = ones(out_mesh(1,3),out_mesh(1,1),out_mesh(1,2))*ndv;
    
end

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
    
    for mm = 1 : 5
        
        model{ww,mm}=reshape(model{ww,mm},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));

        switch flag2

            case 'rem_pad'
                % Temporary line to remove padding cells
                model{ww,mm}(:,[1:arg2(1) (end-arg2(1)+1):end],:) = ndv;
                model{ww,mm}(:,:,[1:arg2(2) (end-arg2(2)+1):end]) = ndv;

        end

        model{ww,mm}=model{ww,mm}-ndv;

        for i=1:mesh{ww}(1,1)

            for j=1:mesh{ww}(1,2)

                columnxy=model{ww,mm}(:,i,j);
                ground = (mesh{ww}(1,3)-nnz(columnxy));

                if ground ~= mesh{ww}(1,3)
                    columnxy(1:ground)=columnxy(ground+1)*ones(ground,1);
                end

                model{ww,mm}(:,i,j)=columnxy;


            end
        end
        model{ww,mm}=model{ww,mm}+ndv;
        
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
    leoutax = 0;
    % Take the top slice of the model
    nullit = reshape(model{ww,1}(1,:,:),mesh{ww}(1,1),mesh{ww}(1,2));
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
            
            for mm = 1 : 5
                numer{mm} = zeros(size(out_model{mm},1),1);
            end
            
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
                    for mm = 1 : 5
                        numer{mm} = numer{mm} + sum( Z_interp{ tile_in(kk) } * (model{ tile_in(kk) , mm }(:,i,j) .* wght) , 2);
                    end
                    denom = denom + sum( Z_interp{ tile_in(kk) } * wght , 2);

                    
                end
                
            end
            
            
            if sum(denom)~=0
            
                for mm = 1 : 5
                    
                   out_model{mm}(:,ii,jj) =  (numer{mm} ./ denom);
                   
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

%% Compute topocheck on final mesh to remove aircells
switch flag1
    
    case 'topofile'
        
        % Load topography file (UBC format)
        topo = read_UBC_topo([work_dir dsep arg1]);      
        
        % Get nodal discretization of final mesh
        xn = [out_mesh(2,1) out_mesh(2,1) + cumsum( out_mesh(3,1:out_mesh(1,1)) ) ];
        yn = [out_mesh(2,2) out_mesh(2,2) + cumsum( out_mesh(4,1:out_mesh(1,2)) ) ];
        zn = [out_mesh(2,3) out_mesh(2,3) - cumsum( out_mesh(5,1:out_mesh(1,3)) ) ];
        
        [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
                
        [nullcell,topocell,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
        
        save([work_dir dsep 'nullcell.dat'],'-ascii','nullcell')
        fprintf('Nullcell has been saved to: \n %s \n for future use\n',work_dir);
        
        for mm = 1 : 2
            
            % Assign air cells
            out_model{mm}(nullcell==0) = ndv;    
        end
        
    case 'nullfile'

        nullcell = load([work_dir dsep arg1]);
        for mm = 1 : 2
            
            % Assign air cells
            out_model{mm}(nullcell==0) = ndv; 
            
        end
        
    otherwise
        
           
end
%% Write back to file in UBC format 

for mm = 1 : 5
    
    out_model{mm}=out_model{mm}(:);
    out_model{mm}(isnan(out_model{mm})) = ndv;
    
end    
    
fld_model = [out_model{3} out_model{4} out_model{5}];
fld_model(fld_model==ndv) = 0;

rem_model = out_model{1};
ind_model = out_model{2};
amp_model = ( sum(fld_model.^2,2) ).^0.5;
amp_model(nullcell==0) = ndv;

save([work_dir dsep 'Merged_M_model.fld'],'-ascii','fld_model');
save([work_dir dsep 'Merged_Rem_model.rem'],'-ascii','rem_model');
save([work_dir dsep 'Merged_Ind_model.ind'],'-ascii','ind_model');
save([work_dir dsep 'Merged_Amp_model.amp'],'-ascii','amp_model');

fprintf('Program TileMerge succesfully ended\n')
% fprintf('Merged model saved: %s %s\n',work_dir,filename);



    