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

clear all
close all
root_dir = pwd;
addpath(root_dir)

%% INPUTS
% Load files
work_dir ='C:\Projects\Research\VoxMerge\center_2it';

cd(work_dir);

model{1}=importdata([work_dir '\tile1_inv06.con']);
meshfile{1}=[work_dir '\tile1_inv06.msh'];

model{2}=importdata([work_dir '\tile2_inv03.con']);
meshfile{2}=[work_dir '\tile2_inv02.msh'];

model{3}=importdata([work_dir '\tile3_inv04.con']);
meshfile{3}=[work_dir '\tile3_inv04.msh'];

model{4}=importdata([work_dir '\tile4_inv03.con']);
meshfile{4}=[work_dir '\tile4_inv03.msh'];

model{5}=importdata([work_dir '\center_2it.con']);
meshfile{5}=[work_dir '\mesh_centerb.msh'];

%% Specify final mesh file 

m_final=size(model,2)+1;
meshfile{m_final}=[work_dir '\merge_mesh_ALL.msh'];

% Load topo model to substract air cells
% topo = load('topo_model.txt');

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=1e-008;

%% Load mesh files for all the tiles
for ii = 1:size(meshfile,2);
    
    mesh{ii} = get_UBCmesh(meshfile{ii});
    
end

%% Create final mesh
% Find the extent of all the tiles
model{m_final}=ones(mesh{m_final}(1,3),mesh{m_final}(1,1),mesh{m_final}(1,2))*ndv;



%Create final mesh
    Ax0=mesh{m_final}(2,1);
    Ay0=mesh{m_final}(2,2);
    nx = mesh{m_final}(1,1);
    ny = mesh{m_final}(1,2);
    dxf=mesh{m_final}(3,1:nx);
    dyf=mesh{m_final}(4,1:ny); 
    [finalx,finaly]=meshgrid((Ax0+cumsum(dxf)-dxf/2),...
                            (Ay0+cumsum(dyf)-dyf/2));

%% Propagate data upward
% This steps deals with conflicting topo between adjacente tiles.
% The finale mesh must be edited using TOPOCHECK to substrat the real topo

for ww=1:(size(model,2))-1
    model{ww}=reshape(model{ww},mesh{ww}(1,3),mesh{ww}(1,1),mesh{ww}(1,2));
    model{ww}=model{ww}-ndv;

    for i=1:mesh{ww}(1,1)
        for j=1:mesh{ww}(1,2)
            temp=model{ww}(:,i,j);
            ground = (mesh{ww}(1,3)-nnz(temp));
            
            if ground ~= mesh{ww}(1,3)
                temp(1:ground)=temp(ground+1)*ones(ground,1);
            end
            
            model{ww}(:,i,j)=temp;

        end
    end
    model{ww}=model{ww}+ndv;
end


    
 
%% Beggining of merging computations
% Iterate over all zones as center node
for kk=1:(size(model,2))-1
    
    
    % Iterate over all neighbours
    for jj=1:(size(model,2))-1 
    if jj==kk
        continue
    end
        
    % Get coordinate grids
    Ax0=mesh{kk}(2,1);
    Ay0=mesh{kk}(2,2);
    nx = mesh{kk}(1,1);
    ny = mesh{kk}(1,2);
    dx=mesh{kk}(3,1:nx);
    dy=mesh{kk}(4,1:ny); 
    [Ax,Ay]=meshgrid((Ax0+cumsum(dx)-dx/2),...
                            (Ay0+cumsum(dy)-dy/2));
    Bx0=mesh{jj}(2,1);
    By0=mesh{jj}(2,2);  
    nx = mesh{jj}(1,1);
    ny = mesh{jj}(1,2);
    dx=mesh{jj}(3,1:nx);
    dy=mesh{jj}(4,1:ny);                    
    [Bx,By]=meshgrid((Bx0+cumsum(dx)-dx/2),...
                            (By0+cumsum(dy)-dy/2));

    
    % Find the overlapping cells
    overlapX_A=zeros(1,size(Ax,2));
    overlapX_B=zeros(1,size(Bx,2));
   
    overlapY_A=zeros(1,size(Ay,1));
    overlapY_B=zeros(1,size(By,1));
    
    % Find the overlaps in X
    for oo=1:size(Ax,2)
        for pp=1:size(Bx,2)
            
            if abs(Ax(1,oo)-Bx(1,pp)) < min(dx)/2
                overlapX_A(oo)=1;
                overlapX_B(pp)=1;
            end
            
        end
    end
    
    %Compute the overlap in X
    mullet=sum(overlapX_A);
    mullet=min([mullet-1,floor(size(Ax,2)/2),floor(size(Bx,2)/2)]);
    
    if mullet==0;
        sprintf('tiles are not overlapping')
        continue
    end
    
    for oo=1:size(Ay,1)
        for pp=1:size(By,1)
            
            if abs(Ay(oo,1)-By(pp,1)) < min(dy)/2
                overlapY_A(oo)=1;
                overlapY_B(pp)=1;
            end
            
        end
    end
    

    %Compute the length to trim in Y    
    bangs=sum(overlapY_A);
    bangs=min([bangs-1,floor(size(Ay,1)/2),floor(size(By,1)/2)]);
        
    % Find centroid
    centroid{1}=floor([size(Ax,2)/2 size(Ay,1)/2]);
    centroid{2}=floor([size(Bx,2)/2 size(By,1)/2]);
    
%     azmth=centroid{1}-centroid{2};
    rampupX_A=ones(1,size(Ax,2));
    rampupX_B=ones(1,size(Bx,2));
    rampupY_A=ones(1,size(Ay,1));
    rampupY_B=ones(1,size(By,1));
    rampX=(1e-6:1:mullet)/mullet;
    rampY=(1e-6:1:bangs)/bangs;
    
    rampupX_A(1:mullet)=rampX;
    rampupX_A(end-mullet+1:1:end)=flipdim(rampX,2);
    rampupX_A=-0.5*cos(-rampupX_A*pi)+0.5;
    
    rampupX_B(1:mullet)=rampX;
    rampupX_B(end-mullet+1:1:end)=flipdim(rampX,2);
    rampupX_B=-0.5*cos(-rampupX_B*pi)+0.5;
    
    rampupY_A(1:bangs)=rampY;
    rampupY_A(end-bangs+1:1:end)=flipdim(rampY,2);
    rampupY_A=-0.5*cos(-rampupY_A*pi)+0.5;
    
    rampupY_B(1:bangs)=rampY;
    rampupY_B(end-bangs+1:1:end)=flipdim(rampY,2);
    rampupY_B=-0.5*cos(-rampupY_B*pi)+0.5;
    
    
        [Z_A,X_A,Y_A]=ndgrid(1:mesh{kk}(1,3),rampupX_A,rampupY_A);
        [Z_B,X_B,Y_B]=ndgrid(1:mesh{jj}(1,3),rampupX_B,rampupY_B);
        
        X_A=X_A.*Y_A;
        Y_A=X_A;
        
        X_B=X_B.*Y_B;
        Y_B=X_B;
        
    %Flip the Y-overlap because UBC counts from the south
%     overlapY_A=flipdim(overlapY_A,2);
%     overlapY_B=flipdim(overlapY_B,2);
    
    model{kk}(:,overlapX_A==1,overlapY_A==1)=...
        (model{kk}(:,overlapX_A==1,overlapY_A==1).*X_A(:,overlapX_A==1,overlapY_A==1)+...
        model{kk}(:,overlapX_A==1,overlapY_A==1).*Y_A(:,overlapX_A==1,overlapY_A==1)+...
        model{jj}(:,overlapX_B==1,overlapY_B==1).*X_B(:,overlapX_B==1,overlapY_B==1)+...
        model{jj}(:,overlapX_B==1,overlapY_B==1).*Y_B(:,overlapX_B==1,overlapY_B==1))./...
        (X_A(:,overlapX_A==1,overlapY_A==1)+Y_A(:,overlapX_A==1,overlapY_A==1)+...
        X_B(:,overlapX_B==1,overlapY_B==1)+Y_B(:,overlapX_B==1,overlapY_B==1));
    
    model{jj}(:,overlapX_B==1,overlapY_B==1)=model{kk}(:,overlapX_A==1,overlapY_A==1);
        

    
    %Update final matrix
    % Find the overlapping cells
    overlapX_A=zeros(1,size(finalx,2));
    overlapY_A=zeros(1,size(finaly,1));  
    
      for oo=1:size(Ax,2)
        for pp=1:size(finalx,2)
            
            if abs(Ax(1,oo) - finalx(1,pp)) < min(dxf)/2;
                overlapX_A(pp)=1;
            end
            
        end
      end
    
    for oo=1:size(Ay,1)
        for pp=1:size(finaly,1)
            
            if abs(Ay(oo,1) - finaly(pp,1)) < min(dyf)/2;
                overlapY_A(pp)=1;
            end
            
        end
    end
    
%     overlapY_A=flipdim(overlapY_A,2);
    model{m_final}(:,overlapX_A==1,overlapY_A==1)=model{kk};

    % Reset the overlapping cells for second model
    overlapX_A=overlapX_A*0;
    overlapY_A=overlapY_A*0;
    
      for oo=1:size(Bx,2)
        for pp=1:size(finalx,2)
            
            if abs(Bx(1,oo) - finalx(1,pp)) < min(dxf)/2;
                overlapX_A(pp)=1;
            end
            
        end
      end
    
    for oo=1:size(By,1)
        for pp=1:size(finaly,1)
            
            if abs(By(oo,1) - finaly(pp,1)) < min(dyf)/2;
                overlapY_A(pp)=1;
            end
            
        end
    end 
%     overlapY_A=flipdim(overlapY_A,2);
    model{m_final}(:,overlapX_A==1,overlapY_A==1)=model{jj};
    end
    
    clear overlapX_A
    clear overlapY_A

    
end


%% Write back to file in UBC format
for ww=(size(model,2))
    
    
    model{ww}=reshape(model{ww},mesh{ww}(1,3)*mesh{ww}(1,2)*mesh{ww}(1,1),1);
%     model{ww}(topo==0)=ndv;
    out=model{ww};
    filename='Merged_tiles_VM_v1.dat';
    save([work_dir '\' filename],'-ascii','out')
end


    