function [nullcell] = topo_2_model(work_dir,meshfile,topofile)
% topo_2_model(meshfile,topofile)
% Create active cell matrix from discretized topography
% Inputs//
% meshfile: UBC format
% topofile: UBC format
%
% Output//
% nullcell: 1D vector, in UBC format, of binary values 1 (active) or
% 0(inactive) depending if below or above topography

%% FOR DEV
% clear all
% close all
% 
% work_dir = 'C:\Local Cloud\Private\dominiquef\Projects\Research\Modelling\Inversion\Synthetic\Block_40m';
% 
% obsfile = 'Synthetic_INDU_TMI_2pc_noise.obs';
% meshfile = 'Mesh_40m.msh';
% model_sus = 'cube.sus';
% % model_azm = 'Horseshoe_azm.sus';
% % model_dip = 'Horseshoe_dip.sus';
% topofile = 'Gaussian.topo';

%% SCRIPT START HERE
% Load mesh size
[meshfile]=get_UBC_mesh([work_dir '\' meshfile]);
nx = meshfile(1,1); %size(X,1);    %number of cell in X
ny = meshfile(1,2); %size(X,2);    %number of cell in Y
nz = meshfile(1,3); %size(X,3);    %number of cell in Z

% mcell = nx*ny*nz;

% Cell size array
dx = meshfile(3,1:nx);
dy = meshfile(4,1:ny);
dz = meshfile(5,1:nz);

% Corner grid
x0 = meshfile(2,1);
y0 = meshfile(2,2);
z0 = meshfile(2,3);

% Create 3D cell node array
xn = [x0 x0 + cumsum(dx)];
yn = [y0 y0 + cumsum(dy)];
zn = [z0 (z0 - cumsum(dz))]; % Compute top of cell only

% Load topography - UBC style
fprintf('Loading Topofile ... please standby\n')

if isempty(topofile)==1
    
    fprintf('Topo file name is empty\n')
    fprintf('Nullcell model of all ones created\n')
    
    % Pre-allocate space
    nullcell = ones(nz,nx,ny);
    
else
    
    % Load header
    fid = fopen([work_dir '\' topofile],'r');
    line = fgets(fid);

    nnodes = str2double(line);
    topo = zeros(nnodes,3);
    for ii = 1 : nnodes


        topo(ii,:) = str2num(fgets(fid));

    end

    fclose(fid);
    fprintf('Completed!\n')

    % Start looking for elevation of topography at each cell
    fprintf('Start computing topography over the mesh ... \n')

    % Pre-allocate space
    nullcell = ones(nz,nx,ny);

    progress = 0;
    tic
    for jj = 1 : ny

        for ii = 1 : nx
            r = zeros(nnodes,4);
            
            % First compute distance of cell to all points in topofile (2D)
            r(:,1) = ( ( xn(ii) - topo(:,1) ).^2 +...
                ( yn(jj) - topo(:,2) ).^2 ) .^ 0.5;

            r(:,2) = ( ( xn(ii+1) - topo(:,1) ).^2 +...
                ( yn(jj) - topo(:,2) ).^2 ) .^ 0.5;
            
            r(:,3) = ( ( xn(ii) - topo(:,1) ).^2 +...
                ( yn(jj+1) - topo(:,2) ).^2 ) .^ 0.5;
            
            r(:,4) = ( ( xn(ii+1) - topo(:,1) ).^2 +...
                ( yn(jj+1) - topo(:,2) ).^2 ) .^ 0.5;
            
            % Compute inverse distance + small number to avoid singularity
            r = 1./(r+1e-8);

            % Sort points by distance and keep n smallest
            [rsort,index] = sort(r,'descend');

            % Compute weighted average
            numer = zeros(1,4);
            denom = zeros(1,4);

            for kk = 1 : 3

                for rr = 1 : 4
                    
                    numer(rr) = numer(rr) + topo(index(kk,rr),3)*rsort(kk,rr);
                    denom(rr) = denom(rr) + rsort(kk,rr);
                    
                end

            end

            ztopo = min(numer./denom);

            % Change values along the column for 0 until reaches topography
            count = 1;
            while zn(count) > ztopo

                nullcell(count,ii,jj) = 0;
                count = count +1;

            end

            % Prompt iteration and time message
            d_iter = floor(ii*jj/(nx*ny)*20);

            if  d_iter > progress

                fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)

                progress = d_iter;        
                tic

            end

        end

    end

end

    nullcell = reshape(nullcell,nx*ny*nz,1);

