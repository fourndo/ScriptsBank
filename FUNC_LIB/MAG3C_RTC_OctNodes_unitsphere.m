function [OctNodes] = MAG3C_RTC_OctNodes_unitsphere(Xn,Yn,Zn,cellID,OctLev)
% Function [OctNodes] = MAG3C_RTC_OctNodes(Xn,Yn,Zn,topocell,ztopo_n,OctLev)
% Takes a tensor mesh defined by the node locations and topographic
% cells. Breaks down the topographic cells for various Octree level and
% returns the nodes of octree active cells
%
% INPUT
% Xn, Yn, Zn: 3D array of node location (X, Y, Z) coordinates
% topocell: ID of the cells in mesh intersected by topography
% ztopo_n : Z-value for topographyon all nodes of the mesh
% OctLev : Octree levels to be computed
%           Level   Octree cells
%           0            1
%           1       	 8 
%           2       	 64 
%           3       	 512 
%
% SUB_FUNCTIONS CALL
% topocheck.m
%
% OUTPUT
% OctNodes: 3D array for the coordinate location of all active octree cells
% Dim-1: Topocell ID
% Dim-2: Octree Level
% Dim-3: [ z11 x11 y11 z12 x12 y12 ... zn1 xn1 yn1 zn2 xn2 yn2 ]

% Create empty array for index
nz = size(Xn,1)-1 ;
nx = size(Xn,2)-1 ;
ny = size(Xn,3)-1 ;

nOctLev = length(OctLev);
mcell   = length(cellID);

% Pre-allocate space
OctNodes = zeros( mcell , nOctLev , 8^OctLev(end) * 6 );

for ii = 1 : mcell
    
    % Get [i,j,k] for current topocell
    [k,i,j] = ind2sub([nz nx ny],cellID(ii));
    
    % Extract node location for current topocell
    xn = [Xn(k,i,j) Xn(k,i+1,j)];
    yn = [Yn(k,i,j) Yn(k,i,j+1)];
    zn = [Zn(k,i,j) Zn(k+1,i,j)];
    
    % Extract topo on nodes for current cell
%     topo = [xn(1) yn(1) ztopo_n(k,i,j);...
%             xn(1) yn(2) ztopo_n(k,i,j+1);...
%             xn(2) yn(1) ztopo_n(k,i+1,j);...
%             xn(2) yn(2) ztopo_n(k,i+1,j+1) ];
    
    for jj = 1 : nOctLev
        
        edgecell = 2^OctLev(jj);% Number of cells on edge
        
        % Generate new array of nodes within current topocell
        if OctLev(jj)==0
            
            oxn = [xn(1) xn(2)];
            oyn = [yn(1) yn(2)];
            ozn = [zn(1) zn(2)];
            
        else
            
            oxn = linspace( xn(1) , xn(2) , edgecell+1 );
            oyn = linspace( yn(1) , yn(2) , edgecell+1 );
            ozn = linspace( zn(1) , zn(2) , edgecell+1 );
            
        end
        
        [OZn,OXn,OYn] = ndgrid(ozn,oxn,oyn);
        
        % Call topocheck to get active octree cells within
        [keep,~] = topocheck_unitsphere(OXn,OYn,OZn);
        
        
        % Save active nodes
        [subk,subi,subj]=ind2sub([edgecell edgecell edgecell],find(keep==1));
        
        if isempty(subk)==0
            
        % Create temporary array for active node (6 nodes per active cell)
        z1n = OZn(subk,1,1); 
        z2n = OZn(subk+1,1,1); 
        x1n = OXn(1,subi,1); 
        x2n = OXn(1,subi+1,1); 
        y1n = OYn(1,1,subj); 
        y2n = OYn(1,1,subj+1); 
        
        temp = [z1n(:) x1n(:) y1n(:) z2n(:) x2n(:) y2n(:)]';
        temp = temp(:);
        
        % Store nodes
        OctNodes(ii,jj,1:length(temp)) = temp+1e-20;
        
        end
        
    end
    
end
    
    