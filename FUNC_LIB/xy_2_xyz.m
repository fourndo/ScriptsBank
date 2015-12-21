function [X,Y,Z] = xy_2_xyz(X,Y, dZ,topofile,n)
% Load X,Y file and assign Z values from topography using a weighted
% averaging on n nearest neigbours. 
% Zavg = sum(1/r1 * z1 + ... + 1/rn * zn) / sum(1/r1 + ... 1/rn)
% Input:
% XYfile    : [X,Y,[O]] coordinate to be referenced. The last column of
% the input file must be the offset [O] specified by the user.
% If omitted = directly on topography. +O above topography
%
% topofile  : Topofile in UBC.topo format
%
% n         : Number of neighbouring nodes to average


%% FOR DEV ONLY
% clear all
% 
% close all
% 
% work_dir = 'C:\Users\Thomasc\Dropbox\EOSC 556b (2013) EM\TKCdata\DIGHEM\Codes\Test\';
% 
% XYfile = [work_dir 'DIGHEMtest.obs'];
% topofile = [work_dir 'CDED_076c05_NAD27.topo'];
% n = 3;




%% SCRIPT STARTS HERE

% Load XY file
ndata = size(X,1);

if size(Y,1) ~= ndata
    
    fprintf('Number of X and Y locations are different, please revise\n')
    
    return
    
end

% Create array where Z values will be added
Z = zeros(ndata,1); 


% Load topography - UBC style
fprintf('Loading Topofile ... please standby\n')

% Load header
fid = fopen(topofile,'r');
line = fgets(fid);

nnodes = str2double(line);
topo = zeros(nnodes,3);
for ii = 1 : nnodes
    
   
    topo(ii,:) = str2num(fgets(fid));
    
end

fclose(fid)
fprintf('Completed!\n')

% Start geo-referencing
fprintf('Start geo-referencing z values ... \n')
progress = 0;
tic
for ii = 1 : ndata
    
    
    % First compute distance to all points in topofile
    r = ( ( X(ii) - topo(:,1) ).^2 +...
        ( Y(ii) - topo(:,2) ).^2 ) .^ 0.5;
    
    % Compute inverse distance + small number to avoid singularity
    r = 1./(r+1e-2);
    
    % Sort points by distance and keep n smallest
    [rsort,index] = sort(r,'descend');
    
    % Compute weighted average
    numer = 0;
    denom = 0;
    
    for jj = 1 : n
        
        numer = numer + topo(index==jj,3)*rsort(jj);
        denom = denom + rsort(jj);
        
    end
    
    Z(ii) = numer/denom + dZ(ii);
        
    % Prompt iteration and time message
    d_iter = floor(ii/ndata*20);
   
    if  d_iter > progress
        
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
  
        progress = d_iter;        
        tic
         
    end
    
end