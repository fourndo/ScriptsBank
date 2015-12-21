function [cntm] = MAG3C_RTC_cntmass(OctLev)
% Function [cntm] = MAG3C_RTC_cntmass(OctLev)
% Takes an octree level matrix and compute the center of
% of mass for the highest octree level.
% Center of mass is computed from the nodes location assuming density of 1:
% [Rx,Ry,Rz,V] = 1/N * sum( r(N) * v(N) ), V = sum(dx(N) * dy(N) * dz(N))
% where: N is the number of nodes, r(N) are the [x,y,z] coordinates of each nodes
% and v(N) is the fraction of volume associated to the ith corner.
% 
%
% INPUT
% OctLev: Octree level information as supplied by sub-function MAG3C_RTC_Octlev.m
% 3D array for the coordinate location of all active octree cells
% Dim-1: Topocell ID
% Dim-2: Octree Level
% Dim-3: [ z11 x11 y11 z12 x12 y12 ... zn1 xn1 yn1 zn2 xn2 yn2 ]
%
% The thrid dimension contains the max and min coordinates for every active cell
% at the highest octree level associated to the ith topocell.
%
% OUTPUT
% cntm: Center of mass information computed from the corners of every active cells
% [Rx;Ry;Rz;V] size(4-by-ncell)

ntcell   = size(OctLev,1);
toplev  = size(OctLev,2);

% Pre-allocate space
cntm = zeros(ntcell,4);

% Cycle through tcells and compute center of mass
for ii = 1 : ntcell
    
    	% Initialize variables
    	Rx = 0;	Ry = 0;	Rz = 0;	V  = 0;
    	
    	nnodes = nnz(OctLev(ii,toplev,:));
    	
	for jj = 1 : 6 : nnodes
	
		Rz = Rz + 4 * OctLev(ii,toplev,jj) + 4 * OctLev(ii,toplev,jj+3);
		Rx = Rx + 4 * OctLev(ii,toplev,jj+1) + 4 * OctLev(ii,toplev,jj+4);
		Ry = Ry + 4 * OctLev(ii,toplev,jj+2) + 4 * OctLev(ii,toplev,jj+5);
		V  = V + abs( ( OctLev(ii,toplev,jj+3) - OctLev(ii,toplev,jj) ) *...
		( OctLev(ii,toplev,jj+4) - OctLev(ii,toplev,jj+1) ) *...
		( OctLev(ii,toplev,jj+5) - OctLev(ii,toplev,jj+2) ) );
		
	end
    	
    	% Multiply by 4 since 8 nodes
    	% Divide by 3 since (x,y,z)/node    	
    	cntm(ii,1) = Rx / (4 * nnodes / 3);
    	cntm(ii,2) = Ry / (4 * nnodes / 3);
    	cntm(ii,3) = Rz / (4 * nnodes / 3);
    	cntm(ii,4) = V;
    	
end
    
