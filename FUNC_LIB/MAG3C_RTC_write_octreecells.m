function MAG3C_RTC_write_octreecells(work_dir,filename,OctLev)
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

%% FOR DEV ONLY%%
fid = fopen([work_dir '\' filename],'w');

%% \\\\\\\\\\\\\\\\%%
% SCRIPT STARTS HERE%
%%\\\\\\\\\\\\\\\\\\%
ntcell   = size(OctLev,1);

% Cycle through tcells and compute center of mass
for jj = 1 : ntcell
    
    	% Initialize variables
    	tz1 = OctLev(jj,end,1:6:end);
        tz2 = OctLev(jj,end,4:6:end);
        tx1 = OctLev(jj,end,2:6:end);
        tx2 = OctLev(jj,end,5:6:end);
        ty1 = OctLev(jj,end,3:6:end);
        ty2 = OctLev(jj,end,6:6:end);
        
        tZm = (tz2+tz1)/2;
        tXm = (tx2+tx1)/2;
        tYm = (ty2+ty1)/2;
    	
        index = (tZm(:)~=0&tXm(:)~=0&tYm(:)~=0);
        XYZ = [tXm(index) tYm(index) tZm(index)];
        for ii = 1 : size(XYZ,1)
            fprintf(fid,'%12.5f %12.5f %12.5f\n',XYZ(ii,:));
        end
    	
end
fclose(fid);

