function [mesh]=get_UBCmesh2D(meshfile)
% Read UBC mesh file and extract parameters
% Works for the condenced version (20 * 3) --> [20 20 20] 
fid=fopen(meshfile,'rt');


% Go through the log file and extract data and the last achieved misfit
         	
line=fgets(fid);
    
% First line: number of cells X 
mesh(1,1) = str2num(line);

line=fgets(fid);
numstring = str2num(line);    
mesh(2,1) =  numstring(1);
dx = abs(numstring(1) - numstring(2))/numstring(3);
xin = numstring(2);
dX = ones(1,numstring(3)) * dx;
xin = numstring(2);% Second line: origin coordinate (X,Y,Z)
    for jj= 1 : mesh(1,1)-1
        line=fgets(fid);
        
        numstring = str2num(line);
        
        dx = abs( numstring(1) - xin ) / numstring(2);
        
        dX = [dX ones(1,numstring(2)) * dx];
        
        xin = numstring(1);
        
    end
    
mesh(3,1:length(dX)) = dX;

line=fgets(fid);

clear dX

% Keep reading file for z-axis
line=fgets(fid);
mesh(1,2) = str2num(line);
line=fgets(fid);
numstring = str2num(line);
mesh(2,2) =  numstring(1);
dx = abs(numstring(1) - numstring(2))/numstring(3);
xin = numstring(2);

dX = ones(1,numstring(3)) * dx;

for jj= 1 : mesh(1,2)-1
        line=fgets(fid);
        
        numstring = str2num(line);
        
        
        dx = abs( numstring(1) - xin ) / numstring(2);
        
        dX = [dX ones(1,numstring(2)) * dx];
        
        xin = numstring(1);
    end
    
    mesh(4,1:length(dX)) = dX;
    % Other lines for the dX, dY ,dZ
    
fclose(fid);    