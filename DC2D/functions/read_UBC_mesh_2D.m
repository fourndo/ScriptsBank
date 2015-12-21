function [xn,zn] = read_UBC_mesh_2D(meshfile)
% Read UBC mesh 2D and output nodal location

xn = [];
zn = [];

fid = fopen(meshfile,'r');

line = fgets(fid);

% Read the first line
nlines = str2num(line);

% Get the first node and iterate down the file
line = fgets(fid);
var = str2num(line);
xn(1) = var(1);

nd = var(3);
dx = (var(2) - xn(end)) / nd ;

xn = [xn xn(end)+cumsum(ones(1,nd)*dx)];
    
for ii = 2 : nlines
    
    line = fgets(fid);
    var = str2num(line);
    
    nd = var(2);
    dx = (var(1) - xn(end)) / nd ;
    
    xn = [xn xn(end)+cumsum(ones(1,nd)*dx)];
        
end

line = fgets(fid);
while isempty( str2num(line) )
    
    line = fgets(fid);
    
end

var = str2num(line);

nlines = str2num(line);

% Get the first node and iterate down the file
line = fgets(fid);
var = str2num(line);
zn(1) = var(1);

nd = var(3);
dz = (var(2) - zn(end)) / nd ;

zn = [zn zn(end)+cumsum(ones(1,nd)*dz)];

for ii = 2 : nlines
    
    line = fgets(fid);
    var = str2num(line);
    
    nd = var(2);
    dz = (var(1) - zn(end)) / nd ;
    
    zn = [zn zn(end)+cumsum(ones(1,nd)*dz)];
        
end

fclose(fid);