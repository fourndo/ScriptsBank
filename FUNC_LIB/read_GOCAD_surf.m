function [vrtx,trgl] = read_GOCAD_surf(tsfile)
% Load GOCAD surface *.ts file and return vertices and triangles
%
% INPUT
% tsfile : Triangulated surface file from OGCAD
% 
% OUTPUT
% vrtx: List of vertices [x1,y1,z1;...;xn,yn,zn]
% trgl: List of vertices used for each triangles [1,2,3;...;]


%% SCRIPT STARTS HERE
%% Load Horizontal dipole experiment
vrtx = [];
trgl = [];
 
fid=fopen(tsfile,'rt');    

line=fgets(fid); %gets next line

while isempty(regexp(line,'TFACE','match'))
    
    line=fgets(fid);
    
end

line=fgets(fid);


while ~isempty(regexp(line,'VRTX','match'))
       
    temp = regexp(line,'\s','split');
    vrtx = [vrtx;[str2num(temp{3}) str2num(temp{4}) str2num(temp{5})]];
    
    line=fgets(fid);
    
end

while ~isempty(regexp(line,'TRGL','match'))
       
    temp = regexp(line,'\s','split');
    trgl = [trgl;[str2num(temp{2}) str2num(temp{3}) str2num(temp{4})]];
    
    line=fgets(fid);
    
end

fclose(fid);

