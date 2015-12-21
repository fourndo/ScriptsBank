% Matlab contours to Gocad line.
clear all
close all

work_dir='C:\Projects\tempo\DEM';

% Load topography grid
Topo=load([work_dir '\' 'CDED_115g_InvDistgrid10m_core.dat']);
Topo2D = Topo(:,3);
Zmax = ceil(max(Topo2D)/10) * 10;
Zmin = floor(min(Topo2D)/10) * 10;
Topo2D = reshape(Topo2D,2011,1598);

% Generate contours
C=contour(Topo2D,Zmin:5:Zmax);
save([work_dir 'C_10m'],'C');
% load C_10m;

% Write contours to file
fid = fopen([work_dir '\' 'Contours_5m_core.dat'],'w');
fprintf(fid,'GOCAD PLine 1\n')
fprintf(fid,'HEADER {\n')
fprintf(fid,'name:Contours_10m_global\n')
fprintf(fid,'}\n')
fprintf(fid,'GOCAD_ORIGINAL_COORDINATE_SYSTEM\n')
fprintf(fid,'NAME Default\n')
fprintf(fid,'AXIS_NAME "X" "Y" "Z"\n')
fprintf(fid,'AXIS_UNIT "m" "m" "m"\n')
fprintf(fid,'ZPOSITIVE Elevation\n')
fprintf(fid,'END_ORIGINAL_COORDINATE_SYSTEM\n')
fprintf(fid,'PROPERTY_CLASS_HEADER Z {\n')
fprintf(fid,'is_z:on\n')
fprintf(fid,'}\n')


count = 1;
node = 1;
X0 = Topo(1,1);
Y0 = Topo(1,2);
while count<=size(C,2)
    fprintf(fid,'ILINE\n');
    nnodes = C(2,count);
    Z = C(1,count);
    count=count+1;
    vrtx = [];
    for ii=1:nnodes
        
        vrtx(ii) = node;
        fprintf(fid,'VRTX %i %12.3f %12.3f %8.3f\n',vrtx(ii),X0+C(2,count)*10-10,Y0+C(1,count)*10-10,Z);
        
        count=count+1;
        node = node+1;
        
    end
    
    for jj = 1:nnodes-1
        
        fprintf(fid,'SEG %i %i\n',vrtx(jj),vrtx(jj+1));
        
    end
    
end

fprintf(fid,'END\n');
fclose(fid);