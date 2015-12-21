function make_UBC_mesh(XYZo,dX,dY,dZ)

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

%Create UBC mesh file
fid1=fopen('UBC_mesh.msh', 'w');
fprintf(fid1, '%i %i %i\n', nX, nY, nZ);
fprintf(fid1, '%i %i %i\n', XYZo(1), XYZo(2), XYZo(3));

for jj=1:nX
    fprintf(fid1, '%4.2f ', dX(jj));    
end

fprintf(fid1,'\n');

for ii=1:nY
           fprintf(fid1,'%4.2f ', dY(ii));
end

fprintf(fid1, '\n');

for kk=1 : nZ
       fprintf(fid1, '%4.2f ', dZ(kk));
end

fclose(fid1);