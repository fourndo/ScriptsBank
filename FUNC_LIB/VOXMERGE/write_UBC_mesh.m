function write_UBC_mesh(work_dir,x0,y0,z0,dx,dy,dz)
% Generate *.msh file in UBC format for giving mesh parameters

nx = length(dx);
ny = length(dy);
nz = length(dz);

%Create UBC mesh file
fid1=fopen([work_dir '\UBC_mesh.msh'], 'w');
fprintf(fid1, '%i %i %i\n', nx, ny, nz);
fprintf(fid1, '%i %i %i\n', x0, y0, z0);

for jj=1:nx
    fprintf(fid1, '%4.2f ', dx(jj));    
end

fprintf(fid1,'\n');

for ii=1:ny
           fprintf(fid1,'%4.2f ', dy(ii));
end

fprintf(fid1, '\n');

for kk=1 : nz
       fprintf(fid1, '%4.2f ', dz(kk));
end

fclose(fid1);