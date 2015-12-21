function write2UBC(z0,vector)
% Function write a UBC mesh file from a 1D vector for Z only

fid=fopen('output.msh','w');
nbcell=length(vector);

%Creates header for number of cells 
fprintf(fid,'%i %i %i\n',4,4,nbcell);

% Set horizontal cells
fprintf(fid,'226790.00 9680990.00 %8.2f\n',z0);
fprintf(fid,'%8.5e %8.5e %8.5e %8.5e\n',500,500,500,500);
fprintf(fid,'%8.5e %8.5e %8.5e %8.5e\n',500,500,500,500);

for ii=1:nbcell
    fprintf(fid,'%8.5e ',vector(ii));
end

fclose(fid);
end