function write_MAG3D_TMI(file_name,B0,BI,BD,MI,MD,obsx,obsy,obsz,TMI,wd)
% function write_MAG3D_TMI(file_name,B0,BI,BD,MI,MD,obsx,obsy,obsz,TMI,wd)
% Write 3-component obs file in UBC format with prescribed values

ndata = length(obsx);

fid = fopen(file_name,'w');

fprintf(fid,'%6.2f %6.2f %8.2f\n',BI,BD,B0);
fprintf(fid,'%6.2f %6.2f %8.2f\n',MI,MD,1);
fprintf(fid,'%i\n',ndata);

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %12.6e %12.6e\n',obsx(ii),obsy(ii),obsz(ii),TMI(ii),wd(ii));
    
end


fclose(fid);
