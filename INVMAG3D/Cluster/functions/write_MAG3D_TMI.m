function write_MAG3D_TMI(out_dir,file_name,H,I,Dazm,obsx,obsy,obsz,TMI,wd)
% function write_MAG3D_TMI(out_dir,file_name,H,I,D,obsx,obsy,obsz,TMI,wd)
% Write 3-component obs file in UBC format with prescribed values

ndata = length(obsx);

fid = fopen([out_dir '\' file_name],'w');

fprintf(fid,'%6.2f %6.2f %8.2f\n',I,Dazm,H);
fprintf(fid,'%6.2f %6.2f %8.2f\n',I,Dazm,1);
fprintf(fid,'%i\n',ndata);

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %12.8e %12.8e\n',obsx(ii),obsy(ii),obsz(ii),TMI(ii),wd(ii));
    
end


fclose(fid);
