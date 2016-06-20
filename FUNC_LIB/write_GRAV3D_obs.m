function write_GRAV3D_obs(file_name,obsx,obsy,obsz,data,wd)
% function write_MAG3D_TMI(out_dir,file_name,H,I,D,obsx,obsy,obsz,TMI,wd)
% Write 3-component obs file in UBC format with prescribed values

ndata = length(obsx);

fid = fopen(file_name,'w');
fprintf(fid,'%i\n',ndata);

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %12.6e %12.6e\n',obsx(ii),obsy(ii),obsz(ii),data(ii),wd(ii));
    
end


fclose(fid);
