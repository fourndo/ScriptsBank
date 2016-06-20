function write_MAG3D_3C(file_name,H,I,D,obsx,obsy,obsz,Bx,By,Bz,wdx,wdy,wdz)
% Write 3-component obs file in UBC format with prescribed values

ndata = length(obsx);

fid = fopen(file_name,'w');

fprintf(fid,'%6.2f %6.2f %8.2f\n',I,D,H);
fprintf(fid,'%6.2f %6.2f %8.2f\n',I,D,0);
fprintf(fid,'%i\n',ndata*3);

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %5.2f %5.2f %12.8e %12.8e\n',obsx(ii),obsy(ii),obsz(ii),0,90,Bx(ii),wdx(ii));
    
end

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %5.2f %5.2f %12.8e %12.8e\n',obsx(ii),obsy(ii),obsz(ii),0,0,By(ii),wdy(ii));
    
end

for ii = 1 : ndata
    
    fprintf(fid,'%12.8e %12.8e %12.8e %5.2f %5.2f %12.8e %12.8e\n',obsx(ii),obsy(ii),obsz(ii),90,0,Bz(ii),wdz(ii));
    
end

fclose(fid);
