function write_h3d_obs(work_dir,filename,x,y,z,tx,times)
% Write to file h3d observation file

mcell = sum(isnan(x)==0);

%% Write data file header
ndv=-99999;

fid = fopen([work_dir '\' filename],'w');

% Header
fprintf(fid,'IGNORE %i\n',ndv);
fprintf(fid,'\n');

fprintf(fid,'N_TRX 1\n');
fprintf(fid,'\n');

%% write transmitter
fprintf(fid,'TRX_ORIG\n');
nnode=size(tx,1);
fprintf(fid,'%i\n',nnode);
for ii = 1 : nnode
    fprintf(fid,'%6.2f %6.2f %6.2f\n',tx(ii,:));
end

fprintf(fid,'\n');


%% write data    
ntime=size(times,1);

%Data head for Ex
fprintf(fid,'N_RECV %i\n',mcell);
fprintf(fid,'N_TIME %i\n',ntime);
fprintf(fid,'\n');

progress = -1;
tic
for ii = 1 : mcell
    
    for jj= 1 : ntime  
        
        fprintf(fid,'%8.2f  %8.2f %8.2f %8e ',x(ii),y(ii),z(ii),times(jj));
        
        for mm=1:18
            fprintf(fid,'%8.0f ',1);
        end
        
        fprintf(fid,'\n');
    end
    
    d_iter = floor(ii/mcell*20);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
    
end
fclose all;