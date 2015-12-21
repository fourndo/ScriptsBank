function write_e3d_obs(filename,freq,tx,rx)
% Function write_e3d_obs(filename,freq,tx,rx)
% Write obsservation file for E3Dinv
% TO FINISH DESCRIPTION...

fid = fopen(filename,'w');

fprintf(fid,'! Exported from write_e3d_obs: D.Fournier\n');
fprintf(fid,'IGNORE NaN\n\n');
fprintf(fid,'N_TRX %i\n', size(tx,1) / 3 );

for ii = 1 : size(tx,1)/3;
    
    for jj = 1 : length(freq)
        
        fprintf(fid,'TRX_LOOP\n');
        fprintf(fid,'%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n',tx(ii,1:end-1));
        fprintf(fid,'FREQUENCY %12.8e\n',freq(jj));
        fprintf(fid,'N_RECV %i\n',size(rx,1));
        
        for kk = 1 : size(rx,1)
            
            fprintf(fid,'%12.8e %12.8e %12.8e ',rx(kk,1:3));
            
%             for ll = 1 : 24
%                 
%                 fprintf(fid,'NaN ');
%                 
%             end
            
            fprintf(fid,'\n');
            
        end
        
    end
    
end

fclose(fid);
