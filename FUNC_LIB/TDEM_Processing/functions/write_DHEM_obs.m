function write_DHEM_obs(work_dir,filename,dobs,tx,recv,t_out,component,pc_err)

ndv = 99999;
% Create obs file
fid = fopen([work_dir '\' filename],'w');

% Create receiver array

%% Write transmiter header
fprintf(fid,'IGNORE %i\n\n',ndv);
fprintf(fid,'N_TRX %i\n\n',1);
fprintf(fid,'TRX_ORIG\n');
fprintf(fid,'%i\n',size(tx,1));

for ii = 1 : size(tx,1)
    
    fprintf(fid,'%12.3f %12.3f %12.3f\n',tx(ii,1),tx(ii,2),tx(ii,3));
    
end

nrecv = length(unique(recv));
fprintf(fid,'\nN_RECV %i\n', nrecv);
fprintf(fid,'N_TIME %i\n',length(t_out));


%% Write data
for ii = 1 : nrecv;
    
    subdata = dobs( recv == ii , : );
    
    for jj = 1 : length(t_out);
               
        index = t_out(jj);
        
        % If no data for a specific receiver/time, post no data
        if isnan(subdata(index,17))==1 || isnan(subdata(index,19))==1 || isnan(subdata(index,21))==1
            
            fprintf(fid,'%12.3f %12.3f %12.3f %3.8e %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n',...
            subdata(index,1),subdata(index,2),subdata(index,3),subdata(index,4),...
            ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,...
            ndv,ndv,ndv,ndv,ndv,ndv);
        
        else
            
            % Assign error on magnitude of field instead of components
            % individually
%             magdB = ( Atlas.dBx(index(ii),time_out(jj))^2 +...
%                 Atlas.dBy(index(ii),time_out(jj))^2 +...
%                 Atlas.dBz(index(ii),time_out(jj))^2 )^ 0.5;
            
            err_dBx = abs(subdata(index,17))*pc_err + pc_err*subdata(index,18);
            err_dBy = abs(subdata(index,19))*pc_err + pc_err*subdata(index,20);
            err_dBz = abs(subdata(index,21))*pc_err + pc_err*subdata(index,22);

            fprintf(fid,'%12.3f %12.3f %12.3f %3.8e %i %i %i %i %i %i %i %i %i %i %i %i %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ',...
                subdata(index,1),subdata(index,2),subdata(index,3),subdata(index,4),...
                ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv,...
                subdata(index,17),err_dBx,subdata(index,19),err_dBy,subdata(index,21),err_dBz);

            % Write dB/dt for input components in coil
%             for kk = 1:3
% 
%                 if kk==1 && sum(coil==kk)~=0
% 
%                     fprintf(fid,'%12.8e %12.8e ',...
%                     Atlas.dBx(index(ii),t_out(jj)),err_dBx);
% 
%                 elseif kk==2 && sum(coil==kk)~=0
% 
%                     fprintf(fid,'%12.8e %12.8e ',...
%                     Atlas.dBy(index(ii),t_out(jj)),err_dBy);
% 
%                 elseif kk==3 && sum(coil==kk)~=0
% 
%                     fprintf(fid,'%12.8e %12.8e ',...
%                     Atlas.dBz(index(ii),t_out(jj)),err_dBz);
% 
%                 else
% 
%                     fprintf(fid,'%i %i ', ndv,ndv);
% 
%                 end
% 
%             end
        
            fprintf(fid,'\n');
            
        end
    
    end
    
        fprintf(fid,'\n');
end

fclose(fid);