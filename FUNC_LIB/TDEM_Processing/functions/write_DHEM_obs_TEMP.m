
clear all
close all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\VTEM\3D\DO27\FWR';
filename = 'VTEM_DO27_dBdz.dat';

load([work_dir '\VTEM_DIGHEM_DO27']);

ndv = 99999;

pc_err = 0.05;

data = DATA_L150 * pi*12^2 * 1e-12;

%% Pre-compute uncertainty floor
uncert_flr = zeros(length(time),1);
for ii = 1 : length(time)
    
    uncert_flr(ii) = std(data(:,ii));
    
end
% Create obs file
fid = fopen([work_dir '\' filename],'w');

%% Write transmiter header and data
fprintf(fid,'IGNORE %i\n\n',ndv);
fprintf(fid,'N_TRX %i\n\n',size(data,1));


for ii = 1 : size(data,1);
    
    fprintf(fid,'TRX_LOOP\n');
    fprintf(fid,'%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n', X_L150(ii) , Y_L150(ii) , 468 , 13 , 0 , 0  );
    
    fprintf(fid,'\nN_RECV %i\n',1);
    fprintf(fid,'N_TIME %i\n',length(time));
    
    for jj = 1 : length(time);
        
        % If no data for a specific receiver/time, post no data
            
        % Assign error on magnitude of field instead of components
        % individually

        err_Hz = abs(data(ii,jj))*pc_err + pc_err*uncert_flr(jj);
%             err_dBy = abs(magdB)*pc_err + pc_err*std_dBy(time_out(jj));
%             err_dBz = abs(magdB)*pc_err + pc_err*std_dBz(time_out(jj));

        fprintf(fid,'%12.3f %12.3f %12.3f %3.8f ',...
            X_L150(ii) , Y_L150(ii) , 468 ,time(jj));

        for kk = 1 : 16 

            fprintf(fid,'%i ',ndv);

        end

        % Write Z-component
%         fprintf(fid,'%12.8e %12.8e ',data(ii,jj), err_Hz);

        % Write dB/dt for input components in coil
%         for kk = 1 : 5 
% 
%             fprintf(fid,'%i ',ndv);
% 
%         end

        fprintf(fid,'%12.8e %12.8e ',data(ii,jj), err_Hz);
        
        fprintf(fid,'\n');
            
           
    end
    
        fprintf(fid,'\n');
end

fclose(fid);