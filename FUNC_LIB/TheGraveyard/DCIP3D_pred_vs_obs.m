% function rem_check_sign(check_file,obs_file)
%
% Load a UBC-DCIP3D "check_sign.txt" file and remove bad lines from the
% corresponding obs file

% Written by: D. Fournier
% Last update: January 3th, 2014

clear all
close all

%% USER INPUT
work_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Modelling\Inversion\BrokenHammer\Inv1';
pred_file = [work_dir '\dcinv3d_01.pre'];
obs_file = [work_dir '\4180_DC_obs.dat'];

% Name for edited output observation file
obs_file_out = [work_dir '\4180_DC_obs_signchecked.dat'];


%% ||| SCRIPT STARTS HERE |||
pred = fopen(pred_file,'r');
% Import check sign data
line = fgets(pred);

% Pre-allocate 2D array
d_pred = zeros(1,15);
count = 0;
while line~=-1
    
    
    while isempty(str2num(line))==1
        
        
        line = fgets(pred);
        
        
    end
    
    % Read transmitter line
    tx = str2num(line);
    nrx = tx(7);
    
    for ii = 1 : nrx
        
        count= count+1;
        line = fgets(pred);
       
        d_pred(count,1:6) = tx(1:6);
        
        d_pred(count,7:end) = str2num(line);
        
        
    end
    
    
    line = fgets(pred);
    
end

fclose(pred);

% Open obsfile and change sign check sign flag
obs_in = fopen(obs_file,'r');
line = fgets(obs_in);

% Write new obs file
obs_out = fopen(obs_file_out,'w');
count = 0;
count_reject = 0;
while line~=-1
    
    
        
    while isempty(str2num(line))==1
        
        fprintf(obs_out,'%s',line);
        line = fgets(obs_in);
        
        
    end
    
    % Read transmitter line
    tx = str2num(line);
    nrx = tx(7);
    
    % Change for 1 receiver line per transmitter
    tx(7) = 1;
    
    for ii = 1 : nrx
        
        count= count+1;
        line = fgets(obs_in);
        
        data = str2num(line);
        
        % If data could not be fitted within 100%
        if abs(data(7)/d_pred(count,13)) > 1e+2 || abs(data(7)/d_pred(count,13))<1e-2
        
            fprintf('Data point #%i has been discarded\n',count);   
            count_reject = count_reject + 1;
            
        
        else
            
            % If data has wrong sign, multiply by -1
            if data(7)/d_pred(count,13) < 0
            
            data(7) = data(7)*-1;
            
            end
            
            % Write transmitter line
            for jj = 1 : length(tx)-1
                fprintf(obs_out,'%15.8e ',tx(jj));
            end
            
            fprintf(obs_out,'%i\n',1);
            
            % Write Receiver line
            for jj = 1 : length(data)
                fprintf(obs_out,'%15.8e ',data(jj));
            end

            fprintf(obs_out,'\n');
 
        end
        
    end
    
    
    line = fgets(obs_in);
    
end
        

fclose(obs_in);
fclose(obs_out);

