function DCIP_Batch_STD(work_directory,argin)
% Extract data and assign weighted standard deviation from DC log file
% Cycles through all the lines and all the folders and extract data from
% the *.log files. Observation files must be the only files starting by dc_
% and ip_ respectively. The reste of the name file will be used to create
% the new observation file with calculate error from the achieved misfit.
% 
% Author: D Fournier
% Last update : April 10th, 2013

home_dir = pwd;
fprintf('***START ASSIGN NEW ERRORS***\n\n')
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>
cd (work_directory)

%% Driver
max_num_lines=30000;

DCline_list=ls;

cond_model = [];

line_range = argin;
% count=0;
switch line_range 
    case 'all'
       range = 1:size(DCline_list,1)-2;
       
    otherwise
        range = str2num(argin);

end

% Cycle through all the DC lines
for oo=range
    
    % Move to Line folder
    cd (DCline_list(oo+2,:))
    

    files_in = ls;
    
    % Extract observation file name for ip_ and dc_ from the current
    % folder.
    for file=1:size(files_in,1)-2;

        look_at = strtrim(files_in(file+2,:));

        if strcmp(look_at(1:3),'dc_')==1 && strcmp(look_at(end-3:end),'.dat')==1
            look_at=strtrim(look_at);
            dc_obs_file = look_at;
        end

        if strcmp(look_at(1:3),'ip_')==1 && strcmp(look_at(end-3:end),'.dat')==1
            look_at=strtrim(look_at);
            ip_obs_file = look_at;
        end
    end
    
    %% Extract data    


    dc_file = importdata([dc_obs_file]);
    dc_data = dc_file.data; 
    if isnan(dc_data(1))==1
        
        dc_data = dc_data(:,2:end);
        
    end


    ip_file = importdata([ip_obs_file]);
    ip_data = ip_file.data; 
    if isnan(ip_data(1))==1
        
        ip_data = ip_data(:,2:end);
        
    end

        
     
    %% Read the dc log file and extract information
    fid=fopen('dcinv2d.log','rt');

    
    
    % Go through the log file and extract data and the last achieved misfit
    for ii=1:max_num_lines         	
    line=fgets(fid); %gets next line 

        if line==-1
            fprintf('File ended at line %i\n',ii);
            break
        end
        
        
        % Extract the conductivity from all the files and ouput it at the
        % end
%         if length(strtrim(line))>=length('Using default reference model')
%             description = strtrim(line);
%             if strcmp(description(1:29),'Using default reference model')==1
%                 cond_model = [cond_model str2num(description(31:end-6))];
%             end
%         end
        
        if isempty(regexp(line,'(reference conductivity model)','match'))==0
            
            temp = regexp(line,'(\d*+[.E-]+\d*)','match');
            cond_model = [cond_model str2double([temp{1} temp{2}])];
            
        end
        
        %Extract the last achieved misfit
        if isempty(regexp(line,'(achieved misfit)','match'))==0
            
            temp = regexp(line,'(\d*+[-.E+]+\d*)','match');
            misfit = str2double([temp{1} temp{2}]);
            
        end
    end    
    % Assign new error on data
    % STD * sqrt(achieved misfit / N)
    if isempty(dc_data)==1
        fprintf(['Could not read dc log file' dc_obs_file '**Please revise**'])
        break
    end
    dc_data(:,end) = dc_data(:,end) * sqrt(misfit/size(dc_data,1));
    
    % Save data to file
    
    wrt2file = fopen([dc_obs_file(1:end-4) '_std.obs'],'w');
    fprintf(wrt2file,'MIRA - DC Observations with weighted error\n');
%     fprintf(wrt2file,'pole-dipole\n');
    for kk = 1:size(dc_data,1)
        
                % Check sign for the pole configuration
                 if ((dc_data(kk,1) >= dc_data(kk,3) && dc_data(kk,3) >= dc_data(kk,4)) || (dc_data(kk,1) <= dc_data(kk,3) && dc_data(kk,3) <= dc_data(kk,4))) && dc_data(kk,5)<0
                    fprintf('Obs found with the wrong negative sign %8.4e ---> %8.4e\n',dc_data(kk,5),abs(dc_data(kk,5)))
                    
                    dc_data(kk,5)=abs(dc_data(kk,5));
                end
                
%                 if ((data(kk,1) <= data(kk,3) && data(kk,3) >= data(kk,4)) || (data(kk,1) >= data(kk,3) && data(kk,3) <= data(kk,4))) && data(kk,5)>0
%                     fprintf('Obs found with the wrong positive sign %8.4e ---> %8.4e\n',data(kk,5),-abs(data(kk,5)))
%                     
%                     data(kk,5)=-abs(data(kk,5));
%                 end
                
                % Write to file
        fprintf(wrt2file,'%8.6f\t%8.6f\t%8.6f\t%8.6f\t%12.8e\t%12.8e\t\n',dc_data(kk,1),dc_data(kk,2),dc_data(kk,3),dc_data(kk,4),dc_data(kk,5),dc_data(kk,6));
    end
    fclose(wrt2file);
    
    clear data
    
%% Repeat the operation for the IP
        %% Read the dc log file and extract information
    fid=fopen('ipinv2d.log','rt');

    
    % Go through the log file and extract data and the last achieved misfit
    for ii=1:max_num_lines         	
    line=fgets(fid); %gets next line 

        if line==-1
            fprintf('File ended at line %i\n',ii);
            break
        end
        
    

        %Extract the last achieved misfit
        %Extract the last achieved misfit
        if isempty(regexp(line,'(achieved misfit)','match'))==0
            
            temp = regexp(line,'(\d*+[-.E+]+\d*)','match');
            misfit = str2double([temp{1} temp{2}]);
            
        end
    end 
    
    % Assign new error on data
    % STD * sqrt(achieved misfit / N)
    if isempty(ip_data)==1
        fprintf(['Could not read ip log file' ip_obs_file '**Please revise**'])
        break
    end
    ip_data(:,end) = ip_data(:,end) * sqrt(misfit/size(ip_data,1));
    
    % Save data to file
    new_IP_obs = [ip_obs_file(1:end-4) '_std.obs'];
    wrt2file = fopen(new_IP_obs,'w');
    fprintf(wrt2file,'MIRA - IP Observations with weighted error\n');
    fprintf(wrt2file,'IPTYPE=1\n');
    for kk = 1:size(ip_data,1)
        fprintf(wrt2file,'%8.6f\t%8.6f\t%8.6f\t%8.6f\t%12.8e\t%12.8e\t\n',ip_data(kk,1),ip_data(kk,2),ip_data(kk,3),ip_data(kk,4),ip_data(kk,5),ip_data(kk,6));
    end
    fclose(wrt2file);
    
%% Back to current folder and repeat...
        
    
    cd ..

end
    cond_model = cond_model(:);
    cd ..
    save('Reference_Cond_models_ALL.dat','-ascii','cond_model')
    fprintf('\nAverage conductivity model used : %8.5e\n',mean(cond_model))
    fprintf('***END OF DCIP_Batch_STD***\n\n')
%     count
cd(home_dir);    
    
    
%         test_freq_match=strcmp(line(1:length(match_freq)),match_freq);