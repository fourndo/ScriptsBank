function DCIP_Batch_STD_FB(work_directory,lines,dtype)
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

line_range = lines;
% count=0;
switch line_range 
    case 'all'
       range = 1:size(DCline_list,1)-2;
       
    otherwise
        range = str2num(lines);

end

% Cycle through all the DC lines
for oo=range
    
    % Move to Line folder
    cd (DCline_list(oo+2,:))
    
    folder_list = ls;
    
    for folder=1:size(folder_list,1)-2;
        
        % Move to first direction
        cd (folder_list(folder+2,:))

        files_in = ls;
        data = [];
        
        % Extract observation file name for ip_ and dc_ from the current
        % folder.
        for file=1:size(files_in,1)-2;
            
            look_at = strtrim(files_in(file+2,:));
            
            if strcmp(look_at(1:3),'dc_')==1 && strcmp(look_at(end-6:end),'edt.dat')==1
                look_at=strtrim(look_at);
                dc_obs_file = look_at(1:end-4);
            end
            
            if strcmp(look_at(1:3),'ip_')==1 && strcmp(look_at(end-6:end),'edt.dat')==1
                look_at=strtrim(look_at);
                ip_obs_file = look_at(1:end-4);
            end
        end
        
        

        switch dtype

            case 'dcip'
                %% Read the dc log file and extract information
                fid=fopen('dcinv2d.log','rt');

                if isempty(fid)

                    fprintf('Could not find a *.log for DC at: %s\n',[DCline_list(oo+2,:) folder_list(folder+2,:)])

                end
                line=fgets(fid); %gets next line

                find_start_data=('pole-dipole');
                find_achieved_misfit=('achieved misfit');


                % Go through the log file and extract data and the last achieved misfit
                while line~=-1         	

                    %Extract data 
                    if length(line)>=length(find_start_data)

                        if strcmp(strtrim(line),'pole-dipole')==1

                            line=fgets(fid);
                            x = str2num(line);

                            counter=1;

                            %Extract data line-by-line
                            while isempty(x)==0

                                data(counter,:)=x; 

                                line=fgets(fid);

                                x = str2num(line);
                                counter=counter+1;
                            end
                            
                        end

                    end

                    % Extract the conductivity from all the files and ouput it at the
                    % end
                    if length(strtrim(line))>=...
                            length('Using default reference model')
                        description = strtrim(line);
                        
                        if strcmp(description(1:29),'Using default reference model')==1
                            cond_model = [cond_model str2num(description(31:end-6))];
                        end
                        
                    end

                    %Extract the last achieved misfit
                    if length(strtrim(line))>=length(find_achieved_misfit)
                        
                        description = strtrim(line);
                        
                        if strcmp(description(1:15),find_achieved_misfit)==1
                            misfit=str2num(description(17:end));

                        end
                        
                    end

                    line=fgets(fid); %gets next line

                end
                
                % Assign new error on data
                % STD * sqrt(achieved misfit / N)
                if isempty(data)==1
                    fprintf(['Could not read dc log file' dc_obs_file '**Please revise**'])
                    break
                end
                
                data(:,end) = data(:,end) * sqrt(misfit/size(data,1));

                % Save data to file

                wrt2file = fopen([dc_obs_file '_std.obs'],'w');
                fprintf(wrt2file,'MIRA - DC Observations with weighted error\n');
                fprintf(wrt2file,'pole-dipole\n');
                
                for kk = 1:size(data,1)

                    % Check sign for the pole configuration
                    if ((data(kk,1) >= data(kk,3) && data(kk,3) >= data(kk,4)) || (data(kk,1) <= data(kk,3) && data(kk,3) <= data(kk,4))) && data(kk,5)<0
                        
                        fprintf('Obs found with the wrong negative sign %8.4e ---> %8.4e\n',data(kk,5),abs(data(kk,5)))
                        data(kk,5)=abs(data(kk,5));
                        
                    end

                    % if ((data(kk,1) <= data(kk,3) && data(kk,3) >= data(kk,4)) || (data(kk,1) >= data(kk,3) && data(kk,3) <= data(kk,4))) && data(kk,5)>0
                    %     fprintf('Obs found with the wrong positive sign %8.4e ---> %8.4e\n',data(kk,5),-abs(data(kk,5)))
                    % 
                    %     data(kk,5)=-abs(data(kk,5));
                    % end

                    % Write to file
                    fprintf(wrt2file,'%i %i %i %i %12.8f %12.8f\n',data(kk,1),data(kk,2),data(kk,3),data(kk,4),data(kk,5),data(kk,6));
                    
                end
                
                fclose(wrt2file);

                fprintf('DC Line: %s DONE\n',[DCline_list(oo+2,:) folder_list(folder+2,:)]);

                clear data

                %% Repeat the operation for the IP
                fid=fopen('ipinv2d.log','rt');

                if isempty(fid)

                    fprintf('Could not find a *.log for IP at: %s\n',[DCline_list(oo+2,:) folder_list(folder+2,:)])

                end
                
                line=fgets(fid); %gets next line 

                find_start_data=('pole-dipole');
                find_achieved_misfit=('achieved misfit');

                data = [];
                % Go through the log file and extract data and the last achieved misfit
                while line~=-1

                    %Extract data
                    if length(line)>=length(find_start_data)
                        if strcmp(strtrim(line),'pole-dipole')==1
                            line=fgets(fid);
                            x = str2num(line);
                            counter=1;
                            while isempty(x)==0

                                data(counter,:)=x;  %Extract data line-by-line

                                line=fgets(fid);

                                x = str2num(line);
                                counter=counter+1;
                            end
                        end

                    end

                    %Extract the last achieved misfit
                    if length(strtrim(line))>=length(find_achieved_misfit)
                        description = strtrim(line);
                        if strcmp(description(1:15),find_achieved_misfit)==1
                            misfit=str2num(description(17:end));

                        end
                    end

                    line=fgets(fid); %gets next line

                end 

                % Assign new error on data
                % STD * sqrt(achieved misfit / N)
                if isempty(data)==1
                    
                    fprintf(['Could not read ip log file' ip_obs_file '**Please revise**'])
                    break
                    
                end
                
                data(:,end) = data(:,end) * sqrt(misfit/size(data,1));

                % Save data to file
                new_IP_obs = [ip_obs_file '_std.obs'];
                wrt2file = fopen(new_IP_obs,'w');
                fprintf(wrt2file,'MIRA - IP Observations with weighted error\n');
                fprintf(wrt2file,'pole-dipole\n');
                for kk = 1:size(data,1)
                    fprintf(wrt2file,'%i %i %i %i %12.8f %12.8f\n',data(kk,1),data(kk,2),data(kk,3),data(kk,4),data(kk,5),data(kk,6));
                end
                fclose(wrt2file);

                fprintf('IP Line: %s DONE\n',[DCline_list(oo+2,:) folder_list(folder+2,:)]);

                % Back to current folder and repeat...
                cd ..

            case 'dc'   
                %% Read the dc log file and extract information
                fid=fopen('dcinv2d.log','rt');

                if isempty(fid)

                    fprintf('Could not find a *.log for DC at: %s\n',[DCline_list(oo+2,:) folder_list(folder+2,:)])

                end
                line=fgets(fid); %gets next line

                find_start_data=('pole-dipole');
                find_achieved_misfit=('achieved misfit');


                % Go through the log file and extract data and the last achieved misfit
                while line~=-1         	

                    %Extract data 
                    if length(line)>=length(find_start_data)

                        if strcmp(strtrim(line),'pole-dipole')==1

                            line=fgets(fid);
                            x = str2num(line);

                            counter=1;

                            %Extract data line-by-line
                            while isempty(x)==0

                                data(counter,:)=x; 

                                line=fgets(fid);

                                x = str2num(line);
                                counter=counter+1;
                            end
                            
                        end

                    end

                    % Extract the conductivity from all the files and ouput it at the
                    % end
                    if length(strtrim(line))>=...
                            length('Using default reference model')
                        description = strtrim(line);
                        
                        if strcmp(description(1:29),'Using default reference model')==1
                            cond_model = [cond_model str2num(description(31:end-6))];
                        end
                        
                    end

                    %Extract the last achieved misfit
                    if length(strtrim(line))>=length(find_achieved_misfit)
                        
                        description = strtrim(line);
                        
                        if strcmp(description(1:15),find_achieved_misfit)==1
                            misfit=str2num(description(17:end));

                        end
                        
                    end

                    line=fgets(fid); %gets next line

                end
                
                % Assign new error on data
                % STD * sqrt(achieved misfit / N)
                if isempty(data)==1
                    fprintf(['Could not read dc log file' dc_obs_file '**Please revise**'])
                    break
                end
                
                data(:,end) = data(:,end) * sqrt(misfit/size(data,1));

                % Save data to file

                wrt2file = fopen([dc_obs_file '_std.obs'],'w');
                fprintf(wrt2file,'MIRA - DC Observations with weighted error\n');
                fprintf(wrt2file,'pole-dipole\n');
                
                for kk = 1:size(data,1)

                    % Check sign for the pole configuration
                    if ((data(kk,1) >= data(kk,3) && data(kk,3) >= data(kk,4)) || (data(kk,1) <= data(kk,3) && data(kk,3) <= data(kk,4))) && data(kk,5)<0
                        
                        fprintf('Obs found with the wrong negative sign %8.4e ---> %8.4e\n',data(kk,5),abs(data(kk,5)))
                        data(kk,5)=abs(data(kk,5));
                        
                    end

                    % if ((data(kk,1) <= data(kk,3) && data(kk,3) >= data(kk,4)) || (data(kk,1) >= data(kk,3) && data(kk,3) <= data(kk,4))) && data(kk,5)>0
                    %     fprintf('Obs found with the wrong positive sign %8.4e ---> %8.4e\n',data(kk,5),-abs(data(kk,5)))
                    % 
                    %     data(kk,5)=-abs(data(kk,5));
                    % end

                    % Write to file
                    fprintf(wrt2file,'%i %i %i %i %12.8f %12.8f\n',data(kk,1),data(kk,2),data(kk,3),data(kk,4),data(kk,5),data(kk,6));
                    
                end
                
                fclose(wrt2file);

                fprintf('DC Line: %s DONE\n',[DCline_list(oo+2,:) folder_list(folder+2,:)]);

                clear data
                
                % Back to current folder and repeat...
                cd ..
        end
    end

                cd ..
end
            
cond_model = cond_model(:);
cd ..
save('Reference_Cond_models_ALL.dat','-ascii','cond_model')
fprintf('\nAverage conductivity model used : %8.5e\n',mean(cond_model))
fprintf('***END OF DCIP_Batch_STD***\n\n')

end
    
    
%         test_freq_match=strcmp(line(1:length(match_freq)),match_freq);