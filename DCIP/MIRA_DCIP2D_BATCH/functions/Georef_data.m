function Georef_data(survey_full,work_dir,argin)
% Add XYZ coordinates to line/stations data

cd (work_dir)
cd ..
home_dir = pwd;
fprintf('***START GEO-REFERENCING***\n\n')
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>

cd (work_dir)


%% Driver
DCline_list=ls;

cond_model = [];

line_range = argin;

switch line_range 
    case 'all'
       range = 1:size(DCline_list,1)-2;
       
    otherwise
        range = str2num(argin);

end

% Cycle through all the DC lines
        wrt2filedc_all = fopen([home_dir '\DC_line_ALL_UTM.obs'],'w');
        fprintf(wrt2filedc_all,'!!MIRA - DCIP3D UTM BATCH\n');
%         fprintf(wrt2file,[dc_obs_file '\n']);
        fprintf(wrt2filedc_all,'\n\n');
        
        wrt2fileip_all = fopen([home_dir '\IP_line_ALL_UTM.obs'],'w');
        fprintf(wrt2fileip_all,'!!MIRA - DCIP3D UTM BATCH\n');
        fprintf(wrt2fileip_all,'IPTYPE=1\n');
        fprintf(wrt2fileip_all,'\n\n');
  
% Sort rows of surveys along line and stations
survey_full = sortrows(survey_full,[5 4]);

for oo=range
    
    % Move to Line folder
    cd (DCline_list(oo+2,:))
    
    

    files_in = ls;
    
        % Extract observation file name for ip_ and dc_ from the current
        % folder.
        for file=1:size(files_in,1)-2;
            
            look_at = strtrim(files_in(file+2,:));
            
            if strcmp(look_at(1:3),'dc_')==1 && strcmp(look_at(end-7:end),'_std.obs')==1
                look_at=strtrim(look_at);
                dc_obs_file = look_at(1:end-4);
            end
            
            if strcmp(look_at(1:3),'ip_')==1 && strcmp(look_at(end-7:end),'_std.obs')==1
                look_at=strtrim(look_at);
                ip_obs_file = look_at(1:end-4);
            end
        end
    
        %% Read the dc log file and extract information
        dc_obs = importdata([dc_obs_file '.obs']);
        ip_obs = importdata([ip_obs_file '.obs']);
        
        dc_data = dc_obs.data;
        ip_data = ip_obs.data;

        %Create new matrix (line,station,X,Y,Z,prop1,...)
%         dc_data_geo=zeros(size(dc_data,1),size(dc_data,2)+3);
%         ip_data_geo=zeros(size(ip_data,1),size(ip_data,2)+3);
        
        
        % Extract line number
        line_name = regexp(dc_obs_file,'[-+]?\d*[NSEW]','match');
        
        % Put negative sign to differentiate between North/South or
        % East/West lines
        if strcmp(line_name{1}(end),'S') || strcmp(line_name{1}(end),'E')
            
            line_num = str2double(line_name{1}(1:end-1)) * -1;
            
        else
            
            line_num = str2double(line_name{1}(1:end-1));
            
        end
        
        survey = survey_full(survey_full(:,5)==line_num,:);
        
        %Densify survey files by interpolating between stations
        %Set dl interval along line
        dl_steps=5; 
        meter_per_foot=1.0;
        counter=1;
        survey_interp=[];
        
        for jj=1:size(survey,1)-1
            interpolate=[];
            
            % Check that we are still on the same line
            % If yes, interpolate stations
            if survey(jj,5) == survey(jj+1,5)    
                %Compute length between points
                L = survey(jj+1,4)-survey(jj,4);
                dX = abs(survey(jj+1,1)-survey(jj,1));
                dY = abs(survey(jj+1,2)-survey(jj,2));
                dZ = survey(jj+1,3)-survey(jj,3);

                if dX==0 && dY==0
                    continue
                end
                azm(counter)=atan(dX/dY);
                counter=counter+1;
                %Compute number of stations to add
                num_stations=abs(round(L/dl_steps));

                dl = L / (num_stations);
                dx = dX / (num_stations);
                dy = dY / (num_stations);
                dz = dZ / (num_stations);

        %         interpolate=survey(jj,:);
                interpolate(1:num_stations,5)=survey(jj,5);
                interpolate(1:num_stations,4)=survey(jj,4) + dl * (0:num_stations-1)' ;
                interpolate(1:num_stations,1)=survey(jj,1) + dx * (0:num_stations-1)' ;
                interpolate(1:num_stations,2)=survey(jj,2) + dy * (0:num_stations-1)' ;
                interpolate(1:num_stations,3)=survey(jj,3) + dz * (0:num_stations-1)' ;

                survey_interp=[survey_interp;interpolate];
                clear interpolate
            end
            
        end
        
        % Merge last row
        survey_interp=[survey_interp;survey(end,:)];
        


        X=survey(:,1);
        Y=survey(:,2);
        Z=survey(:,3);


        station_start = min(min(min(dc_data(:,1:4))),min(min(ip_data(:,1:4))));   %Smallest station measured at
        station_in = min(survey(:,4));
        Xin=X(1);
        Yin=Y(1);
        xdir=(X(end)-X(1)+1e-20)/abs(X(end)-X(1)+1e-20);
        ydir=(Y(end)-Y(1)+1e-20)/abs(Y(end)-Y(1)+1e-20);

        % Extrapolate survey line beyond first survey station
        while station_in > station_start
            station_in = station_in - dl_steps;
            Xin = Xin - xdir * dl_steps * meter_per_foot * abs(sin(mean(azm)));
            Yin = Yin - ydir * dl_steps * meter_per_foot * abs(cos(mean(azm)));

            interpolate(1,5)= line_num;
            interpolate(1,4)= station_in;
            interpolate(1,1)= Xin;
            interpolate(1,2)= Yin;
            interpolate(1,3) = mean(Z) ;

            survey_interp=[survey_interp;interpolate];       
        end


        % Extrapolate survey line beyond last survey stations
        station_in=max(survey(:,4));
        station_stop=max(max(max(dc_data(:,1:4))),max(max(ip_data(:,1:4))));
        Xin=X(end);
        Yin=Y(end);  
        xdir=(X(end)-X(1)+1e-20)/abs(X(end)-X(1)+1e-20);
        ydir=(Y(end)-Y(1)+1e-20)/abs(Y(end)-Y(1)+1e-20);
        while station_in < station_stop
            station_in = station_in + dl_steps;
            Xin = Xin + xdir * dl_steps * meter_per_foot * abs(sin(mean(azm)));
            Yin = Yin + ydir * dl_steps * meter_per_foot * abs(cos(mean(azm)));

            interpolate(1,5)= line_num;
            interpolate(1,4)= station_in;
            interpolate(1,1)= Xin;
            interpolate(1,2)= Yin;
            interpolate(1,3) = mean(Z) ;

            survey_interp=[survey_interp;interpolate];       
        end


        clear azm Z

%             end
        
        dc_data=sortrows(dc_data,1);
        survey = survey_interp;
        
        
        %Assign coordinates to DC data
        %line,station -> line,station,X,Y,Z
        wrt2file = fopen([dc_obs_file '_UTM.obs'],'w');
        fprintf(wrt2file,'!!MIRA - DCIP3D UTM BATCH\n');
%         fprintf(wrt2file,[dc_obs_file '\n']);
        fprintf(wrt2file,'\n\n');
        
        wrt2file_2D = fopen([dc_obs_file '_UTM2D.obs'],'w');
        fprintf(wrt2file_2D,'!!MIRA - DCIP2D UTM BATCH\n');
        fprintf(wrt2file_2D,'pole-dipole');
        fprintf(wrt2file_2D,'\n');
        
        count =1;
        while count<size(dc_data,1)
            
            % Get all the dipoles from the same transmiter
            num_pole = sum (dc_data(:,1)==dc_data(count,1) & dc_data(:,2)==dc_data(count,2));
            
            % Find coordinate of transmiter
            tx1 = survey(dc_data(count,1)==survey(:,4),1);
            ty1 = survey(dc_data(count,1)==survey(:,4),2);
            tx2 = survey(dc_data(count,2)==survey(:,4),1);
            ty2 = survey(dc_data(count,2)==survey(:,4),2);
            fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %i\n',tx1,ty1,tx2,ty2,num_pole);
            fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %i\n',tx1,ty1,tx2,ty2,num_pole);
            
            for kk = 1:num_pole
                
                % Find coordinate of receivers
                rx1 = survey(dc_data(count,3)==survey(:,4),1);
                ry1 = survey(dc_data(count,3)==survey(:,4),2);
                
                rx2 = survey(dc_data(count,4)==survey(:,4),1);
                ry2 = survey(dc_data(count,4)==survey(:,4),2);
                
                % Write 3D format
                fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',rx1,ry1,rx2,ry2,dc_data(count,5),dc_data(count,6));
                fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',rx1,ry1,rx2,ry2,dc_data(count,5),dc_data(count,6));

                % Write 2D format for DOI
                fprintf(wrt2file_2D,'%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\n',tx1,tx2,rx1,rx2,dc_data(count,5),dc_data(count,6));
                
                count = count+1;
            end
            fprintf(wrt2file,'\n');
            fprintf(wrt2filedc_all,'\n');
        end
        fclose(wrt2file);
        fclose(wrt2file_2D);
        
        ip_data=sortrows(ip_data,1);
        %Assign coordinates to IP data
        %line,station -> line,station,X,Y,Z
        wrt2file = fopen([ip_obs_file '_UTM.obs'],'w');
        fprintf(wrt2file,'!!MIRA - DCIP3D UTM BATCH\n');
        fprintf(wrt2file,'IPTYPE=1\n');
        fprintf(wrt2file,'\n');
        
        wrt2file_2D = fopen([ip_obs_file '_UTM2D.obs'],'w');
        fprintf(wrt2file_2D,'!!MIRA - DCIP2D UTM BATCH\n');
        fprintf(wrt2file_2D,'pole-dipole');
        fprintf(wrt2file_2D,'\n');
        
        count =1;
        while count<size(ip_data,1)
            
            % Get all the dipoles from the same transmiter
            num_pole = sum (ip_data(:,1)==ip_data(count,1) & ip_data(:,2)==ip_data(count,2));
            
            % Find coordinate of transmiter
            tx1 = survey(ip_data(count,1)==survey(:,4),1);
            ty1 = survey(ip_data(count,1)==survey(:,4),2);
            tx2 = survey(ip_data(count,2)==survey(:,4),1);
            ty2 = survey(ip_data(count,2)==survey(:,4),2);
            fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %i\n',tx1,ty1,tx2,ty2,num_pole);
            fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %i\n',tx1,ty1,tx2,ty2,num_pole);
            for kk = 1:num_pole
                
                % Find coordinate of receivers
                rx1 = survey(ip_data(count,3)==survey(:,4),1);
                ry1 = survey(ip_data(count,3)==survey(:,4),2);
                
                rx2 = survey(ip_data(count,4)==survey(:,4),1);
                ry2 = survey(ip_data(count,4)==survey(:,4),2);
                
                % Write to file 3D
                fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',rx1,ry1,rx2,ry2,ip_data(count,5),ip_data(count,6));
                fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',rx1,ry1,rx2,ry2,ip_data(count,5),ip_data(count,6));
                
               % Write 2D format for DOI
                fprintf(wrt2file_2D,'%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\n',tx1,tx2,rx1,rx2,ip_data(count,5),ip_data(count,6));
               
                count = count+1;
            end
            
            fprintf(wrt2file,'\n');
            fprintf(wrt2fileip_all,'\n');
        end
        fclose(wrt2file);
        fclose(wrt2file_2D);
        
            

        cd ..
end

fclose(wrt2filedc_all);
fclose(wrt2fileip_all);   
% fclose all;

fprintf('**DONE GEO-REFERENCING**\n\n')