% function Georef_data_FB(survey_full,work_dir,out_dir,units)
% Add XYZ coordinates to line/stations data
close all
clear all

out_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\Titan\3D_Offline';
work_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\From_Client\2013_12_20_Titan\Results\3D_OfflineTx\Raw';
survey =load('C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\Titan\4180_survey_Titan_ROT40.dat');
argin ='all';
units ='meter';

fprintf('***START GEO-REFERENCING***\n\n')

%% Driver

%Densify survey files by interpolating between stations
%Set dl interval along line
dl_steps=1;

cond_model = [];

line_range = argin;

switch units
    case 'foot'
        
        meter_per_foot=3.280839;
                
    case 'meter'
        
        meter_per_foot=1.0;
                
    otherwise
        
        sprintf('Scripts takes units of "foot" or "meter"\nRevise inputs')
        return
        
end

% Write DC data to file
wrt2filedc_all = fopen([out_dir '\DC_line_3DOffline_UTM.obs'],'w');
fprintf(wrt2filedc_all,'!!MIRA - DCIP3D UTM BATCH\n');
fprintf(wrt2filedc_all,'\n\n');

wrt2fileip_all = fopen([out_dir '\IP_line_3DOffline_UTM.obs'],'w');
fprintf(wrt2fileip_all,'!!MIRA - DCIP3D UTM BATCH\n');
fprintf(wrt2fileip_all,'IPTYPE=1\n');
fprintf(wrt2fileip_all,'\n\n');
        
        
cd (work_dir);

files_in = ls;

% Extract observation file name for ip_ and dc_ from the current
% folder.
for file=1:size(files_in,1)-2;

    dc_obs_file= strtrim(files_in(file+2,:));


    %% Extract data from CSV file
    % Extract data from the file
    fid=fopen(dc_obs_file,'rt');
    line = fgets(fid);
    line = fgets(fid);
    
    countbck = 0;
    countfwr = 0;
    dataDCIP = zeros(1,8);
    
    % Write result to file
    wrt2filedc = fopen([out_dir '\dc_' dc_obs_file(1:end-4) '_UTM.obs'],'w');
    fprintf(wrt2filedc,'!!MIRA - DCIP3D UTM BATCH\n');
    fprintf(wrt2filedc,'\n\n');

    % Write result to file
    wrt2fileip = fopen([out_dir '\ip_' dc_obs_file(1:end-4) '_UTM.obs'],'w');
    fprintf(wrt2fileip,'!!MIRA - DCIP3D UTM BATCH\n');
    fprintf(wrt2fileip,'\n\n');

    while length(line)>1;
        
        data = regexp(line,',','split');
        
        % Tx Station, Line
        dataDCIP(1) = str2double(data{4});
        dataDCIP(2) = str2double(data{5});
        
        % P1 Station, Line
        dataDCIP(3) = str2double(data{8});
        dataDCIP(4) = str2double(data{9});
        
        % P2 Station, Line
        dataDCIP(5) = str2double(data{10});
        dataDCIP(6) = str2double(data{11});
        
        % DC and IP data
        dataDCIP(7) = str2double(data{12});
        dataDCIP(8) = str2double(data{14});
        
        data_UTM = zeros(1,6);
        % Find closest points in survey file and interpolated if needed
        for ii = 1 : 2 : 5
           
            line = survey( survey(:,5) == dataDCIP(ii+1) , :) ;
            stn = line(:,4) == dataDCIP(ii);
            
            % If survey file has the exact line and station
            if sum(stn)==1
                
                data_UTM(ii) = line(stn,1);
                data_UTM(ii+1) = line(stn,2);
                
            % If station is between two survey points    
            else
                
                % Compute distance to stations
                dstn = line(:,4) - dataDCIP(ii);
                line = [line abs(dstn)];
                line = sortrows(line,6);
                
                % Find two closest stations
                STN1 = line(1,:);
                STN2 = line(2,:);
                
                % Compute azimuth and interpolate
                dx = (STN1(1) - STN2(1));
                dy = (STN1(2) - STN2(2));
                               
                azm = atan( dx / dy );
                % Rotate azimuth in the right quadrant
                if  ( sign(dx)>0 && sign(dy)<0 ) || ( sign(dx)<0 && sign(dy)<0 )
                    
                    azm = pi + azm;
                    
                elseif sign(dx)<0 && sign(dy)>0
                    
                    azm = 2*pi + azm;
                    
                end
                
                dl = STN2(4) - dataDCIP(ii);
                
                data_UTM(ii) = STN2(1) + dl * meter_per_foot * sin(azm);
                data_UTM(ii+1) = STN2(2) + dl * meter_per_foot * cos(azm);
                
            end               
            
        end
        

        % DC DATA
        % Write Tx location
        fprintf(wrt2filedc,'%15.8e %15.8e %15.8e %15.8e %i\n',data_UTM(1),data_UTM(2),data_UTM(1),data_UTM(2),1);
        fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %i\n',data_UTM(1),data_UTM(2),data_UTM(1),data_UTM(2),1);


        % Write Rx data
        fprintf(wrt2filedc,'%15.8e %15.8e %15.8e %15.8e %15.8e\n\n',data_UTM(3),data_UTM(4),data_UTM(5),data_UTM(6),dataDCIP(7));
        fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %15.8e\n\n',data_UTM(3),data_UTM(4),data_UTM(5),data_UTM(6),dataDCIP(7));

        % IP DATA
        % Write Tx location
        fprintf(wrt2fileip,'%15.8e %15.8e %15.8e %15.8e %i\n',data_UTM(1),data_UTM(2),data_UTM(1),data_UTM(2),1);
        fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %i\n',data_UTM(1),data_UTM(2),data_UTM(1),data_UTM(2),1);

        % Write Rx data
        fprintf(wrt2fileip,'%15.8e %15.8e %15.8e %15.8e %15.8e\n\n',data_UTM(3),data_UTM(4),data_UTM(5),data_UTM(6),dataDCIP(8));
        fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %15.8e\n\n',data_UTM(3),data_UTM(4),data_UTM(5),data_UTM(6),dataDCIP(8));
        
        line = fgets(fid);
        
    end

    fclose(wrt2fileip);
    fclose(wrt2filedc);
    
end



fclose(wrt2filedc_all);
fclose(wrt2fileip_all);

% for jj=1:size(survey,1)-1
%     interpolate=[];
% 
%     % Check that we are still on the same line
%     % If yes, interpolate stations
%     if survey(jj,5) == survey(jj+1,5)    
%         %Compute length between points
%         L = survey(jj+1,4)-survey(jj,4);
%         dX = abs(survey(jj+1,2)-survey(jj,2));
%         dY = abs(survey(jj+1,1)-survey(jj,1));
%         dZ = survey(jj+1,3)-survey(jj,3);
% 
%         if dX==0 && dY==0
%             continue
%         end
%         azm(counter)=atan(dX/dY);
%         counter=counter+1;
%         %Compute number of stations to add
%         num_stations=round(L/dl_steps);
% 
%         dl = L / (num_stations);
%         dx = dX / (num_stations);
%         dy = dY / (num_stations);
%         dz = dZ / (num_stations);
% 
% %         interpolate=survey(jj,:);
%         interpolate(1:num_stations,5)=survey(jj,5);
%         interpolate(1:num_stations,4)=survey(jj,4) + dl * (0:num_stations-1)' ;
%         interpolate(1:num_stations,2)=survey(jj,2) + dx * (0:num_stations-1)' ;
%         interpolate(1:num_stations,1)=survey(jj,1) + dy * (0:num_stations-1)' ;
%         interpolate(1:num_stations,3)=survey(jj,3) + dz * (0:num_stations-1)' ;
% 
%         survey_interp=[survey_interp;interpolate];
%         clear interpolate
%     end
% 
% end
% survey_interp=[survey_interp;survey(end,:)];
% 
% %             if survey(jj,1) ~= survey(jj+1,1) || jj==(size(survey,1)-1)
% %                 if jj==(size(survey,1)-1)
% %                 interpolate=survey(jj+1,:);    
% %                 else
% %                 interpolate=survey(jj,:);
% %                 end
% %                 survey_interp=[survey_interp;interpolate];
% %                 azimuth = mean(azm);
% 
% %                 line_num = survey(jj,1);
% 
%         X=survey(:,2);
%         Y=survey(:,1);
%         Z=survey(:,3);
% 
% 
%         station_start = min(min(min(dc_data(:,1:4))),min(min(ip_data(:,1:4))));   %Smallest station measured at
%         station_in = min(survey(:,4));
%         Xin=X(1);
%         Yin=Y(1);
%         xdir=(X(end)-X(1)+1e-20)/abs(X(end)-X(1)+1e-20);
%         ydir=(Y(end)-Y(1)+1e-20)/abs(Y(end)-Y(1)+1e-20);
% 
%         % Extrapolate survey line beyond first survey station
%         while station_in > station_start
%             station_in = station_in - dl_steps;
%             Xin = Xin - xdir * dl_steps * meter_per_foot * abs(sin(mean(azm)));
%             Yin = Yin - ydir * dl_steps * meter_per_foot * abs(cos(mean(azm)));
% 
%             interpolate(1,5)= line_num;
%             interpolate(1,4)= station_in;
%             interpolate(1,2)= Xin;
%             interpolate(1,1)= Yin;
%             interpolate(1,3) = mean(Z) ;
% 
%             survey_interp=[survey_interp;interpolate];       
%         end
% 
% 
%         % Extrapolate survey line beyond last survey stations
%         station_in=max(survey(:,4));
%         station_stop=max(max(max(dc_data(:,1:4))),max(max(ip_data(:,1:4))));
%         Xin=X(end);
%         Yin=Y(end);  
%         xdir=(X(end)-X(1)+1e-20)/abs(X(end)-X(1)+1e-20);
%         ydir=(Y(end)-Y(1)+1e-20)/abs(Y(end)-Y(1)+1e-20);
%         while station_in < station_stop
%             station_in = station_in + dl_steps;
%             Xin = Xin + xdir * dl_steps * meter_per_foot * abs(sin(mean(azm)));
%             Yin = Yin + ydir * dl_steps * meter_per_foot * abs(cos(mean(azm)));
% 
%             interpolate(1,5)= line_num;
%             interpolate(1,4)= station_in;
%             interpolate(1,2)= Xin;
%             interpolate(1,1)= Yin;
%             interpolate(1,3) = mean(Z) ;
% 
%             survey_interp=[survey_interp;interpolate];       
%         end
% 
% 
%         clear azm Z
% 
% %             end
% 
% dc_data=sortrows(dc_data,1);
% survey = survey_interp;
% 
% 
% %Assign coordinates to DC data
% %line,station -> line,station,X,Y,Z
% wrt2file = fopen([dc_obs_file '_UTM.obs'],'w');
% fprintf(wrt2file,'!!MIRA - DCIP3D UTM BATCH\n');
% %         fprintf(wrt2file,[dc_obs_file '\n']);
% fprintf(wrt2file,'\n\n');
% 
% wrt2file_2D = fopen([dc_obs_file '_UTM2D.obs'],'w');
% fprintf(wrt2file_2D,'!!MIRA - DCIP2D UTM BATCH\n');
% fprintf(wrt2file_2D,'pole-dipole');
% fprintf(wrt2file_2D,'\n');
% 
% count =1;
% while count<size(dc_data,1)
% 
%     % Get all the dipoles from the same transmiter
%     num_pole = sum (dc_data(:,1)==dc_data(count,1));
% 
%     % Find coordinate of transmiter
%     Y = survey(dc_data(count,1)==survey(:,4),1);
%     X = survey(dc_data(count,1)==survey(:,4),2);
%     fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %i\n',X,Y,X,Y,num_pole);
%     fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %i\n',X,Y,X,Y,num_pole);
% 
%     for kk = 1:num_pole
% 
%         % Find coordinate of receivers
%         Y1 = survey(dc_data(count,3)==survey(:,4),1);
%         X1 = survey(dc_data(count,3)==survey(:,4),2);
% 
%         Y2 = survey(dc_data(count,4)==survey(:,4),1);
%         X2 = survey(dc_data(count,4)==survey(:,4),2);
% 
%         % Write 3D format
%         fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',X1,Y1,X2,Y2,dc_data(count,5),dc_data(count,6));
%         fprintf(wrt2filedc_all,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',X1,Y1,X2,Y2,dc_data(count,5),dc_data(count,6));
% 
%         % Write 2D format for DOI
%         fprintf(wrt2file_2D,'%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\n',X,X,X1,X2,dc_data(count,5),dc_data(count,6));
% 
%         count = count+1;
%     end
%     fprintf(wrt2file,'\n');
%     fprintf(wrt2filedc_all,'\n');
% end
% fclose(wrt2file);
% fclose(wrt2file_2D);
% 
% ip_data=sortrows(ip_data,1);
% %Assign coordinates to IP data
% %line,station -> line,station,X,Y,Z
% wrt2file = fopen([ip_obs_file '_UTM.obs'],'w');
% fprintf(wrt2file,'!!MIRA - DCIP3D UTM BATCH\n');
% fprintf(wrt2file,'IPTYPE=1\n');
% fprintf(wrt2file,'\n');
% 
% wrt2file_2D = fopen([ip_obs_file '_UTM2D.obs'],'w');
% fprintf(wrt2file_2D,'!!MIRA - DCIP2D UTM BATCH\n');
% fprintf(wrt2file_2D,'pole-dipole');
% fprintf(wrt2file_2D,'\n');
% 
% count =1;
% while count<size(ip_data,1)
% 
%     % Get all the dipoles from the same transmiter
%     num_pole = sum (ip_data(:,1)==ip_data(count,1));
% 
%     % Find coordinate of transmiter
%     Y = survey(ip_data(count,1)==survey(:,4),1);
%     X = survey(ip_data(count,1)==survey(:,4),2);
%     fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %i\n',X,Y,X,Y,num_pole);
%     fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %i\n',X,Y,X,Y,num_pole);
%     for kk = 1:num_pole
% 
%         % Find coordinate of receivers
%         Y1 = survey(ip_data(count,3)==survey(:,4),1);
%         X1 = survey(ip_data(count,3)==survey(:,4),2);
% 
%         Y2 = survey(ip_data(count,4)==survey(:,4),1);
%         X2 = survey(ip_data(count,4)==survey(:,4),2);
% 
%         % Write to file 3D
%         fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',X1,Y1,X2,Y2,ip_data(count,5),ip_data(count,6));
%         fprintf(wrt2fileip_all,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',X1,Y1,X2,Y2,ip_data(count,5),ip_data(count,6));
% 
%        % Write 2D format for DOI
%         fprintf(wrt2file_2D,'%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\n',X,X,X1,X2,ip_data(count,5),ip_data(count,6));
% 
%         count = count+1;
%     end
%     fprintf(wrt2file,'\n');
% 
% 
% end
% 
% end
% fclose(wrt2file);
% fclose(wrt2file_2D);

fprintf('**DONE GEO-REFERENCING**\n\n')