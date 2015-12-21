% Add XYZ coordinates to line/stations data
clear all
close all

data=importdata('3893_Goldstrike_inter.dat');

survey=importdata('Goldstrike_merged.dat');

%Create new matrix (line,station,X,Y,Z,prop1,...)
data_geo=zeros(size(data,1),size(data,2)+3);
survey_interp=[];

%Densify survey files by interpolating between stations
%Set dl interval along line
dl_steps=50; 
meter_per_foot=0.3048;
counter=1;
for jj=1:size(survey,1)-1
    interpolate=[];

    if survey(jj,1) == survey(jj+1,1)    
        %Compute length between points
        L = survey(jj+1,2)-survey(jj,2);
        dX = survey(jj+1,3)-survey(jj,3);
        dY = survey(jj+1,4)-survey(jj,4);
        dZ = survey(jj+1,5)-survey(jj,5);
        
        azm(counter)=atan(dX/dY);
        counter=counter+1;
        %Compute number of stations to add
        num_stations=round(L/dl_steps);
        
        dl = L / (num_stations);
        dx = dX / (num_stations);
        dy = dY / (num_stations);
        dz = dZ / (num_stations);
        
%         interpolate=survey(jj,:);
        interpolate(1:num_stations,1)=survey(jj,1);
        interpolate(1:num_stations,2)=survey(jj,2) + dl * (0:num_stations-1)' ;
        interpolate(1:num_stations,3)=survey(jj,3) + dx * (0:num_stations-1)' ;
        interpolate(1:num_stations,4)=survey(jj,4) + dy * (0:num_stations-1)' ;
        interpolate(1:num_stations,5)=survey(jj,5) + dz * (0:num_stations-1)' ;
        
        survey_interp=[survey_interp;interpolate];
        clear interpolate
    end
    if survey(jj,1) ~= survey(jj+1,1) || jj==(size(survey,1)-1)
        if jj==(size(survey,1)-1)
        interpolate=survey(jj+1,:);    
        else
        interpolate=survey(jj,:);
        end
        survey_interp=[survey_interp;interpolate];
        azimuth = mean(azm);
        
        line_num = survey(jj,1);

        X=survey(survey(:,1)==line_num,3);
        Y=survey(survey(:,1)==line_num,4);
        Z=survey(survey(:,1)==line_num,5);
        
        
        station_start = min(data(data(:,1)==line_num,2));   %Smallest station measured at
        station_in = min(survey(survey(:,1)==line_num,2));
        Xin=X(1);
        Yin=Y(1);
        xdir=(X(end)-X(1))/abs(X(end)-X(1));
        ydir=(Y(end)-Y(1))/abs(Y(end)-Y(1));
        % Extrapolate survey line beyond first survey station
        while station_in > station_start
            station_in = station_in - dl_steps;
            Xin = Xin - xdir * dl_steps * meter_per_foot * abs(sin(azimuth));
            Yin = Yin - ydir * dl_steps * meter_per_foot * abs(cos(azimuth));
            
            interpolate(1,1)= line_num;
            interpolate(1,2)= station_in;
            interpolate(1,3)= Xin;
            interpolate(1,4)= Yin;
            interpolate(1,5) = mean(Z) ;
            
            survey_interp=[survey_interp;interpolate];       
        end
        
        
        % Extrapolate survey line beyond last survey stations
        station_in=survey(jj,2);
        station_stop=max(data(data(:,1)==line_num,2));
        Xin=X(end);
        Yin=Y(end);  
        xdir=(X(end)-X(1))/abs(X(end)-X(1));
        ydir=(Y(end)-Y(1))/abs(Y(end)-Y(1));
        while station_in < station_stop
            station_in = station_in + dl_steps;
            Xin = Xin + xdir * dl_steps * meter_per_foot * abs(sin(azimuth));
            Yin = Yin + ydir * dl_steps * meter_per_foot * abs(cos(azimuth));
            
            interpolate(1,1)= line_num;
            interpolate(1,2)= station_in;
            interpolate(1,3)= Xin;
            interpolate(1,4)= Yin;
            interpolate(1,5) = mean(Z) ;
            
            survey_interp=[survey_interp;interpolate];       
        end
        
        counter=1;
        clear azm Z
            
    end
end    

%Assign coordinates to data
%line,station -> line,station,X,Y,Z
for ii=1:size(data,1)
    
    data_geo(ii,1:2)=data(ii,1:2);
    data_geo(ii,6:14)=data(ii,3:11);
    
    [line_rows,xxx]=find(survey_interp==data(ii,1));
    [station,xxx]=find(survey_interp(line_rows,2)==data(ii,2));
    
    if isempty(station)==1
        continue
    end
    
    data_geo(ii,3:5)=survey_interp(line_rows(station),3:5);
    
end

figure; hist(data_geo(:,6),1e+6);set(gca,'xscale','log') 
title('\bfHistogram of measured frequencies')
ylabel('\bfN datum')
xlabel('\bflog(freq)')

save('Goldstrike_georef_inter.dat','data_geo','-ascii')