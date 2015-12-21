% function [Atlas] = DAD2Atlas(work_dir,DHEM_raw)
% Assign an xyz coordinate to downhole data.

clear all
close all

%% USER INPUT
% Interpolation distance along hole 
dlmin = 1; % meter
floor_DC = 1e-3;
floor_IP = 1e-6;
pct = 0.05;
f2m = 0.3048;
% Set Survey and Collar file
work_dir= 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\BrokenHammer';
data_file = 'BHIPDB_BrokenHammer.csv';
fid2 = fopen([work_dir '\' 'BH ddh surveys.csv'],'r');
fid1 = fopen([work_dir '\' 'BH ddh collars_v2.csv'],'r');

% Set transmitter location
Tx(1,:) = [497367.000000 5178976.000000 403.232758];
Tx(2,:) = [497441.000000 5178716.000000 410.825378];
Tx(3,:) = [497037.000000 5178377.000000 398.071350];
Tx(4,:) = [496738.000000 5178932.000000 397.051147];  
    
%/////////////////////%
%% SCRIPT START HERE %%
%/////////////////////%

%% Load Collar and Survey information
line = fgets(fid1); % Skip header
line = fgets(fid1); % Read first line
count = 1;
while (line)~=-1
    
    data = regexp(line,'[\t ,]','split');
    
    Atlas.holeid{count,1} = data{1};
    Atlas.holeid{count,2} = count;
    Atlas.holeid{count,3} = 0; % Flag used for importing survey info
    
    % Import collar information for hole [ East, North, Z, Length ]
    Atlas.collar{count,1}(1) = str2double(data{2});
    Atlas.collar{count,1}(2) = str2double(data{3});
    Atlas.collar{count,1}(3) = str2double(data{4});
%     Atlas.collar{count,1}(4) = str2double(data{5});
    
    %Next line
    count = count + 1;
    line = fgets(fid1);
end
  

% Load survey information
line = fgets(fid2);
line = fgets(fid2);
count = 1;
while (line)~=-1
    
    data = regexp(line,'[\t ,]','split');
    holeid = [];
    
    % Find hole collar id#, and check if already has a survey
    for ii = 1:size(Atlas.holeid,1)
        
        if strcmp(data{1},Atlas.holeid{ii,1}) == 1;
            
            holeid = Atlas.holeid{ii,2};
            
            if Atlas.holeid{ii,3}==0 % Beggining of survey
                
                countid = 1;
                Atlas.holeid{ii,3} = 1; % Change flag
                
            else
                
                countid = size(Atlas.survey{holeid,1},1) + 1;
                
            end
            
            break
            
        end
        
    end
    
    if isempty(holeid)
        
        fprintf('Program did not find a holeid match at line %i\n',count);
        line = fgets(fid2);
        continue
        
    else
    
        % Create survey array [ DEPTH, AZIMUTH, DIP ]
        Atlas.survey{holeid,1}(countid,1) = str2double(data{2});
        Atlas.survey{holeid,1}(countid,2) = str2double(data{3});
        Atlas.survey{holeid,1}(countid,3) = str2double(data{4});
    
    end
    
    %Next line
    count = count + 1;
    line = fgets(fid2);
    
end

% Write collar file for GOCAD import
wrt2coll = fopen([work_dir '\Collar_interp.dat'],'w');
fprintf(wrt2coll,'HOLEID X Y Z DEPTH\n');

wrt2surv = fopen([work_dir '\Survey_interp.dat'],'w');
    fprintf(wrt2surv,'HOLEID DEPTH AZIMUTH DIP\n');
    
%% Interpolate between survey points
count = 0;
for ii = 1:size(Atlas.survey,1)
    
    
    XYZ0 = Atlas.collar{ii}(1:3);
    depth = 0;
        
    survey = Atlas.survey{ii}(:,1:3);
    
    fprintf(wrt2coll,'%s %12.8e %12.8e %12.8e %12.8e\n',...
        Atlas.holeid{ii,1},Atlas.collar{ii}(1),...
        Atlas.collar{ii}(2),Atlas.collar{ii}(3),...
        Atlas.survey{ii}(end,1));
    
    % Set number of interpolated stations
    nstn = floor(Atlas.survey{ii}(end,1) / dlmin)  ; 
%     count = count + nstn ;
    
    % Pre-allocate space for interpolated survey
    survey_xyz = zeros(nstn,7);
    
    
    % For first station, replicate first survey   
    survey_xyz(1,1) = survey(1,1);
    survey_xyz(1,2) = survey(1,2);
    survey_xyz(1,3) = survey(1,3);
%     survey_xyz(1:nstn,4) = 0:(nstn-1);
%     count = count+nstn;
    
    % For all other stations, do a linear interpolation between
    count = 1;
    for jj = 2:size(survey,1)
        
%         fprintf(fid2,'%s %12.4f %12.4f %12.4f\n',...
%         char(DHEM_raw.name{ii}),survey(jj,4),survey(jj,1),...
%         -survey(jj,2));
    
%         nstn = survey(jj,3);
        nstn = floor( (survey(jj,1) - survey(jj-1,1) ) / dlmin) ;
        indx = 0:(nstn); indx = indx(:);
        
        % Interpolate azimuth, takes care of angles crossing 360 and 180
        if survey(jj,2)-survey(jj-1,2) > 180
            
            temp = survey(jj-1,2) -...
            (360 - survey(jj,2) + survey(jj-1,2))/nstn * indx;
        
            temp(temp<=0) = temp(temp<=0) +360;
            
            survey_xyz(count:count+nstn,2) = temp;
            
        elseif survey(jj,2)-survey(jj-1,2) < -180
            
            temp = survey(jj-1,2) +...
            (survey(jj,2) + 360 - survey(jj-1,2))/nstn * indx;
        
            temp(temp>=360) = temp(temp>=360) - 360;
            
            survey_xyz(count:count+nstn,2) = temp;
            
        else
            
            survey_xyz(count:count+nstn,2) = survey(jj-1,2) +...
                (survey(jj,2)-survey(jj-1,2))/nstn * indx;
            
        end
        
        % Interpolate depth
        survey_xyz(count:count+nstn,1) = survey(jj-1,1) +...
            (survey(jj,1)-survey(jj-1,1))/nstn * indx;
        
        % Interpolate dip
        survey_xyz(count:count+nstn,3) = survey(jj-1,3) +...
        (survey(jj,3)-survey(jj-1,3))/nstn * indx;

        % Add last survey station
        survey_xyz(count+nstn+1,1:3) = survey(end,:);
        
        count = count+nstn;
                
    end
    
    %Assign x,y,z values to stations   
    survey_xyz(1,5:7) = Atlas.collar{ii}(1:3);
   

    
    
    
   for jj = 2:size(survey_xyz,1)
       
       dl = survey_xyz(jj,1) - survey_xyz(jj-1,1);
        
       azm = survey_xyz(jj,2);
       dip = survey_xyz(jj,3);
       
       dx = dl * cosd(dip) * sind(azm);
       dy = dl * cosd(dip) * cosd(azm);
       dz = dl * sind(dip);
       
       survey_xyz(jj,4) = jj-1;
       
       survey_xyz(jj,5:7) = survey_xyz(jj-1,5:7) + [dx dy dz];
%        fprintf(fid3,'%12.4f %12.4f %12.4f\n',survey_xyz(jj-1,5:7));
       fprintf(wrt2surv,'%s %12.8e %12.8e %12.8e\n',...
        Atlas.holeid{ii,1},survey_xyz(jj,1)-survey_xyz(1,1),...
        azm,dip);   
   end
       
   Atlas.survey{ii,2} = survey_xyz;
   
end

fclose(wrt2coll);
fclose(wrt2surv);
%% Write station location to ASCII file
% wrt2file = fopen([work_dir '\' 'STn_loc.dat'],'w');
% fprintf(wrt2file,'X Y Z STN HoleID\n');
% for jj = 1 : size(Atlas.survey,1)
% 
%         for kk = 1 : ( size(Atlas.survey{jj,2},1))
%             
%             fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %s\n',...
%             Atlas.survey{jj,2}(kk,5),...
%             Atlas.survey{jj,2}(kk,6),...
%             Atlas.survey{jj,2}(kk,7),...
%             kk,Atlas.holeid{jj,1});
%             
%             
%         end        
%         
% end

%% Load data and re-assign HOLE-STATION to XYZ

fid = fopen([work_dir '\' data_file],'r');

line = fgets(fid); % Read first line

while strcmp(line(1),'/')==1;
    
    line = fgets(fid);
    
end

count = 1;
holename = ' ';
P = zeros(1,4);
data_DC_UTM = zeros(1,14);
data_IP_UTM = zeros(1,14);
while (line)~=-1
    
    temp = regexp(line,'[\t ,]','split');
    
    % Get code
    array_code = str2double(temp{2});
    
    % Get pole stations (C1,C2,P1,P2) and modify according to array code
    P(count,1) = str2double(temp{3});
    P(count,2) = str2double(temp{4});
    P(count,3) = str2double(temp{5});
    P(count,4) = str2double(temp{6});
    
    if array_code==3 || array_code==4 || array_code==5 || array_code==6
        
        P(count,1) = P(count,3)+P(count,1);
        P(count,2) = P(count,1);
        
    end
        
    % Find the hole-ID
    % Find hole collar id# if different from previous
    if strcmp(holename,temp{1})==0
        for ii = 1:size(Atlas.holeid,1)

            if strcmp(temp{1},Atlas.holeid{ii,1}) == 1;

                holeid = Atlas.holeid{ii,2};

            end

        end
        
    end
    
    % Find station for each pole location
    for jj = 1 : 4
    
        % If a remote transmitter, assign value
        % 1: 10000 , 2:20000 , 3:30000 , 4:40000
        if P(count,jj)==10000

            data_DC_UTM(count,3*(jj-1)+(1:3)) = Tx(1,1:3);

        elseif P(count,jj)==20000

            data_DC_UTM(count,3*(jj-1)+(1:3)) = Tx(2,1:3);

        elseif P(count,jj)==30000

            data_DC_UTM(count,3*(jj-1)+(1:3)) = Tx(3,1:3);

        elseif P(count,jj)==40000

            data_DC_UTM(count,3*(jj-1)+(1:3)) = Tx(4,1:3);
            
        else
        
        % Else find the station from interpolated survey
        index = Atlas.survey{holeid,2}(:,4) == round(P(count,jj)*f2m);
        
        % Skip data if station mismatch - too bad
        if sum(index)~=1
            
            break
            
        end
        
        data_DC_UTM(count,3*(jj-1)+(1:3)) = Atlas.survey{holeid,2}(index,5:7);
        data_IP_UTM(count,3*(jj-1)+(1:3)) = Atlas.survey{holeid,2}(index,5:7);
        end
        % Assign data Vp/I
        data_DC_UTM(count,13) = str2double(temp{10});
        
        % Secondary Potential  SP/I
        data_IP_UTM(count,13) = str2double(temp{11}) / str2double(temp{8}); 
    
        % Assign error
        data_DC_UTM(count,14) = pct * abs(data_DC_UTM(count,13)) + floor_DC;
        data_IP_UTM(count,14) = pct * abs(data_IP_UTM(count,13)) + floor_IP;
        
    end
        
    % Skip data if station mismatch - too bad
    if sum(index)~=1

        fprintf('No station found in survey for hole %s : Station:%i\n',Atlas.holeid{holeid},P(jj));
        line = fgets(fid);
        continue

    end
        
        line = fgets(fid);
        
        count = count+1;
        
end
        
    
    
%% Create DCIP3D Obs.loc file
wrt2fileDC = fopen([work_dir '\' '4180_DC_obs.dat'],'w');
fprintf(wrt2fileDC,'!!MIRA - DC  Data for DHIP\n');
%         fprintf(wrt2file,[dc_obs_file '\n']);
fprintf(wrt2fileDC,'\n\n');

wrt2fileIP = fopen([work_dir '\' '4180_IP_obs.dat'],'w');
fprintf(wrt2fileIP,'!!MIRA - IP Data (Secondary Potential for DHIP\n');
%         fprintf(wrt2file,[dc_obs_file '\n']);
fprintf(wrt2fileIP,'IPTYPE=2\n\n');

wrt2GOCAD = fopen([work_dir '\' '4180_GOCAD_pseudoXYZ.dat'],'w');
fprintf(wrt2GOCAD,'X Y Z Vp StdDev\n');


for ii = 1 : size(data_DC_UTM,1)
    
    % Compute half way point for pseudo plot in GOCAD
    rxhalfXYZ = ( data_DC_UTM(ii,4:6) + data_DC_UTM(ii,7:9))/2;
    pseudoXYZ = ( rxhalfXYZ +  data_DC_UTM(ii,1:3) ) / 2;
    
    fprintf(wrt2GOCAD,'%12.2f %12.2f %12.2f %12.2f  %12.2f\n',...
        rxhalfXYZ(1),rxhalfXYZ(2),rxhalfXYZ(3),data_DC_UTM(ii,13),data_DC_UTM(ii,14));
    
    
    % Write location + data to file
    for jj = 1 : size(data_DC_UTM,2)
    
        fprintf(wrt2fileDC,'%15.8e ',data_DC_UTM(ii,jj));
        fprintf(wrt2fileIP,'%15.8e ',data_DC_UTM(ii,jj));
        
        if jj == 6
            
        fprintf(wrt2fileDC,'%i\n',1);
        fprintf(wrt2fileIP,'%i\n',1);
        
        end
        
    end
    
    fprintf(wrt2fileDC,'\n');
    fprintf(wrt2fileIP,'\n');
    
end

fclose(wrt2fileIP);
fclose(wrt2fileDC);
fclose(wrt2GOCAD);
