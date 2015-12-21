% function [Atlas] = DAD2Atlas(work_dir,DHEM_raw)
% Assign an xyz coordinate to downhole data.

clear all
close all

%% USER INPUT
% Interpolation distance along hole 
dlmin = 1; % meter

% Set Survey and Collar file
work_dir= 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP';
fid2 = fopen([work_dir '\' 'BH_ddh_surveys.dat'],'r');
fid1 = fopen([work_dir '\' 'BH_ddh_collars.dat'],'r');

% Set transmitter location
Tx = [497367.000000 5178976.000000 403.232758;
        497441.000000 5178716.000000 410.825378;
        497037.000000 5178377.000000 398.071350;
        496738.000000 5178932.000000 397.051147];
    
    
%/////////////////////%
%% SCRIPT START HERE %%
%/////////////////////%

%% Load Collar and Survey information
line = fgets(fid1); % Skip header
line = fgets(fid1); % Read first line
count = 1;
while (line)~=-1
    
    data = regexp(line,'[\t ]','split');
    
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
    
    data = regexp(line,'[\t ]','split');
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
    nstn = floor(Atlas.survey{ii}(end,1) / dlmin) +1 ; 
    count = count + nstn ;
    
    % Pre-allocate space for interpolated survey
    survey_xyz = zeros(nstn,7);
    
    
    % For first station, replicate first survey   
    survey_xyz(1,1) = survey(1,1);
    survey_xyz(1,2) = survey(1,2);
    survey_xyz(1,3) = survey(1,3);
%     survey_xyz(1:nstn,4) = 0:(nstn-1);
%     count = count+nstn;
    
    % For all other stations, do a linear interpolation between
    for jj = 2:size(survey,1)
        
%         fprintf(fid2,'%s %12.4f %12.4f %12.4f\n',...
%         char(DHEM_raw.name{ii}),survey(jj,4),survey(jj,1),...
%         -survey(jj,2));
    
%         nstn = survey(jj,3);
        indx = 0:(nstn-1); indx = indx(:);
        
        % Interpolate azimuth, takes care of angles crossing 360 and 180
        if survey(jj,2)-survey(jj-1,2) > 180
            
            temp = survey(jj-1,2) -...
            (360 - survey(jj,2) + survey(jj-1,2))/nstn * indx;
        
            temp(temp<=0) = temp(temp<=0) +360;
            
            survey_xyz(1:nstn,2) = temp;
            
        elseif survey(jj,2)-survey(jj-1,2) < -180
            
            temp = survey(jj-1,2) +...
            (survey(jj,2) + 360 - survey(jj-1,2))/nstn * indx;
        
            temp(temp>=360) = temp(temp>=360) - 360;
            
            survey_xyz(1:nstn,2) = temp;
            
        else
            
            survey_xyz(1:nstn,2) = survey(jj-1,2) +...
                (survey(jj,2)-survey(jj-1,2))/nstn * indx;
            
        end
        
        % Interpolate depth
        survey_xyz(1:nstn,1) = survey(jj-1,1) +...
            (survey(jj,1)-survey(jj-1,1))/nstn * indx;
        
        % Interpolate dip
        survey_xyz(1:nstn,3) = survey(jj-1,3) +...
        (survey(jj,3)-survey(jj-1,3))/nstn * indx;

        % Add last survey station
        survey_xyz(nstn+1,1:3) = survey(end,:);
            
                
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
wrt2file = fopen([work_dir '\' 'STn_loc.dat'],'w');
fprintf(wrt2file,'X Y Z STN HoleID\n');
for jj = 1 : size(Atlas.survey,1)

        for kk = 1 : ( size(Atlas.survey{jj,2},1))
            
            fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %s\n',...
            Atlas.survey{jj,2}(kk,5),...
            Atlas.survey{jj,2}(kk,6),...
            Atlas.survey{jj,2}(kk,7),...
            kk,Atlas.holeid{jj,1});
            
            
        end        
        
end
%% Create DCIP3D Obs.loc file
wrt2file = fopen([work_dir '\' 'Obs_loc.dat'],'w');
fprintf(wrt2file,'!!MIRA - Obs loc for Forward modeling DCIP3D\n');
%         fprintf(wrt2file,[dc_obs_file '\n']);
fprintf(wrt2file,'\n\n');

for ii = 1 : size(Tx,1)/2
    
    txid = 2^(ii);
    % Write transmiter coordinates
    fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %i\n',...
        Tx(txid-1,1),Tx(txid-1,2),Tx(txid-1,3),Tx(txid,1),Tx(txid,2),Tx(txid,3),count);

    for jj = 1 : size(Atlas.survey,1)

        for kk = 1 : ( size(Atlas.survey{jj,2},1)-1)
            
            fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',...
            Atlas.survey{jj,2}(kk,5),...
            Atlas.survey{jj,2}(kk,6),...
            Atlas.survey{jj,2}(kk,7),...
            Atlas.survey{jj,2}(kk+1,5),...
            Atlas.survey{jj,2}(kk+1,6),...
            Atlas.survey{jj,2}(kk+1,7));
        
        end


    end
    
    fprintf(wrt2file,'\n');
    
end

fclose(wrt2file);