function [DHEM_UTM] = DAD2UTM(work_dir,DHEM_raw)
% Assign an xyz coordinate to downhole data.

% For DEV only
% work_dir= 'C:\Projects\UBC_Lalor_DHEM\Data';
% load ([work_dir '\' 'Lalor_DHEM_RAW']);

DHEM_UTM = DHEM_raw;

% Write to file for Gocad view
fid1 = fopen([work_dir '\' 'Lalor_DHEM_Collar.dat'],'w');
fprintf(fid1,'HOLE-ID X Y Z LENGTH\n');

fid2 = fopen([work_dir '\' 'Lalor_DHEM_Survey.dat'],'w');
fprintf(fid1,'HOLE-ID DISTANCE AZIMUTH DIP\n');

fid3 = fopen([work_dir '\' 'Lalor_DHEM_XYZ.dat'],'w');
for ii = 1:size(DHEM_raw.collar,1)
    
    
    XYZ0 = DHEM_raw.collar{ii}(:);
    
    fprintf(fid1,'%s %12.4f %12.4f %12.4f %12.4f\n',...
        char(DHEM_raw.name{ii}),DHEM_raw.collar{ii}(1),DHEM_raw.collar{ii}(2),...
        DHEM_raw.collar{ii}(3),DHEM_raw.survey{ii}(end,4));
    
    survey = DHEM_raw.survey{ii}(:,1:4);
    
    % Pre-allocate space for interpolated survey
    survey_xyz = zeros(round(survey(end,4))-round(survey(1,4))+1,7);
    count = 0;
    
    % For first station, replicate first survey   
    nstn = survey(1,3);
    survey_xyz((1:nstn),1) = survey(1,1);
    survey_xyz((1:nstn),2) = survey(1,2);
    survey_xyz((1:nstn),3) = 1;
    survey_xyz(1:nstn,4) = 0:(nstn-1);
    count = count+nstn;
    
    % For all other stations, do a linear interpolation between
    for jj = 2:size(survey,1)
        
        fprintf(fid2,'%s %12.4f %12.4f %12.4f\n',...
        char(DHEM_raw.name{ii}),survey(jj,4),survey(jj,1),...
        -survey(jj,2));
    
        nstn = survey(jj,3);
        indx = 0:(nstn-1); indx = indx(:);
        
        % Deal with vec
        if survey(jj,1)-survey(jj-1,1) > 180
            
            temp = survey(jj-1,1) -...
            (360 - survey(jj,1) + survey(jj-1,1))/nstn * indx;
        
            temp(temp<=0) = temp(temp<=0) +360;
            
            survey_xyz(count+(1:nstn),1) = temp;
            
        elseif survey(jj,1)-survey(jj-1,1) < -180
            
            temp = survey(jj-1,1) +...
            (survey(jj,1) + 360 - survey(jj-1,1))/nstn * indx;
        
            temp(temp>=360) = temp(temp>=360) - 360;
            
            survey_xyz(count+(1:nstn),1) = temp;
            
        else
            
            survey_xyz(count+(1:nstn),1) = survey(jj-1,1) +...
                (survey(jj,1)-survey(jj-1,1))/nstn * indx;
            
        end
        survey_xyz(count+(1:nstn),2) = survey(jj-1,2) +...
            (survey(jj,2)-survey(jj-1,2))/nstn * indx;

        survey_xyz(count+(1:nstn),3) = 1;

        survey_xyz(count+(1:nstn),4) = survey(jj-1,4) + indx;

        count = count+nstn;
            
                
    end
    
    %Assign x,y,z values to stations
    
   survey_xyz(1,5:7) = XYZ0;
   
   for jj = 2:size(survey_xyz,1)
       

        
       azm = survey_xyz(jj,1);
       dip = survey_xyz(jj,2);
       
       dx = cosd(dip) * sind(azm);
       dy = cosd(dip) * cosd(azm);
       dz = -sind(dip);
       
       survey_xyz(jj,5:7) = survey_xyz(jj-1,5:7) + [dx dy dz];
       fprintf(fid3,'%12.4f %12.4f %12.4f\n',survey_xyz(jj-1,5:7));
           
   end
       
   DHEM_UTM.survey{ii} = survey_xyz;
   
end

fclose(fid1);
fclose(fid2);
fclose(fid3);

% Merge XYZ and rotate data: Create new array vectors
Rz = @(x)   [cosd(x) sind(x) 0;
            -sind(x) cosd(x) 0;
            0 0 1];
Rx = @(x)   [1 0 0;
            0 cosd(x) sind(x);
            0 -sind(x) cosd(x)]; 
                
Ry = @(x)   [cosd(x) 0 sind(x);
            0 1 0;
            -sind(x) 0 cosd(x)];

% fidx = fopen([work_dir '\' 'Lalor_DHEM_dbxdt.dat'],'w');
% fidy = fopen([work_dir '\' 'Lalor_DHEM_dbydt.dat'],'w');
% fidz = fopen([work_dir '\' 'Lalor_DHEM_dbzdt.dat'],'w');
counter = 0;        
for hh = 1 : size(DHEM_UTM.data,1)
    
    count = 1;
    
    for dd = 1 : size(DHEM_UTM.data{hh,1}{1},1)
        
        target = DHEM_UTM.data{hh,1}{1}(dd,1)==DHEM_UTM.survey{hh}(:,4);
        
        if sum(target)==0 || DHEM_UTM.data{hh,1}{1}(dd,1)~=DHEM_UTM.data{hh,1}{2}(dd,1)
            
            continue
            
        end
        
        % Find Z component, ignore data if no match with XY
        target_z = DHEM_UTM.data{hh,1}{1}(dd,1)==DHEM_UTM.data{hh,1}{3}(:,1);
        
        if sum(target_z) == 0
            
            continue             
                     
        elseif sum(target_z) > 1
            
            target_z = find(target_z==1);
            target_z = target_z(1);
            
        end
        
        counter = counter+1;        
        % Build new data vector with three-components stacked
%         DHEM_UTM.data{hh,2}{1}(count,1) = DHEM_UTM.survey{hh}(target,1);
%         DHEM_UTM.data{hh,2}{1}(count,2) = DHEM_UTM.survey{hh}(target,2);
        
        DHEM_UTM.data{hh,2}{1}(count,:) = DHEM_UTM.data{hh,1}{1}(dd,:);
        DHEM_UTM.data{hh,2}{1}(count+1,:) = DHEM_UTM.data{hh,1}{2}(dd,:);
        DHEM_UTM.data{hh,2}{1}(count+2,:) = DHEM_UTM.data{hh,1}{3}(target_z,:);
        
%         vector(count,1) = DHEM_UTM.survey{hh}(target,5);
%         vector(count,2) = DHEM_UTM.survey{hh}(target,6);
%         vector(count,3) = DHEM_UTM.survey{hh}(target,7);
        
        % Extract azimuth and dip for station
        % Substract 90d from azimuth for y reference
        % Substract 90 from dip for vertical z-vector
        
        azm = DHEM_UTM.survey{hh}(target,1)+90 ;
        dip = 90-(DHEM_UTM.survey{hh}(target,2));      
                  
        % Create UTM, <i,j,k> for rotated data and normalized by current
        % and convert to T/s
        DHEM_UTM.data{hh,2}{2}(count:count+2,:)=DHEM_UTM.data{hh,2}{1}(count:count+2,:);
        DHEM_UTM.data{hh,2}{2}(count:count+2,2:end) = Rz(azm)*Ry(dip) *...
                            DHEM_UTM.data{hh,2}{1}(count:count+2,2:end) / DHEM_UTM.tx{hh,2}(1) * 1e-9;

%         % Add heather first
%         if hh == 1
%             
%             fprintf(fidx,'HID DEPTH ');
%             fprintf(fidy,'HID DEPTH ');
%             fprintf(fidz,'HID DEPTH ');
%         
%             for ii = 1:size(DHEM_UTM.data{hh,2}{2},2)
%             
%             fprintf(fidx,'Time%i ',ii);
%             fprintf(fidy,'Time%i ',ii);
%             fprintf(fidz,'Time%i ',ii);
%         
%             end
%             
%         end
%         
%         fprintf(fidx,'\n%s ',DHEM_UTM.name{hh}{:});
%         fprintf(fidy,'\n%s ',DHEM_UTM.name{hh}{:});
%         fprintf(fidz,'\n%s ',DHEM_UTM.name{hh}{:});
%         
%         for ii = 1:size(DHEM_UTM.data{hh,2}{2},2)
%             
%             fprintf(fidx,'%12.5f ',DHEM_UTM.data{hh,2}{2}(count,ii));
%             fprintf(fidy,'%12.5f ',DHEM_UTM.data{hh,2}{2}(count+1,ii));
%             fprintf(fidz,'%12.5f ',DHEM_UTM.data{hh,2}{2}(count+2,ii));
%         
%         end
        
        
        count = count + 3;
    end
    
%     figure(hh+1);hold off
    
end

DHEM_UTM.nholes = hh;
DHEM_UTM.nstations = counter;
