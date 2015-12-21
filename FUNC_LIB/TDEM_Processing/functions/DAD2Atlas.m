function [Atlas] = DAD2Atlas(work_dir,DHEM_raw)
% Assign an xyz coordinate to downhole data.

% For DEV only
% work_dir= 'C:\Projects\UBC_Lalor_DHEM\Data';
% load ([work_dir '\' 'Lalor_DHEM_RAW']);

% DHEM_UTM = DHEM_raw;
Atlas.survey = DHEM_raw.survey;
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
    
    survey_xyz = [survey_xyz;[survey(end,1:2) 1 survey(end,4) zeros(1,3)]];
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
       
   Atlas.survey{ii} = survey_xyz;
   
end

fclose(fid1);
fclose(fid2);
fclose(fid3);

% Merge XYZ and rotate data: Create new array vectors
%% PEM data in U,V,A converted to UTM
% A: Z-component along hole pointing towards collar
% U: X-component perpendicular to Z along vertical plane
% V: Y-component perp to AU, get it from RHandRule A x U
Rz = @(x)   [cosd(x) -sind(x) 0;
            sind(x) cosd(x) 0;
            0 0 1];
% Rx = @(x)   [1 0 0;
%             0 cosd(x) sind(x);
%             0 -sind(x) cosd(x)]; 
                    
Ry = @(x)   [cosd(x) 0 sind(x);
            0 1 0;
            -sind(x) 0 cosd(x)];

% rz = @(x) [cosd(x) -sind(x);
%             sind(x) cosd(x)] * [0;1];
%         
% Rv = @(x,ux,uy) [cosd(x) + ux^2*(1-cosd(x)) uy^2*ux^2*(1-cosd(x)) uy^2*sind(x);
%                  uy^2*ux^2*(1-cosd(x)) cosd(x) + uy^2*(1-cosd(x)) -ux^2*sind(x);
%                  -uy^2*sind(x) ux^2*sind(x) cosd(x)];
 %%
 
count = 1;

% Compute number total of stations
nstn = 0;
ntimes = 0;
for ii = 1 : size(DHEM_raw.data,1)
    
    nstn = nstn + size(DHEM_raw.data{ii}{1},1);
    ntimes = max([ntimes size(DHEM_raw.data{ii}{1},2)]);
    
end

% Memory allocation
Atlas.dBx = zeros(nstn,ntimes);
Atlas.dBy = zeros(nstn,ntimes);
Atlas.dBz = zeros(nstn,ntimes);
Atlas.DXYZ = zeros(nstn,4);
Atlas.Htag = zeros(nstn,1);

fid = fopen([work_dir '\' 'Lalor_Obs.loc'],'w');

for hh = 1 : size(DHEM_raw.data,1)    

    
    for dd = 1 : size(DHEM_raw.data{hh,1}{1},1)
        
        
        countdat = size(DHEM_raw.data{hh,1}{1},2)-1;
        
        target = DHEM_raw.data{hh,1}{1}(dd,1)==Atlas.survey{hh}(:,4);
        
        if sum(target)==0 || DHEM_raw.data{hh,1}{1}(dd,1)~=DHEM_raw.data{hh,1}{2}(dd,1)
            
            continue
            
        end
        
        % Build new data vector with three-components stacked
        temp = zeros(3,size(DHEM_raw.data{hh,1}{1},2)-1);
        
        temp(1,:) = DHEM_raw.data{hh,1}{1}(dd,2:end);
        temp(2,:) = DHEM_raw.data{hh,1}{2}(dd,2:end);
        
        % Find corresponding Z component,
        % if no match then linearly interpolate data
        target_z = DHEM_raw.data{hh,1}{1}(dd,1)==DHEM_raw.data{hh,1}{3}(:,1);
        
        if sum(target_z) == 0
            
            % Find closest stations
            dx = abs(DHEM_raw.data{hh,1}{1}(dd,1) - DHEM_raw.data{hh,1}{3}(:,1));
            
            [~,index] = sort(dx);
            temp(3,:) = mean( DHEM_raw.data{hh,1}{3}(index,2:end , 1) );
            
                     
        elseif sum(target_z) > 1
            
            target_z = find(target_z==1);
            target_z = target_z(1);
            temp(3,:) = DHEM_raw.data{hh,1}{3}(target_z,2:end);
            
        else
            
            temp(3,:) = DHEM_raw.data{hh,1}{3}(target_z,2:end);
            
        end
                        
        
        % Extract azimuth and dip for station
        % Convert from geographic to cartesian angle        
%         azm = mod(450-Atlas.survey{hh}(target,1),360) ;
        dip = mod(450-Atlas.survey{hh}(target,2),360);      
         azm = mod(450-Atlas.survey{hh}(target,1),360);
         
%          v_vec = rz(azm);
         
        % Create UTM, <i,j,k> for rotated data and normalized by current
        % and convert to T/s

%         Fx = eye(3);Fx(1,1)=-1;
%         Fy = eye(3);Fy(2,2)=-1;
        
        % Rotate data and convert units
%         temp =  Rz(azm) * (Rv(dip,v_vec(1),v_vec(2))* temp) / DHEM_raw.tx{hh,2}(1) * 1e-9;
        temp =  Ry(dip)* (Rz(-azm) *temp) / DHEM_raw.tx{hh,2}(1) * 1e-9;
        
        
        Atlas.dBx(count,1:countdat) = temp(1,:);
        Atlas.dBy(count,1:countdat) = temp(2,:);
        Atlas.dBz(count,1:countdat) = -temp(3,:);
        
        DXYZ = Atlas.survey{hh}(Atlas.survey{hh}(:,4) == DHEM_raw.data{hh,1}{1}(dd,1),4:7);
        
        Atlas.id{count,1} = DHEM_raw.name{hh};
        Atlas.DXYZ(count,:) = DXYZ;
        Atlas.Htag(count) = hh;
        
        fprintf(fid,'%12.8f %12.8f %12.8f %12.8f \n',DXYZ(2),DXYZ(3),DXYZ(4),hh);
        
        count = count + 1;
    end
    
%     figure(hh+1);hold off
    
end
fclose(fid);
% Collapse Atlas for no-data stations
nstn = nnz(Atlas.dBx(:,1));

Atlas.dBx = Atlas.dBx(1:nstn,:);
Atlas.dBy = Atlas.dBy(1:nstn,:);
Atlas.dBz = Atlas.dBz(1:nstn,:);
Atlas.DXYZ = Atlas.DXYZ(1:nstn,:);
Atlas.Htag = Atlas.Htag(1:nstn);