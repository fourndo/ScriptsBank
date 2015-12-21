function hole = Read_PEM(work_dir)

% Batch read PEM file
% Input: Working directory with subfolders containing PEM files
%
% Output: Structure array containing
% 1- name  : String array of hole names
% 2- collar: XYZ coordinated of holes at surface
% 3- survey: Along hole Azimuth, Dip and Distance
% 4- tx    : XYZ for transmitter location
% 5- data  : <STN#> dBx/dt for each time channels
%            <STN#> dBy/dt for each time channels
%            <STN#> dBz/dt for each time channels
% 6- time  : Bounds for time channels (n-1 channels)

% For DEV only
% work_dir='C:\Projects\UBC_Lalor_DHEM\Data\data_from_client\Phase2';
folders = ls(work_dir);

% Extract observation file in each folder
hole.name{1} = cellstr('empty');
hole.collar = [];
hole.survey = [];
hole.tx = [];
hole.data = [];
hole.time = [];
counth = 0;
for ii=3:size(folders,1);
    
    files = ls([work_dir '\' folders(ii,:)]);
    
    for jj = 3:size(files,1);
        
        file = strtrim(files(jj,:));
        
        % Skip file if not a *.PEM file
        if strcmp(file(end-2:end),'PEM')==0
            
            continue
            
        end
        
        temp = cellstr(strtrim(file(1:end-6)));
        
        fid=fopen([work_dir '\' folders(ii,:) '\' file],'rt');
        line=fgets(fid);
        new_id = 0;    
        % Extract hole name, collar and survey if first time encountered
        for kk = 1:length(hole.name)
            
            if strcmp(temp,hole.name{kk,1})==1
                
                new_id = 1;
                
            end
            
        end
        
        if new_id==0
            
            counth = counth+1;
            hole.name{counth,1} = temp;

            % Extract peak current
            while isempty(strfind(line,'<CUR>'))==1
                line=fgets(fid);
            end
            
            data = regexp(line,'\s+','split');
            hole.tx{counth,2} = str2num(data{2});
            
            % Move to beggining of Loop-coordinates
            while isempty(strfind(line,'Transmitter Loop'))==1
                line=fgets(fid);
            end
            
            % Extract Loop-coordinates
            line = fgets(fid);
            count = 1;
            while isempty(strfind(line,'Hole/Profile'))==1
               
                data = regexp(line,'\s+','split');
                hole.tx{counth,1}(count,1) = str2num(data{2});
                hole.tx{counth,1}(count,2) = str2num(data{3});
                hole.tx{counth,1}(count,3) = str2num(data{4});
                
                line = fgets(fid);
                count = count+1;
                
            end
            
            % Extract Collar and Survey
            line = fgets(fid);
            data = regexp(line,'\s+','split');
            hole.collar{counth,1}(1) = str2num(data{2});
            hole.collar{counth,1}(2) = str2num(data{3});
            hole.collar{counth,1}(3) = str2num(data{4});
            
            line = fgets(fid);
            count = 1;
            while isempty(strfind(line,'<P'))==0
                
                
                data = regexp(line,'\s+','split');
                hole.survey{counth,1}(count,1) = str2num(data{2});
                hole.survey{counth,1}(count,2) = str2num(data{3});
                hole.survey{counth,1}(count,3) = str2num(data{4});
                hole.survey{counth,1}(count,4) = str2num(data{6});
                
                line = fgets(fid);
                count = count+1;
                
            end
            
            % Extract Time channels
            line = fgets(fid);
            while isempty(strfind(line(1),'#'))==1
                line=fgets(fid);
            end
            line = fgets(fid);
            line = fgets(fid);
            time = [];  %Initialize empty string for time channels
            count = 1;
            while isempty(strfind(line,'$'))==1
                
                
                temp = str2num(line);

                time = [time temp];

                line=fgets(fid);
                
            end
            
            hole.time{counth,1} = time;
            
        end
        % Move to beggining of Data
        while isempty(strfind(line(1),'$'))==1
            line=fgets(fid);
        end

        % Read until the end of file

        % Count the number of receivers  
        countx = 0;
        county = 0;
        countz = 0;
        line=fgets(fid);

        while line~=-1
            data = [];

            header = regexp(line,'\s+','split');                    
%                 hole.data{counth}{comp} = 
            data(1) = str2num(header{1});

            if isempty(strfind(line,'XR'))==0

                comp = 1;

            elseif isempty(strfind(line,'YR'))==0

                comp = 2;

            elseif isempty(strfind(line,'ZR'))==0

                comp = 3;

             end

            % not sure of header info for now
            line=fgets(fid);

            % Extract data
            line=fgets(fid);
            while length(line) ~= 1

                temp = str2num(line);

                data = [data temp];

                line=fgets(fid);

            end

            % Check for which component and extract station along hole
            % Data format:  {X} <STN1> Time1, 2 , 3 ,... 
            %                   <STN2> Time1, 2 , 3 ,...
            %                   <....>
            %
            %               {Y} <STN> Time1, 2 , 3 ,... 
            %               {Z} <STN> Time1, 2 , 3 ,... 

            if comp==1

                countx = countx+1;
                hole.data{counth,1}{comp}(countx,:) = data;

            elseif comp==2

                county = county+1;
                hole.data{counth,1}{comp}(county,:) = data;

            elseif comp==3

                countz = countz+1;
                hole.data{counth,1}{comp}(countz,:) = data;

            end

            line=fgets(fid);

        end
            
        fclose(fid);
        
        
    end
    

    
end
 
end

