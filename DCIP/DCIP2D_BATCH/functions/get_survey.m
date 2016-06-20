function [survey] = get_survey(work_directory,argin)
% Extract UTM station coordinates from survey files. 
% Read files with format:
% Line 20N,NAD83 z17T,,,,,,,,,

% Latitude,Longtitude,Northing UTM,Easting UTM,Elevation,Grid Northing,Grid Easting,Elevation,description,,
% 01,"46ø43'23.85""N",""80ø51'59.29""W",5174424,510204,316,5174424,510204,316,100L20N
% 02,"46ø43'23.82""N",""80ø51'58.07""W",5174423,510230,324,5174423,510230,324,125L20N
% ...
%
% Output:
% X_UTM Y_UTM Elevation Station Line
%  ...
%
% Author: D Fournier
% Last update : December 28th, 2013

home_dir = pwd;
fprintf('***START READING SURVEY FILES***\n')
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>
cd (work_directory)

%% Driver
Sline_list=ls;

line_range = argin;

survey = zeros(1,5);
switch line_range 
    case 'all'
       range = 1:size(Sline_list,1)-2;
       
    otherwise
        range = str2num(argin);

end

    count = 0;
% Cycle through all the survey lines
for oo=range

    
    survey_file = Sline_list(oo+2,:);
    
    fprintf('Extracting survey information from file: %s\n',survey_file);
    
    %% Read the survey file and extract information
    fid=fopen(survey_file,'rt');
    line=fgets(fid); %gets next line

    %Start looking at line 4
    for ii = 1:3
        
        line=fgets(fid);
        
    end
    
    
    while line~=-1         	         
              
        count=count+1;


        %Extract data
        data = regexp(line,',','split');

        survey(count,1:3)= [str2num(data{5}) str2num(data{4}) str2num(data{6})];

        % Extrat station number and line number
        line_ID = regexp(strtrim(data{10}),'L','split');
        stn_num = str2num(line_ID{1});
        
        line_num = str2num(line_ID{2}(1:end-1));
        
        % Put negative sign to differentiate between North/South or
        % East/West lines
        if strcmp(line_ID{2}(end),'S') || strcmp(line_ID{2}(end),'W')
            
            line_num = line_num * -1;
            
        end
        
        survey(count,4:5) = [stn_num line_num];

        line=fgets(fid);

    end
    
    

end

fclose(fid);
fprintf('***END OF get_survey***\n');
cd(home_dir);    
    
    
%         test_freq_match=strcmp(line(1:length(match_freq)),match_freq);