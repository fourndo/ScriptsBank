function [survey] = get_survey_v2(survey_file)
% Extract UTM station coordinates from survey files. 
% Read files with format:
% /      X        Y        STATION    filt       n=1       n=2       n=3       n=4
%/
%     line          0
%  497960.9   5178861      -975     10443      6297  
%  ...
%
% Output:
% X_UTM Y_UTM Elevation Station Line
%  ...
%
% Author: D Fournier
% Last update : December 21th, 2013

%% FOR DEV ONLY
% close all
% clear all
% 
% cd 'C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Data\Processing'
% survey_file ='IP-UTM.XYZ'

home_dir = pwd;
fprintf('***START READING SURVEY FILES***\n')
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>

    fid=fopen(survey_file,'rt');
    
    line=fgets(fid); %gets next line 
    
    % Move to beggining of data
    while isempty(strfind(strtrim(line),'Line'))==1
        
        line=fgets(fid);
        
    end
    
    count = 1;
    % Go through the log file and extract data and the last achieved misfit
    while line~=-1         	

        % If beggining of line, get line #       
        if isempty(strfind(strtrim(line),'Line'))==0 
            
           line_num = str2double(regexp(line,'[-+]?\d*','match'));
           
        elseif isempty(strfind(strtrim(line),'/'))==1
             %Extract data
            XYStn = str2num(line);
            survey(count,1:3)= [XYStn(1) XYStn(2) 1000];
            survey(count,4:5) = [XYStn(3) line_num];
            
            count = count + 1;
        end

        line=fgets(fid);
        
    end


fprintf('***END OF get_survey***\n')