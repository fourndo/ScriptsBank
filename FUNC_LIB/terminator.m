% function terminator(work_directory,argin)
% Search and destroy files
% Author: D Fournier
% Last Update: January 3th, 2014

clear all
close all

% home_dir = pwd;
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>
cd ('C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\Titan\Processed_UBC2D')

%% Driver
DCline_list=ls;

range = 1:size(DCline_list,1)-2;

% Cycle through all the DC lines
for oo=range
    
    
    cd (DCline_list(oo+2,:))
    
    folder_list=ls;
    
    % Cycle through all the observations within each DC line
    for ff = 1:size(folder_list,1)-2
        
        cd (folder_list(ff+2,:))

        % Find the files to delete observation file
        file_list = ls;

        for ii = 1:size(file_list,1)-2;

            look_at = strtrim(file_list(ii+2,:));

            if (isempty(strfind(look_at,'.dat'))==0 || isempty(strfind(look_at,'.obs'))==0 ) && isempty(strfind(look_at,'edt'))==1

            dos(['del ' look_at]);

            end


        end

        cd ..
        
    end
    
    cd ..
    
end