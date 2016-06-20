% function CSV_2_UBC_FB(input_dir,outdir)
% Extract files from the output of TQIP and create a file structure for the
% 2D Inversions. User need to specify the INPUT directory of the lines and
% an Ouput directory for the file structures. THe program cycle through all
% the folders and extract the data by seperating the forward and backward
% lines.
% 
% Author: D Fournier
% Last update : December 24th, 2013 


%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>> for DEV only
clear all
close all

% INPUT
input_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\From_Client\2013_12_20_Titan\Results';
% OUTPUT
outdir ='C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\Titan\Raw';


%% Driver
dos (['rmdir  ' outdir ' /S /Q']);
dos (['mkdir  ' outdir]);

cd(input_dir)

% Extract DC and IP lines from CSV files from all folders
DC_lines = ls;


nlines=size(DC_lines,1)-3;

% Cycle through all the file in the INPUT folder
for oo=1:nlines
   
    cd(DC_lines(oo+3,:));
    line_name = strtrim(DC_lines(oo+3,:));
    
    files_in = ls;
    
    for ff = 1 : size(files_in,1)-2;
        
        lookat = strtrim(files_in(ff+2,:));
        
        % Search for csv file
        if strcmp(lookat(end-2:end),'csv')==1;
            
            datafile = lookat;
            break
            
        end
        
    end
    
    % Create UBC files for DC and IP, forward and backward
    DC_file = fopen([outdir '\dc_' line_name '.obs'],'w');
    fprintf(DC_file,'UBC-GIF DC Model\npole-dipole\n');
    
    IP_file = fopen([outdir '\ip_' line_name '.obs'],'w');
    fprintf(IP_file,'UBC-GIF IP Model\npole-dipole\n');
    
    % Extract data from the file
    fid=fopen(datafile,'rt');
    line = fgets(fid);
    line = fgets(fid);
    
    count = 0;

    DC = zeros(1,5);
    IP = zeros(1,5);

    while length(line)>1;
        
        data = regexp(line,',','split');
        C1 = str2num(data{4});
        P1 = str2num(data{8});
        P2 = str2num(data{10});
        Vp = str2num(data{12});
        Ip = str2num(data{14});
        
            
        count = count + 1;

        DC(count,:) = [C1 C1 P1 P2 Vp];
        IP(count,:) = [C1 C1 P1 P2 Ip];

        line = fgets(fid);
        
    end
    
    % Sort data and write to file
    DC = sortrows(DC,1);
    IP = sortrows(IP,1);

    for ii = 1 : count

        fprintf(DC_file,'%i %i %i %i %8.5e \n',DC(ii,1),DC(ii,2),DC(ii,3),DC(ii,4),DC(ii,5));

        fprintf(IP_file,'%i %i %i %i %8.5e \n',IP(ii,1),IP(ii,2),IP(ii,3),IP(ii,4),IP(ii,5));

    end

    fclose(DC_file); fclose(IP_file);
        
    
    cd ..
        
end            

    
    

