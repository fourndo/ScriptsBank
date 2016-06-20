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
    DCfwr_file = fopen([outdir '\dc_' line_name '_forward.dat'],'w');
    fprintf(DCfwr_file,'UBC-GIF DC Model\npole-dipole\n');
    
    DCbck_file = fopen([outdir '\dc_' line_name '_back.dat'],'w');
    fprintf(DCbck_file,'UBC-GIF DC Model\npole-dipole\n');
    
    IPfwr_file = fopen([outdir '\ip_' line_name '_forward.dat'],'w');
    fprintf(IPfwr_file,'UBC-GIF IP Model\npole-dipole\n');
    
    IPbck_file = fopen([outdir '\ip_' line_name '_back.dat'],'w');
    fprintf(IPbck_file,'UBC-GIF IP Model\npole-dipole\n');
    
    % Extract data from the file
    fid=fopen(datafile,'rt');
    line = fgets(fid);
    line = fgets(fid);
    
    countbck = 0;
    countfwr = 0;
    DCfwr = zeros(1,5);
    DCbck = zeros(1,5);
    IPfwr = zeros(1,5);
    IPbck = zeros(1,5);
    while length(line)>1;
        
        data = regexp(line,',','split');
        C1 = str2num(data{4});
        P1 = str2num(data{8});
        P2 = str2num(data{10});
        Vp = str2num(data{12});
        Ip = str2num(data{14});
        
        % Split data into forward and backward file
        if C1 > P1 && C1 > P2
            
            % Check if sign is good
            if C1 > P2 && Vp > 0
                
                Vp = Vp*-1;
                
            end
            
            countfwr = countfwr + 1;
            
            DCfwr(countfwr,:) = [C1 C1 P1 P2 Vp];
            IPfwr(countfwr,:) = [C1 C1 P1 P2 Ip];
            
            
            
        elseif C1 < P1 && C1 < P2
            
            % Check if sign is good
            if C1 < P2 && Vp < 0
                
                Vp = Vp*-1;
                
            end
            
            countbck = countbck + 1;
            
            DCbck(countbck,:) = [C1 C1 P1 P2 Vp];
            IPbck(countbck,:) = [C1 C1 P1 P2 Ip];
            
        end
        

        
        line = fgets(fid);
        
    end
    
    % Sort data and write to file
    DCbck = sortrows(DCbck,1);
    IPbck = sortrows(IPbck,1);
    DCfwr = sortrows(DCfwr,1);
    IPfwr = sortrows(IPfwr,1);

    for ii = 1 : countfwr

        fprintf(DCfwr_file,'%i %i %i %i %8.5e \n',DCfwr(ii,1),DCfwr(ii,2),DCfwr(ii,3),DCfwr(ii,4),DCfwr(ii,5));

        fprintf(IPfwr_file,'%i %i %i %i %8.5e \n',IPfwr(ii,1),IPfwr(ii,2),IPfwr(ii,3),IPfwr(ii,4),IPfwr(ii,5));

    end

    for ii = 1 : countbck

        fprintf(DCbck_file,'%i %i %i %i %8.5e \n',DCbck(ii,1),DCbck(ii,2),DCbck(ii,3),DCbck(ii,4),DCbck(ii,5));

        fprintf(IPbck_file,'%i %i %i %i %8.5e \n',IPbck(ii,1),IPbck(ii,2),IPbck(ii,3),IPbck(ii,4),IPbck(ii,5));

    end

    fclose(DCbck_file); fclose(IPbck_file); fclose(DCfwr_file); fclose(IPfwr_file);
        
    
    cd ..
        
end            

    
    

