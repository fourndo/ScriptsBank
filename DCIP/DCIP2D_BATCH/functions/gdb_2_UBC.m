function gdb_2_UBC(input_dir,output_dir,DC_file,IP_file)
% Extract files from the output of TQIP and create a file structure for the
% 2D Inversions. User need to specify the INPUT directory of the lines and
% an Ouput directory for the file structures. THe program cycle through all
% the folders and extract the data by seperating the forward and backward
% lines.
% 
% Author: D Fournier
% Last update : April 15th, 2013


%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>> for DEV only
% INPUT
% input_dir = 'C:\Projects\4081_Teck_HVC_DCIP3D\Data\Processing';
% DC_file = 'TECK_Data_gdb_DC.dat'; IP_File ='TECK_Data_gdb_IP.dat';
% OUTPUT
% output_dir ='C:\Projects\4081_Teck_HVC_DCIP3D\Data\Processing\MIRA_Inv2D';


%% Driver
dos (['rmdir  ' output_dir ' /S /Q']);
dos (['mkdir  ' output_dir]);

cd(input_dir)

% Extract DC and IP lines from ASCII files
DC_data = load(DC_file);
jj = 1;
ii =1;

Line.data = zeros(1,5);

while ii <= size(DC_data,1)
    
    Line.ID = num2str(abs(DC_data(ii,1)));
    Line.data(jj,:) = DC_data(ii,2:end);
    
    if ii < size(DC_data,1)
        
        if DC_data(ii,1) ~= DC_data(ii+1,1)
        
            fid=fopen([output_dir '\dc_L' Line.ID '_Teck.dat'],'w');
            fprintf(fid,'UBC-GIF DC Model\npole-dipole\n');
            
            for kk = 1:size(Line.data,1)
                
                % Check the sign of data
                if (Line.data(kk,2) > Line.data(kk,3)) && (Line.data(kk,3) > Line.data(kk,4)) && (Line.data(kk,5) < 0)
                    
                    Line.data(kk,5) = Line.data(kk,5) * -1;
                    
                elseif (Line.data(kk,2) > Line.data(kk,3)) && (Line.data(kk,3) < Line.data(kk,4)) && (Line.data(kk,5) > 0)
                    
                    Line.data(kk,5) = Line.data(kk,5) * -1;

                end
                
                fprintf(fid,'%8.5f %8.5f %8.5f %8.5f %8.5f\n',Line.data(kk,:));
            end
            
            fclose(fid);
            jj =0;
            Line.data = zeros(1,5);
            
        end
        
    end
    
    % Takes care of last line
    if ii == size(DC_data,1)
        
        fid=fopen([output_dir '\dc_L' Line.ID '_Teck.dat'],'w');
        fprintf(fid,'UBC-GIF DC Model\npole-dipole');

        for kk = 1:size(Line.data,1)
            
             % Check the sign of data
            if (Line.data(kk,2) > Line.data(kk,3)) && (Line.data(kk,3) > Line.data(kk,4)) && (Line.data(kk,5) < 0)

                Line.data(kk,5) = Line.data(kk,5) * -1;

            elseif (Line.data(kk,2) > Line.data(kk,3)) && (Line.data(kk,3) < Line.data(kk,4)) && (Line.data(kk,5) > 0)

                Line.data(kk,5) = Line.data(kk,5) * -1;

            end
             
            fprintf(fid,'%8.5f %8.5f %8.5f %8.5f %8.5f\n',Line.data(kk,:));
        end

        fclose(fid);
        jj =0;
        Line.data = zeros(1,5);
 
    end
    
    jj = jj + 1;

    ii = ii + 1;
    
end
    
    
IP_data = load(IP_file);
jj = 1;
ii =1;

Line.data = zeros(1,5);

while ii <= size(IP_data,1)
    
    Line.ID = num2str(abs(IP_data(ii,1)));
    Line.data(jj,:) = IP_data(ii,2:end);
    
    if ii < size(IP_data,1)
        
        if IP_data(ii,1) ~= IP_data(ii+1,1)
        
            fid=fopen([output_dir '\ip_L' Line.ID '_Teck.dat'],'w');
            fprintf(fid,'UBC-GIF IP Model\npole-dipole\n');
            
            for kk = 1:size(Line.data,1)
                
%                 % Check the sign of data
%                 if (Line.data(kk,2) > Line.data(kk,3)) && (Line.data(kk,3) > Line.data(kk,4)) && (Line.data(kk,5) < 0)
%                     
%                     Line.data(kk,5) = Line.data(kk,5) * -1;
%                     
%                 end
                
                fprintf(fid,'%8.5f %8.5f %8.5f %8.5f %8.5f\n',Line.data(kk,:));
            end
            
            fclose(fid);
            jj =0;
            Line.data = zeros(1,5);
            
        end
        
    end
    
    if ii == size(IP_data,1)
           
        fid=fopen([output_dir '\ip_L' Line.ID '_Teck.dat'],'w');
        fprintf(fid,'UBC-GIF DC Model\npole-dipole');

        for kk = 1:size(Line.data,1)
            fprintf(fid,'%8.5f %8.5f %8.5f %8.5f %8.5f\n',Line.data(kk,:));
        end

        fclose(fid);
        jj =0;
        Line.data = zeros(1,5);
 
    end
    
    jj = jj + 1;

    ii = ii + 1;
    
end

