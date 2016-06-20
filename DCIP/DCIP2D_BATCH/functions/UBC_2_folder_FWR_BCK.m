function UBC_2_folder_FWR_BCK(work_directory,output_folder,argin)
% Extract files from the output of TQIP and create a file structure for the
% 2D Inversions. User need to specify the INPUT directory of the lines and
% an Ouput directory for the file structures. THe program cycle through all
% the folders and extract the data by seperating the forward and backward
% lines.
% 
% Author: D Fournier
% Last update : April 15th, 2013


%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>
% INPUT
% input_folder = 'C:\Projects\4001_Wallbridge_SkynnerLake_DCIP\Processing\Extracted_lines\TQIPdb_output';

% OUTPUT
% output_folder = 'C:\Projects\4001_Wallbridge_SkynnerLake_DCIP\Processing\MIRA_Inv2D\Inv_2D';

%% Driver

home_dir = pwd;

cd(work_directory)

DCline_list=ls;

nb_DC_lines=size(DCline_list,1)-2;

line_range = argin;

switch line_range 
    case 'all'
       range = 1:nb_DC_lines;
       
    otherwise
        range = str2num(argin);

end

for DCline=range
    

    file_name =  strtrim(DCline_list(DCline+2,:));
    
%     if strcmp(file_name(8:9),'-W')==1 || strcmp(file_name(8:9),'-E')
%         
%         line_name = ([file_name(4:9)]);
%     
%     else
%         line_name = ([file_name(4:7)]);
%         
%     end

    line_name = regexp(file_name,'_','split');
    line_name = line_name{2};
    
    
    if isempty(strfind(file_name,'back'))==0
                
        out_dir = [output_folder '\' line_name '\Back'];

        %         dos (['copy ' line_name ' ' bck_dir]);
    
    elseif isempty(strfind(file_name,'forward'))==0
        
        out_dir = [output_folder '\' line_name '\Forward'];
        
%         dos (['copy ' line_name ' ' fwr_dir]);
    
    end
    
    dos (['mkdir ' out_dir]);
        
    if strcmp(file_name(1:3),'dc_')==1
        
        dc_file = importdata(file_name);
        dc_obs = dc_file.data;
        dc_file_name = (['dc_' line_name]) ;
                
        fob=fopen([out_dir '\' dc_file_name '_edt.dat'],'w');
        fprintf(fob,'UBC-GIF DC Model\npole-dipole');
         for q=1:size(dc_obs,1)
        %      if ((dc_obs(q,1) > dc_obs(q,3)))
                fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6e ',dc_obs(q,1), dc_obs(q,2), dc_obs(q,3), dc_obs(q,4), dc_obs(q,5));
        %      end
         end
        fclose('all');


    end

    if strcmp(file_name(1:3),'ip_')==1
        ip_file = importdata(file_name);
        ip_obs = ip_file.data; ip_obs(:,5) = ip_obs(:,5) /1000;
        ip_file_name = (['ip_' line_name]);
        
        
        fob=fopen([out_dir '\' ip_file_name '_edt.dat'],'w');
        fprintf(fob,'UBC-GIF IP Model\npole-dipole');
         for q=1:size(ip_obs,1)
        %     if ((ip_obs(q,1) > ip_obs(q,3)))
                fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6e ',ip_obs(q,1), ip_obs(q,2), ip_obs(q,3), ip_obs(q,4), ip_obs(q,5));
        %      end
         end
         
        fclose('all');
        
    end

    


   
   
end
    
cd(home_dir);
