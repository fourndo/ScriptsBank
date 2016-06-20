function TQIP_2_folder(work_directory,output_folder,dc_pct,dc_flr,ip_pct,ip_flr,argin)
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
    

    line_name =  strtrim(DCline_list(DCline+2,:));
    line_name = (['L_' line_name(6:end)]);
    % Move to first direction
    cd (DCline_list(DCline+2,:))

    files_in = ls;
    
        % Extract observation file name for ip_ and dc_ from the current
        % folder.
        for file=1:size(files_in,1)-2;
            
            look_at = strtrim(files_in(file+2,:));
            
            if strcmp(look_at(end-5:end),'vp.obs')==1
                dc_file = importdata([look_at]);
                dc_obs = dc_file.data; 
                % Assign uncertainties
                dc_obs(:,6) = abs(dc_obs(:,5))*dc_pct + dc_flr;
                dc_file_name = (['dc_' line_name]) ;
            end
            
            if strcmp(look_at(end-5:end),'ip.obs')==1
                ip_file = importdata([look_at]);
                ip_obs = ip_file.data; 
                % Convert from mV/V to V/V
                ip_obs(:,5) = ip_obs(:,5) /1000;
                ip_obs(:,6) = abs(ip_obs(:,5))*ip_pct + ip_flr;
                % Assign uncertainties
                
                ip_file_name = (['ip_' line_name]);
            end
        end
    

%% Create new directory in the OUPUT folder and split the data in forward 
% and backward lines

    dos (['mkdir ' output_folder '\' line_name]);
fob=fopen([output_folder '\' line_name '\' dc_file_name '_edt.dat'],'w');
fprintf(fob,'UBC-GIF DC Model');
 for q=1:size(dc_obs,1)
     
     % Check if electrodes don't overlap
     if dc_obs(q,1)==dc_obs(q,3) || dc_obs(q,1)==dc_obs(q,4) || dc_obs(q,2)==dc_obs(q,3) || dc_obs(q,2)==dc_obs(q,4)
         
         fprintf('Obs removed - Electrodes overlap\n');
         continue
         
     end
%      if ((dc_obs(q,1) > dc_obs(q,3)))
        fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6e\t%8.6e ',dc_obs(q,1), dc_obs(q,2), dc_obs(q,3), dc_obs(q,4), dc_obs(q,5), dc_obs(q,6));
%      end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen([output_folder '\' line_name '\' ip_file_name '_edt.dat'],'w');
fprintf(fob,'UBC-GIF IP Model\nIPTYPE=1');
 for q=1:size(ip_obs,1)
     
        if ip_obs(q,1)==ip_obs(q,3) || ip_obs(q,1)==ip_obs(q,4) || ip_obs(q,2)==ip_obs(q,3) || ip_obs(q,2)==ip_obs(q,4)

            fprintf('Obs removed - Electrodes overlap\n');
            continue

        end
     
%     if ((ip_obs(q,1) > ip_obs(q,3)))
        fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6e\t%8.6e ',ip_obs(q,1), ip_obs(q,2), ip_obs(q,3), ip_obs(q,4), ip_obs(q,5), ip_obs(q,6));
%      end
 end
fclose('all');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dos (['mkdir ' output_folder '\' line_name '\Forward']);
% fob=fopen([output_folder '\' line_name '\Forward' '\dc_' line_name '_Forward_edt.dat'],'w');
% fprintf(fob,'UBC-GIF DC Model\npole-dipole');
%  for q=1:size(dc_obs,1)
%      if ((dc_obs(q,1) < dc_obs(q,3)))
%         fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',dc_obs(q,1), dc_obs(q,2), dc_obs(q,3), dc_obs(q,4), dc_obs(q,5));
%      end
%  end
% fclose('all');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fob=fopen([output_folder '\' line_name '\Forward' '\ip_' line_name '_Forward_edt.dat'],'w');
% fprintf(fob,'UBC-GIF IP Model\npole-dipole');
%  for q=1:size(ip_obs,1)
%     if ((ip_obs(q,1) < ip_obs(q,3)))
%         fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',ip_obs(q,1), ip_obs(q,2), ip_obs(q,3), ip_obs(q,4), ip_obs(q,5));
%      end
%  end
% fclose('all');

cd ..
   
end
    
cd(home_dir);
