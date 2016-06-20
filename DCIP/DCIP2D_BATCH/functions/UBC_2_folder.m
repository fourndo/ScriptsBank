function UBC_2_folder(work_directory,output_folder,argin)
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
    

        line_name = ([file_name(4:7)]);
        line_name = regexprep(line_name,'[_]','');

    
    
    
%     if strcmp(file_name((end-7):(end-4)),'back')==1
        
    dos (['mkdir ' output_folder '\' line_name]);
    
    dos (['copy ' file_name ' ' output_folder '\' line_name '/Y']);
    
%     elseif strcmp(file_name((end-10):(end-4)),'forward')==1
%         
%     dos (['mkdir ' output_folder '\' line_name '\Forward']);
%     
%     dos (['copy ' file_name ' ' output_folder '\' line_name '\Forward']);
%     
%     end
     

   
   
end
    
cd(home_dir);
