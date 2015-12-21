function IP_Batch_Inv_FB(work_directory,argin)
% IP Inversion script for all DCIP lines
% >>>>>>(Must run the DC_Batch_Inv first)
% Cycles through all the folders, create the UBC input file and perform the
% inversion. Each folder must contain ONE observation file with the prefixe 
% dc_ and extension .dat (e.g. ip_#####.dat ).
% Observation files must be at the location of the inverted conductivity 
% model. (Not a concern if DC_Batch_Inv has been runned beforehand) 
% 
% INPUTS
% User only has to change the first line of this script to point to the
% HEAD directory. (e.g. HEAD -> all DCIP lines in seperate folders -> Forward and
% Backward observations in seperate sub-folders.)
%
%    Default UBC parameters used:
%    0 50  ! niter, irest
%    'NULL  ! chifact
%    obs_file name 
%    dcinv2d.con  ! conductivity file
%    dcinv2d.msh  ! mesh
%    NULL  ! topography
%    NULL  ! initial model
%    NULL  ! reference model
%    NULL  ! alpha
%    NULL  ! w.dat
%
% Author: D Fournier
% Last Update: April 10th, 2013

home_dir = pwd;
%% CHANGE DIRECTORY HERE AND RUN >>>>>>>>>>
cd (work_directory)

%% Driver
DCline_list=ls;
line_range = argin;

switch line_range 
    case 'all'
       range = 1:size(DCline_list,1)-2;
       
    otherwise
        range = str2num(argin);

end

% Cycle through all the DC lines
for oo=range
    
    
    cd (DCline_list(oo+2,:))
    
    folder_list=ls;
    
    % Cycle through all the observations within each DC line
    for ff = 1:size(folder_list,1)-2
        
        cd (folder_list(ff+2,:))

        % Find the IP observation file
        file_list = ls;
        
        ip_obs_file={};

            for ii = 1:size(file_list,1)-2;

                look_at = strtrim(file_list(ii+2,:));
            
                if strcmp(look_at(1:3),'ip_')==1 && strcmp(look_at(end-6:end),'edt.dat')==1
                    
                ip_obs_file=look_at;
                
                end


            end
            
            if isempty(ip_obs_file)==1

                fprintf(['Program could not find the IP file for line ' DCline_list(oo+2,:) folder_list(ff+2,:) '\n'])
                fprintf('Make sure that the file has the right format (ip_XXXX.dat)\n')
                cd ..
                break

            end
        
        % Write UBC control file
        fid=fopen('IP.inp','w');

        fprintf(fid,'0 50  ! niter, irest\n');
        fprintf(fid,'NULL  ! chifact\n');
        fprintf(fid,[(ip_obs_file) '\n']);
        fprintf(fid,'dcinv2d.con  ! conductivity file\n');
        fprintf(fid,'dcinv2d.msh  ! mesh\n');
        fprintf(fid,'NULL  ! topography\n');
        fprintf(fid,'NULL  ! initial model\n');
        fprintf(fid,'NULL  ! reference model\n');
        fprintf(fid,'NULL  ! alpha\n');
        fprintf(fid,'NULL  ! w.dat\n');

        fclose(fid);

        dos ('ipinv2d IP.inp')
        cd ..
        
    end
    
    cd ..
    
end

cd(home_dir);