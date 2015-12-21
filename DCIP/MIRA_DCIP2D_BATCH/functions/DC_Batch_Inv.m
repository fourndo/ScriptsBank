function DC_Batch_Inv(work_directory,argin)
% DC Inversion 
% Batch script for all DCIP lines
% Cycles through all the folders, create the UBC input file and perform the
% inversion. Each folder must contain ONE observation file with the prefixe 
% dc_ and extension .dat (e.g. ip_#####.dat ).
% 
% INPUTS
% User only has to change the first line of this script to point to the
% HEAD directory. e.g. HEAD -> all DCIP lines in seperate folders -> Forward and
% backward observations in seperate sub-folders.
%
%    Default UBC parameters used:
%    0 60  ! niter, irest
%    'NULL  ! chifact
%    obs_file name 
%    NULL  ! mesh
%    NULL  ! topography
%    NULL  ! initial model
%    NULL  ! reference model
%    NULL  ! alpha
%    NULL  ! w.dat
%
% Author: D Fournier
% Last Update: April 10th, 2013

home_dir = pwd;
%% CHANGE DIRECTORY AND RUN >>>>>>>>>>
cd(work_directory);

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
    

        % Find the IP observation file
        file_list = ls;
        
        dc_obs_file={};

            for ii = 1:size(file_list,1)-2;

                look_at = strtrim(file_list(ii+2,:));
            
                if strcmp(look_at(1:3),'dc_')==1 && strcmp(look_at(end-3:end),'.dat')==1
                    
                dc_obs_file=look_at;
                
                end


            end
            
        if isempty(dc_obs_file)==1
            
            fprintf(['Program could not find the DC file for line ' DCline_list(oo+2,:) '\n'])
            fprintf('Make sure that the file has the right format (dc_XXXX.dat)\n')
            cd ..
            break
            
        end
        
        % Write UBC control file
        fid=fopen('DC.inp','w');

%         fprintf(fid,'0 60  ! niter, irest\n');
%         fprintf(fid,'NULL  ! chifact\n');
%         fprintf(fid,[(dc_obs_file) '\n']);
%         fprintf(fid,'NULL  ! mesh\n');
%         fprintf(fid,'NULL  ! topography\n');
%         fprintf(fid,'NULL  ! initial model\n');
%         fprintf(fid,'NULL  ! reference model\n');
%         fprintf(fid,'NULL  ! alpha\n');
%         fprintf(fid,'NULL  ! w.dat\n');

        fprintf(fid,['OBS LOC_X ' dc_obs_file '\n']);
        fprintf(fid,'TOPO DEFAULT\n');
        fprintf(fid,'MESH DEFAULT\n');
        fprintf(fid,'ALPHA VALUE 1.e-3 1.0 1.0\n');
        fprintf(fid,'CHIFACT DEFAULT\n');
        fprintf(fid,'INIT_MOD DEFAULT\n');
        fprintf(fid,'REF_MOD DEFAULT\n');
        fprintf(fid,'INVMODE SVD\n');

        fclose(fid);

        dos ('dcinv2d DC.inp')
        cd ..
        
    
    
end

cd(home_dir);