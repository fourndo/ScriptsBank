% EM3DTM_batch
% Function prepares the directory and files required for the inversion of
% time domain EM survey lines.
clear all
close all

root_dir = pwd;
addpath functions

%% USER INPUTS
% Directory where the obs files are
data_dir = 'C:\LC\Private\dominiquef\Projects\4234_Anglo_Mosku_3DTM\Block14_data';

% Local Directory to send the line folders
work_dir = 'C:\LC\Private\dominiquef\Projects\4234_Anglo_Mosku_3DTM\Block14_Inv1';

% Copy in data_dir and remote_dir
topo_file = '../Sakatti_topo.topo';
wave_file = 'wave_2p5ms.txt';

% Remote directory to send the inversion
home_dir = '/home/dominiquef';
% remote_dir_MIRA = '/home/scottn/scott/AngloSakatti/3DOctree/Block8/Inv_dir';
remote_dir_AWS  = '../';

password = '3KzUjk3ghdGUZfZm9e8q';
ID_server   = 'dominiquef@173.231.120.107';
AWS_server = 'dominiquef@ec2-54-90-94-65.compute-1.amazonaws.com';

%% Octree mesh specs
% Min cell size [dx dy dz]
dx = [25 25 25];

% Maximum expansion distance [X Y Z Air]
expan_x = [3000 2500 1500 2300];

% Depth core
core = 400;

% Number of partition
parts = 1;

%% Inversion specs
% Initial model
mod_in = 0.01;

% Reference model
mod_ref = 0.01;

% beta_max  beta_min  beta_factor
beta = [1e-2  1.e-60  0.1];

% [ alpha_s  alpha_x  alpha_y  alpha_z]
alpha = [0.001 50 500 50];

% Target chi factor
target = 1e-4;

% Number of iteration per beta
iter_per_beta = 3;

%% SCRIPT STARTS HERE
%% Move observation files to directory
data_2_folder(data_dir,work_dir,'_','all')


%% Run the octree meshing scipt
run_octree_mesh(work_dir,data_dir,topo_file,dx,expan_x,core,parts,'all');


%% Write the input file for EM3DTD inversion
write_em3dtd_input(work_dir,remote_dir_AWS,wave_file,mod_in,mod_ref,beta,alpha,target,iter_per_beta,'all');

%% Compress inversion files and transfer to remote server
cd(work_dir)
% cd ..
% % system(['7za a ' work_dir '.7z ' work_dir]) 
% 
% fidread = fopen('Work_list.dat','r');
% count = 1;
% 
% line = strtrim(fgets(fidread));
% while line~=-1
%     
%     
%     
%     
%      
% 
%     fprintf('Start Transfering Inversion folder to the cluster: %s:%s\n', ID_server,home_dir)
% 
%     fid = fopen('Run.sh','w');
%     fprintf(fid,'#!/bin/bash\n');
%     fprintf(fid,['cd ' remote_dir_AWS '/' line '/ && qsub ' line '.pbs \n']);
%     fclose(fid);
% 
% %     system(['pscp -pw ' password ' Run.sh ' ID_server ':' home_dir])
% %     system(['plink -pw ' password ' ' ID_server ' bash Run.sh']); 
% 
% 
%     system(['pscp Run.sh ' AWS_server ':' home_dir])
%     system(['plink ' AWS_server ' bash Run.sh']);
% 
%     line = strtrim(fgets(fidread));
%     count = count + 1;
%     
% end
% fclose(fidread);

%% Restart inversion.
% Transfert mesh file and latest inversion to new output directory
cd(work_dir)
fid = fopen('Work_list.dat','r');
lines = strtrim(fgets(fid));

while lines~=-1
    

    
    
        cd ([data_dir lines])
    

        % Find the Observation and topo file
        file_list = ls;
        
        con_file    = [];
        
        for ii = 1:size(file_list,1)-2;
            
            look_at = strtrim(file_list(ii+2,:));
            
            if strcmp(look_at(end-3:end),'.con')==1

                con_file = look_at;
                iter = str2double(look_at(5:6));
                break
                
            end
            
        end
        

    
    if isempty(con_file)==0
        
        system(['move ' data_dir lines '\' con_file ' ' work_dir '\' lines]);
    
    end
    system(['move ' data_dir lines '\3D_core_mesh.txt ' work_dir '\' lines]);
    system(['move ' data_dir lines '\3D_mesh.txt ' work_dir '\' lines]);
    system(['move ' data_dir lines '\octree_mesh_000001.txt ' work_dir '\' lines]);
    system(['move ' data_dir lines '\octree_mesh_big.txt ' work_dir '\' lines]);
    system(['move ' data_dir lines '\active_cells.txt ' work_dir '\' lines]);
    lines = strtrim(fgets(fid));
    
end