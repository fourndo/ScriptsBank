function [home_dir,UBC_dir,iter_start,chifact,pvec,qvec,lvec,cool_beta,iter_max] = read_input(input_file)
% [home_dir,UBC,dir,iter_start,chifact,pvec,qvec,lvec,cool_beta,iter_max] = read_input(input_file);
% Read input file for compact inversion
% 
% INPUT//
% File and directory
% 
% OUTPUT//
% home_dir  : Directory on local drive containing *.mtx, *.msh, *.topo files
% UBC_dir   : Directory for the uncontrained UBC results (models & *.log)
% iter_start: Starting UBC model #
% chifact   : Target chi-factor
% pvec      : Value(s) for the norm on gradient in range [0 2]
% qvec      : Value(s) for the norm on model in range [0 2]
% lvec      : Value(s) for the scaling between gradient and model norm in range [0 2]
% // If multiple values are used for any of p, q or l, the algorithm will
% invert for every combinaison. Results will be stored in seperate folders.
% cool_beta : Cooling schedule for the trade-off parameter (default 0.5)
% iter_max  : Maximum number of iteration

fid=fopen(input_file,'rt');

% Go through the input file and extract parameters
    
% First line: home_directory 
line=fgets(fid);
delm=regexp(line,'!');
home_dir = strtrim(line(1:delm-1));

% Second line: UBC_results directory 
line=fgets(fid);
delm=regexp(line,'!');
UBC_dir = strtrim(line(1:delm-1));

% Third line: Starting UBC model 
line=fgets(fid);
delm=regexp(line,'!');
iter_start = strtrim(line(1:delm-1));

% Fourth line: Target chi factor 
line=fgets(fid);
delm=regexp(line,'!');
chifact = str2num(strtrim(line(1:delm-1)));

% Fifth line: lp norm on gradient 
line=fgets(fid);
delm=regexp(line,'!');
pvec = str2num(strtrim(line(1:delm-1)));

% Sixth line: lq norm on gradient 
line=fgets(fid);
delm=regexp(line,'!');
qvec = str2num(strtrim(line(1:delm-1)));

% Seventh line: lambda scaling 
line=fgets(fid);
delm=regexp(line,'!');
lvec = str2num(strtrim(line(1:delm-1)));

% Eighth line: lp norm on gradient 
line=fgets(fid);
delm=regexp(line,'!');
cool_beta = str2num(strtrim(line(1:delm-1)));

% Nineth line: lp norm on gradient 
line=fgets(fid);
delm=regexp(line,'!');
iter_max = str2num(strtrim(line(1:delm-1)));