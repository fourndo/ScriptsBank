function [ndata,beta_in,L_scale,phid,phi,meshfile,obs_file,logfile]=read_input


% addpath ([home_dir '\functions']);  
file_list=ls;

mesh = [];
logfile = [];
topo = [];
topo_model = [];
topo_check = [];
obs_file = [];
phid = [];
phi = [];

% Extract file names
for ii = 1:size(file_list,1)-2;

    look_at = strtrim(file_list(ii+2,:));

        if strcmp(look_at,'maginv3d.log')==1

            logfile =look_at;

        elseif strcmp(look_at(end-2:end),'msh')==1

            meshfile = look_at;

        elseif strcmp(look_at(end-3:end),'topo')==1

            topo =  look_at;

        elseif strcmp(look_at,'topo_model.txt')==1

            topo_check =  look_at;

        elseif strcmp(look_at(end-2:end),'obs')==1

            obs_file =  look_at;

        end

end


    if isempty(topo_check)==1

    % Run a topo check if not done already
    fprintf('Computing Topocheck - might take few minutes')
    dos (['topocheck ' meshfile ' ' topo])
    end




%% Read information from log file
fid=fopen(logfile,'rt');


max_num_lines = 30000;
% Go through the log file and extract data and the last achieved misfit
for ii=1:max_num_lines         	
line=fgets(fid); %gets next line 

    if line==-1
        fprintf('File ended at line %i\n',ii);
        fprintf('Did not find the information needed - review log file\n')
        break
    end

    if length(strtrim(line))>=length('# of data:')
        description = strtrim(line);
        if strcmp(description(1:10),'# of data:')==1
            ndata = str2num(description(11:end));
        end
    end

    if length(strtrim(line))>=length('multiplier:')
        description = strtrim(line);
        if strcmp(description(1:11),'multiplier:')==1
            beta_in = str2num(description(12:end));
        end
    end

    % Extract alpha values
    if length(strtrim(line))>=length('Le, Ln, Lz:')
        description = strtrim(line);
        if strcmp(description(1:11),'Le, Ln, Lz:')==1
            L_scale = str2num(description(12:end));
        end
    end

    if length(strtrim(line))>=length('data misfit:')
        description = strtrim(line);
        if strcmp(description(1:12),'data misfit:')==1
            phid = str2num(description(13:end));
        end
    end

    if length(strtrim(line))>=length('total objective:')
        description = strtrim(line);
        if strcmp(description(1:16),'total objective:')==1
            phi = str2num(description(17:end));
        end
    end

    if isempty(phid)==0 && isempty(phi)==0
        break
    end

end
fclose(fid);


        