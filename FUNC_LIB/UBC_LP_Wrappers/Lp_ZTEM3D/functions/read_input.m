function [meshfile,obsfile,topofile]=read_input


% addpath ([home_dir '\functions']);  
file_list=ls;

meshfile = [];
topofile = [];
topo_check = [];
obsfile = [];


% Extract file names
for ii = 1:size(file_list,1)-2;

    look_at = strtrim(file_list(ii+2,:));

        if strcmp(look_at(end-2:end),'msh')==1

            meshfile = look_at;

        elseif length(look_at)>3 && strcmp(look_at(end-3:end),'topo')==1

            topofile =  look_at;

        elseif strcmp(look_at,'topo_model.txt')==1

            topo_check =  look_at;

        elseif strcmp(look_at(end-2:end),'obs')==1

            obsfile =  look_at;

        end

end


    if isempty(topo_check)==1

    % Run a topo check if not done already
    fprintf('Computing Topocheck - might take few minutes')
    dos (['topocheck ' meshfile ' ' topofile])
    end







        