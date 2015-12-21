function [meshfile,obsfile]=fetch_input(home_dir)

cd(home_dir)

% addpath ([home_dir '\functions']);  
file_list=ls;

meshfile = [];
obsfile = [];

% Extract file names
for ii = 1:size(file_list,1)-2;

    look_at = strtrim(file_list(ii+2,:));

        if strcmp(look_at(end-2:end),'msh')==1

            meshfile = look_at;

        elseif strcmp(look_at(end-2:end),'obs')==1

            obsfile =  look_at;

        end

end



        