function [ndata,phid,meshfile,obsfile,data,logfile,topofile,topo_check]=read_input


% addpath ([home_dir '\functions']);  
file_list=ls;

meshfile = [];
logfile = [];
topofile = [];
topo_check = [];
obsfile = [];
phid = [];

% Extract file names
for ii = 1:size(file_list,1)-2;

    look_at = strtrim(file_list(ii+2,:));

        if strcmp(look_at,'dcinv2d.log')==1

            logfile =look_at;

        elseif strcmp(look_at(end-2:end),'msh')==1

            meshfile = look_at;

        elseif length(look_at)>3 && strcmp(look_at(end-3:end),'topo')==1

            topofile =  look_at;

        elseif strcmp(look_at,'topo_model.txt')==1

            topo_check =  look_at;

        elseif strcmp(look_at(end-2:end),'obs')==1

            obsfile =  look_at;

        end

end

fid=fopen(obsfile,'rt');
line = fgets(fid);
line = fgets(fid);
line = fgets(fid);
count=1;
while line~=-1 
    data(count,:) = str2num(line);
    count = count +1;
    line = fgets(fid);
end
ndata = count-1;
fclose(fid);



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


    if length(strtrim(line))>=length('achieved misfit=')
        description = strtrim(line);
        if strcmp(description(1:16),'achieved misfit=')==1
            phid = str2num(description(17:end));
        end
    end

    if isempty(phid)==0
        break
    end

end
fclose(fid);


        