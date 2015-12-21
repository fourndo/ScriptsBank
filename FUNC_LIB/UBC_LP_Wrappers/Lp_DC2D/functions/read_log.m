function [phid,iter]=read_log(logfile)
% Opens UBC log file and extract misfit and iteration number.

fid=fopen(logfile,'rt');

max_num_lines = 30000;
% Go through the log file and extract data and the last achieved misfit
for ii=1:max_num_lines         	
line=fgets(fid); %gets next line 

    if line==-1
        break
    end

    if length(strtrim(line))>=length('achieved misfit=')
        description = strtrim(line);
        if strcmp(description(1:16),'achieved misfit=')==1
            phid = str2num(description(17:end));
        end
    end

    if length(strtrim(line))>=length('Iteration  ')
        description = strtrim(line);
        if strcmp(description(1:11),'Iteration  ')==1
            iter = str2num(description(12:end));

        end
    end

end


fclose(fid);
